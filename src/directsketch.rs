use anyhow::{Result, Context, anyhow, bail};
use camino::Utf8PathBuf as PathBuf;
use csv::Reader;
use regex::Regex;
use reqwest::Client;
use std::{collections::HashSet, fs};
use std::fs::{File, create_dir_all};
use csv::Writer;
use niffler::get_reader;
use needletail::parse_fastx_reader;
use std::io::Cursor;
use std::hash::Hash;
use std::hash::Hasher;
use sourmash::signature::Signature;
use sourmash::cmd::ComputeParameters;
use tokio::task;


#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Params {
    pub ksize: u32,
    pub track_abundance: bool,
    pub num: u32,
    pub scaled: u64,
    pub seed: u32,
    pub is_protein: bool,
    pub is_dna: bool,
}

impl Hash for Params {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.ksize.hash(state);
        self.track_abundance.hash(state);
        self.num.hash(state);
        self.scaled.hash(state);
        self.seed.hash(state);
        self.is_protein.hash(state);
        self.is_dna.hash(state);
    }
}

fn parse_params_str(params_strs: String) -> Result<Vec<Params>, String> {
    let mut unique_params: std::collections::HashSet<Params> = std::collections::HashSet::new();

    // split params_strs by _ and iterate over each param
    for p_str in params_strs.split('_').collect::<Vec<&str>>().iter() {
        let items: Vec<&str> = p_str.split(',').collect();

        let mut ksizes = Vec::new();
        let mut track_abundance = false;
        let mut num = 0;
        let mut scaled = 1000;
        let mut seed = 42;
        let mut is_protein = false;
        let mut is_dna = true;

        for item in items.iter() {
            match *item {
                _ if item.starts_with("k=") => {
                    let k_value = item[2..]
                        .parse()
                        .map_err(|_| format!("cannot parse k='{}' as a number", &item[2..]))?;
                    ksizes.push(k_value);
                }
                "abund" => track_abundance = true,
                "noabund" => track_abundance = false,
                _ if item.starts_with("num=") => {
                    num = item[4..]
                        .parse()
                        .map_err(|_| format!("cannot parse num='{}' as a number", &item[4..]))?;
                }
                _ if item.starts_with("scaled=") => {
                    scaled = item[7..]
                        .parse()
                        .map_err(|_| format!("cannot parse scaled='{}' as a number", &item[7..]))?;
                }
                _ if item.starts_with("seed=") => {
                    seed = item[5..]
                        .parse()
                        .map_err(|_| format!("cannot parse seed='{}' as a number", &item[5..]))?;
                }
                "protein" => {
                    is_protein = true;
                    is_dna = false;
                }
                "dna" => {
                    is_protein = false;
                    is_dna = true;
                }
                _ => return Err(format!("unknown component '{}' in params string", item)),
            }
        }

        for &k in &ksizes {
            let param = Params {
                ksize: k,
                track_abundance,
                num,
                scaled,
                seed,
                is_protein,
                is_dna,
            };
            unique_params.insert(param);
        }
    }

    Ok(unique_params.into_iter().collect())
}

fn build_siginfo(params: &[Params], moltype: &str) -> Vec<Signature> {
    let mut sigs = Vec::new();

    for param in params.iter().cloned() {
        match moltype {
            // if dna, only build dna sigs. if protein, only build protein sigs
            "dna" | "DNA" if !param.is_dna => continue,
            "protein" if !param.is_protein => continue,
            _ => (),
        }

        // Adjust ksize value based on the is_protein flag
        let adjusted_ksize = if param.is_protein {
            param.ksize * 3
        } else {
            param.ksize
        };

        let cp = ComputeParameters::builder()
            .ksizes(vec![adjusted_ksize])
            .scaled(param.scaled)
            .protein(param.is_protein)
            .dna(param.is_dna)
            .num_hashes(param.num)
            .track_abundance(param.track_abundance)
            .build();

        let sig = Signature::from_params(&cp);
        sigs.push(sig);
    }

    sigs
}

async fn fetch_genbank_filename(client: &Client, accession: &str) -> Result<(String, String)> {
    let (db, acc) = accession.trim().split_once('_').ok_or_else(|| anyhow!("Invalid accession format"))?;
    let (number, _) = acc.split_once('.').unwrap_or((acc, "1"));
    let number_path = number.chars().collect::<Vec<_>>().chunks(3).map(|chunk| chunk.iter().collect::<String>()).collect::<Vec<_>>().join("/");

    let base_url = format!("https://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}", db, number_path);
    let directory_response = client.get(&base_url).send().await?;
    if !directory_response.status().is_success() {
        return Err(anyhow!("Failed to open genome directory: HTTP {}", directory_response.status()));
    }
    let text = directory_response.text().await?;
    let link_regex = Regex::new(r#"<a href="([^"]*)""#)?;

    for cap in link_regex.captures_iter(&text) {
        let name = &cap[1];
        // Check if name ends with '/', and remove it if so
        let clean_name = if name.ends_with('/') { name.strip_suffix('/').unwrap() } else { name };
        if clean_name.starts_with(db) && clean_name.split('_').nth(1).map_or(false, |x| x.starts_with(number)) {
            // Formulate the correct URL by ensuring no double slashes
            return Ok((format!("{}/{}", base_url, clean_name), clean_name.into()));
        }
    }

    Err(anyhow!("No matching genome found for accession {}", accession))
}

// download and return data directly instead of saving to file
async fn download_with_retry(client: &Client, url: &str, retry_count: u32) -> Result<Vec<u8>> {
    let mut attempts = retry_count;
    while attempts > 0 {
        let response = client.get(url).send().await;
        match response {
            Ok(resp) if resp.status().is_success() => {
                let data = resp.bytes().await.context("Failed to read bytes from response")?;
                return Ok(data.to_vec());  // Return the downloaded data as Vec<u8>
            },
            _ => {
                eprintln!("Failed to download file: {}. Retrying...", url);
                attempts -= 1;
            }
        }
    }

    Err(anyhow!("Failed to download file after {} retries: {}", retry_count, url))
}

enum FileType {
    Genomic,
    Protein,
    AssemblyReport,
    Checksum,
}

impl FileType {
    fn suffix(&self) -> &'static str {
        match self {
            FileType::Genomic => "_genomic.fna.gz",
            FileType::Protein => "_protein.faa.gz",
            FileType::AssemblyReport => "_assembly_report.txt",
            FileType::Checksum => "md5checksums.txt",
        }
    }

    fn filename(&self, accession: &str) -> String {
        match self {
            FileType::Checksum => format!("{}_{}", accession, self.suffix()),
            _ => format!("{}{}", accession, self.suffix()),
        }
    }

    fn url(&self, base_url: &str, full_name: &str) -> String {
        format!("{}/{}{}", base_url, full_name, self.suffix())
    }
}

async fn decompress_and_parse_data(compressed_data: Vec<u8>) -> Result<()> {
    task::block_in_place(|| {
        let cursor = Cursor::new(compressed_data);

        let mut fastx_reader = parse_fastx_reader(cursor)
            .context("Failed to parse FASTA/FASTQ data")?;

        while let Some(record) = fastx_reader.next() {
            let record = record.context("Failed to read record")?;
            // println!("Record ID: {}", std::str::from_utf8(record.id())?);
            // println!("Sequence: {}", std::str::from_utf8(record.seq())?);
            // Process each record
        }

        Ok(())
    })
}

async fn process_accession(client: &Client, accession: &str, location: &PathBuf, retry: Option<u32>, keep_fastas: bool) -> Result<()> {
    let retry_count = retry.unwrap_or(3);  // Default retry count

    let (base_url, full_name) = fetch_genbank_filename(client, accession).await?;

    // Combine all file types into a single vector
    let file_types = vec![
        FileType::Genomic,
        FileType::Protein,
        FileType::AssemblyReport,
        FileType::Checksum, // Including standalone files like checksums here
    ];

    for file_type in &file_types {
        let url = file_type.url(&base_url, &full_name); 
        let data = download_with_retry(client, &url, retry_count).await?;

        if keep_fastas {
            let file_name = file_type.filename(accession);
            let path = location.join(&file_name);
            fs::write(&path, &data).context("Failed to write data to file")?;
        }
        match file_type {
            FileType::Genomic => {
                decompress_and_parse_data(data).await?;
            },
            // FileType::Protein => {
            //     sketch_protein(&data).await?;
            // },
            _ => {}  // Do nothing for other file types
        }
    }
    // to download md5sum. To do -- download this + check md5sum on the fly?
    // for filename in standalone {
    //     let url = format!("{}/{}", base_url, filename);
    //     let data = download_with_retry(client, &url, retry_count).await?;
    //     let file_name = format!("{}_{}", accession, filename);  // Generate file name using the directory name and suffix
    //     let path = location.join(&file_name);  // Create the full path for the file
    //     download_with_retry(client, &url, retry_count).await?;
    // }

    Ok(())
}

#[tokio::main]
pub async fn download_and_sketch(
    input_csv: String,
    param_str: String,
    failed_csv: String,
    retry_times: u32,
    fasta_location: String,
    keep_fastas: bool,
   ) -> Result<(), anyhow::Error> {

    let download_path = PathBuf::from(fasta_location);
    if !download_path.exists() {
        create_dir_all(&download_path)?;
    }

    // Open the file containing the accessions synchronously
    let file = File::open(input_csv)?;
    let mut rdr = Reader::from_reader(file);

    // Initialize a HashSet to store unique accessions
    let mut accessions = HashSet::new();

    // Read all accessions into the HashSet to remove duplicates
    for result in rdr.deserialize::<(String,)>() {
        let record = result?;
        let accession = record.0.trim();
        if !accession.is_empty() {
            accessions.insert(accession.to_string());
        }
    }

     // parse param string into params_vec, print error if fail
     let param_result = parse_params_str(param_str);
     let params_vec = match param_result {
         Ok(params) => params,
         Err(e) => {
             eprintln!("Error parsing params string: {}", e);
             bail!("Failed to parse params string");
         }
     };

    // Initialize the HTTP client
    let client = Client::new();
    let mut failed_writer = Writer::from_path(failed_csv)?;
    failed_writer.write_record(&["accession", "url"])?;

    // Process each unique accession in the HashSet
    for accession in &accessions {
        match process_accession(&client, accession, &download_path, Some(retry_times), keep_fastas).await {
            Ok(_) => println!("Successfully processed accession: {}", accession),
            Err(e) => {
                let err_message = e.to_string();
                let parts: Vec<&str> = err_message.split("retries: ").collect();
                let failed_url = parts.get(1).unwrap_or(&"Unknown URL").trim();
                
                failed_writer.write_record(&[accession, failed_url])?;
                eprintln!("Failed to process accession: {}. Error: {}", accession, err_message);
            }
        }
    }


    Ok(())
}

// #[tokio::main]
// pub async fn download_accessions(
//     input_csv: String,
//     failed_csv: String,
//     retry_times: u32,
//     location: String,
// }