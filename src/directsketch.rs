use anyhow::{anyhow, bail, Context, Result};
use camino::Utf8PathBuf as PathBuf;
use csv::Writer;
use regex::Regex;
use reqwest::Client;
use std::fs::{self, create_dir_all}; //File
use std::path::Path;
// use niffler::get_reader;
use needletail::parse_fastx_reader;
use std::io::Cursor;
use tokio::task;

use crate::utils::{build_siginfo, load_accession_info, parse_params_str}; //, sigwriter, Params, ZipMessage};

enum GenBankFileType {
    Genomic,
    Protein,
    AssemblyReport,
    Checksum,
}

impl GenBankFileType {
    fn suffix(&self) -> &'static str {
        match self {
            GenBankFileType::Genomic => "_genomic.fna.gz",
            GenBankFileType::Protein => "_protein.faa.gz",
            GenBankFileType::AssemblyReport => "_assembly_report.txt",
            GenBankFileType::Checksum => "md5checksums.txt",
        }
    }

    fn filename(&self, accession: &str) -> String {
        match self {
            GenBankFileType::Checksum => format!("{}_{}", accession, self.suffix()),
            _ => format!("{}{}", accession, self.suffix()),
        }
    }

    fn url(&self, base_url: &str, full_name: &str) -> String {
        format!("{}/{}{}", base_url, full_name, self.suffix())
    }
}

async fn fetch_genbank_filename(client: &Client, accession: &str) -> Result<(String, String)> {
    let (db, acc) = accession
        .trim()
        .split_once('_')
        .ok_or_else(|| anyhow!("Invalid accession format"))?;
    let (number, _) = acc.split_once('.').unwrap_or((acc, "1"));
    let number_path = number
        .chars()
        .collect::<Vec<_>>()
        .chunks(3)
        .map(|chunk| chunk.iter().collect::<String>())
        .collect::<Vec<_>>()
        .join("/");

    let base_url = format!(
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}",
        db, number_path
    );
    let directory_response = client.get(&base_url).send().await?;
    if !directory_response.status().is_success() {
        return Err(anyhow!(
            "Failed to open genome directory: HTTP {}",
            directory_response.status()
        ));
    }
    let text = directory_response.text().await?;
    let link_regex = Regex::new(r#"<a href="([^"]*)""#)?;

    for cap in link_regex.captures_iter(&text) {
        let name = &cap[1];
        // Check if name ends with '/', and remove it if so
        let clean_name = if name.ends_with('/') {
            name.strip_suffix('/').unwrap()
        } else {
            name
        };
        if clean_name.starts_with(db)
            && clean_name
                .split('_')
                .nth(1)
                .map_or(false, |x| x.starts_with(number))
        {
            // Formulate the correct URL by ensuring no double slashes
            return Ok((format!("{}/{}", base_url, clean_name), clean_name.into()));
        }
    }

    Err(anyhow!(
        "No matching genome found for accession {}",
        accession
    ))
}

// download and return data directly instead of saving to file
async fn download_with_retry(client: &Client, url: &str, retry_count: u32) -> Result<Vec<u8>> {
    let mut attempts = retry_count;
    while attempts > 0 {
        let response = client.get(url).send().await;
        match response {
            Ok(resp) if resp.status().is_success() => {
                let data = resp
                    .bytes()
                    .await
                    .context("Failed to read bytes from response")?;
                return Ok(data.to_vec()); // Return the downloaded data as Vec<u8>
            }
            _ => {
                eprintln!("Failed to download file: {}. Retrying...", url);
                attempts -= 1;
            }
        }
    }

    Err(anyhow!(
        "Failed to download file after {} retries: {}",
        retry_count,
        url
    ))
}

async fn decompress_and_parse_data(
    _accession: &str,
    _name: &str,
    compressed_data: Vec<u8>,
) -> Result<()> {
    task::block_in_place(|| {
        let cursor = Cursor::new(compressed_data);

        let mut fastx_reader =
            parse_fastx_reader(cursor).context("Failed to parse FASTA/FASTQ data")?;

        while let Some(record) = fastx_reader.next() {
            let record = record.context("Failed to read record")?;
            // Process each record
            println!("Record ID: {}", std::str::from_utf8(&record.id())?);
            println!("Sequence: {}", std::str::from_utf8(&record.seq())?);
        }

        Ok(())
    })
}

async fn process_accession(
    client: &Client,
    accession: String,
    name: String,
    location: &PathBuf,
    retry: Option<u32>,
    keep_fastas: bool,
) -> Result<()> {
    let retry_count = retry.unwrap_or(3); // Default retry count

    let (base_url, full_name) = fetch_genbank_filename(client, accession.as_str()).await?;

    // Combine all file types into a single vector
    let file_types = vec![
        GenBankFileType::Genomic,
        GenBankFileType::Protein,
        // GenBankFileType::AssemblyReport,
        // GenBankFileType::Checksum, // Including standalone files like checksums here
    ];

    for file_type in &file_types {
        let url = file_type.url(&base_url, &full_name);
        let data = download_with_retry(client, &url, retry_count).await?;

        if keep_fastas {
            let file_name = file_type.filename(&accession);
            let path = location.join(&file_name);
            fs::write(&path, &data).context("Failed to write data to file")?;
        }
        match file_type {
            // also pass in hashfunction to determine dna vs prot sketch
            GenBankFileType::Genomic => {
                decompress_and_parse_data(accession.as_str(), name.as_str(), data).await?;
            }
            GenBankFileType::Protein => {
                decompress_and_parse_data(accession.as_str(), name.as_str(), data).await?;
            }
            _ => {} // Do nothing for other file types
        }
    }

    Ok(())
}

#[tokio::main]
pub async fn download_and_sketch(
    input_csv: String,
    output_sigs: String,
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

    // if sig output doesn't end in zip, bail
    if Path::new(&output_sigs)
        .extension()
        .map_or(true, |ext| ext != "zip")
    {
        bail!("Output must be a zip file.");
    }

    // Open the file containing the accessions synchronously
    let (accession_info, _n_accs) = load_accession_info(input_csv)?;

    // parse param string into params_vec, print error if fail
    let param_result = parse_params_str(param_str);
    let params_vec = match param_result {
        Ok(params) => params,
        Err(e) => {
            eprintln!("Error parsing params string: {}", e);
            bail!("Failed to parse params string");
        }
    };
    let _dna_sig_templates = build_siginfo(&params_vec, "DNA");
    let _prot_sig_templates = build_siginfo(&params_vec, "protein");

    // if no sigs to build, skip this iteration
    // if dna_sig_templates.is_empty() && prot_sig_templates.is_empty() {
    //     if not keep_fastas:
    //         bail!("No signatures to build. Not storin")
    // }

    // Initialize the HTTP client
    let client = Client::new();
    let mut failed_writer = Writer::from_path(failed_csv)?;
    failed_writer.write_record(&["accession", "url"])?;

    // Process each accession
    for accinfo in &accession_info {
        match process_accession(
            &client,
            accinfo.accession.clone(),
            accinfo.name.clone(),
            &download_path,
            Some(retry_times),
            keep_fastas,
        )
        .await
        {
            Ok(_) => println!("Successfully processed accession: {}", &accinfo.accession),
            Err(e) => {
                let err_message = e.to_string();
                let parts: Vec<&str> = err_message.split("retries: ").collect();
                let failed_url = parts.get(1).unwrap_or(&"Unknown URL").trim();

                failed_writer.write_record(&[&accinfo.accession, &failed_url.to_string()])?;
                eprintln!(
                    "Failed to process accession: {}. Error: {}",
                    &accinfo.accession, err_message
                );
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
