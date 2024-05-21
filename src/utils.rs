use anyhow::{anyhow, Result};
use reqwest::Url;
use sourmash::cmd::ComputeParameters;
use sourmash::signature::Signature;
use std::fmt;
use std::hash::Hash;
use std::hash::Hasher;

#[derive(Clone)]
pub enum InputMolType {
    Dna,
    Protein,
}

impl InputMolType {}

impl fmt::Display for InputMolType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            InputMolType::Dna => write!(f, "dna"),
            InputMolType::Protein => write!(f, "protein"),
        }
    }
}

impl std::str::FromStr for InputMolType {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "dna" => Ok(InputMolType::Dna),
            "protein" => Ok(InputMolType::Protein),
            _ => Err(()),
        }
    }
}

#[allow(dead_code)]
pub enum GenBankFileType {
    Genomic,
    Protein,
    AssemblyReport,
    Checksum,
}

impl GenBankFileType {
    pub fn suffix(&self) -> &'static str {
        match self {
            GenBankFileType::Genomic => "_genomic.fna.gz",
            GenBankFileType::Protein => "_protein.faa.gz",
            GenBankFileType::AssemblyReport => "_assembly_report.txt",
            GenBankFileType::Checksum => "md5checksums.txt",
        }
    }

    //use for checksums
    pub fn server_filename(&self, full_name: &str) -> String {
        format!("{}{}", full_name, self.suffix())
    }

    pub fn filename_to_write(&self, accession: &str) -> String {
        match self {
            GenBankFileType::Checksum => format!("{}_{}", accession, self.suffix()),
            _ => format!("{}{}", accession, self.suffix()),
        }
    }

    pub fn url(&self, base_url: &Url, full_name: &str) -> Url {
        match self {
            GenBankFileType::Checksum => base_url
                .join(&format!("{}/{}", full_name, self.suffix()))
                .unwrap(),
            _ => base_url
                .join(&format!("{}/{}{}", full_name, full_name, self.suffix()))
                .unwrap(),
        }
    }

    pub fn moltype(&self) -> String {
        match self {
            GenBankFileType::Genomic => "DNA".to_string(),
            GenBankFileType::Protein => "protein".to_string(),
            _ => "".to_string(),
        }
    }
}
#[allow(dead_code)]
#[derive(Clone)]
pub struct AccessionData {
    pub accession: String,
    pub name: String,
    pub moltype: InputMolType,
    pub url: reqwest::Url,
    pub expected_md5sum: Option<String>,
    pub download_filename: Option<String>, // need to require this if --keep-fastas are used
}

#[derive(Clone)]
pub struct GBAssemblyData {
    pub accession: String,
    pub name: String,
    pub url: Option<reqwest::Url>,
}

pub fn load_gbassembly_info(input_csv: String) -> Result<(Vec<GBAssemblyData>, usize)> {
    let mut results = Vec::new();
    let mut row_count = 0;
    let mut processed_rows = std::collections::HashSet::new();
    let mut duplicate_count = 0;
    let mut url_count = 0; // Counter for entries with URL
                           // to do - maybe use HashSet for accessions too to avoid incomplete dupes
    let mut rdr = csv::Reader::from_path(input_csv)?;

    // Check column names
    let header = rdr.headers()?;
    let expected_header = vec!["accession", "name", "ftp_path"];
    if header != expected_header {
        return Err(anyhow!(
            "Invalid column names in CSV file. Columns should be: {:?}",
            expected_header
        ));
    }

    for result in rdr.records() {
        let record = result?;
        let row_string = record.iter().collect::<Vec<_>>().join(",");
        if processed_rows.contains(&row_string) {
            duplicate_count += 1;
            continue;
        }
        processed_rows.insert(row_string.clone());
        row_count += 1;

        // require acc, name
        let acc = record
            .get(0)
            .ok_or_else(|| anyhow!("Missing 'accession' field"))?
            .to_string();
        let name = record
            .get(1)
            .ok_or_else(|| anyhow!("Missing 'name' field"))?
            .to_string();
        // optionally get url
        let url = record.get(2).and_then(|s| {
            if s.is_empty() {
                None
            } else {
                let trimmed_s = s.trim_end_matches('/');
                reqwest::Url::parse(trimmed_s).map_err(|_| ()).ok()
            }
        });

        if url.is_some() {
            url_count += 1;
        }
        // store accession data
        results.push(GBAssemblyData {
            accession: acc.to_string(),
            name: name.to_string(),
            url,
        });
    }

    // Print warning if there were duplicated rows.
    if duplicate_count > 0 {
        println!("Warning: {} duplicated rows were skipped.", duplicate_count);
    }
    println!(
        "Loaded {} rows (including {} rows with valid URL).",
        row_count, url_count
    );

    Ok((results, row_count))
}

pub fn load_accession_info(
    input_csv: String,
    keep_fasta: bool,
) -> Result<(Vec<AccessionData>, usize)> {
    let mut results = Vec::new();
    let mut row_count = 0;
    let mut processed_rows = std::collections::HashSet::new();
    let mut duplicate_count = 0;
    let mut md5sum_count = 0; // Counter for entries with MD5sum
                              // to do - maybe use HashSet for accessions too to avoid incomplete dupes
    let mut rdr = csv::Reader::from_path(input_csv)?;

    // Check column names
    let header = rdr.headers()?;
    let expected_header = vec![
        "accession",
        "name",
        "moltype",
        "md5sum",
        "download_filename",
        "url",
    ];
    if header != expected_header {
        return Err(anyhow!(
            "Invalid column names in CSV file. Columns should be: {:?}",
            expected_header
        ));
    }

    for result in rdr.records() {
        let record = result?;
        let row_string = record.iter().collect::<Vec<_>>().join(",");
        if processed_rows.contains(&row_string) {
            duplicate_count += 1;
            continue;
        }
        processed_rows.insert(row_string.clone());
        row_count += 1;

        // require acc, name
        let acc = record
            .get(0)
            .ok_or_else(|| anyhow!("Missing 'accession' field"))?
            .to_string();
        let name = record
            .get(1)
            .ok_or_else(|| anyhow!("Missing 'name' field"))?
            .to_string();
        let moltype = record
            .get(2)
            .ok_or_else(|| anyhow!("Missing 'moltype' field"))?
            .parse::<InputMolType>()
            .map_err(|_| anyhow!("Invalid 'moltype' value"))?;
        let expected_md5sum = record.get(3).map(|s| s.to_string());
        let download_filename = record.get(4).map(|s| s.to_string());
        if keep_fasta && download_filename.is_none() {
            return Err(anyhow!("Missing 'download_filename' field"));
        }
        let url = record
            .get(5)
            .ok_or_else(|| anyhow!("Missing 'url' field"))?
            .split(',')
            .filter_map(|s| {
                if s.starts_with("http://") || s.starts_with("https://") || s.starts_with("ftp://")
                {
                    reqwest::Url::parse(s).ok()
                } else {
                    None
                }
            })
            .next()
            .ok_or_else(|| anyhow!("Invalid 'url' value"))?;
        // count entries with url and md5sum
        if expected_md5sum.is_some() {
            md5sum_count += 1;
        }
        // store accession data
        results.push(AccessionData {
            accession: acc,
            name,
            moltype,
            url,
            expected_md5sum,
            download_filename,
        });
    }

    // Print warning if there were duplicated rows.
    if duplicate_count > 0 {
        println!("Warning: {} duplicated rows were skipped.", duplicate_count);
    }
    println!(
        "Loaded {} rows (including {} rows with MD5sum).",
        row_count, md5sum_count
    );

    Ok((results, row_count))
}

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

pub fn parse_params_str(params_strs: String) -> Result<Vec<Params>, String> {
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

pub fn build_siginfo(params: &[Params], moltype: &str) -> Vec<Signature> {
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
