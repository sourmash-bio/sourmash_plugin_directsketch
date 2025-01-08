use anyhow::{anyhow, Result};
use reqwest::Url;
use sourmash::collection::Collection;
use std::collections::HashMap;
use std::fmt;

pub mod buildutils;
use crate::utils::buildutils::{BuildManifest, BuildRecord};

#[derive(Clone, PartialEq)]
pub enum InputMolType {
    Dna,
    Protein,
}

impl InputMolType {}

impl fmt::Display for InputMolType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            InputMolType::Dna => write!(f, "DNA"),
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
    pub url_info: Vec<UrlInfo>,
    pub download_filename: Option<String>, // Need to require this if --keep-fastas are used
}

#[derive(Clone)]
pub struct UrlInfo {
    pub url: reqwest::Url,
    pub md5sum: Option<String>,
    pub range: Option<(usize, usize)>,
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
        "range",
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

        // Parse URLs
        let url_field = record
            .get(5)
            .ok_or_else(|| anyhow!("Missing 'url' field"))?;

        let urls: Vec<reqwest::Url> = url_field
            .split(';')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .filter_map(|s| reqwest::Url::parse(s).ok())
            .collect();

        if urls.is_empty() {
            return Err(anyhow!("No valid URLs found in 'url' field"));
        }

        // Parse MD5sums (optional)
        let md5sum_field = record.get(3).unwrap_or("");
        let md5sums: Vec<Option<String>> = if !md5sum_field.trim().is_empty() {
            md5sum_field
                .split(';')
                .map(|s| Some(s.trim().to_string()))
                .collect()
        } else {
            vec![None; urls.len()]
        };

        if md5sums.len() != urls.len() {
            return Err(anyhow!(
                "Number of MD5 sums ({}) does not match the number of URLs ({}) for accession '{}'",
                md5sums.len(),
                urls.len(),
                acc
            ));
        }

        // Check if there are any MD5 sums in this row
        if md5sums.iter().any(|md5| md5.is_some()) {
            md5sum_count += 1;
        }

        // Parse ranges (optional)
        let range_field = record.get(6).unwrap_or("");
        let ranges: Vec<Option<(usize, usize)>> = if !range_field.trim().is_empty() {
            range_field
                .split(';')
                .map(|s| {
                    let s = s.trim(); // trim whitespace
                    let parts: Vec<&str> = s.split('-').collect();
                    if parts.len() == 2 {
                        if let (Ok(start), Ok(end)) = (parts[0].parse(), parts[1].parse()) {
                            if start < end {
                                return Some((start, end));
                            }
                        }
                    }
                    None
                })
                .collect()
        } else {
            vec![None; urls.len()]
        };

        // Combine URLs, MD5 sums, and ranges into UrlInfo
        let url_info: Vec<UrlInfo> = urls
            .into_iter()
            .zip(md5sums)
            .zip(ranges)
            .map(|((url, md5sum), range)| UrlInfo { url, md5sum, range })
            .collect();

        let download_filename = record.get(4).map(|s| s.to_string());
        if keep_fasta && download_filename.is_none() {
            return Err(anyhow!("Missing 'download_filename' field"));
        }
        // store accession data
        results.push(AccessionData {
            accession: acc,
            name,
            moltype,
            url_info,
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

// this should be replaced with branchwater's MultiCollection when it's ready
#[derive(Clone)]
pub struct MultiCollection {
    collections: Vec<Collection>,
}

impl MultiCollection {
    pub fn new(collections: Vec<Collection>) -> Self {
        Self { collections }
    }

    pub fn is_empty(&self) -> bool {
        self.collections.is_empty()
    }

    pub fn build_recordsmap(&self) -> HashMap<String, BuildManifest> {
        let mut records_map = HashMap::new();
        // Iterate over all collections in MultiCollection
        for collection in &self.collections {
            // Iterate over all records in the current collection
            for (_, record) in collection.iter() {
                // Get the record's name or fasta filename
                let record_name = record.name().clone();

                // Create template buildrecord from this record
                let build_record = BuildRecord::from_record(record);

                // If the name is already in the HashMap, extend the existing HashSet
                // Otherwise, create a new BuildManifest and insert the BuildRecord
                records_map
                    .entry(record_name)
                    .or_insert_with(BuildManifest::default) // Create a new HashSet if the key doesn't exist
                    .add_record(build_record); // add buildrecord to buildmanifest
            }
        }

        records_map
    }
}
