use anyhow::{anyhow, Result};
use reqwest::Url;
use sourmash::collection::Collection;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt;

pub mod buildutils;
use crate::utils::buildutils::BuildParams;

#[derive(Clone)]
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

    pub fn buildparams_hashmap(&self) -> HashMap<String, HashSet<u64>> {
        let mut name_params_map = HashMap::new();

        // Iterate over all collections in MultiCollection
        for collection in &self.collections {
            // Iterate over all records in the current collection
            for (_, record) in collection.iter() {
                // Get the record's name or fasta filename
                let record_name = record.name().clone();

                // Calculate the hash of the Params for the current record
                let params_hash = BuildParams::from_record(record).calculate_hash();

                // If the name is already in the HashMap, extend the existing HashSet
                // Otherwise, create a new HashSet and insert the hashed Params
                name_params_map
                    .entry(record_name)
                    .or_insert_with(HashSet::new) // Create a new HashSet if the key doesn't exist
                    .insert(params_hash); // Insert the hashed Params into the HashSet
            }
        }

        name_params_map
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use camino::Utf8PathBuf;
    #[test]
    fn test_buildparams_hashmap() {
        // read in zipfiles to build a MultiCollection
        // load signature + build record
        let mut filename = Utf8PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("tests/test-data/GCA_000961135.2.sig.zip");
        let path = filename.clone();

        let mut collections = Vec::new();
        let coll: Collection = Collection::from_zipfile(&path).unwrap();
        collections.push(coll);
        let mc = MultiCollection::new(collections);

        //  Call build_params_hashmap
        let name_params_map = mc.buildparams_hashmap();

        // Check that the HashMap contains the correct names
        assert_eq!(
            name_params_map.len(),
            1,
            "There should be 1 unique names in the map"
        );

        let mut hashed_params = Vec::new();
        for (name, params_set) in name_params_map.iter() {
            eprintln!("Name: {}", name);
            for param_hash in params_set {
                eprintln!("  Param Hash: {}", param_hash);
                hashed_params.push(param_hash);
            }
        }
        let expected_params1 = BuildParams {
            ksize: 31,
            track_abundance: true,
            ..Default::default()
        };

        let expected_params2 = BuildParams {
            ksize: 21,
            track_abundance: true,
            ..Default::default()
        };

        let expected_hash1 = expected_params1.calculate_hash();
        let expected_hash2 = expected_params2.calculate_hash();

        assert!(
            hashed_params.contains(&&expected_hash1),
            "Expected hash1 should be in the hashed_params"
        );
        assert!(
            hashed_params.contains(&&expected_hash2),
            "Expected hash2 should be in the hashed_params"
        );
    }
}
