use anyhow::{anyhow, Result};
use csv::ReaderBuilder;
use reqwest::Url;
use serde::Deserialize;
use sourmash::collection::Collection;
use std::collections::{HashMap, HashSet};
use std::fmt;
use tokio::io::AsyncWriteExt;

pub mod buildutils;
use crate::utils::buildutils::{BuildManifest, BuildRecord};

#[derive(Clone, PartialEq, Debug, Deserialize)]
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
#[derive(PartialEq, Clone, Eq, Hash)]
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
    #[allow(dead_code)]
    pub fn server_filename(&self, full_name: &str) -> String {
        format!("{}{}", full_name, self.suffix())
    }

    pub fn filename_to_write(&self, accession: &str) -> String {
        match self {
            GenBankFileType::Checksum => format!("{}_{}", accession, self.suffix()),
            _ => format!("{}{}", accession, self.suffix()),
        }
    }

    #[allow(dead_code)]
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

#[derive(Debug, Deserialize, Clone)]
pub struct AccessionData {
    pub accession: String,
    pub name: String,
    pub moltype: InputMolType,

    #[serde(deserialize_with = "deserialize_urlinfo")]
    pub url_info: Vec<UrlInfo>,

    #[serde(default)] // Allows missing `download_filename`
    pub download_filename: Option<String>,
}

#[derive(Debug, Clone)]
pub struct UrlInfo {
    pub url: Url,
    pub md5sum: Option<String>,
    pub range: Option<(usize, usize)>,
}

#[derive(Debug, Deserialize, Clone)]
pub struct GBAssemblyData {
    pub name: String,

    #[serde(deserialize_with = "deserialize_accession_with_range")]
    pub accession_info: Vec<GBAssemblyInfo>,
}

#[derive(Debug, Clone)]
pub struct GBAssemblyInfo {
    pub accession: String,
    pub range: Option<(usize, usize)>,
}

// Custom deserialization for `accession` and `range`
fn deserialize_accession_with_range<'de, D>(
    deserializer: D,
) -> Result<Vec<GBAssemblyInfo>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    #[derive(Deserialize)]
    struct RawRow {
        accession: String,
        #[serde(default)] // Allows missing "range" column
        range: Option<String>,
    }

    let raw: RawRow = Deserialize::deserialize(deserializer)?;

    let accessions: Vec<String> = raw
        .accession
        .split(';')
        .map(|s| s.trim().to_string())
        .collect();

    let ranges = match raw.range {
        Some(r) => parse_ranges(&r, accessions.len()).map_err(serde::de::Error::custom)?,
        None => vec![None; accessions.len()], // Default to None for each accession if no range
    };

    if accessions.len() != ranges.len() {
        return Err(serde::de::Error::custom(format!(
            "Mismatch between accessions ({}) and ranges ({})",
            accessions.len(),
            ranges.len()
        )));
    }

    let result = accessions
        .into_iter()
        .zip(ranges.into_iter())
        .map(|(acc, range)| GBAssemblyInfo {
            accession: acc,
            range,
        })
        .collect();

    Ok(result)
}

pub fn load_gbassembly_info<P: AsRef<std::path::Path>>(
    input_csv: P,
) -> Result<(Vec<GBAssemblyData>, usize)> {
    let required_columns = &["accession", "name"]; // Only enforce required columns exist

    let file = std::fs::File::open(&input_csv).map_err(|e| {
        anyhow!(
            "Failed to open CSV file '{}': {}",
            input_csv.as_ref().display(),
            e
        )
    })?;
    let mut rdr = ReaderBuilder::new().has_headers(true).from_reader(file);

    // Validate headers as we read
    let headers = rdr
        .headers()
        .map_err(|e| anyhow!("Failed to read CSV headers: {}", e))?;

    for &column in required_columns {
        if !headers.iter().any(|h| h == column) {
            return Err(anyhow!("Missing required column: '{}' in CSV file", column));
        }
    }

    let mut results = Vec::new();
    let mut row_count = 0;
    let mut processed_rows = HashSet::new();
    let mut duplicate_count = 0;

    for result in rdr.deserialize() {
        let record: GBAssemblyData =
            result.map_err(|e| anyhow!("Failed to deserialize CSV row: {}", e))?;

        let row_string = format!("{:?}", record);
        if !processed_rows.insert(row_string) {
            duplicate_count += 1;
            continue;
        }
        row_count += 1;
        results.push(record);
    }

    if duplicate_count > 0 {
        println!("Warning: {} duplicated rows were skipped.", duplicate_count);
    }
    println!("Loaded {} rows", row_count);

    Ok((results, row_count))
}

// Custom deserialization for `url` and `md5sum`
fn deserialize_urlinfo<'de, D>(deserializer: D) -> Result<Vec<UrlInfo>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    #[derive(Deserialize)]
    struct RawRow {
        url: String,
        #[serde(default)]
        md5sum: Option<String>,
        #[serde(default)]
        range: Option<String>,
    }

    let raw: RawRow = Deserialize::deserialize(deserializer)?;
    let urls: Vec<&str> = raw.url.split(';').map(|s| s.trim()).collect();

    let md5sums: Vec<Option<String>> = match raw.md5sum {
        Some(ref md5s) => md5s
            .split(';')
            .map(|s| {
                let trimmed = s.trim();
                if trimmed.is_empty() {
                    None
                } else {
                    Some(trimmed.to_string())
                }
            })
            .collect(),
        None => vec![None; urls.len()], // Default to None if no MD5s provided
    };

    let ranges = match raw.range {
        Some(ref r) => parse_ranges(r, urls.len()).map_err(serde::de::Error::custom)?,
        None => vec![None; urls.len()], // Default to None for each URL if no range
    };

    if urls.len() != md5sums.len() || urls.len() != ranges.len() {
        return Err(serde::de::Error::custom(format!(
            "Mismatch between URLs ({}) and MD5 sums ({}) or ranges ({})",
            urls.len(),
            md5sums.len(),
            ranges.len()
        )));
    }

    let mut url_infos = Vec::new();
    for ((url, md5), range) in urls
        .into_iter()
        .zip(md5sums.into_iter())
        .zip(ranges.into_iter())
    {
        let parsed_url = Url::parse(url)
            .map_err(|e| serde::de::Error::custom(format!("Invalid URL '{}': {}", url, e)))?;
        url_infos.push(UrlInfo {
            url: parsed_url,
            md5sum: md5,
            range,
        });
    }

    Ok(url_infos)
}

fn parse_ranges(
    range_field: &str,
    expected_num_ranges: usize,
) -> Result<Vec<Option<(usize, usize)>>, String> {
    if range_field.trim().is_empty() {
        // Return a vector of None for each expected range
        return Ok(vec![None; expected_num_ranges]);
    }

    let ranges: Vec<&str> = range_field.split(';').collect();

    // Check if the number of ranges matches expected_num_ranges
    if ranges.len() != expected_num_ranges {
        return Err(format!(
            "Number of ranges ({}) does not match expected number of ranges ({})",
            ranges.len(),
            expected_num_ranges
        ));
    }

    ranges
        .into_iter()
        .map(|s| {
            let s = s.trim(); // Trim whitespace
            if s.is_empty() {
                return Ok(None); // Treat empty range as None
            }
            let parts: Vec<&str> = s.split('-').collect();
            if parts.len() == 2 {
                let start = parts[0]
                    .parse::<usize>()
                    .map_err(|_| format!("Invalid start value in range: {}", s))?;
                let end = parts[1]
                    .parse::<usize>()
                    .map_err(|_| format!("Invalid end value in range: {}", s))?;
                if start < end {
                    Ok(Some((start, end))) // Return Some for valid ranges
                } else {
                    Err(format!(
                        "Start value must be less than end value in range: {}",
                        s
                    ))
                }
            } else {
                Err(format!("Invalid range format: {}", s))
            }
        })
        .collect()
}

pub fn load_accession_info<P: AsRef<std::path::Path>>(
    input_csv: P,
    keep_fastas: bool,
) -> Result<(Vec<AccessionData>, usize)> {
    let mut required_columns = vec!["accession", "name", "moltype", "url"]; // Enforce required columns

    // If `keep_fastas` is true, add `"download_filename"` to required columns
    if keep_fastas {
        required_columns.push("download_filename");
    }

    let file = std::fs::File::open(&input_csv).map_err(|e| {
        anyhow!(
            "Failed to open CSV file '{}': {}",
            input_csv.as_ref().display(),
            e
        )
    })?;
    let mut rdr = ReaderBuilder::new().has_headers(true).from_reader(file);

    // Validate headers as we read
    let headers = rdr
        .headers()
        .map_err(|e| anyhow!("Failed to read CSV headers: {}", e))?;

    for &column in &required_columns {
        if !headers.iter().any(|h| h == column) {
            return Err(anyhow!("Missing required column: '{}' in CSV file", column));
        }
    }

    let mut results = Vec::new();
    let mut row_count = 0;
    let mut processed_rows = HashSet::new();
    let mut duplicate_count = 0;

    for result in rdr.deserialize() {
        let record: AccessionData =
            result.map_err(|e| anyhow!("Failed to deserialize CSV row: {}", e))?;

        // If `keep_fastas` is true, ensure `download_filename` is present and not None
        if keep_fastas && record.download_filename.is_none() {
            return Err(anyhow!(
                "Missing required value in 'download_filename' column for accession '{}'",
                record.accession
            ));
        }

        let row_string = format!("{:?}", record);
        if !processed_rows.insert(row_string) {
            duplicate_count += 1;
            continue;
        }
        row_count += 1;
        results.push(record);
    }

    if duplicate_count > 0 {
        println!("Warning: {} duplicated rows were skipped.", duplicate_count);
    }
    println!("Loaded {} rows", row_count);

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

#[derive(Clone)]
pub struct FailedDownload {
    accession: String,
    name: String,
    moltype: String,
    md5sum: String,
    download_filename: String,
    url: String,
    range: String,
}

impl FailedDownload {
    /// Build a `FailedDownload` from `GBAssemblyData` with detailed information
    pub fn from_gbassembly(
        accession: String,
        name: String,
        moltype: String,
        md5sum: Option<String>,            // Single MD5 checksum
        download_filename: Option<String>, // Download filename
        url: Option<reqwest::Url>,         // URL for the file
        range: Option<(usize, usize)>,     // Optional range for the download
    ) -> Self {
        Self {
            accession,
            name,
            moltype,
            md5sum: md5sum.unwrap_or_default(),
            download_filename: download_filename.unwrap_or_default(),
            url: url.map(|u| u.to_string()).unwrap_or_default(),
            range: range
                .map(|(start, end)| format!("{}-{}", start, end))
                .unwrap_or_default(), // Format range or use ""
        }
    }

    fn parse_to_separated_string<T, F>(url_info: &[UrlInfo], mut extractor: F) -> String
    where
        F: FnMut(&UrlInfo) -> Option<T>,
        T: ToString,
    {
        let results: Vec<String> = url_info
            .iter()
            .map(|info| extractor(info).map_or("".to_string(), |v| v.to_string())) // Map `None` to empty string
            .collect();

        if results.iter().all(|entry| entry.is_empty()) {
            "".to_string() // If all entries are empty, return `""`
        } else {
            results.join(";") // Otherwise, join with `;`
        }
    }

    /// Build a `FailedDownload` from `AccessionData`
    pub fn from_accession_data(acc_data: &AccessionData) -> Self {
        Self {
            accession: acc_data.accession.clone(),
            name: acc_data.name.clone(),
            moltype: acc_data.moltype.to_string(),
            md5sum: Self::parse_to_separated_string(&acc_data.url_info, |info| info.md5sum.clone()),
            download_filename: acc_data.download_filename.clone().unwrap_or_default(),
            url: Self::parse_to_separated_string(&acc_data.url_info, |info| {
                Some(info.url.to_string())
            }),
            range: Self::parse_to_separated_string(&acc_data.url_info, |info| {
                info.range.map(|(start, end)| format!("{}-{}", start, end))
            }),
        }
    }

    pub fn to_csv_record(&self) -> String {
        format!(
            "{},{},{},{},{},{},{}\n",
            self.accession,
            self.name,
            self.moltype,
            self.md5sum,
            self.download_filename,
            self.url,
            self.range,
        )
    }

    pub fn csv_header() -> &'static str {
        "accession,name,moltype,md5sum,download_filename,url,range\n"
    }

    /// Write a `FailedDownload` to a CSV writer
    pub async fn to_writer<W: tokio::io::AsyncWrite + Unpin>(
        &self,
        writer: &mut W,
    ) -> Result<(), std::io::Error> {
        writer.write_all(self.to_csv_record().as_bytes()).await
    }
}

pub struct FailedChecksum {
    accession: String,
    name: String,
    moltype: String,
    md5sum_url: Option<Url>,
    download_filename: Option<String>,
    url: Option<Url>,
    expected_md5sum: Option<String>,
    reason: String,
}

impl FailedChecksum {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        accession: String,
        name: String,
        moltype: String,
        md5sum_url: Option<reqwest::Url>,
        download_filename: Option<String>,
        url: Option<reqwest::Url>,
        expected_md5sum: Option<String>,
        reason: String,
    ) -> Self {
        Self {
            accession,
            name,
            moltype,
            md5sum_url,
            download_filename,
            url,
            expected_md5sum,
            reason,
        }
    }

    /// Convert a `FailedChecksum` to a CSV-formatted string
    pub fn to_csv_record(&self) -> String {
        let md5sum_url_str = self
            .md5sum_url
            .as_ref()
            .map(|u| u.to_string())
            .unwrap_or_default();

        let url_str = self.url.as_ref().map(|u| u.to_string()).unwrap_or_default();

        format!(
            "{},{},{},{},{},{},{},{}\n",
            self.accession,
            self.name,
            self.moltype,
            md5sum_url_str,
            self.download_filename.clone().unwrap_or_default(),
            url_str,
            self.expected_md5sum.clone().unwrap_or_default(),
            self.reason,
        )
    }

    /// Get the CSV header for a `FailedChecksum`
    pub fn csv_header() -> &'static str {
        "accession,name,moltype,md5sum_url,download_filename,url,expected_md5sum,reason\n"
    }

    /// Write a `FailedChecksum` to a CSV writer
    pub async fn to_writer<W: tokio::io::AsyncWrite + Unpin>(
        &self,
        writer: &mut W,
    ) -> Result<(), std::io::Error> {
        writer.write_all(self.to_csv_record().as_bytes()).await
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use reqwest::Url;

    #[test]
    fn test_parse_urls_valid_urls() {
        let url_field = Some("http://example.com; https://example.org");
        let result = parse_urls(url_field).unwrap();

        assert_eq!(
            result,
            vec![
                Url::parse("http://example.com").unwrap(),
                Url::parse("https://example.org").unwrap()
            ]
        );
    }

    #[test]
    fn test_parse_urls_with_whitespace() {
        let url_field = Some("   http://example.com   ;   https://example.org   ");
        let result = parse_urls(url_field).unwrap();

        assert_eq!(
            result,
            vec![
                Url::parse("http://example.com").unwrap(),
                Url::parse("https://example.org").unwrap()
            ]
        );
    }

    #[test]
    fn test_parse_urls_with_empty_entries() {
        let url_field = Some("http://example.com;;https://example.org");
        let result = parse_urls(url_field);

        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "Empty URL entry found in 'url' field"
        );
    }

    #[test]
    fn test_parse_urls_invalid_url() {
        let url_field = Some("http://example.com; invalid-url; https://example.org");
        let result = parse_urls(url_field);

        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "Invalid URL 'invalid-url': relative URL without a base"
        );
    }

    #[test]
    fn test_parse_urls_empty_field() {
        let url_field = Some("");
        let result = parse_urls(url_field);

        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "Empty URL entry found in 'url' field"
        );
    }

    #[test]
    fn test_parse_urls_missing_field() {
        let url_field = None;
        let result = parse_urls(url_field);

        assert!(result.is_err());
        assert_eq!(result.unwrap_err().to_string(), "Missing 'url' field");
    }

    #[test]
    fn test_parse_urls_all_invalid() {
        let url_field = Some("invalid-url; still-not-a-url");
        let result = parse_urls(url_field);

        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "Invalid URL 'invalid-url': relative URL without a base"
        );
    }

    #[test]
    fn test_parse_ranges_valid() {
        let range_field = "1-10;20-30;40-50";
        let expected_num_ranges = 3;
        let result = parse_ranges(range_field, expected_num_ranges);

        assert!(result.is_ok());
        assert_eq!(
            result.unwrap(),
            vec![Some((1, 10)), Some((20, 30)), Some((40, 50))]
        );
    }

    #[test]
    fn test_parse_ranges_empty_field() {
        let range_field = "   ";
        let expected_num_ranges = 3;
        let result = parse_ranges(range_field, expected_num_ranges);

        assert!(result.is_ok());
        assert_eq!(result.unwrap(), vec![None, None, None]);
    }

    #[test]
    fn test_parse_ranges_neg_start() {
        let range_field = "1-10;-20-30";
        let expected_num_ranges = 2;
        let result = parse_ranges(range_field, expected_num_ranges);

        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Invalid range format: -20-30");
    }

    #[test]
    fn test_parse_ranges_invalid_start() {
        let range_field = "1-10;bar-30";
        let expected_num_ranges = 2;
        let result = parse_ranges(range_field, expected_num_ranges);

        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Invalid start value in range: bar-30");
    }

    #[test]
    fn test_parse_ranges_invalid_end() {
        let range_field = "1-10;20-bar";
        let expected_num_ranges = 2;
        let result = parse_ranges(range_field, expected_num_ranges);

        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Invalid end value in range: 20-bar");
    }

    #[test]
    fn test_parse_ranges_start_not_less_than_end() {
        let range_field = "30-10";
        let expected_num_ranges = 1;
        let result = parse_ranges(range_field, expected_num_ranges);

        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "Start value must be less than end value in range: 30-10"
        );
    }

    #[test]
    fn test_parse_ranges_extra_ranges() {
        let range_field = "1-10;20-30;40-50";
        let expected_num_ranges = 5; // Expecting more ranges than provided
        let result = parse_ranges(range_field, expected_num_ranges);

        assert!(result.is_err()); // Now expecting an error
        assert_eq!(
            result.unwrap_err(),
            "Number of ranges (3) does not match expected number of ranges (5)"
        );
    }

    #[test]
    fn test_parse_ranges_fewer_ranges() {
        let range_field = "1-10;20-30";
        let expected_num_ranges = 3; // Expecting more ranges than provided
        let result = parse_ranges(range_field, expected_num_ranges);

        assert!(result.is_err()); // Now expecting an error
        assert_eq!(
            result.unwrap_err(),
            "Number of ranges (2) does not match expected number of ranges (3)"
        );
    }

    #[test]
    fn test_parse_ranges_with_empty_values() {
        let range_field = "1-10;;20-30";
        let expected_num_ranges = 3;
        let result = parse_ranges(range_field, expected_num_ranges).unwrap();

        assert_eq!(result, vec![Some((1, 10)), None, Some((20, 30))]);
    }

    #[test]
    fn test_parse_md5sums_valid() {
        let md5sum_field = "abcd1234;efgh5678;ijkl9012";
        let expected_num_urls = 3;
        let accession = "ACC123";

        let result = parse_md5sums(md5sum_field, expected_num_urls, accession).unwrap();
        assert_eq!(
            result,
            (
                vec![
                    Some("abcd1234".to_string()),
                    Some("efgh5678".to_string()),
                    Some("ijkl9012".to_string())
                ],
                3
            )
        );
    }

    #[test]
    fn test_parse_md5sums_empty_field() {
        let md5sum_field = "";
        let expected_num_urls = 2;
        let accession = "ACC123";

        let result = parse_md5sums(md5sum_field, expected_num_urls, accession).unwrap();
        assert_eq!(result, (vec![None, None], 0));
    }

    #[test]
    fn test_parse_md5sums_mismatched_count_more_md5s() {
        let md5sum_field = "abcd1234;efgh5678;ijkl9012";
        let expected_num_urls = 2; // Fewer URLs than MD5 sums
        let accession = "ACC123";

        let result = parse_md5sums(md5sum_field, expected_num_urls, accession);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "Number of MD5 sums (3) does not match the number of URLs (2) for accession 'ACC123'"
        );
    }

    #[test]
    fn test_parse_md5sums_mismatched_count_fewer_md5s() {
        let md5sum_field = "abcd1234;efgh5678";
        let expected_num_urls = 3; // More URLs than MD5 sums
        let accession = "ACC123";

        let result = parse_md5sums(md5sum_field, expected_num_urls, accession);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "Number of MD5 sums (2) does not match the number of URLs (3) for accession 'ACC123'"
        );
    }

    #[test]
    fn test_parse_md5sums_with_whitespace() {
        let md5sum_field = "  abcd1234  ; efgh5678 ;  ijkl9012 ";
        let expected_num_urls = 3;
        let accession = "ACC123";

        let result = parse_md5sums(md5sum_field, expected_num_urls, accession).unwrap();
        assert_eq!(
            result,
            (
                vec![
                    Some("abcd1234".to_string()),
                    Some("efgh5678".to_string()),
                    Some("ijkl9012".to_string())
                ],
                3
            )
        );
    }

    #[test]
    fn test_parse_md5sums_some_empty_entries() {
        let md5sum_field = "abcd1234;;ijkl9012";
        let expected_num_urls = 3;
        let accession = "ACC123";

        let result = parse_md5sums(md5sum_field, expected_num_urls, accession).unwrap();
        assert_eq!(
            result,
            (
                vec![
                    Some("abcd1234".to_string()),
                    None, // Empty MD5 sum
                    Some("ijkl9012".to_string())
                ],
                2 // Only count non-empty Some values
            )
        );
    }

    #[test]
    fn test_failed_download_from_gbassembly_valid() {
        let accession = "ACC123".to_string();
        let name = "Sample Name".to_string();
        let moltype = "DNA".to_string();
        let md5sum = Some("abcd1234".to_string());
        let download_filename = Some("file.fasta".to_string());
        let url = Some(Url::parse("http://example.com/file.fasta").unwrap());
        let range = Some((10, 20));

        let failed_download = FailedDownload::from_gbassembly(
            accession.clone(),
            name.clone(),
            moltype.clone(),
            md5sum.clone(),
            download_filename.clone(),
            url.clone(),
            range.clone(),
        );

        assert_eq!(failed_download.accession, accession);
        assert_eq!(failed_download.name, name);
        assert_eq!(failed_download.moltype, moltype);
        assert_eq!(failed_download.md5sum, "abcd1234");
        assert_eq!(failed_download.download_filename, "file.fasta");
        assert_eq!(failed_download.url, "http://example.com/file.fasta");
        assert_eq!(failed_download.range, "10-20");
    }

    #[test]
    fn test_failed_download_from_gbassembly_defaults() {
        let accession = "ACC123".to_string();
        let name = "Sample Name".to_string();
        let moltype = "DNA".to_string();

        let failed_download = FailedDownload::from_gbassembly(
            accession.clone(),
            name.clone(),
            moltype.clone(),
            None, // No MD5 checksum
            None, // No filename
            None, // No URL
            None, // No range
        );

        assert_eq!(failed_download.accession, accession);
        assert_eq!(failed_download.name, name);
        assert_eq!(failed_download.moltype, moltype);
        assert_eq!(failed_download.md5sum, "");
        assert_eq!(failed_download.download_filename, "");
        assert_eq!(failed_download.url, "");
        assert_eq!(failed_download.range, "");
    }

    #[test]
    fn test_failed_download_from_accession_data() {
        let url_info = vec![
            UrlInfo {
                url: Url::parse("http://example.com/file1").unwrap(),
                md5sum: Some("abcd1234".to_string()),
                range: Some((10, 20)),
            },
            UrlInfo {
                url: Url::parse("http://example.com/file2").unwrap(),
                md5sum: None,
                range: Some((30, 40)),
            },
        ];

        let acc_data = AccessionData {
            accession: "ACC123".to_string(),
            name: "Sample Name".to_string(),
            moltype: InputMolType::Dna,
            url_info,
            download_filename: Some("file.fasta".to_string()),
        };

        let failed_download = FailedDownload::from_accession_data(&acc_data);

        assert_eq!(failed_download.accession, "ACC123");
        assert_eq!(failed_download.name, "Sample Name");
        assert_eq!(failed_download.moltype, "DNA");
        assert_eq!(failed_download.md5sum, "abcd1234;");
        assert_eq!(failed_download.download_filename, "file.fasta");
        assert_eq!(
            failed_download.url,
            "http://example.com/file1;http://example.com/file2"
        );
        assert_eq!(failed_download.range, "10-20;30-40");
    }

    #[test]
    fn test_parse_to_separated_string() {
        let url_info = vec![
            UrlInfo {
                url: Url::parse("http://example.com/file1").unwrap(),
                md5sum: Some("abcd1234".to_string()),
                range: Some((10, 20)),
            },
            UrlInfo {
                url: Url::parse("http://example.com/file2").unwrap(),
                md5sum: None,
                range: Some((30, 40)),
            },
        ];

        let md5sum_result =
            FailedDownload::parse_to_separated_string(&url_info, |info| info.md5sum.clone());
        assert_eq!(md5sum_result, "abcd1234;");

        let url_result =
            FailedDownload::parse_to_separated_string(&url_info, |info| Some(info.url.to_string()));
        assert_eq!(
            url_result,
            "http://example.com/file1;http://example.com/file2"
        );

        let range_result = FailedDownload::parse_to_separated_string(&url_info, |info| {
            info.range.map(|(start, end)| format!("{}-{}", start, end))
        });
        assert_eq!(range_result, "10-20;30-40");
    }

    #[test]
    fn test_parse_to_separated_string_2() {
        let url_info = vec![
            UrlInfo {
                url: Url::parse("http://example.com/file1").unwrap(),
                md5sum: Some("abcd1234".to_string()),
                range: Some((10, 20)),
            },
            UrlInfo {
                url: Url::parse("http://example.org/file2").unwrap(),
                md5sum: Some("efgh5678".to_string()),
                range: Some((30, 40)),
            },
            UrlInfo {
                url: Url::parse("http://example.net/file3").unwrap(),
                md5sum: Some("ijkl9012".to_string()),
                range: Some((50, 60)),
            },
        ];

        let md5sum_result =
            FailedDownload::parse_to_separated_string(&url_info, |info| info.md5sum.clone());
        assert_eq!(md5sum_result, "abcd1234;efgh5678;ijkl9012");

        let url_result =
            FailedDownload::parse_to_separated_string(&url_info, |info| Some(info.url.to_string()));
        assert_eq!(
            url_result,
            "http://example.com/file1;http://example.org/file2;http://example.net/file3"
        );

        let range_result = FailedDownload::parse_to_separated_string(&url_info, |info| {
            info.range.map(|(start, end)| format!("{}:{}", start, end))
        });
        assert_eq!(range_result, "10:20;30:40;50:60");
    }
}
