use anyhow::{anyhow, Result};
use reqwest::Url;
use sourmash::collection::Collection;
use std::collections::HashMap;
use std::fmt;
use tokio::io::AsyncWriteExt;

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
#[derive(PartialEq)]
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

fn parse_urls(url_field: Option<&str>) -> Result<Vec<reqwest::Url>, anyhow::Error> {
    let url_field = url_field.ok_or_else(|| anyhow!("Missing 'url' field"))?;

    let mut urls = Vec::new();

    for s in url_field.split(';').map(|s| s.trim()) {
        if s.is_empty() {
            return Err(anyhow!("Empty URL entry found in 'url' field"));
        }

        let parsed_url =
            reqwest::Url::parse(s).map_err(|e| anyhow!("Invalid URL '{}': {}", s, e))?;
        urls.push(parsed_url);
    }

    if urls.is_empty() {
        return Err(anyhow!("No valid URLs found in 'url' field"));
    }

    Ok(urls)
}

fn parse_md5sums(
    md5sum_field: &str,
    expected_num_urls: usize,
    accession: &str,
) -> Result<(Vec<Option<String>>, usize), anyhow::Error> {
    if md5sum_field.trim().is_empty() {
        // Return a vector of None for each expected URL and a count of 0
        return Ok((vec![None; expected_num_urls], 0));
    }

    let md5sums: Vec<Option<String>> = md5sum_field
        .split(';')
        .map(|s| {
            let trimmed = s.trim();
            if trimmed.is_empty() {
                None
            } else {
                Some(trimmed.to_string())
            }
        })
        .collect();

    // Validate the number of MD5 sums matches the expected number of URLs
    if md5sums.len() != expected_num_urls {
        return Err(anyhow::anyhow!(
            "Number of MD5 sums ({}) does not match the number of URLs ({}) for accession '{}'",
            md5sums.len(),
            expected_num_urls,
            accession
        ));
    }

    // Count the number of non-None MD5 sums
    let count = md5sums.iter().filter(|md5| md5.is_some()).count();

    Ok((md5sums, count))
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
        let url_result = parse_urls(record.get(5));
        let urls = match url_result {
            Ok(urls) => {
                if urls.is_empty() {
                    return Err(anyhow!("No valid URLs found in 'url' field"));
                }
                urls
            }
            Err(e) => return Err(e), // Propagate the error if parsing fails
        };

        // Parse MD5sums (optional)
        let (md5sums, md5sum_count_in_row) =
            parse_md5sums(record.get(3).unwrap_or(""), urls.len(), &acc)?;
        // Update the overall MD5 sum count
        md5sum_count += md5sum_count_in_row;

        // Parse ranges (optional)
        let range_field = record.get(6).unwrap_or("");
        let ranges = parse_ranges(range_field, urls.len()).map_err(|e| anyhow!("{}", e))?;

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
