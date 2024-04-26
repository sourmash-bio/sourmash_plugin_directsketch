use anyhow::{Result, Context, anyhow};
// use std::path::Path;
use camino::Utf8PathBuf as PathBuf;
use csv::Reader;
use regex::Regex;
use reqwest::Client;
use std::{collections::HashSet, fs};
use std::fs::{File, create_dir_all};
use csv::Writer;

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


async fn download_with_retry(client: &Client, url: &str, file_name: PathBuf, retry_count: u32) -> Result<()> {
    let mut attempts = retry_count;
    while attempts > 0 {
        let response = client.get(url).send().await;
        match response {
            Ok(resp) if resp.status().is_success() => {
                let data = resp.bytes().await.context("Failed to read bytes from response")?;
                fs::write(file_name, &data).context("Failed to write data to file")?;
                return Ok(());
            },
            _ => {
                eprintln!("Failed to download file: {}. Retrying...", url);
                attempts -= 1;
            }
        }
    }

    Err(anyhow!("Failed to download file after {} retries: {}", retry_count, url))
}


async fn process_accession(client: &Client, accession: &str, location: &PathBuf, retry: Option<u32>) -> Result<()> {
    let retry_count = retry.unwrap_or(3);  // Default retry count

    let (base_url, full_name) = fetch_genbank_filename(client, accession).await?;

    
    let suffixes = vec!["_genomic.fna.gz", "_protein.faa.gz", "_assembly_report.txt"];
    let standalone = vec!["md5checksums.txt"];

    for suffix in suffixes.iter() {
        let url = format!("{}/{}{}", base_url, full_name, suffix);  // Correctly format the URL for each file type
        let file_name = format!("{}{}", accession, suffix);  // Generate file name using the directory name and suffix
        let path = location.join(&file_name);  // Create the full path for the file
        download_with_retry(client, &url, path, retry_count).await?;
    }
    
    // download standalone files (mostly md5checksums.txt)
    for filename in standalone {
        let url = format!("{}/{}", base_url, filename);
        let file_name = format!("{}_{}", accession, filename);  // Generate file name using the directory name and suffix
        let path = location.join(&file_name);  // Create the full path for the file
        download_with_retry(client, &url, path, retry_count).await?;
    }

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
) -> Result<(), Box<dyn std::error::Error>> {

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

    // Initialize the HTTP client
    let client = Client::new();
    let mut failed_writer = Writer::from_path(failed_csv)?;
    failed_writer.write_record(&["accession", "url"])?;

    // Process each unique accession in the HashSet
    for accession in &accessions {
        match process_accession(&client, accession, &download_path, Some(retry_times)).await {
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