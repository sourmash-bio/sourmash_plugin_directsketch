use anyhow::{anyhow, bail, Context, Error, Result};
use async_zip::base::write::ZipFileWriter;
use camino::Utf8PathBuf as PathBuf;
use regex::Regex;
use reqwest::Client;
use sourmash::collection::Collection;
use std::cmp::max;
use std::collections::{HashMap, HashSet};
use std::fs::{self, create_dir_all};
use std::panic;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use tokio::fs::File;
use tokio::io::{AsyncWriteExt, BufWriter};
use tokio::sync::Semaphore;
use tokio_util::compat::Compat;

use pyo3::prelude::*;

use crate::utils::{
    load_accession_info, load_gbassembly_info, parse_params_str, AccessionData, BuildCollection,
    BuildManifest, GBAssemblyData, GenBankFileType, InputMolType, MultiBuildCollection,
    MultiCollection,
};
use reqwest::Url;

async fn find_genome_directory(
    client: &Client,
    db: &str,
    number_path: &str,
    accession: &str,
    acc_number: &str,
    version: &str,
) -> Result<(Url, String)> {
    let base_url = format!(
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}",
        db, number_path
    );
    let directory_response = client.get(&base_url).send().await?;

    if !directory_response.status().is_success() {
        return Err(anyhow!(
            "Failed to open genome directory: HTTP {}, {}",
            directory_response.status(),
            directory_response
                .status()
                .canonical_reason()
                .unwrap_or("Unknown reason")
        ));
    }
    let text = directory_response.text().await?;
    let link_regex = Regex::new(r#"<a href="([^"]+)/""#)?;

    for cap in link_regex.captures_iter(&text) {
        let name = &cap[1];

        // Use acc numerical identifier and version as expected pattern
        let expected_pattern = format!("{}.{}", acc_number, version);

        // Check if directory name contains contains the right acc number and version
        if name.contains(&expected_pattern) {
            let url = format!("{}/{}", base_url, name).parse()?;
            return Ok((url, name.to_string()));
        }
    }
    Err(anyhow!(
        "No matching directory found for accession {}",
        accession
    ))
}

async fn fetch_genbank_filename(
    client: &Client,
    accession: &str,
    url: Option<Url>,
) -> Result<(Url, String)> {
    let (db, acc) = accession
        .trim()
        .split_once('_')
        .ok_or_else(|| anyhow!("Invalid accession format"))?;
    let (acc_number, version) = acc.split_once('.').unwrap_or((acc, "1"));
    let number_path = acc_number
        .chars()
        .collect::<Vec<_>>()
        .chunks(3)
        .map(|chunk| chunk.iter().collect::<String>())
        .collect::<Vec<_>>()
        .join("/");

    let (url, name) = if let Some(url) = url {
        let url_parts: Vec<&str> = url
            .path_segments()
            .ok_or_else(|| anyhow!("Failed to extract path segments from URL"))?
            .collect();
        let name = url_parts
            .last()
            .ok_or_else(|| anyhow!("Failed to extract name from URL"))?
            .to_string();
        (url, name)
    } else {
        find_genome_directory(client, db, &number_path, accession, acc_number, version).await?
    };

    Ok((url, name))
}

async fn download_and_parse_md5(client: &Client, url: &Url) -> Result<HashMap<String, String>> {
    let response = client
        .get(url.clone())
        .send()
        .await
        .context("Failed to send request")?;

    // Check if the status code is 200 OK
    if !response.status().is_success() {
        return Err(anyhow!(
            "Failed to download MD5 checksum file from URL {}: HTTP status {}",
            url,
            response.status()
        ));
    }

    let content = response
        .text()
        .await
        .context("Failed to download MD5 checksum file")?;

    let mut checksums = HashMap::new();

    // Iterate over each line in the checksum file
    for line in content.lines() {
        let parts: Vec<&str> = line.splitn(2, ' ').collect();
        if parts.len() == 2 {
            // Trim the filename to remove any leading " ./" if present
            let filename = parts[1].trim_start_matches(" ./"); // remove any ' ', '.', '/' from front
            checksums.insert(filename.to_string(), parts[0].to_string());
        } else {
            return Err(anyhow!(
                "Invalid checksum line format in URL '{}': '{}'",
                url,
                line
            ));
        }
    }

    Ok(checksums)
}

// download and return data directly instead of saving to file
async fn download_with_retry(
    client: &Client,
    url: &Url,
    expected_md5: Option<&str>,
    retry_count: u32,
) -> Result<Vec<u8>> {
    let mut attempts = retry_count;
    let mut last_error: Option<anyhow::Error> = None;

    while attempts > 0 {
        let response = client.get(url.clone()).send().await;

        match response {
            Ok(resp) if resp.status().is_success() => {
                let data = resp
                    .bytes()
                    .await
                    .context("Failed to read bytes from response")?;

                if let Some(md5) = expected_md5 {
                    let computed_hash = format!("{:x}", md5::compute(&data));
                    if computed_hash == md5 {
                        return Ok(data.to_vec());
                    } else {
                        last_error = Some(anyhow!(
                            "MD5 hash does not match. Expected: '{}'; Found: '{}'",
                            md5,
                            computed_hash
                        ));
                    }
                } else {
                    return Ok(data.to_vec()); // If no expected MD5 is provided, just return the data
                }
            }
            Ok(resp) => {
                last_error = Some(anyhow!(
                    "Server error status code {}: {}. Retrying...",
                    resp.status(),
                    url
                ));
            }
            Err(e) => {
                last_error = Some(anyhow!("Failed to download file: {}. Error: {}.", url, e));
            }
        }

        attempts -= 1;
        if attempts == 0 {
            break;
        }
    }
    Err(last_error.unwrap_or_else(|| {
        anyhow!(
            "Failed to download file after {} retries: {}",
            url,
            retry_count
        )
    }))
}

pub struct FailedDownload {
    accession: String,
    name: String,
    moltype: String,
    md5sum: Option<String>,
    download_filename: Option<String>,
    url: Option<Url>,
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

#[allow(clippy::too_many_arguments)]
async fn dl_sketch_assembly_accession(
    client: &Client,
    accinfo: GBAssemblyData,
    location: &PathBuf,
    retry: Option<u32>,
    keep_fastas: bool,
    dna_sigs: &mut BuildCollection,
    prot_sigs: &mut BuildCollection,
    genomes_only: bool,
    proteomes_only: bool,
    download_only: bool,
) -> Result<(
    MultiBuildCollection,
    Vec<FailedDownload>,
    Vec<FailedChecksum>,
)> {
    let retry_count = retry.unwrap_or(3); // Default retry count
    let mut built_sigs = MultiBuildCollection::new();
    let mut download_failures = Vec::<FailedDownload>::new();
    let mut checksum_failures = Vec::<FailedChecksum>::new();

    let name = accinfo.name;
    let accession = accinfo.accession;

    // keep track of any accessions for which we fail to find URLs
    let (base_url, full_name) =
        match fetch_genbank_filename(client, accession.as_str(), accinfo.url).await {
            Ok(result) => result,
            Err(_err) => {
                // Add accession to failed downloads with each moltype
                if !proteomes_only {
                    let failed_download_dna = FailedDownload {
                        accession: accession.clone(),
                        name: name.clone(),
                        moltype: "dna".to_string(),
                        md5sum: None,
                        download_filename: None,
                        url: None,
                    };
                    download_failures.push(failed_download_dna);
                }
                if !genomes_only {
                    let failed_download_protein = FailedDownload {
                        accession: accession.clone(),
                        name: name.clone(),
                        moltype: "protein".to_string(),
                        md5sum: None,
                        download_filename: None,
                        url: None,
                    };
                    download_failures.push(failed_download_protein);
                }

                return Ok((built_sigs, download_failures, checksum_failures));
            }
        };
    let md5sum_url = GenBankFileType::Checksum.url(&base_url, &full_name);

    let mut file_types = vec![
        GenBankFileType::Genomic,
        GenBankFileType::Protein,
        // GenBankFileType::AssemblyReport,
    ];
    if genomes_only {
        file_types = vec![GenBankFileType::Genomic];
    } else if proteomes_only {
        file_types = vec![GenBankFileType::Protein];
    }

    let checksums = match download_and_parse_md5(client, &md5sum_url).await {
        Ok(cs) => cs,
        Err(err) => {
            // capture the error message as a string
            let error_message = err.to_string();
            // if we can't download/parse the md5sum file, write to checksum failures file to allow manual troubleshooting
            for file_type in &file_types {
                // get filename, filetype info to facilitate downstream
                let url = file_type.url(&base_url, &full_name);
                let file_name = file_type.filename_to_write(&accession);
                let failed_checksum_download: FailedChecksum = FailedChecksum {
                    accession: accession.clone(),
                    name: name.clone(),
                    moltype: file_type.moltype(),
                    md5sum_url: Some(md5sum_url.clone()),
                    download_filename: Some(file_name),
                    url: Some(url),
                    expected_md5sum: None,
                    reason: error_message.clone(), // write full error message
                };
                checksum_failures.push(failed_checksum_download);
            }
            // return early from function b/c we can't check any checksums
            return Ok((built_sigs, download_failures, checksum_failures));
        }
    };

    for file_type in &file_types {
        let url = file_type.url(&base_url, &full_name);
        let expected_md5 = checksums.get(&file_type.server_filename(&full_name));
        let file_name = file_type.filename_to_write(&accession);
        let data =
            match download_with_retry(client, &url, expected_md5.map(|x| x.as_str()), retry_count)
                .await
            {
                Ok(data) => data,
                Err(e) => {
                    let error_message = e.to_string();
                    // did we have a checksum error or a download error?
                    // here --> keep track of accession errors + filetype
                    if error_message.contains("MD5 hash does not match") {
                        let checksum_mismatch: FailedChecksum = FailedChecksum {
                            accession: accession.clone(),
                            name: name.clone(),
                            moltype: file_type.moltype(),
                            md5sum_url: Some(md5sum_url.clone()),
                            download_filename: Some(file_name.clone()),
                            url: Some(url.clone()),
                            expected_md5sum: expected_md5.cloned(),
                            reason: error_message.clone(),
                        };
                        checksum_failures.push(checksum_mismatch);
                    } else {
                        let failed_download = FailedDownload {
                            accession: accession.clone(),
                            name: name.clone(),
                            moltype: file_type.moltype(),
                            md5sum: expected_md5.map(|x| x.to_string()),
                            download_filename: Some(file_name),
                            url: Some(url),
                        };
                        download_failures.push(failed_download);
                    }
                    continue;
                }
            };

        if keep_fastas {
            let path = location.join(&file_name);
            fs::write(&path, &data).context("Failed to write data to file")?;
        }
        if !download_only {
            // sketch data
            match file_type {
                GenBankFileType::Genomic => {
                    dna_sigs.build_sigs_from_data(data, "dna", name.clone(), file_name.clone())?;
                    built_sigs.add_collection(dna_sigs);
                }
                GenBankFileType::Protein => {
                    prot_sigs.build_sigs_from_data(
                        data,
                        "protein",
                        name.clone(),
                        file_name.clone(),
                    )?;
                    built_sigs.add_collection(prot_sigs);
                }
                _ => {} // Do nothing for other file types
            };
        }
    }

    Ok((built_sigs, download_failures, checksum_failures))
}

#[allow(clippy::too_many_arguments)]
async fn dl_sketch_url(
    client: &Client,
    accinfo: AccessionData,
    location: &PathBuf,
    retry: Option<u32>,
    _keep_fastas: bool,
    dna_sigs: &mut BuildCollection,
    prot_sigs: &mut BuildCollection,
    _genomes_only: bool,
    _proteomes_only: bool,
    download_only: bool,
) -> Result<(
    MultiBuildCollection,
    Vec<FailedDownload>,
    Vec<FailedChecksum>,
)> {
    let retry_count = retry.unwrap_or(3); // Default retry count
    let mut built_sigs = MultiBuildCollection::new();
    let mut download_failures = Vec::<FailedDownload>::new();
    let mut checksum_failures = Vec::<FailedChecksum>::new();

    let name = accinfo.name;
    let accession = accinfo.accession;
    let url = accinfo.url;
    let expected_md5 = accinfo.expected_md5sum;
    let download_filename = accinfo.download_filename;
    let moltype = accinfo.moltype;

    match download_with_retry(client, &url, expected_md5.as_deref(), retry_count).await {
        Ok(data) => {
            // check keep_fastas instead??
            if let Some(ref download_filename) = download_filename {
                let path = location.join(download_filename);
                fs::write(path, &data).context("Failed to write data to file")?;
            }
            if !download_only {
                let filename = download_filename.clone().unwrap_or("".to_string());
                // sketch data
                match moltype {
                    InputMolType::Dna => {
                        dna_sigs.build_sigs_from_data(
                            data,
                            "dna",
                            name.clone(),
                            filename.clone(),
                        )?;
                        built_sigs.add_collection(dna_sigs);
                    }
                    InputMolType::Protein => {
                        prot_sigs.build_sigs_from_data(
                            data,
                            "protein",
                            name.clone(),
                            filename.clone(),
                        )?;
                        built_sigs.add_collection(prot_sigs);
                    }
                };
            }
        }
        Err(err) => {
            let error_message = err.to_string();
            // did we have a checksum error or a download error?
            // here --> keep track of accession errors + filetype
            if error_message.contains("MD5 hash does not match") {
                let checksum_mismatch: FailedChecksum = FailedChecksum {
                    accession: accession.clone(),
                    name: name.clone(),
                    moltype: moltype.to_string(),
                    md5sum_url: None,
                    download_filename,
                    url: Some(url.clone()),
                    expected_md5sum: expected_md5.clone(),
                    reason: error_message.clone(),
                };
                checksum_failures.push(checksum_mismatch);
            } else {
                let failed_download = FailedDownload {
                    accession: accession.clone(),
                    name: name.clone(),
                    moltype: moltype.to_string(),
                    md5sum: expected_md5.map(|x| x.to_string()),
                    download_filename,
                    url: Some(url),
                };
                download_failures.push(failed_download);
            }
        }
    }

    Ok((built_sigs, download_failures, checksum_failures))
}

// Load existing batch files into MultiCollection, skipping corrupt files
async fn load_existing_zip_batches(outpath: &PathBuf) -> Result<(MultiCollection, usize)> {
    // Remove the .zip extension to get the base name
    let outpath_base = outpath.with_extension("");

    // Regex to match the exact zip filename and its batches (e.g., "outpath.zip", "outpath.1.zip", "outpath.2.zip", etc.)
    let zip_file_pattern = Regex::new(&format!(
        r"^{}(?:\.(\d+))?\.zip$",
        regex::escape(outpath_base.file_name().unwrap())
    ))
    .unwrap();

    // Initialize a vector to store valid collections
    let mut collections = Vec::new();
    let mut highest_batch = 0; // Track the highest batch number

    // Read the directory containing the outpath
    let dir = outpath_base
        .parent()
        .ok_or_else(|| anyhow::anyhow!("Could not get parent directory"))?;
    let mut dir_entries = tokio::fs::read_dir(dir).await?;

    // Scan through all files in the directory
    while let Some(entry) = dir_entries.next_entry().await? {
        let entry_path: PathBuf = entry.path().try_into()?;

        if let Some(file_name) = entry_path.file_name() {
            // Check if the file matches the base zip file or any batched zip file (outpath.zip, outpath.1.zip, etc.)
            if let Some(captures) = zip_file_pattern.captures(file_name) {
                // Wrap the `from_zipfile` call in `catch_unwind` to prevent panic propagation
                let result = panic::catch_unwind(|| Collection::from_zipfile(&entry_path));
                match result {
                    Ok(Ok(collection)) => {
                        // Successfully loaded the collection, push to `collections`
                        collections.push(collection);

                        // Extract the batch number (if it exists) and update the highest_batch
                        if let Some(batch_str) = captures.get(1) {
                            if let Ok(batch_num) = batch_str.as_str().parse::<usize>() {
                                highest_batch = max(highest_batch, batch_num);
                            }
                        }
                    }
                    Ok(Err(e)) => {
                        // Handle the case where `from_zipfile` returned an error
                        eprintln!(
                            "Warning: Failed to load zip file '{}'. Error: {:?}",
                            entry_path, e
                        );
                        continue; // Skip the file and continue
                    }
                    Err(_) => {
                        // The code inside `from_zipfile` panicked
                        eprintln!("Warning: Invalid zip file '{}'; skipping.", entry_path);
                        continue; // Skip the file and continue
                    }
                }
            }
        }
    }
    // Return the loaded MultiCollection and the max batch index, even if no collections were found
    Ok((MultiCollection::new(collections), highest_batch))
}

/// create zip file depending on batch size and index.
async fn create_or_get_zip_file(
    outpath: &PathBuf,
    batch_size: usize,
    batch_index: usize,
) -> Result<ZipFileWriter<Compat<File>>, anyhow::Error> {
    let batch_outpath = if batch_size == 0 {
        // If batch size is zero, use provided outpath (contains .zip extension)
        outpath.clone()
    } else {
        // Otherwise, modify outpath to include the batch index
        let outpath_base = outpath.with_extension(""); // remove .zip extension
        outpath_base.with_file_name(format!(
            "{}.{}.zip",
            outpath_base.file_stem().unwrap(),
            batch_index
        ))
    };
    let file = File::create(&batch_outpath)
        .await
        .with_context(|| format!("Failed to create file: {:?}", batch_outpath))?;

    Ok(ZipFileWriter::with_tokio(file))
}

pub fn zipwriter_handle(
    mut recv_sigs: tokio::sync::mpsc::Receiver<MultiBuildCollection>,
    output_sigs: Option<String>,
    batch_size: usize,      // Tunable batch size
    mut batch_index: usize, // starting batch index
    error_sender: tokio::sync::mpsc::Sender<anyhow::Error>,
) -> tokio::task::JoinHandle<()> {
    tokio::spawn(async move {
        let mut md5sum_occurrences = HashMap::new();
        let mut zip_manifest = BuildManifest::new();
        let mut wrote_sigs = false;
        let mut file_count = 0; // Count of files in the current batch
        let mut zip_writer = None;

        if let Some(outpath) = output_sigs {
            let outpath: PathBuf = outpath.into();

            while let Some(mut multibuildcoll) = recv_sigs.recv().await {
                if zip_writer.is_none() {
                    // create zip file if needed
                    zip_writer =
                        match create_or_get_zip_file(&outpath, batch_size, batch_index).await {
                            Ok(writer) => Some(writer),
                            Err(e) => {
                                let _ = error_sender.send(e).await;
                                return;
                            }
                        };
                }

                if let Some(zip_writer) = zip_writer.as_mut() {
                    // write all sigs from sigcoll. Note that this method updates each record's internal location
                    for sigcoll in &mut multibuildcoll.collections {
                        match sigcoll
                            .async_write_sigs_to_zip(zip_writer, &mut md5sum_occurrences)
                            .await
                        {
                            Ok(_) => {
                                file_count += sigcoll.size();
                                wrote_sigs = true;
                            }
                            Err(e) => {
                                let error = e.context("Error processing signature");
                                if error_sender.send(error).await.is_err() {
                                    return;
                                }
                            }
                        }
                        // Add all records from sigcoll manifest
                        zip_manifest.extend_from_manifest(&sigcoll.manifest);
                    }
                }

                // if batch size is non-zero and is reached, close the current zip
                if batch_size > 0 && file_count >= batch_size {
                    eprintln!("writing batch {}", batch_index);
                    if let Some(mut zip_writer) = zip_writer.take() {
                        if let Err(e) = zip_manifest
                            .async_write_manifest_to_zip(&mut zip_writer)
                            .await
                        {
                            let _ = error_sender.send(e).await;
                        }
                        if let Err(e) = zip_writer.close().await {
                            let error = anyhow::Error::new(e).context("Failed to close ZIP file");
                            let _ = error_sender.send(error).await;
                            return;
                        }
                    }
                    // Start a new batch
                    batch_index += 1;
                    file_count = 0;
                    zip_manifest.clear();
                    zip_writer = None; // reset zip_writer so a new zip will be created when needed
                }
            }

            if file_count > 0 {
                // write the final manifest
                if let Some(mut zip_writer) = zip_writer.take() {
                    if let Err(e) = zip_manifest
                        .async_write_manifest_to_zip(&mut zip_writer)
                        .await
                    {
                        let _ = error_sender.send(e).await;
                    }

                    // close final zip file
                    if let Err(e) = zip_writer.close().await {
                        let error = anyhow::Error::new(e).context("Failed to close ZIP file");
                        let _ = error_sender.send(error).await;
                        return;
                    }
                }
            }
            if !wrote_sigs {
                // If no signatures were written at all
                let error = anyhow::Error::new(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    "No signatures written",
                ));
                let _ = error_sender.send(error).await;
            }
        }
    })
}

pub fn failures_handle(
    failed_csv: String,
    mut recv_failed: tokio::sync::mpsc::Receiver<FailedDownload>,
    error_sender: tokio::sync::mpsc::Sender<Error>,
) -> tokio::task::JoinHandle<()> {
    tokio::spawn(async move {
        match File::create(&failed_csv).await {
            Ok(file) => {
                let mut writer = BufWriter::new(file);

                // Attempt to write CSV headers
                if let Err(e) = writer
                    .write_all(b"accession,name,moltype,md5sum,download_filename,url\n")
                    .await
                {
                    let error = Error::new(e).context("Failed to write headers");
                    let _ = error_sender.send(error).await;
                    return; // Exit the task early after reporting the error
                }

                while let Some(FailedDownload {
                    accession,
                    name,
                    md5sum,
                    download_filename,
                    url,
                    moltype,
                }) = recv_failed.recv().await
                {
                    let record = format!(
                        "{},{},{},{},{},{}\n",
                        accession,
                        name,
                        moltype,
                        md5sum.unwrap_or("".to_string()),
                        download_filename.unwrap_or("".to_string()),
                        url.map(|u| u.to_string()).unwrap_or("".to_string())
                    );
                    // Attempt to write each record
                    if let Err(e) = writer.write_all(record.as_bytes()).await {
                        let error = Error::new(e).context("Failed to write record");
                        let _ = error_sender.send(error).await;
                        continue; // Optionally continue to try to write next records
                    }
                }

                // Attempt to flush the writer
                if let Err(e) = writer.flush().await {
                    let error = Error::new(e).context("Failed to flush writer");
                    let _ = error_sender.send(error).await;
                }
            }
            Err(e) => {
                let error = Error::new(e).context("Failed to create file");
                let _ = error_sender.send(error).await;
            }
        }
        drop(error_sender);
    })
}

pub fn checksum_failures_handle(
    checksum_failed_csv: String,
    mut recv_failed: tokio::sync::mpsc::Receiver<FailedChecksum>,
    error_sender: tokio::sync::mpsc::Sender<Error>, // Additional parameter for error channel
) -> tokio::task::JoinHandle<()> {
    tokio::spawn(async move {
        match File::create(&checksum_failed_csv).await {
            Ok(file) => {
                let mut writer = BufWriter::new(file);

                // Attempt to write CSV headers
                if let Err(e) = writer
                    .write_all(b"accession,name,moltype,md5sum_url,download_filename,url,expected_md5sum,reason\n")
                    .await
                {
                    let error = Error::new(e).context("Failed to write headers");
                    let _ = error_sender.send(error).await;
                    return; // Exit the task early after reporting the error
                }

                while let Some(FailedChecksum {
                    accession,
                    name,
                    moltype,
                    md5sum_url,
                    download_filename,
                    url,
                    expected_md5sum,
                    reason,
                }) = recv_failed.recv().await
                {
                    let record = format!(
                        "{},{},{},{},{},{},{},{}\n",
                        accession,
                        name,
                        moltype,
                        md5sum_url.map(|u| u.to_string()).unwrap_or("".to_string()),
                        download_filename.unwrap_or("".to_string()),
                        url.map(|u| u.to_string()).unwrap_or("".to_string()),
                        expected_md5sum.unwrap_or("".to_string()),
                        reason,
                    );
                    // Attempt to write each record
                    if let Err(e) = writer.write_all(record.as_bytes()).await {
                        let error = Error::new(e).context("Failed to write failed checksum record");
                        let _ = error_sender.send(error).await;
                        continue; // continue to try to write next records
                    }
                }

                // Attempt to flush the writer
                if let Err(e) = writer.flush().await {
                    let error = Error::new(e).context("Failed to flush failed checksum writer");
                    let _ = error_sender.send(error).await;
                }
            }
            Err(e) => {
                let error = Error::new(e).context("Failed to create failed checksum file");
                let _ = error_sender.send(error).await;
            }
        }
        drop(error_sender);
    })
}

pub fn error_handler(
    mut recv_errors: tokio::sync::mpsc::Receiver<anyhow::Error>,
    error_flag: Arc<AtomicBool>,
) -> tokio::task::JoinHandle<()> {
    tokio::spawn(async move {
        while let Some(error) = recv_errors.recv().await {
            eprintln!("Error: {}", error);
            if error.to_string().contains("No signatures written") {
                error_flag.store(true, Ordering::SeqCst);
                break;
            }
        }
    })
}

#[tokio::main]
#[allow(clippy::too_many_arguments)]
pub async fn gbsketch(
    py: Python,
    input_csv: String,
    param_str: String,
    failed_csv: String,
    failed_checksums_csv: String,
    retry_times: u32,
    fasta_location: String,
    keep_fastas: bool,
    genomes_only: bool,
    proteomes_only: bool,
    download_only: bool,
    batch_size: u32,
    output_sigs: Option<String>,
) -> Result<(), anyhow::Error> {
    let batch_size = batch_size as usize;
    let mut batch_index = 1;
    let mut name_params_map: HashMap<String, HashSet<u64>> = HashMap::new();
    let mut filter = false;
    if let Some(ref output_sigs) = output_sigs {
        // Create outpath from output_sigs
        let outpath = PathBuf::from(output_sigs);

        // Check if the extension is "zip"
        if outpath.extension().map_or(true, |ext| ext != "zip") {
            bail!("Output must be a zip file.");
        }
        // find and read any existing sigs
        let (existing_sigs, max_existing_batch_index) = load_existing_zip_batches(&outpath).await?;
        // Check if there are any existing batches to process
        if !existing_sigs.is_empty() {
            name_params_map = existing_sigs.buildparams_hashmap();

            batch_index = max_existing_batch_index + 1;
            eprintln!(
                "Found {} existing valid zip batch(es). Starting new sig writing at batch {}",
                max_existing_batch_index, batch_index
            );
            filter = true;
        } else {
            // No existing batches, skipping signature filtering
            eprintln!("No valid existing signature batches found; building all signatures.");
        }
    }

    // set up fasta download path
    let download_path = PathBuf::from(fasta_location);
    if !download_path.exists() {
        create_dir_all(&download_path)?;
    }
    // create channels. buffer size here is 4 b/c we can do 3 downloads simultaneously
    let (send_sigs, recv_sigs) = tokio::sync::mpsc::channel::<MultiBuildCollection>(4);
    let (send_failed, recv_failed) = tokio::sync::mpsc::channel::<FailedDownload>(4);
    let (send_failed_checksums, recv_failed_checksum) =
        tokio::sync::mpsc::channel::<FailedChecksum>(4);
    // Error channel for handling task errors
    let (error_sender, error_receiver) = tokio::sync::mpsc::channel::<anyhow::Error>(1);

    // Set up collector/writing tasks
    let mut handles = Vec::new();

    let sig_handle = zipwriter_handle(
        recv_sigs,
        output_sigs,
        batch_size,
        batch_index,
        error_sender.clone(),
    );

    let failures_handle = failures_handle(failed_csv, recv_failed, error_sender.clone());
    let checksum_failures_handle = checksum_failures_handle(
        failed_checksums_csv,
        recv_failed_checksum,
        error_sender.clone(),
    );
    let critical_error_flag = Arc::new(AtomicBool::new(false));
    let error_handle = error_handler(error_receiver, critical_error_flag.clone());
    handles.push(sig_handle);
    handles.push(failures_handle);
    handles.push(error_handle);
    handles.push(checksum_failures_handle);

    // Worker tasks
    let semaphore = Arc::new(Semaphore::new(3)); // Limiting concurrent downloads
    let client = Arc::new(Client::new());

    // Open the file containing the accessions synchronously
    let (accession_info, n_accs) = load_gbassembly_info(input_csv)?;
    if n_accs == 0 {
        bail!("No accessions to download and sketch.")
    }

    // parse param string into params_vec, print error if fail
    let param_result = parse_params_str(param_str);
    let params_vec = match param_result {
        Ok(params) => params,
        Err(e) => {
            bail!("Failed to parse params string: {}", e);
        }
    };
    let dna_template_collection = BuildCollection::from_buildparams(&params_vec, "DNA");
    // prot will build protein, dayhoff, hp
    let prot_template_collection = BuildCollection::from_buildparams(&params_vec, "protein");

    let mut genomes_only = genomes_only;
    let mut proteomes_only = proteomes_only;

    // Check if dna_sig_templates is empty and not keep_fastas
    if dna_template_collection.manifest.is_empty() && !keep_fastas {
        eprintln!("No DNA signature templates provided, and --keep-fasta is not set.");
        proteomes_only = true;
    }
    // Check if protein_sig_templates is empty and not keep_fastas
    if prot_template_collection.manifest.is_empty() && !keep_fastas {
        eprintln!("No protein signature templates provided, and --keep-fasta is not set.");
        genomes_only = true;
    }
    if genomes_only {
        if !download_only {
            eprintln!("Downloading and sketching genomes only.");
        } else {
            eprintln!("Downloading genomes only.");
        }
    } else if proteomes_only {
        if !download_only {
            eprintln!("Downloading and sketching proteomes only.");
        } else {
            eprintln!("Downloading proteomes only.");
        }
    }

    // report every 1 percent (or every 1, whichever is larger)
    let reporting_threshold = std::cmp::max(n_accs / 100, 1);

    for (i, accinfo) in accession_info.into_iter().enumerate() {
        py.check_signals()?; // If interrupted, return an Err automatically

        let mut dna_sigs = dna_template_collection.clone();
        let mut prot_sigs = prot_template_collection.clone();

        // filter template sigs based on existing sigs
        if filter {
            if let Some(existing_paramset) = name_params_map.get(&accinfo.name) {
                // If the key exists, filter template sigs
                dna_sigs.filter(existing_paramset);
                prot_sigs.filter(existing_paramset);
            }
        }

        // clone remaining utilities
        let semaphore_clone = Arc::clone(&semaphore);
        let client_clone = Arc::clone(&client);
        let send_sigs = send_sigs.clone();
        let send_failed = send_failed.clone();
        let checksum_send_failed = send_failed_checksums.clone();
        let download_path_clone = download_path.clone(); // Clone the path for each task
        let send_errors = error_sender.clone();

        tokio::spawn(async move {
            let _permit = semaphore_clone.acquire().await;
            // progress report when the permit is available and processing begins
            if (i + 1) % reporting_threshold == 0 {
                let percent_processed = (((i + 1) as f64 / n_accs as f64) * 100.0).round();
                println!(
                    "Starting accession {}/{} ({}%)",
                    (i + 1),
                    n_accs,
                    percent_processed
                );
            }
            // Perform download and sketch
            let result = dl_sketch_assembly_accession(
                &client_clone,
                accinfo.clone(),
                &download_path_clone,
                Some(retry_times),
                keep_fastas,
                &mut dna_sigs,
                &mut prot_sigs,
                genomes_only,
                proteomes_only,
                download_only,
            )
            .await;
            match result {
                Ok((coll, failed_downloads, failed_checksums)) => {
                    if !coll.is_empty() {
                        if let Err(e) = send_sigs.send(coll).await {
                            eprintln!("Failed to send signatures: {}", e);
                            let _ = send_errors.send(e.into()).await; // Send the error through the channel
                        }
                    }
                    for fail in failed_downloads {
                        if let Err(e) = send_failed.send(fail).await {
                            eprintln!("Failed to send failed download info: {}", e);
                            let _ = send_errors.send(e.into()).await; // Send the error through the channel
                        }
                    }
                    for fail in failed_checksums {
                        if let Err(e) = checksum_send_failed.send(fail).await {
                            eprintln!("Failed to send failed checksum info: {}", e);
                            let _ = send_errors.send(e.into()).await; // Send the error through the channel
                        }
                    }
                }
                Err(e) => {
                    let _ = send_errors.send(e).await;
                }
            }
            drop(send_errors);
        });
    }
    // drop senders as we're done sending data
    drop(send_sigs);
    drop(send_failed);
    drop(send_failed_checksums);
    drop(error_sender);
    // Wait for all tasks to complete
    for handle in handles {
        if let Err(e) = handle.await {
            eprintln!("Handle join error: {}.", e);
        }
    }
    // since the only critical error is not having written any sigs
    // check this here at end. Bail if we wrote expected sigs but wrote none.
    if critical_error_flag.load(Ordering::SeqCst) & !download_only {
        bail!("No signatures written, exiting.");
    }

    Ok(())
}

#[tokio::main]
#[allow(clippy::too_many_arguments)]
pub async fn urlsketch(
    py: Python,
    input_csv: String,
    param_str: String,
    failed_csv: String,
    retry_times: u32,
    fasta_location: String,
    keep_fastas: bool,
    download_only: bool,
    batch_size: u32,
    output_sigs: Option<String>,
    failed_checksums_csv: Option<String>,
) -> Result<(), anyhow::Error> {
    let batch_size = batch_size as usize;
    let mut batch_index = 1;
    let mut name_params_map: HashMap<String, HashSet<u64>> = HashMap::new();
    let mut filter = false;
    if let Some(ref output_sigs) = output_sigs {
        // Create outpath from output_sigs
        let outpath = PathBuf::from(output_sigs);

        // Check if the extension is "zip"
        if outpath.extension().map_or(true, |ext| ext != "zip") {
            bail!("Output must be a zip file.");
        }
        // find and read any existing sigs
        let (existing_sigs, max_existing_batch_index) = load_existing_zip_batches(&outpath).await?;
        // Check if there are any existing batches to process
        if !existing_sigs.is_empty() {
            name_params_map = existing_sigs.buildparams_hashmap();

            batch_index = max_existing_batch_index + 1;
            eprintln!(
                "Found {} existing zip batches. Starting new sig writing at batch {}",
                max_existing_batch_index, batch_index
            );
            filter = true;
        } else {
            // No existing batches, skipping signature filtering
            eprintln!("No existing signature batches found; building all signatures.");
        }
    }

    // set up fasta download path
    let download_path = PathBuf::from(fasta_location);
    if !download_path.exists() {
        create_dir_all(&download_path)?;
    }

    // create channels. buffer size here is 4 b/c we can do 3 downloads simultaneously
    let (send_sigs, recv_sigs) = tokio::sync::mpsc::channel::<MultiBuildCollection>(4);
    let (send_failed, recv_failed) = tokio::sync::mpsc::channel::<FailedDownload>(4);
    let (send_failed_checksums, recv_failed_checksum) =
        tokio::sync::mpsc::channel::<FailedChecksum>(4);
    // Error channel for handling task errors
    let (error_sender, error_receiver) = tokio::sync::mpsc::channel::<anyhow::Error>(1);

    // Set up collector/writing tasks
    let mut handles = Vec::new();

    let sig_handle = zipwriter_handle(
        recv_sigs,
        output_sigs,
        batch_size,
        batch_index,
        error_sender.clone(),
    );

    let failures_handle = failures_handle(failed_csv, recv_failed, error_sender.clone());

    let mut write_failed_checksums = false;
    if let Some(ref failed_checksums) = failed_checksums_csv {
        let checksum_failures_handle = checksum_failures_handle(
            failed_checksums.clone(),
            recv_failed_checksum,
            error_sender.clone(),
        );
        write_failed_checksums = true;
        handles.push(checksum_failures_handle);
    }

    let critical_error_flag = Arc::new(AtomicBool::new(false));
    let error_handle = error_handler(error_receiver, critical_error_flag.clone());
    handles.push(sig_handle);
    handles.push(failures_handle);
    handles.push(error_handle);

    // Worker tasks
    let semaphore = Arc::new(Semaphore::new(3)); // Limiting concurrent downloads
    let client = Arc::new(Client::new());

    // Open the file containing the accessions synchronously
    let (accession_info, n_accs) = load_accession_info(input_csv, keep_fastas)?;
    if n_accs == 0 {
        bail!("No accessions to download and sketch.")
    }

    // parse param string into params_vec, print error if fail
    let param_result = parse_params_str(param_str);
    let params_vec = match param_result {
        Ok(params) => params,
        Err(e) => {
            bail!("Failed to parse params string: {}", e);
        }
    };
    let dna_template_collection = BuildCollection::from_buildparams(&params_vec, "DNA");
    let prot_template_collection = BuildCollection::from_buildparams(&params_vec, "protein");

    let mut genomes_only = false;
    let mut proteomes_only = false;

    // Check if dna_sig_templates is empty and not keep_fastas
    if dna_template_collection.manifest.is_empty() && !keep_fastas {
        eprintln!("No DNA signature templates provided, and --keep-fastas is not set.");
        proteomes_only = true;
    }
    // Check if protein_sig_templates is empty and not keep_fastas
    if prot_template_collection.manifest.is_empty() && !keep_fastas {
        eprintln!("No protein signature templates provided, and --keep-fastas is not set.");
        genomes_only = true;
    }
    if genomes_only {
        if !download_only {
            eprintln!("Downloading and sketching genomes only.");
        } else {
            eprintln!("Downloading genomes only.");
        }
    } else if proteomes_only {
        if !download_only {
            eprintln!("Downloading and sketching proteomes only.");
        } else {
            eprintln!("Downloading proteomes only.");
        }
    }

    // report every 1 percent (or every 1, whichever is larger)
    let reporting_threshold = std::cmp::max(n_accs / 100, 1);

    for (i, accinfo) in accession_info.into_iter().enumerate() {
        py.check_signals()?; // If interrupted, return an Err automatically
        let mut dna_sigs = dna_template_collection.clone();
        let mut prot_sigs = prot_template_collection.clone();

        // filter template sigs based on existing sigs
        if filter {
            if let Some(existing_paramset) = name_params_map.get(&accinfo.name) {
                // If the key exists, filter template sigs
                dna_sigs.filter(existing_paramset);
                prot_sigs.filter(existing_paramset);
            }
        }

        let semaphore_clone = Arc::clone(&semaphore);
        let client_clone = Arc::clone(&client);
        let send_sigs = send_sigs.clone();
        let send_failed = send_failed.clone();
        let checksum_send_failed = send_failed_checksums.clone();
        let download_path_clone = download_path.clone(); // Clone the path for each task
        let send_errors = error_sender.clone();

        tokio::spawn(async move {
            let _permit = semaphore_clone.acquire().await;
            // progress report when the permit is available and processing begins
            if (i + 1) % reporting_threshold == 0 {
                let percent_processed = (((i + 1) as f64 / n_accs as f64) * 100.0).round();
                println!(
                    "Starting accession {}/{} ({}%)",
                    (i + 1),
                    n_accs,
                    percent_processed
                );
            }
            // Perform download and sketch
            let result = dl_sketch_url(
                &client_clone,
                accinfo.clone(),
                &download_path_clone,
                Some(retry_times),
                keep_fastas,
                &mut dna_sigs,
                &mut prot_sigs,
                genomes_only,
                proteomes_only,
                download_only,
            )
            .await;
            match result {
                Ok((sigs, failed_downloads, failed_checksums)) => {
                    if !sigs.is_empty() {
                        if let Err(e) = send_sigs.send(sigs).await {
                            eprintln!("Failed to send signatures: {}", e);
                            let _ = send_errors.send(e.into()).await; // Send the error through the channel
                        }
                    }
                    for fail in failed_downloads {
                        if let Err(e) = send_failed.send(fail).await {
                            eprintln!("Failed to send failed download info: {}", e);
                            let _ = send_errors.send(e.into()).await; // Send the error through the channel
                        }
                    }
                    if write_failed_checksums {
                        for fail in failed_checksums {
                            if let Err(e) = checksum_send_failed.send(fail).await {
                                eprintln!("Failed to send failed checksum info: {}", e);
                                let _ = send_errors.send(e.into()).await; // Send the error through the channel
                            }
                        }
                    } else {
                        // if we don't have a failed checksum file, convert to failed downloads + write there
                        for fail in failed_checksums {
                            let dl_fail: FailedDownload = FailedDownload {
                                accession: fail.accession,
                                name: fail.name,
                                moltype: fail.moltype,
                                md5sum: fail.expected_md5sum,
                                download_filename: fail.download_filename,
                                url: fail.url,
                            };
                            if let Err(e) = send_failed.send(dl_fail).await {
                                eprintln!("Failed to send failed download info: {}", e);
                                let _ = send_errors.send(e.into()).await; // Send the error through the channel
                            }
                        }
                    }
                }
                Err(e) => {
                    let _ = send_errors.send(e).await;
                }
            }
            drop(send_errors);
        });
    }
    // drop senders as we're done sending data
    drop(send_sigs);
    drop(send_failed);
    drop(error_sender);
    drop(send_failed_checksums);
    // Wait for all tasks to complete
    for handle in handles {
        if let Err(e) = handle.await {
            eprintln!("Handle join error: {}.", e);
        }
    }
    // since the only critical error is not having written any sigs
    // check this here at end. Bail if we wrote expected sigs but wrote none.
    if critical_error_flag.load(Ordering::SeqCst) & !download_only {
        bail!("No signatures written, exiting.");
    }

    Ok(())
}
