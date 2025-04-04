use anyhow::{anyhow, bail, Context, Error, Result};
use async_zip::base::write::ZipFileWriter;
use camino::Utf8Path;
use camino::Utf8PathBuf as PathBuf;
use futures::stream::{self, StreamExt};
use needletail::parse_fastx_reader;
use needletail::parser::SequenceRecord;
use regex::Regex;
use reqwest::{header::HeaderMap, Client};
use serde_json::json;
use sourmash::collection::Collection;
use std::cmp::max;
use std::collections::{HashMap, HashSet};
use std::fs::create_dir_all;
use std::io::{Cursor, Read, Write};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Arc;
use tokio::fs::{File, OpenOptions};
use tokio::io::{duplex, AsyncWriteExt, BufWriter};
use tokio::sync::mpsc::Sender;
use tokio::sync::{oneshot, Semaphore};
use tokio::task::spawn_blocking;
use tokio_util::compat::Compat;
use tokio_util::io::SyncIoBridge;
use tokio_util::sync::CancellationToken;
use zip::read::ZipArchive;

use pyo3::prelude::*;

use crate::utils::{
    load_accession_info, load_gbsketch_info, AccessionData, FailedChecksum, InputMolType,
    MultiCollection, TempFastaFile, ToCsvRow,
};

use crate::utils::buildutils::{BuildCollection, BuildManifest, MultiSelect, MultiSelection};
use reqwest::Url;

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

/// Processes FASTX records from the stream. Writes FASTA entries if a file is given
/// and adds sequences to signatures in the provided `BuildCollection`.
/// Note, if internal checksum (e.g. CRC32) or decompression issues occur, it will return an error.
pub fn process_fastx_from_reader<R: std::io::Read + Send + 'static>(
    reader: R,
    mut file: Option<std::fs::File>,
    mut sigs: BuildCollection,
    moltype: &str,
    name: String,
    filename: String,
    range: Option<(usize, usize)>,
) -> Result<BuildCollection> {
    let mut fastx_reader = parse_fastx_reader(reader)
        .context("Failed to parse FASTX data (possibly CRC32 or decompression issue)")?;

    while let Some(record) = fastx_reader.next() {
        let record = record.context("Failed to read FASTX record")?;
        let subseq = extract_range_from_record(&record, range)?;

        // Optionally write the sequence in FASTA format
        if let Some(ref mut file) = file {
            let fasta_entry = format!(
                ">{}\n{}\n",
                String::from_utf8_lossy(record.id()),
                String::from_utf8_lossy(&subseq)
            );
            file.write_all(fasta_entry.as_bytes())
                .context("Failed to write FASTA entry")?;
        }

        // Add the sequence to all matching records/signatures. Will not impact empty collection.
        sigs.add_sequence(moltype, &subseq)?;
    }

    // Update the sig/manifest info once per call. Will not impact empty collection.
    sigs.update_info(name, filename);
    Ok(sigs)
}

#[allow(clippy::too_many_arguments)]
pub async fn stream_and_process_with_retry(
    client: &Client,
    accinfo: AccessionData,
    location: &PathBuf,
    retry: Option<u32>,
    keep_fastas: bool,
    mut sigs: BuildCollection,
    _download_only: bool,
    write_checksum_fail: bool,
    cancel_token: CancellationToken,
) -> Result<(BuildCollection, Vec<AccessionData>, Vec<FailedChecksum>)> {
    if cancel_token.is_cancelled() {
        return Err(anyhow!("Aborting due to cancellation"));
    }

    let retry_count = retry.unwrap_or(3);
    let mut download_failures = Vec::new();
    let mut checksum_failures = Vec::new();
    let empty_coll = BuildCollection::new();

    let name = accinfo.name.clone();
    let accession = accinfo.accession.clone();
    let download_filename = accinfo.download_filename.clone().unwrap_or_default();
    let moltype = accinfo.moltype.clone();
    let filename = download_filename.clone();
    let merged_sample = accinfo.url_info.len() > 1;

    let mut temp_fasta: Option<TempFastaFile> = None;
    let mut file: Option<std::fs::File> = if keep_fastas {
        let (fasta_file, handle) = TempFastaFile::new(location, &download_filename)?;
        temp_fasta = Some(fasta_file);
        Some(handle)
    } else {
        None
    };

    for uinfo in &accinfo.url_info {
        if cancel_token.is_cancelled() {
            return Err(anyhow!("Cancelled before starting new URL"));
        }
        let url = &uinfo.url;
        let expected_md5 = &uinfo.md5sum;
        let range = uinfo.range;

        let mut attempts_left = retry_count;
        let mut last_error: Option<anyhow::Error> = None;

        while attempts_left > 0 {
            if cancel_token.is_cancelled() {
                return Err(anyhow!("Cancelled before starting new URL"));
            }
            let response = match client.get(url.clone()).send().await {
                Ok(resp) if resp.status().is_success() => resp,
                Ok(resp) => {
                    last_error = Some(anyhow!("HTTP error {} from {}", resp.status(), url));
                    attempts_left -= 1;
                    continue;
                }
                Err(e) => {
                    last_error = Some(anyhow!(e).context("Failed to send request"));
                    attempts_left -= 1;
                    continue;
                }
            };

            // create tokio duplex stream for reading and writing
            let (mut writer, reader) = duplex(1024 * 1024);
            let (tx, rx) = oneshot::channel();

            let moltype_cloned = moltype.as_str();
            let name_cloned = name.clone();
            let filename_cloned = filename.clone();
            let file_taken = file.take();
            let sigs_taken = std::mem::take(&mut sigs);

            let reader = SyncIoBridge::new(reader);
            spawn_blocking(move || {
                let result = process_fastx_from_reader(
                    reader,
                    file_taken,
                    sigs_taken,
                    moltype_cloned,
                    name_cloned,
                    filename_cloned,
                    range,
                );
                let _ = tx.send(result);
            });

            let mut stream = response.bytes_stream();

            while let Some(chunk) = stream.next().await {
                if cancel_token.is_cancelled() {
                    return Err(anyhow!("Cancelled during response streaming"));
                }
                if let Err(e) = pyo3::Python::with_gil(|py| py.check_signals()) {
                    // exit if a Python interrupt signal is caught
                    cancel_token.cancel();
                    return Err(anyhow::anyhow!(
                        "Caught Python interrupt: {}, shutting down gracefully",
                        e
                    ));
                }
                match chunk {
                    Ok(bytes) => {
                        // "write" bytes into the in-memory stream
                        if let Err(e) = writer.write_all(&bytes).await {
                            last_error = Some(anyhow!(e).context("Failed to stream FASTX data"));
                            break;
                        }
                    }
                    Err(e) => {
                        last_error = Some(anyhow!(e).context("Failed to read response chunk"));
                        break;
                    }
                }
            }

            writer.shutdown().await.ok();

            if cancel_token.is_cancelled() {
                return Err(anyhow!("Cancelled before starting new URL"));
            }

            match rx.await {
                Ok(Ok(updated_sigs)) => {
                    sigs = updated_sigs;
                    if let Some(fasta) = temp_fasta.take() {
                        fasta.finalize()?;
                    }
                    break;
                }
                Ok(Err(e)) => {
                    last_error = Some(anyhow!(e).context("Failed to process FASTX data"));
                    attempts_left -= 1;
                    if let Some(fasta) = &temp_fasta {
                        fasta.cleanup().await;
                    }
                }
                Err(_) => {
                    last_error = Some(anyhow!("Processing task failed or was cancelled"));
                    if let Some(fasta) = &temp_fasta {
                        fasta.cleanup().await;
                    }
                    attempts_left -= 1;
                }
            }
        }

        // If we exhausted retries for this URL
        if attempts_left == 0 {
            let err_msg = last_error
                .as_ref()
                .map(ToString::to_string)
                .unwrap_or_else(|| "Unknown failure".to_string());

            if err_msg.contains("CRC") || err_msg.contains("decompression") {
                if write_checksum_fail {
                    checksum_failures.push(FailedChecksum::new(
                        accession.clone(),
                        name.clone(),
                        moltype.to_string(),
                        None,
                        Some(download_filename.clone()),
                        Some(url.clone()),
                        expected_md5.clone(),
                        err_msg.clone(),
                    ));
                }

                // Mark the whole merged sample as failed, in addition to writing the checksum failure
                if merged_sample {
                    download_failures.push(accinfo);
                }
            } else {
                download_failures.push(accinfo);
            }

            // If any download completely fails, we should stop processing further URLs (rather than getting incomplete merged sigs/files)
            return Ok((empty_coll, download_failures, checksum_failures));
        }
    }

    sigs.update_info(name, filename);

    Ok((sigs, download_failures, checksum_failures))
}

pub async fn download_parse_dehydrated_zip_with_retry(
    client: &Client,
    url: Url,
    headers: HeaderMap,
    params: serde_json::Value,
    retry_count: u32,
) -> Result<
    HashMap<String, String>, // Direct download links
> {
    let mut attempts = retry_count;
    let mut last_error: Option<anyhow::Error> = None;

    while attempts > 0 {
        attempts -= 1;

        let response = client
            .post(url.clone())
            .headers(headers.clone())
            .json(&params)
            .send()
            .await;

        match response {
            Ok(resp) if resp.status().is_success() => match resp.bytes().await {
                Ok(bytes) => match parse_dehydrated_files(&bytes) {
                    Ok(links) => {
                        if links.is_empty() {
                            last_error = Some(anyhow!(
                                        "Parsed ZIP archive successfully, but no download links were found. Are your accessions valid?"
                                    ));
                        } else {
                            return Ok(links);
                        }
                    }
                    Err(e) => {
                        last_error = Some(e.into());
                    }
                },
                Err(e) => {
                    last_error = Some(e.into());
                }
            },
            Ok(resp) => {
                last_error = Some(anyhow!(
                    "Server error status code {}: {}. Retrying...",
                    resp.status(),
                    url
                ));
            }
            Err(e) => {
                last_error = Some(e.into());
            }
        }
        if attempts > 0 {
            eprintln!("Retrying... ({} attempts left)", attempts);
        }
    }

    Err(anyhow!(
        "Failed to download and parse dehydrated ZIP file after {} attempts: {}\n  Last error: {}",
        retry_count,
        url,
        last_error
            .map(|e| e.to_string())
            .unwrap_or_else(|| "Unknown".to_string())
    ))
}

pub fn parse_dehydrated_files(zip_data: &[u8]) -> Result<HashMap<String, String>> // Direct download links
{
    let reader = Cursor::new(zip_data);
    let mut zip = ZipArchive::new(reader).context("Failed to read ZIP archive")?;

    let mut download_links = HashMap::new();

    for i in 0..zip.len() {
        let mut file = zip
            .by_index(i)
            .context("Failed to access file in ZIP archive")?;
        let file_name = file.name().to_string();

        if file_name.ends_with("fetch.txt") {
            let mut content = Vec::new();
            file.read_to_end(&mut content)
                .context(format!("Failed to read contents of {}", file_name))?;
            download_links = parse_fetch_txt(&content)?;
        }
    }

    Ok(download_links)
}

fn parse_fetch_txt(content: &[u8]) -> Result<HashMap<String, String>> {
    let mut download_links = HashMap::new();

    let content_str =
        std::str::from_utf8(content).context("Failed to convert fetch file to UTF-8")?;

    for line in content_str.lines() {
        let columns: Vec<&str> = line.split_whitespace().collect();
        if columns.len() >= 3 {
            let url = columns[0].to_string();
            let file_path = columns[2].trim_start();

            // Extract accession from the file path and normalize the filename
            if let Some((_, rest)) = file_path.split_once("data/") {
                if let Some((accession, filename)) = rest.split_once('/') {
                    let normalized_filename = if filename == "protein.faa" {
                        format!("{}_protein.faa.gz", accession)
                    } else if filename.ends_with("genomic.fna") {
                        format!("{}_genomic.fna.gz", accession)
                    } else {
                        continue;
                    };

                    download_links.insert(normalized_filename, url);
                }
            }
        }
    }

    Ok(download_links)
}

/// Extracts the specified range from sequences in data and writes to the file in FASTA format.
async fn process_and_write_range(
    data: &[u8],
    file: &mut File,
    range: Option<(usize, usize)>,
) -> Result<()> {
    let cursor = Cursor::new(data);
    let mut fastx_reader =
        needletail::parse_fastx_reader(cursor).context("Failed to parse FASTA/FASTQ data")?;

    while let Some(record) = fastx_reader.next() {
        let record = record.context("Failed to read record")?;
        let sequence_to_write = extract_range_from_record(&record, range)
            .context("Failed to extract range from record")?;

        // Use the `id` and `seq` fields directly to construct the FASTA entry
        let fasta_entry = format!(
            ">{}\n{}\n",
            String::from_utf8_lossy(record.id()),
            String::from_utf8_lossy(&sequence_to_write)
        );

        // Write the FASTA entry to the file
        file.write_all(fasta_entry.as_bytes())
            .await
            .context("Failed to write FASTA entry to file")?;
    }

    Ok(())
}

/// Extracts a range from a `SequenceRecord`. Returns the specified sequence slice as a `Vec<u8>`.
fn extract_range_from_record(
    record: &SequenceRecord,
    range: Option<(usize, usize)>,
) -> Result<Vec<u8>> {
    let full_sequence = record.seq();
    if let Some((start, end)) = range {
        let adjusted_start = start.saturating_sub(1); // Adjust for 1-based indexing
        if adjusted_start >= end || end > full_sequence.len() {
            return Err(anyhow::anyhow!(
                "Invalid range: start={}, end={}, sequence length={}",
                start,
                end,
                full_sequence.len()
            ));
        }
        Ok(full_sequence[adjusted_start..end].to_vec())
    } else {
        Ok(full_sequence.to_vec())
    }
}

/// Opens a file for writing, creating necessary directories and truncating it if it exists.
/// Returns an `Option<File>` if a filename is provided, or `None` if the filename is `None`.
async fn open_file_for_writing(
    location: &PathBuf,
    filename: Option<&String>,
) -> Result<Option<File>> {
    if let Some(download_filename) = filename {
        let path = location.join(download_filename);

        // Create subdirectories if needed
        if let Some(parent) = path.parent() {
            create_dir_all(parent).with_context(|| {
                format!(
                    "Failed to create directories for download filename path {}",
                    &path
                )
            })?;
        }

        // Open the file in write mode (truncate if it exists)
        let file = OpenOptions::new()
            .create(true) // Create the file if it doesn't exist
            .write(true) // Enable write mode
            .truncate(true) // Clear existing content
            .open(&path)
            .await
            .with_context(|| format!("Failed to open file at {}", path))?;
        Ok(Some(file))
    } else {
        Ok(None)
    }
}

#[allow(clippy::too_many_arguments)]
async fn dl_sketch_url(
    client: &Client,
    accinfo: AccessionData,
    location: &PathBuf,
    retry: Option<u32>,
    keep_fastas: bool,
    mut sigs: BuildCollection,
    download_only: bool,
    write_checksum_fail: bool,
) -> Result<(BuildCollection, Vec<AccessionData>, Vec<FailedChecksum>)> {
    let retry_count = retry.unwrap_or(3); // Default retry count
    let empty_coll = BuildCollection::new();
    let mut download_failures = Vec::<AccessionData>::new();
    let mut checksum_failures = Vec::<FailedChecksum>::new();

    let name = accinfo.name.clone();
    let accession = accinfo.accession.clone();
    let download_filename = &accinfo.download_filename;
    let filename = download_filename.clone().unwrap_or("".to_string());
    let moltype = &accinfo.moltype;

    let mut file: Option<File> = if keep_fastas {
        open_file_for_writing(location, download_filename.as_ref()).await?
    } else {
        None
    };

    // are we merging files?
    let merged_sample: bool = accinfo.url_info.len() > 1;
    for uinfo in &accinfo.url_info {
        let url = &uinfo.url;
        let expected_md5 = &uinfo.md5sum;
        let range = uinfo.range;
        match download_with_retry(client, url, expected_md5.as_deref(), retry_count).await {
            Ok(data) => {
                // Write to file if keep_fastas is true and a file is open
                // note, if multiple urls are provided, this will append to the same file
                if let Some(file) = file.as_mut() {
                    if range.is_some() {
                        process_and_write_range(&data, file, range)
                            .await
                            .context("Failed to process and write range to file")?;
                    } else {
                        // Write the entire data if no range is provided
                        file.write_all(&data)
                            .await
                            .context("Failed to write data to file")?;
                    }
                }

                if !download_only {
                    // sketch data

                    match moltype {
                        InputMolType::Dna => {
                            sigs.build_sigs_from_data(
                                data,
                                "DNA",
                                name.clone(),
                                filename.clone(),
                                range,
                            )?;
                        }
                        InputMolType::Protein => {
                            sigs.build_sigs_from_data(
                                data,
                                "protein",
                                name.clone(),
                                filename.clone(),
                                range,
                            )?;
                        }
                    };
                }
            }
            Err(err) => {
                let error_message = err.to_string();
                // did we have a checksum error or a download error?
                // here --> keep track of accession errors + filetype
                if error_message.contains("MD5 hash does not match") && write_checksum_fail {
                    let checksum_mismatch: FailedChecksum = FailedChecksum::new(
                        accession.clone(),
                        name.clone(),
                        moltype.to_string(),
                        None,
                        download_filename.clone(),
                        Some(url.clone()),
                        expected_md5.clone(),
                        error_message.clone(),
                    );
                    checksum_failures.push(checksum_mismatch);
                    // if this is a merged sample, the checksum failure is only for one part of it.
                    // also write a download failure, which is the full entry.
                    // The checksum failures file is mostly for debugging, while the failure csv
                    // can be used to re-run urlsketch.
                    if merged_sample {
                        download_failures.push(accinfo);
                    }
                } else {
                    download_failures.push(accinfo);
                }
                // Clear signatures and return immediately on failure
                return Ok((empty_coll, download_failures, checksum_failures));
            }
        }
    }

    // Update signature info
    sigs.update_info(name, filename);

    Ok((sigs, download_failures, checksum_failures))
}

fn get_current_directory() -> Result<PathBuf> {
    let current_dir =
        std::env::current_dir().context("Failed to retrieve the current working directory")?;
    PathBuf::try_from(current_dir)
        .map_err(|_| anyhow::anyhow!("Current directory is not valid UTF-8"))
}

// Load existing batch files into MultiCollection, skipping corrupt files
async fn load_existing_zip_batches(outpath: &PathBuf) -> Result<(MultiCollection, usize)> {
    // Remove the .zip extension to get the base name
    let filename = outpath.file_name().unwrap();
    let (base, suffix) = if filename.ends_with(".sig.zip") {
        (filename.strip_suffix(".sig.zip").unwrap(), "sig.zip")
    } else if filename.ends_with(".zip") {
        (filename.strip_suffix(".zip").unwrap(), "zip")
    } else {
        return Err(anyhow::anyhow!(
            "Output path must end with .zip or .sig.zip"
        ));
    };

    // Regex to match the exact zip filename and its batches (e.g., "outpath.zip" --> "outpath.N.zip"; "outpath.sig.zip" --> "outpath.N.sig.zip", etc.)
    // Also match .incomplete so we can delete these.
    let zip_file_pattern = Regex::new(&format!(
        r"^{}(?:\.(\d+))?\.{}(?:\.incomplete)?$",
        regex::escape(base),
        regex::escape(suffix)
    ))?;

    // Initialize a vector to store valid collections
    let mut collections = Vec::new();
    let mut highest_batch = 0; // Track the highest batch number

    // find parent dir (or use current dir)
    let dir = outpath
        .parent()
        .filter(|p| !p.as_os_str().is_empty()) // Ensure the parent is not empty
        .map(|p| p.to_path_buf()) // Use the parent if it's valid
        .or_else(|| get_current_directory().ok()) // Fallback to current directory if no valid parent
        .ok_or_else(|| anyhow::anyhow!("Failed to determine a valid directory"))?;

    if !dir.exists() {
        return Err(anyhow::anyhow!(
            "Directory for output zipfile does not exist: {}",
            dir
        ));
    }

    let mut dir_entries = tokio::fs::read_dir(dir).await?;

    // get absolute path for outpath
    let current_dir = std::env::current_dir().context("Failed to retrieve current directory")?;
    let outpath_absolute = outpath
        .parent()
        .filter(|parent| parent.as_std_path().exists())
        .map(|_| outpath.clone())
        .unwrap_or_else(|| {
            PathBuf::from_path_buf(current_dir.join(outpath.as_std_path()))
                .expect("Failed to convert to Utf8PathBuf")
        });

    // Scan through all files in the directory
    while let Some(entry) = dir_entries.next_entry().await? {
        let entry_path: PathBuf = entry.path().try_into()?;
        // Skip the `outpath` itself to loading, as we just overwrite this file for now (not append)
        // TO DO: if we can append to the original output file, we can include this and then just add new signatures

        // get absolute path of entry for comparison with outpath_absolute
        let current_dir =
            std::env::current_dir().context("Failed to retrieve current directory")?;
        let entry_absolute = entry_path
            .parent()
            .filter(|parent| parent.as_std_path().exists())
            .map(|_| entry_path.clone())
            .unwrap_or_else(|| {
                PathBuf::from_path_buf(current_dir.join(entry_path.as_std_path()))
                    .expect("Failed to convert to Utf8PathBuf")
            });

        // For now, skip the `outpath` itself to avoid loading, since we will just overwrite it anyway.
        if entry_absolute == outpath_absolute {
            eprintln!("Skipping the original output file: {}", entry_absolute);
            continue;
        }

        if let Some(file_name) = entry_path.file_name() {
            // Check if the file matches the base zip file or any batched zip file (outpath.zip, outpath.1.zip, etc.)
            if let Some(captures) = zip_file_pattern.captures(file_name) {
                // remove any incomplete batches
                if file_name.ends_with(".incomplete") {
                    eprintln!("Removing incomplete zip file: '{}'", file_name);
                    if let Err(e) = tokio::fs::remove_file(&entry_path).await {
                        eprintln!(
                            "Warning: Failed to remove incomplete file '{}': {:?}",
                            file_name, e
                        );
                    }
                    continue;
                }

                // try to load a valid zip file
                Collection::from_zipfile(&entry_path)
                    .map(|collection| {
                        collections.push(collection);
                        eprintln!("loaded existing collection from {}", &entry_path);

                        // Extract batch number if it exists
                        if let Some(batch_str) = captures.get(1) {
                            if let Ok(batch_num) = batch_str.as_str().parse::<usize>() {
                                highest_batch = max(highest_batch, batch_num);
                            }
                        }
                    })
                    .unwrap_or_else(|e| {
                        eprintln!(
                            "Warning: Failed to load zip file '{}'; skipping. Zipfile Error: {:?}",
                            entry_path, e
                        );
                    });
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
) -> Result<(ZipFileWriter<Compat<File>>, PathBuf, PathBuf)> {
    let batch_outpath = if batch_size == 0 {
        // If batch size is zero, use provided outpath (contains .zip extension)
        outpath.clone()
    } else {
        // Otherwise, modify outpath to include the batch index
        let filename = outpath.file_name().unwrap();

        let new_filename = if filename.ends_with(".sig.zip") {
            let base = filename.strip_suffix(".sig.zip").unwrap();
            format!("{}.{}.sig.zip", base, batch_index)
        } else if filename.ends_with(".zip") {
            let base = filename.strip_suffix(".zip").unwrap();
            format!("{}.{}.zip", base, batch_index)
        } else {
            // This shouldn't happen due to earlier checks, but handle it just in case
            return Err(anyhow::anyhow!("Output path must end with .zip"));
        };

        outpath
            .parent()
            .unwrap_or(camino::Utf8Path::new(""))
            .join(new_filename)
    };
    let temp_outpath = batch_outpath.with_extension(format!(
        "{}.incomplete",
        batch_outpath.extension().unwrap_or_default()
    ));

    let file = File::create(&temp_outpath)
        .await
        .with_context(|| format!("Failed to create file: {:?}", temp_outpath))?;

    Ok((ZipFileWriter::with_tokio(file), temp_outpath, batch_outpath))
}

pub fn zipwriter_handle(
    mut recv_sigs: tokio::sync::mpsc::Receiver<BuildCollection>,
    output_sigs: Option<String>,
    batch_size: usize,      // Tunable batch size
    mut batch_index: usize, // starting batch index
    error_sender: tokio::sync::mpsc::Sender<anyhow::Error>,
) -> tokio::task::JoinHandle<()> {
    tokio::spawn(async move {
        let mut md5sum_occurrences = HashMap::new();
        let mut zip_manifest = BuildManifest::new();
        let mut wrote_sigs = false;
        let mut count = 0; // count the number of accessions/downloaded files sketched
        let mut zip_writer = None;
        let mut current_temp_path: Option<PathBuf> = None;
        let mut current_batch_path: Option<PathBuf> = None;

        if let Some(outpath) = output_sigs {
            let outpath: PathBuf = outpath.into();

            while let Some(mut buildcoll) = recv_sigs.recv().await {
                if zip_writer.is_none() {
                    // create zip file if needed
                    zip_writer =
                        match create_or_get_zip_file(&outpath, batch_size, batch_index).await {
                            Ok((writer, temp_outpath, final_outpath)) => {
                                current_temp_path = Some(temp_outpath);
                                current_batch_path = Some(final_outpath);
                                Some(writer)
                            }
                            Err(e) => {
                                let _ = error_sender.send(e).await;
                                return;
                            }
                        };
                }

                if let Some(zip_writer) = zip_writer.as_mut() {
                    // write all sigs from sigcoll. Note that this method updates each record's internal location
                    match buildcoll
                        .async_write_sigs_to_zip(zip_writer, &mut md5sum_occurrences)
                        .await
                    {
                        Ok(true) => {
                            wrote_sigs = true;
                            // Add all records from buildcoll manifest
                            zip_manifest.extend_from_manifest(&buildcoll.manifest);
                            count += 1
                        }
                        Ok(false) => {
                            continue; // No new signatures written, skip to next
                        }
                        Err(e) => {
                            let error = e.context("Error processing signature");
                            if error_sender.send(error).await.is_err() {
                                return;
                            }
                        }
                    }
                }

                // if batch size is non-zero and is reached, close the current zip
                if batch_size > 0 && count >= batch_size {
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
                    // rename temp to final
                    if let Some(temp_path) = current_temp_path.take() {
                        let final_path = current_batch_path.as_ref().unwrap();
                        if let Err(e) = tokio::fs::rename(&temp_path, final_path).await {
                            let _ = error_sender
                                .send(anyhow::anyhow!(
                                    "Failed to rename {:?} to {:?}: {}",
                                    temp_path,
                                    final_path,
                                    e
                                ))
                                .await;
                            return;
                        }
                    }
                    eprintln!(
                        "finished batch {}: wrote to '{}'",
                        batch_index,
                        current_batch_path.as_ref().unwrap()
                    );
                    // Start a new batch
                    batch_index += 1;
                    count = 0;
                    zip_manifest.clear();
                    zip_writer = None; // reset zip_writer so a new zip will be created when needed
                }
            }

            if count > 0 {
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
                    // rename temp to final
                    if let Some(temp_path) = current_temp_path.take() {
                        let final_path = current_batch_path.as_ref().unwrap();
                        if let Err(e) = tokio::fs::rename(&temp_path, final_path).await {
                            let _ = error_sender
                                .send(anyhow::anyhow!(
                                    "Failed to rename {:?} to {:?}: {}",
                                    temp_path,
                                    final_path,
                                    e
                                ))
                                .await;
                            return;
                        }
                    }
                    // notify about final batch
                    eprintln!(
                        "finished batch {}: wrote to '{}'",
                        batch_index,
                        current_batch_path.as_ref().unwrap()
                    );
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

// sync CSV writer
pub fn write_csv_to_path<T: ToCsvRow>(records: &[T], output_path: &Utf8Path) -> Result<()> {
    if let Some(parent) = output_path.parent() {
        create_dir_all(parent.as_std_path())
            .with_context(|| format!("Failed to create parent directory: {}", parent))?;
    }

    let mut file = std::fs::OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open(output_path.as_std_path())
        .with_context(|| format!("Failed to open file at {}", output_path))?;

    file.write_all(T::csv_header().as_bytes())?;
    for record in records {
        file.write_all(record.csv_record().as_bytes())?;
    }

    Ok(())
}

pub fn csv_writer_handle<T: ToCsvRow + Send + 'static>(
    output_path: String,
    mut receiver: tokio::sync::mpsc::Receiver<T>,
    error_sender: tokio::sync::mpsc::Sender<Error>,
) -> tokio::task::JoinHandle<()> {
    tokio::spawn(async move {
        match File::create(&output_path).await {
            Ok(file) => {
                let mut writer = BufWriter::new(file);

                if let Err(e) = writer.write_all(T::csv_header().as_bytes()).await {
                    let error = Error::new(e).context("Failed to write headers");
                    let _ = error_sender.send(error).await;
                    return;
                }

                while let Some(record) = receiver.recv().await {
                    if let Err(e) = writer.write_all(record.csv_record().as_bytes()).await {
                        let error = Error::new(e).context("Failed to write record");
                        let _ = error_sender.send(error).await;
                    }
                }

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

#[allow(clippy::too_many_arguments)]
pub async fn process_accession_stream(
    acc_data: Vec<AccessionData>,
    concurrency_limit: usize,
    download_counter: Arc<AtomicUsize>,
    existing_records_map: Arc<HashMap<String, BuildManifest>>,
    sig_templates: Arc<BuildCollection>,
    client: Arc<Client>,
    download_path: PathBuf,
    retry_times: u32,
    keep_fastas: bool,
    download_only: bool,
    filter: bool,
    send_sigs: Arc<Sender<BuildCollection>>,
    send_failed: Arc<Sender<AccessionData>>,
    send_failed_checksums: Arc<Sender<FailedChecksum>>,
    error_sender: Arc<Sender<anyhow::Error>>,
    verbose: bool,
    cancel_token: CancellationToken,
) -> Result<()> {
    // get total number of downloads (not acc)
    let n_downloads = acc_data.len();
    // set reporting threshold to every 1 percent (or every 1, whichever is larger)
    let reporting_threshold = std::cmp::max(n_downloads / 100, 1);

    stream::iter(acc_data)
        .for_each_concurrent(Some(concurrency_limit), move |accinfo| {
            let client = client.clone();
            let send_sigs = send_sigs.clone();
            let send_failed = send_failed.clone();
            let send_failed_checksums = send_failed_checksums.clone();
            let error_sender = error_sender.clone();
            let download_path = download_path.clone();
            let download_counter = download_counter.clone();
            let existing_records_map = existing_records_map.clone();
            let sig_templates = sig_templates.clone();
            let cancel_token = cancel_token.clone();

            async move {
                if cancel_token.is_cancelled() {
                    return;
                }
                if let Err(e) = pyo3::Python::with_gil(|py| py.check_signals()) {
                    cancel_token.cancel();
                    eprintln!("Caught Python interrupt: {}", e);
                    return;
                }

                let mut sigs = (*sig_templates).clone();
                // filter template sigs based on existing sigs
                if filter {
                    if let Some(existing_manifest) = existing_records_map.get(&accinfo.name) {
                        // If the key exists, filter template sigs
                        sigs.filter_by_manifest(existing_manifest);
                    }
                }

                let download_index = download_counter.fetch_add(1, Ordering::Relaxed) + 1;

                // progress report when processing begins for a download
                if verbose || download_index % reporting_threshold == 0 {
                    let percent = ((download_index as f64 / n_downloads as f64) * 100.0).round();
                    println!(
                        "Starting download {}/{} ({}%) - accession: '{}', moltype: {}",
                        download_index, n_downloads, percent, accinfo.accession, accinfo.moltype
                    );
                }

                match stream_and_process_with_retry(
                    &client,
                    accinfo,
                    &download_path,
                    Some(retry_times),
                    keep_fastas,
                    sigs,
                    download_only,
                    true,
                    cancel_token.clone(),
                )
                .await
                {
                    Ok((coll, failed_downloads, failed_checksums)) => {
                        if !coll.is_empty() {
                            if let Err(e) = send_sigs.send(coll).await {
                                let _ = error_sender.send(e.into()).await;
                            }
                        }
                        for fail in failed_downloads {
                            let _ = send_failed.send(fail).await;
                        }
                        for fail in failed_checksums {
                            let _ = send_failed_checksums.send(fail).await;
                        }
                    }
                    Err(e) => {
                        let _ = error_sender.send(e).await;
                    }
                }
            }
        })
        .await;

    Ok(())
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
    concurrency_limit: usize,
    api_key: String,
    verbose: bool,
    write_urlsketch_csv: bool,
    output_sigs: Option<String>,
) -> Result<(), anyhow::Error> {
    let batch_size = batch_size as usize;
    let mut batch_index = 1;
    let mut existing_records_map: HashMap<String, BuildManifest> = HashMap::new();
    let mut filter = false;
    // if writing sigs, prepare output and look for existing sig batches
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
            existing_records_map = existing_sigs.build_recordsmap();

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
    let (send_sigs, recv_sigs) = tokio::sync::mpsc::channel::<BuildCollection>(4);
    let (send_failed, recv_failed) = tokio::sync::mpsc::channel::<AccessionData>(4);
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

    let failures_handle = csv_writer_handle(failed_csv, recv_failed, error_sender.clone());
    let checksum_failures_handle = csv_writer_handle(
        failed_checksums_csv,
        recv_failed_checksum,
        error_sender.clone(),
    );
    // error tracker for no sigs written at all
    let critical_error_flag = Arc::new(AtomicBool::new(false));

    // tracker for keyboard interrupt
    let cancel_token = CancellationToken::new();

    let error_handle = error_handler(error_receiver, critical_error_flag.clone());
    handles.push(sig_handle);
    handles.push(failures_handle);
    handles.push(error_handle);
    handles.push(checksum_failures_handle);

    // Worker tasks
    let client = Arc::new(Client::new());

    let mut sig_templates = BuildCollection::new();
    let mut genomes_only = genomes_only;
    let mut proteomes_only = proteomes_only;

    // build sketch params
    if download_only {
        if genomes_only {
            eprintln!("Downloading genomes only.");
        } else if proteomes_only {
            eprintln!("Downloading proteomes only.");
        }
    } else {
        let sig_template_result = BuildCollection::from_param_str(param_str.as_str());
        sig_templates = match sig_template_result {
            Ok(sig_templates) => sig_templates,
            Err(e) => {
                bail!("Failed to parse params string: {}", e);
            }
        };
        // Check if we have dna signature templates and not keep_fastas
        if sig_templates.anydna_size()? == 0 && !keep_fastas {
            eprintln!("No DNA signature templates provided, and --keep-fasta is not set.");
            proteomes_only = true;
        }
        // Check if we have protein signature templates not keep_fastas
        if sig_templates.anyprotein_size()? == 0 && !keep_fastas {
            eprintln!("No protein signature templates provided, and --keep-fasta is not set.");
            genomes_only = true;
        }
        if genomes_only {
            // select only templates built from DNA input
            let multiselection = MultiSelection::from_input_moltype("DNA")?;
            sig_templates.select(&multiselection)?;
            eprintln!("Downloading and sketching genomes only.");
        } else if proteomes_only {
            // select only templates built from protein input
            let multiselection = MultiSelection::from_input_moltype("protein")?;
            sig_templates.select(&multiselection)?;
            eprintln!("Downloading and sketching proteomes only.");
        }
        if sig_templates.is_empty() && !download_only {
            bail!("No signatures to build.")
        }
    }

    // gbsketch specific : set up the HTTP client and headers
    let mut headers = HeaderMap::new();
    headers.insert("Content-Type", "application/json".parse()?);
    let post_url = Url::parse("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/download")?;

    let mut annotation_types = vec!["SEQUENCE_REPORT"];

    let mut include_genome_files = true;
    let mut include_protein_files = true;
    if genomes_only {
        annotation_types.push("GENOME_FASTA");
        include_protein_files = false; // only genomic files
    } else if proteomes_only {
        annotation_types.push("PROT_FASTA");
        include_genome_files = false; // only protein files
    } else {
        annotation_types.push("GENOME_FASTA");
        annotation_types.push("PROT_FASTA");
    }

    // get accessions for retrieval
    // Open the file containing the accessions synchronously
    let (mut accession_info, n_accs) = load_gbsketch_info(
        input_csv.clone(),
        include_genome_files,
        include_protein_files,
    )?;
    if n_accs == 0 {
        bail!("No accessions to download and sketch.")
    }
    let accession_set: HashSet<String> =
        accession_info.iter().map(|a| a.accession.clone()).collect();

    let accessions: Vec<String> = accession_set.into_iter().collect();

    let post_params = json!({
        "accessions": accessions,
        "include_annotation_type": annotation_types,
        "hydrated": "DATA_REPORT_ONLY",
        "api_key": if !api_key.is_empty() { Some(api_key) } else { None }
    });

    // dl + parse dehydrated file
    py.check_signals()?; // If interrupted, return an Err automatically
    let download_links = match download_parse_dehydrated_zip_with_retry(
        &client,
        post_url,
        headers,
        post_params,
        retry_times,
    )
    .await
    {
        Ok(links) => links,
        Err(e) => {
            bail!(
                "Failed to retrieve dehydrated download ZIP. Are your accessions valid? Exiting.\n Reason:\n{}",
                e
            );
        }
    };
    eprintln!("Successfully downloaded and parsed dehydrated zipfile. Now processing accessions.");

    let mut download_failures = Vec::<AccessionData>::new();

    // Now attach URLs to each accession info
    for accinfo in &mut accession_info {
        py.check_signals()?; // If interrupted, return an Err
        let dl_filename = accinfo
            .download_filename
            .as_ref()
            .expect("download_filename should always be set");

        match download_links.get(dl_filename) {
            Some(url_str) => match Url::parse(url_str) {
                Ok(url) => {
                    accinfo.attach_url(url);
                }
                Err(e) => {
                    if verbose {
                        eprintln!(
                            "Failed to parse URL for file {}: {}. Writing to download_failures.",
                            dl_filename, e
                        );
                    }
                    download_failures.push(accinfo.clone());
                }
            },
            None => {
                if verbose {
                    eprintln!(
                        "No download URL found for file: {}. Writing to download_failures.",
                        dl_filename
                    );
                }
                download_failures.push(accinfo.clone());
            }
        }
        // if we failed to get a download link, write to download failures
        for fail in &download_failures {
            if let Err(e) = send_failed.send(fail.clone()).await {
                eprintln!("Failed to send missing download failure info: {}", e);
                let _ = error_sender.send(e.into()).await;
            }
        }
    }

    // only retain AccessionData entries that have valid URLs
    accession_info.retain(|a| !a.url_info.is_empty());

    // optionally write urlsketch csv with all available download links
    if write_urlsketch_csv {
        let urlsketch_path = PathBuf::from(format!("{}.urlsketch.csv", input_csv));
        if let Err(e) = write_csv_to_path(&accession_info, &urlsketch_path) {
            eprintln!("Failed to write urlsketch csv with download links: {:#}", e);
        } else {
            eprintln!(
                "Wrote urlsketch csv with download links to '{}'",
                urlsketch_path
            );
        }
    }

    // these now need to be Arc to be shared across async tasks
    let download_counter = Arc::new(AtomicUsize::new(0));
    let existing_records_map = Arc::new(existing_records_map);
    let sig_templates = Arc::new(sig_templates);

    let send_sigs_arc = Arc::new(send_sigs);
    let send_failed_arc = Arc::new(send_failed);
    let send_failed_checksums_arc = Arc::new(send_failed_checksums);
    let error_sender_arc = Arc::new(error_sender);

    process_accession_stream(
        accession_info,
        concurrency_limit,
        download_counter,
        existing_records_map,
        sig_templates,
        client,
        download_path,
        retry_times,
        keep_fastas,
        download_only,
        filter,
        send_sigs_arc,
        send_failed_arc,
        send_failed_checksums_arc,
        error_sender_arc,
        verbose,
        cancel_token.clone(), // pass the cancellation token to allow for graceful shutdown),
    )
    .await?;

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
    genomes_only: bool,
    proteomes_only: bool,
    download_only: bool,
    batch_size: u32,
    n_permits: usize,
    force: bool,
    verbose: bool,
    output_sigs: Option<String>,
    failed_checksums_csv: Option<String>,
) -> Result<(), anyhow::Error> {
    let batch_size = batch_size as usize;
    let mut batch_index = 1;
    let mut existing_recordsmap: HashMap<String, BuildManifest> = HashMap::new();
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
            existing_recordsmap = existing_sigs.build_recordsmap();

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
    let (send_sigs, recv_sigs) = tokio::sync::mpsc::channel::<BuildCollection>(4);
    let (send_failed, recv_failed) = tokio::sync::mpsc::channel::<AccessionData>(4);
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

    let failures_handle = csv_writer_handle(failed_csv, recv_failed, error_sender.clone());

    let mut write_failed_checksums = false;
    if let Some(failed_checksums) = failed_checksums_csv {
        let checksum_failures_handle =
            csv_writer_handle(failed_checksums, recv_failed_checksum, error_sender.clone());

        write_failed_checksums = true;
        handles.push(checksum_failures_handle);
    }

    let critical_error_flag = Arc::new(AtomicBool::new(false));
    let error_handle = error_handler(error_receiver, critical_error_flag.clone());
    handles.push(sig_handle);
    handles.push(failures_handle);
    handles.push(error_handle);

    // Worker tasks
    let semaphore = Arc::new(Semaphore::new(n_permits)); // Limiting concurrent downloads
    let client = Arc::new(Client::new());

    // Open the file containing the accessions synchronously
    let (accession_info, n_accs) = load_accession_info(input_csv, keep_fastas, force)?;
    if n_accs == 0 {
        bail!("No accessions to download and sketch.")
    }

    let mut sig_templates = BuildCollection::new();
    let mut genomes_only = genomes_only;
    let mut proteomes_only = proteomes_only;
    let dna_multiselection = MultiSelection::from_input_moltype("dna")?;
    let protein_multiselection = MultiSelection::from_input_moltype("protein")?;

    if download_only {
        if genomes_only {
            eprintln!("Downloading genomes only.");
        } else if proteomes_only {
            eprintln!("Downloading proteomes only.");
        }
    } else {
        let sig_template_result = BuildCollection::from_param_str(param_str.as_str());
        sig_templates = match sig_template_result {
            Ok(sig_templates) => sig_templates,
            Err(e) => {
                bail!("Failed to parse params string: {}", e);
            }
        };
        // Check if we have dna signature templates and not keep_fastas
        if sig_templates.anydna_size()? == 0 && !keep_fastas {
            eprintln!("No DNA signature templates provided, and --keep-fasta is not set.");
            proteomes_only = true;
        }
        // Check if we have protein signature templates not keep_fastas
        if sig_templates.anyprotein_size()? == 0 && !keep_fastas {
            eprintln!("No protein signature templates provided, and --keep-fasta is not set.");
            genomes_only = true;
        }
        if genomes_only {
            // select only DNA templates
            sig_templates.select(&dna_multiselection)?;
            eprintln!("Downloading and sketching genomes only.");
        } else if proteomes_only {
            // select only protein templates
            sig_templates.select(&protein_multiselection)?;
            eprintln!("Downloading and sketching proteomes only.");
        }
        if sig_templates.is_empty() && !download_only {
            bail!("No signatures to build.")
        }
    }

    // report every 1 percent (or every 1, whichever is larger)
    let reporting_threshold = std::cmp::max(n_accs / 100, 1);

    for (i, accinfo) in accession_info.into_iter().enumerate() {
        py.check_signals()?; // If interrupted, return an Err automatically
        let mut sigs = sig_templates.clone();

        // filter template sigs based on existing sigs
        if filter {
            if let Some(existing_manifest) = existing_recordsmap.get(&accinfo.name) {
                // If the key exists, filter template sigs
                sigs.filter_by_manifest(existing_manifest);
            }
        }
        // eliminate sigs that won't be added to based on moltype
        // this assumes no translation --> modify as needed if adding that.
        if accinfo.moltype == InputMolType::Dna {
            sigs.select(&dna_multiselection)?;
        } else {
            sigs.select(&protein_multiselection)?;
        }
        if sigs.is_empty() && !download_only {
            continue;
        }

        // eliminate sigs that won't be added to based on moltype
        // this assumes no translation --> modify as needed if adding that.
        if accinfo.moltype == InputMolType::Dna {
            sigs.select(&dna_multiselection)?;
        } else {
            sigs.select(&protein_multiselection)?;
        }
        if sigs.is_empty() && !download_only {
            continue;
        }

        let semaphore_clone = Arc::clone(&semaphore);
        let client_clone = Arc::clone(&client);
        let send_sigs = send_sigs.clone();
        let send_failed = send_failed.clone();
        let checksum_send_failed = send_failed_checksums.clone();
        let download_path_clone = download_path.clone(); // Clone the path for each task
        let send_errors = error_sender.clone();
        let write_checksum_fail = write_failed_checksums;

        tokio::spawn(async move {
            let _permit = semaphore_clone.acquire().await;
            // progress report when the permit is available and processing begins
            if verbose || (i + 1) % reporting_threshold == 0 {
                let percent_processed = (((i + 1) as f64 / n_accs as f64) * 100.0).round();
                println!(
                    "Starting accession {}/{} ({}%) - moltype: {}",
                    (i + 1),
                    n_accs,
                    percent_processed,
                    accinfo.moltype
                );
            }
            // Perform download and sketch
            let result = dl_sketch_url(
                &client_clone,
                accinfo.clone(),
                &download_path_clone,
                Some(retry_times),
                keep_fastas,
                sigs,
                download_only,
                write_checksum_fail,
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
                    // if write_failed_checksums {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::buildutils::BuildRecord;
    use camino::Utf8PathBuf;

    #[test]
    fn test_buildrecordsmap() {
        // read in zipfiles to build a MultiCollection
        let mut filename = Utf8PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("tests/test-data/GCA_000961135.2.sig.zip");
        let path = filename.clone();

        let mut collections = Vec::new();
        let coll: Collection = Collection::from_zipfile(&path).unwrap();
        collections.push(coll);
        let mc = MultiCollection::new(collections);

        // build expected buildmanifest
        let mut refbmf = BuildManifest::new();
        let mut rec1 = BuildRecord::default_dna();
        rec1.set_with_abundance(true);
        refbmf.add_record(rec1);

        //  Call build_recordsmap
        let name_params_map = mc.build_recordsmap();

        // Check that the recordsmap contains the correct names
        assert_eq!(
            name_params_map.len(),
            1,
            "There should be 1 unique names in the map"
        );

        for (name, buildmanifest) in name_params_map.iter() {
            eprintln!("Name: {}", name);
            assert_eq!(
                "GCA_000961135.2 Candidatus Aramenus sulfurataquae isolate AZ1-454",
                name
            );
            assert_eq!(buildmanifest.size(), 2); // should be two records
                                                 // check that we can filter out a record (k=31, abund)
            let filtered = buildmanifest.filter_manifest(&refbmf);
            assert_eq!(filtered.size(), 1)
        }
    }
}
