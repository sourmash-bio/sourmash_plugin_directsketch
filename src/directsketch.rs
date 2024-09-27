use anyhow::{anyhow, bail, Context, Error, Result};
use async_zip::base::write::ZipFileWriter;
use async_zip::{Compression, ZipDateTime, ZipEntryBuilder};
use camino::Utf8PathBuf as PathBuf;
use chrono::Utc;
use futures::future::join;
use needletail::parse_fastx_reader;
use regex::Regex;
use reqwest::Client;
use std::collections::HashMap;
use std::fs::{self, create_dir_all};
use std::io::Cursor;
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use tokio::fs::File;
use tokio::io::{AsyncWriteExt, BufWriter};
use tokio::sync::Semaphore;
use tokio_util::compat::Compat;

use pyo3::prelude::*;

use sourmash::collection::Collection;
use sourmash::manifest::{Manifest, Record};
use sourmash::signature::Signature;

use crate::utils::{
    build_siginfo, build_template_collection, load_accession_info, load_gbassembly_info,
    parse_params_str, AccessionData, GBAssemblyData, GenBankFileType, InputMolType,
    TemplateCollection,
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
                "Invalid checksum line format in URL {}: {}",
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
                            "MD5 hash does not match. Expected: {}, Found: {}",
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

async fn sketch_data(
    name: String,
    filename: String,
    compressed_data: Vec<u8>,
    mut coll: TemplateCollection,
    moltype: String,
) -> Result<TemplateCollection> {
    tokio::task::spawn_blocking(move || {
        let cursor = Cursor::new(compressed_data);
        let mut fastx_reader =
            parse_fastx_reader(cursor).context("Failed to parse FASTA/FASTQ data")?;

        let mut set_name = false;
        while let Some(record) = fastx_reader.next() {
            let record = record.context("Failed to read record")?;
            coll.iter_mut().for_each(|(_tr, sig)| {
                if !set_name {
                    sig.set_name(&name);
                    sig.set_filename(&filename);
                };
                if moltype == "protein" {
                    sig.add_protein(&record.seq())
                        .expect("Failed to add protein");
                } else {
                    sig.add_sequence(&record.seq(), true)
                        .expect("Failed to add sequence");
                    // if not force, panics with 'N' in dna sequence
                }
            });
            if !set_name {
                set_name = true;
            }
        }
        Ok(coll)
    })
    .await?
}
pub struct FailedDownload {
    accession: String,
    name: String,
    moltype: String,
    md5sum: Option<String>,
    download_filename: Option<String>,
    url: Option<Url>,
}

#[allow(clippy::too_many_arguments)]
async fn dl_sketch_assembly_accession(
    client: &Client,
    accinfo: GBAssemblyData,
    location: &PathBuf,
    retry: Option<u32>,
    keep_fastas: bool,
    dna_sigs: TemplateCollection,
    prot_sigs: TemplateCollection,
    genomes_only: bool,
    proteomes_only: bool,
    download_only: bool,
) -> Result<(TemplateCollection, Vec<FailedDownload>)> {
    let retry_count = retry.unwrap_or(3); // Default retry count
                                          // let mut sigs = Vec::<Signature>::new();
    let mut sig_collection = TemplateCollection::new();
    let mut failed = Vec::<FailedDownload>::new();

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
                    failed.push(failed_download_dna);
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
                    failed.push(failed_download_protein);
                }

                return Ok((sig_collection, failed));
            }
        };
    let md5sum_url = GenBankFileType::Checksum.url(&base_url, &full_name);

    let checksums = match download_and_parse_md5(client, &md5sum_url).await {
        Ok(cs) => cs,
        Err(e) => {
            return Err(e);
        }
    };

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

    for file_type in &file_types {
        let url = file_type.url(&base_url, &full_name);
        let expected_md5 = checksums.get(&file_type.server_filename(&full_name));
        let file_name = file_type.filename_to_write(&accession);
        let data =
            match download_with_retry(client, &url, expected_md5.map(|x| x.as_str()), retry_count)
                .await
            {
                Ok(data) => data,
                Err(_err) => {
                    // here --> keep track of accession errors + filetype
                    let failed_download = FailedDownload {
                        accession: accession.clone(),
                        name: name.clone(),
                        moltype: file_type.moltype(),
                        md5sum: expected_md5.map(|x| x.to_string()),
                        download_filename: Some(file_name),
                        url: Some(url),
                    };
                    failed.push(failed_download);
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
                GenBankFileType::Genomic => sig_collection.extend(
                    sketch_data(
                        name.clone(),
                        file_name.clone(),
                        data,
                        dna_sigs.clone(),
                        "dna".to_string(),
                    )
                    .await?,
                ),
                GenBankFileType::Protein => {
                    sig_collection.extend(
                        sketch_data(
                            name.clone(),
                            file_name.clone(),
                            data,
                            prot_sigs.clone(),
                            "protein".to_string(),
                        )
                        .await?,
                    );
                }
                _ => {} // Do nothing for other file types
            };
        }
    }

    Ok((sig_collection, failed))
}

#[allow(clippy::too_many_arguments)]
async fn dl_sketch_url(
    client: &Client,
    accinfo: AccessionData,
    location: &PathBuf,
    retry: Option<u32>,
    _keep_fastas: bool,
    dna_sigs: Vec<Signature>,
    prot_sigs: Vec<Signature>,
    _genomes_only: bool,
    _proteomes_only: bool,
    download_only: bool,
) -> Result<(Vec<Signature>, Vec<FailedDownload>)> {
    let retry_count = retry.unwrap_or(3); // Default retry count
    let mut sigs = Vec::<Signature>::new();
    let mut failed = Vec::<FailedDownload>::new();

    let name = accinfo.name;
    let accession = accinfo.accession;
    let url = accinfo.url;
    let expected_md5 = accinfo.expected_md5sum;
    let download_filename = accinfo.download_filename;
    let moltype = accinfo.moltype;

    match download_with_retry(client, &url, expected_md5.as_deref(), retry_count)
        .await
        .ok()
    {
        Some(data) => {
            // check keep_fastas instead??
            if let Some(ref download_filename) = download_filename {
                let path = location.join(download_filename);
                fs::write(path, &data).context("Failed to write data to file")?;
            }
            if !download_only {
                let filename = download_filename.clone().unwrap_or("".to_string());
                // sketch data
                match moltype {
                    InputMolType::Dna => sigs.extend(
                        sketch_data(
                            name.clone(),
                            filename.clone(),
                            data,
                            dna_sigs.clone(),
                            "dna".to_string(),
                        )
                        .await?,
                    ),
                    InputMolType::Protein => {
                        sigs.extend(
                            sketch_data(
                                name.clone(),
                                filename.clone(),
                                data,
                                prot_sigs.clone(),
                                "protein".to_string(),
                            )
                            .await?,
                        );
                    }
                };
            }
        }
        None => {
            let failed_download = FailedDownload {
                accession: accession.clone(),
                name: name.clone(),
                moltype: moltype.to_string(),
                md5sum: expected_md5.map(|x| x.to_string()),
                download_filename,
                url: Some(url),
            };
            failed.push(failed_download);
        }
    }

    Ok((sigs, failed))
}

async fn write_sig_to_zip(
    sig: &Signature,
    md5sum_occurrences: &mut HashMap<String, usize>,
    manifest_rows: &mut Vec<Record>,
    zip_writer: &mut ZipFileWriter<Compat<tokio::fs::File>>,
) -> Result<()> {
    let md5sum_str = sig.md5sum();
    let count = md5sum_occurrences.entry(md5sum_str.clone()).or_insert(0);
    *count += 1;

    let sig_filename = if *count > 1 {
        format!("signatures/{}_{}.sig.gz", md5sum_str, count)
    } else {
        format!("signatures/{}.sig.gz", md5sum_str)
    };

    let records: Vec<Record> = Record::from_sig(sig, &sig_filename);
    manifest_rows.extend(records);

    let wrapped_sig = vec![sig.clone()];
    let json_bytes = serde_json::to_vec(&wrapped_sig)
        .map_err(|e| anyhow!("Error serializing signature: {}", e))?;

    let gzipped_buffer = {
        let mut buffer = std::io::Cursor::new(Vec::new());
        {
            let mut gz_writer = niffler::get_writer(
                Box::new(&mut buffer),
                niffler::compression::Format::Gzip,
                niffler::compression::Level::Nine,
            )?;
            //     .map_err(|e| anyhow!("Error creating gzip writer: {}", e))?;
            gz_writer.write_all(&json_bytes)?;
            //         .map_err(|e| anyhow!("Error writing gzip data: {}", e))?;
        }

        buffer.into_inner()
    };

    let now = Utc::now();
    let builder = ZipEntryBuilder::new(sig_filename.into(), Compression::Stored)
        .last_modification_date(ZipDateTime::from_chrono(&now))
        .unix_permissions(0o644);
    zip_writer
        .write_entry_whole(builder, &gzipped_buffer)
        .await
        .map_err(|e| anyhow!("Error writing zip entry: {}", e))
}

async fn write_sig_to_tmpdir(
    sig: &Signature,
    md5sum_occurrences: &mut HashMap<String, usize>,
    manifest_rows: &mut Vec<Record>,
    tmpdir: &PathBuf,
) -> Result<()> {
    let md5sum_str = sig.md5sum();
    let count = md5sum_occurrences.entry(md5sum_str.clone()).or_insert(0);
    *count += 1;

    let sig_filename = if *count > 1 {
        format!(
            "{}/signatures/{}_{}.sig.gz",
            tmpdir.as_str(),
            md5sum_str,
            count
        )
    } else {
        format!("{}/signatures/{}.sig.gz", tmpdir.as_str(), md5sum_str)
    };

    let records: Vec<Record> = Record::from_sig(sig, &sig_filename);
    manifest_rows.extend(records);

    let wrapped_sig = vec![sig.clone()];
    let json_bytes = serde_json::to_vec(&wrapped_sig)
        .map_err(|e| anyhow!("Error serializing signature: {}", e))?;

    let gzipped_buffer = {
        let mut buffer = std::io::Cursor::new(Vec::new());
        {
            let mut gz_writer = niffler::get_writer(
                Box::new(&mut buffer),
                niffler::compression::Format::Gzip,
                niffler::compression::Level::Nine,
            )?;
            gz_writer.write_all(&json_bytes)?;
        }
        buffer.into_inner()
    };

    // write sig
    let sig_file_path = Path::new(&sig_filename);
    let mut sig_file = File::create(sig_file_path).await?;
    sig_file.write_all(&gzipped_buffer).await?;

    Ok(())
}

async fn write_manifest_to_zip(
    zip_writer: &mut ZipFileWriter<Compat<tokio::fs::File>>,
    manifest_rows: &Vec<Record>,
) -> Result<()> {
    let manifest_filename = "SOURMASH-MANIFEST.csv".to_string();

    // convert manifest_rows to a manifest and serialize it to a buffer
    let manifest: Manifest = manifest_rows.clone().into();
    let mut manifest_buffer = Vec::new();

    manifest
        .to_writer(&mut manifest_buffer)
        .context("Failed to serialize manifest")?;

    // create the ZipEntryBuilder with the current time and permissions
    let now = Utc::now();
    let builder = ZipEntryBuilder::new(manifest_filename.into(), Compression::Stored)
        .last_modification_date(ZipDateTime::from_chrono(&now))
        .unix_permissions(0o644);

    // srite the manifest entry to the zip, handle potential errors
    zip_writer
        .write_entry_whole(builder, &manifest_buffer)
        .await
        .context("Failed to write manifest to ZIP")?;

    Ok(())
}

pub fn zipwriter_handle(
    mut recv_sigs: tokio::sync::mpsc::Receiver<TemplateCollection>,
    output_sigs: Option<String>,
    batch_size: usize, // Tunable batch size
    error_sender: tokio::sync::mpsc::Sender<anyhow::Error>,
) -> tokio::task::JoinHandle<()> {
    tokio::spawn(async move {
        let mut md5sum_occurrences = HashMap::new();
        let mut manifest_rows = Vec::new();
        let mut wrote_sigs = false;
        let mut file_count = 0; // Count of files in the current batch
        let mut batch_index = 1; // Index to name zip files

        if let Some(outpath) = output_sigs {
            let outpath_base: PathBuf = outpath.into();

            // Function to create a new zip file for each batch
            let create_zip_file = |batch_index: usize| {
                let outpath_base_ref = &outpath_base; // Borrow outpath_base
                let err_send = error_sender.clone();
                async move {
                    let batch_outpath = outpath_base_ref.with_file_name(format!(
                        "{}_batch_{}.zip",
                        outpath_base_ref.file_stem().unwrap(),
                        batch_index
                    ));

                    let file = match File::create(&batch_outpath).await {
                        Ok(file) => file,
                        Err(e) => {
                            let error = anyhow::Error::new(e).context("Failed to create file");
                            let _ = err_send.send(error).await;
                            return None;
                        }
                    };
                    Some(ZipFileWriter::with_tokio(file))
                }
            };

            let mut zip_writer = match create_zip_file(batch_index).await {
                Some(writer) => writer,
                None => return,
            };

            while let Some(sigcoll) = recv_sigs.recv().await {
                for (_r, sig) in sigcoll {
                    match write_sig_to_zip(
                        &sig,
                        &mut md5sum_occurrences,
                        &mut manifest_rows,
                        &mut zip_writer,
                    )
                    .await
                    {
                        Ok(_) => {
                            wrote_sigs = true;
                            file_count += 1; // Increment file count
                        }
                        Err(e) => {
                            let error = e.context("Error processing signature");
                            if error_sender.send(error).await.is_err() {
                                return;
                            }
                        }
                    }

                    // If batch size is reached, close the current ZIP and start a new one
                    if file_count >= batch_size {
                        if let Err(e) = write_manifest_to_zip(&mut zip_writer, &manifest_rows).await
                        {
                            let _ = error_sender.send(e).await;
                        }

                        if let Err(e) = zip_writer.close().await {
                            let error = anyhow::Error::new(e).context("Failed to close ZIP file");
                            let _ = error_sender.send(error).await;
                            return;
                        }

                        // Start a new batch
                        batch_index += 1;
                        file_count = 0;
                        manifest_rows.clear(); // Clear manifest for the new batch
                        zip_writer = match create_zip_file(batch_index).await {
                            Some(writer) => writer,
                            None => return,
                        };
                    }
                }
            }

            // Finalize the last batch if there are remaining signatures
            if wrote_sigs {
                if let Err(e) = write_manifest_to_zip(&mut zip_writer, &manifest_rows).await {
                    let _ = error_sender.send(e).await;
                }

                if let Err(e) = zip_writer.close().await {
                    let error = anyhow::Error::new(e).context("Failed to close ZIP file");
                    let _ = error_sender.send(error).await;
                    return;
                }
            } else {
                let error = anyhow::Error::new(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    "No signatures written",
                ));
                let _ = error_sender.send(error).await;
            }
        }
    })
}

// pub fn zipwriter_handle(
//     mut recv_sigs: tokio::sync::mpsc::Receiver<Vec<Signature>>,
//     output_sigs: Option<String>,
//     error_sender: tokio::sync::mpsc::Sender<anyhow::Error>,
// ) -> tokio::task::JoinHandle<()> {
//     tokio::spawn(async move {
//         let mut md5sum_occurrences = HashMap::new();
//         let mut manifest_rows = Vec::new();
//         let mut wrote_sigs = false;
//         if let Some(outpath) = output_sigs {
//             let outpath: PathBuf = outpath.into();

//             let file = match File::create(&outpath).await {
//                 Ok(file) => file,
//                 Err(e) => {
//                     let error =
//                         anyhow::Error::new(e).context("Failed to create file at specified path");
//                     let _ = error_sender.send(error).await; // Send the error through the channel
//                     return; // Simply exit the task as error handling is managed elsewhere
//                 }
//             };
//             let mut zip_writer = ZipFileWriter::with_tokio(file);

//             while let Some(sigs) = recv_sigs.recv().await {
//                 for sig in sigs {
//                     match write_sig_to_zip(
//                         &sig,
//                         &mut md5sum_occurrences,
//                         &mut manifest_rows,
//                         &mut zip_writer,
//                     )
//                     .await
//                     {
//                         Ok(_) => wrote_sigs = true,
//                         Err(e) => {
//                             let error = e.context("Error processing signature");
//                             if (error_sender.send(error).await).is_err() {
//                                 return; // Exit on failure to send error
//                             }
//                         }
//                     }
//                 }
//             }

//             if wrote_sigs {
//                 if let Err(e) = write_manifest_to_zip(&mut zip_writer, &manifest_rows).await {
//                     let _ = error_sender.send(e).await;
//                     // should i return here if err?
//                 }

//                 if let Err(e) = zip_writer.close().await {
//                     let error = anyhow::Error::new(e).context("Failed to close ZIP file");
//                     let _ = error_sender.send(error).await;
//                     return;
//                 }
//             } else {
//                 let error = anyhow::Error::new(std::io::Error::new(
//                     std::io::ErrorKind::Other,
//                     "No signatures written",
//                 ));
//                 let _ = error_sender.send(error).await; // Send error about no signatures written
//             }
//         }
//         drop(error_sender); // should be dropped automatically as it goes out of scope, but just in case
//     })
// }

/// instead of this, can just read from zipfiles (manifests make this easier!)
pub async fn signature_from_path(path: &PathBuf) -> Result<Vec<Signature>, Error> {
    let path = path.clone(); // Clone the path to move into the blocking thread

    // Use `spawn_blocking` to call `Signature::from_path` (synchronous)
    let sigs = tokio::task::spawn_blocking(move || {
        Signature::from_path(&path).map_err(|e| anyhow!("Error reading signatures: {}", e))
    })
    .await??;

    Ok(sigs)
}

async fn copy_tmpdir_to_zip(tmpdir: &PathBuf, zipfile_path: String) -> Result<()> {
    let zippath: PathBuf = zipfile_path.into();

    // Create the zipfile + propagate errs
    let file = File::create(&zippath).await.with_context(|| {
        format!(
            "Failed to create file at specified path: {}",
            zippath.as_str()
        )
    })?;

    // use block for zip_writer to ensure it gets dropped/finalized before cleaning up tmpdir
    {
        let mut zip_writer = ZipFileWriter::with_tokio(file);

        // find all signature paths
        let sigdir = tmpdir.join("signatures");
        let mut read_dir = tokio::fs::read_dir(&sigdir)
            .await
            .with_context(|| format!("Failed to read signature directory: {}", sigdir.as_str()))?;

        let mut dir_manifest_rows: Vec<Record> = Vec::new();
        let mut md5sum_occurrences: HashMap<String, usize> = HashMap::new();

        // Read + write sigs: iterate over directory entries asynchronously
        while let Some(entry) = read_dir.next_entry().await? {
            let path = entry.path();
            let utf8_path = PathBuf::from_path_buf(path.clone())
                .map_err(|_| anyhow::anyhow!("Invalid UTF-8 path: {:?}", path))?;

            // load sigs from each path. Should just be a single sig, it's just returned in a vector.
            let sigs = signature_from_path(&utf8_path)
                .await
                .with_context(|| format!("Failed to load signature from path: {:?}", utf8_path))?;

            // We could just take out the first (and only) sig here. But iteration is maybe more flexible?
            // Allows you to include sigs made outside of directsketch (that may have mult sigs per file)...
            for sig in &sigs {
                write_sig_to_zip(
                    sig,
                    &mut md5sum_occurrences,
                    &mut dir_manifest_rows,
                    &mut zip_writer,
                )
                .await
                .with_context(|| {
                    format!("Failed to write signature to zip for path: {:?}", utf8_path)
                })?;
            }
        }

        // Write the final manifest to the zip
        write_manifest_to_zip(&mut zip_writer, &dir_manifest_rows)
            .await
            .context("Failed to handle manifest writing process")?;
    } // zip_writer is dropped/finalized upon exiting this block

    // Clean up the tmpdir after zipping
    tokio::fs::remove_dir_all(tmpdir)
        .await
        .with_context(|| format!("Failed to remove tmpdir: {}", tmpdir.as_str()))?;

    Ok(())
}

pub async fn read_from_tmpdir(tmpdir: &PathBuf) -> Result<Manifest> {
    let sigdir = tmpdir.join("signatures");

    // Read directory entries asynchronously
    let mut dir_entries = tokio::fs::read_dir(&sigdir)
        .await
        .with_context(|| format!("Failed to read directory entries from {}", sigdir))?;

    // Collect all paths from the directory and convert them to Utf8PathBuf
    let mut paths = Vec::new();
    while let Some(entry) = dir_entries.next_entry().await? {
        let path = entry.path();

        // Try to convert PathBuf to Utf8PathBuf
        let utf8_path = PathBuf::from_path_buf(path)
            .map_err(|_| anyhow::anyhow!("Invalid UTF-8 path in directory entry"))?;
        paths.push(utf8_path);
    }

    // Load collection from paths (synchronous, so use blocking task)
    let collection = tokio::task::spawn_blocking(move || {
        Collection::from_paths(&paths).map_err(|e| {
            anyhow::Error::new(e).context("Failed to load paths as sourmash collection")
        })
    })
    .await?
    .map_err(|e| anyhow::anyhow!("Failed to execute Collection::from_paths: {}", e))?;

    // Return the collection's manifest; could avoid clone by managing the lifetime differently.
    Ok(collection.manifest().clone())
}

pub fn tmpdir_sigwriter_handle(
    mut recv_sigs: tokio::sync::mpsc::Receiver<Vec<Signature>>,
    tmpdir: String,
    output_zip: Option<String>,
    error_sender: tokio::sync::mpsc::Sender<anyhow::Error>,
) -> tokio::task::JoinHandle<()> {
    tokio::spawn(async move {
        let mut md5sum_occurrences = HashMap::new();
        let mut manifest_rows = Vec::new();
        let mut wrote_sigs = false;
        let tmppath: PathBuf = tmpdir.into();

        let outpath_dir = Path::new(&tmppath).parent().unwrap();
        // will return Ok if dir already exists
        if let Err(e) = tokio::fs::create_dir_all(outpath_dir).await {
            let error =
                anyhow::Error::new(e).context("Failed to create directory at specified path");
            let _ = error_sender.send(error).await; // send the error through the channel
            return;
        }

        while let Some(sigs) = recv_sigs.recv().await {
            for sig in sigs {
                match write_sig_to_tmpdir(
                    &sig,
                    &mut md5sum_occurrences,
                    &mut manifest_rows,
                    &tmppath,
                )
                .await
                {
                    Ok(_) => wrote_sigs = true,
                    Err(e) => {
                        let error = e.context("Error processing signature");
                        if (error_sender.send(error).await).is_err() {
                            return; // Exit on failure to send error
                        }
                    }
                }
            }
        }
        if wrote_sigs {
            // if output zipfile specified, copy to zip and write manifest
            if let Some(zippath) = output_zip {
                // copy sigs from tmpdir into zipfile
                match copy_tmpdir_to_zip(&tmppath, zippath).await {
                    Ok(()) => {}
                    Err(e) => {
                        let _ = error_sender.send(e).await;
                    }
                }
            }
        } else {
            let error = anyhow::Error::new(std::io::Error::new(
                std::io::ErrorKind::Other,
                "No signatures written",
            ));
            let _ = error_sender.send(error).await; // Send error about no signatures written
        }
        drop(error_sender); // should be dropped automatically as it goes out of scope, but just in case
    })
}

pub fn failures_handle(
    failed_csv: String,
    mut recv_failed: tokio::sync::mpsc::Receiver<FailedDownload>,
    error_sender: tokio::sync::mpsc::Sender<Error>, // Additional parameter for error channel
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
    retry_times: u32,
    fasta_location: String,
    keep_fastas: bool,
    genomes_only: bool,
    proteomes_only: bool,
    download_only: bool,
    output_sigs: Option<String>,
    tmpdir: Option<String>,
) -> Result<(), anyhow::Error> {
    // if sig output provided but doesn't end in zip, bail
    let batch_size = 1000;
    if let Some(ref output_sigs) = output_sigs {
        if Path::new(&output_sigs)
            .extension()
            .map_or(true, |ext| ext != "zip")
        {
            bail!("Output must be a zip file.");
        }
    }
    // set up fasta download path
    let download_path = PathBuf::from(fasta_location);
    if !download_path.exists() {
        create_dir_all(&download_path)?;
    }
    // create channels. buffer size here is 4 b/c we can do 3 downloads simultaneously
    let (send_sigs, recv_sigs) = tokio::sync::mpsc::channel::<Vec<Signature>>(4);
    let (send_failed, recv_failed) = tokio::sync::mpsc::channel::<FailedDownload>(4);
    // Error channel for handling task errors
    let (error_sender, error_receiver) = tokio::sync::mpsc::channel::<anyhow::Error>(1);

    // if tmpdir is provided, read from it or create it.
    // Initialize an optional Manifest to hold existing signatures
    let mut existing_sigs: Option<Manifest> = None;

    // If tmpdir is provided, try to read from it or create it if it doesn't exist.
    if let Some(ref td) = tmpdir {
        let outdir: PathBuf = td.into();

        if outdir.exists() {
            match read_from_tmpdir(&outdir).await {
                Ok(manifest) => {
                    existing_sigs = Some(manifest);
                }
                Err(e) => {
                    let error = e.context("Failed to read from existing tmpdir");
                    let _ = error_sender.send(error).await;
                }
            }
        } else {
            // If the directory does not exist, create it
            if let Err(e) = tokio::fs::create_dir_all(&outdir).await {
                let error =
                    anyhow::Error::new(e).context("Failed to create directory at specified path");
                let _ = error_sender.send(error).await;
            }
        }
    }

    // Set up collector/writing tasks
    let mut handles = Vec::new();

    // write to tmpdir OR directly to zip
    // let sig_handle = match tmpdir {
    //     Some(ref td) => {
    //         tmpdir_sigwriter_handle(recv_sigs, td.clone(), output_sigs, error_sender.clone())
    //     }
    //     None => zipwriter_handle(recv_sigs, output_sigs, batch_size, error_sender.clone()),
    // };
    let sig_handle = zipwriter_handle(recv_sigs, output_sigs, batch_size, error_sender.clone());

    let failures_handle = failures_handle(failed_csv, recv_failed, error_sender.clone());
    let critical_error_flag = Arc::new(AtomicBool::new(false));
    let error_handle = error_handler(error_receiver, critical_error_flag.clone());
    handles.push(sig_handle);
    handles.push(failures_handle);
    handles.push(error_handle);

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
    // let dna_sig_templates = build_siginfo(&params_vec, "DNA");
    let dna_template_collection = build_template_collection(&params_vec, "DNA");
    // let dna_sig_collection: Collection = Collection::from_sigs(dna_sig_templates.clone())?;
    // let dna_template_manifest = dna_template_collection.manifest();
    // prot will build protein, dayhoff, hp
    // let prot_sig_templates = build_siginfo(&params_vec, "protein");
    let prot_template_collection = build_template_collection(&params_vec, "protein");
    // let prot_sig_collection: Collection = Collection::from_sigs(prot_sig_templates.clone())?;
    // let prot_sig_manifest = prot_sig_collection.manifest();

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
        let semaphore_clone = Arc::clone(&semaphore);
        let client_clone = Arc::clone(&client);
        let send_sigs = send_sigs.clone();
        let send_failed = send_failed.clone();
        let download_path_clone = download_path.clone(); // Clone the path for each task
        let send_errors = error_sender.clone();
        let dna_sigs = dna_template_collection.clone();
        let prot_sigs = prot_template_collection.clone();
        // clone existing sig manifest
        let e_siginfo = existing_sigs.clone();

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
                dna_sigs,
                prot_sigs,
                genomes_only,
                proteomes_only,
                download_only,
            )
            .await;
            match result {
                Ok((sigs, failed_downloads)) => {
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
    output_sigs: Option<String>,
    tmpdir: Option<String>,
) -> Result<(), anyhow::Error> {
    // if sig output provided but doesn't end in zip, bail
    if let Some(ref output_sigs) = output_sigs {
        if Path::new(&output_sigs)
            .extension()
            .map_or(true, |ext| ext != "zip")
        {
            bail!("Output must be a zip file.");
        }
    }
    // set up fasta download path
    let download_path = PathBuf::from(fasta_location);
    if !download_path.exists() {
        create_dir_all(&download_path)?;
    }

    // create channels. buffer size here is 4 b/c we can do 3 downloads simultaneously
    let (send_sigs, recv_sigs) = tokio::sync::mpsc::channel::<TemplateCollection>(4);
    let (send_failed, recv_failed) = tokio::sync::mpsc::channel::<FailedDownload>(4);
    // Error channel for handling task errors
    let (error_sender, error_receiver) = tokio::sync::mpsc::channel::<anyhow::Error>(1);

    // if tmpdir is provided, try to create it
    if let Some(ref td) = tmpdir {
        let outdir: PathBuf = td.into();
        if let Err(e) = tokio::fs::create_dir_all(outdir).await {
            let error =
                anyhow::Error::new(e).context("Failed to create directory at specified path");
            let _ = error_sender.send(error).await; // Send the error through the channel
        }
    }

    // Set up collector/writing tasks
    let mut handles = Vec::new();

    // write to tmpdir OR directly to zip
    let batch_size = 1000;
    // let sig_handle = match tmpdir {
    //     Some(ref td) => {
    //         tmpdir_sigwriter_handle(recv_sigs, td.clone(), output_sigs, error_sender.clone())
    //     }
    //     None => zipwriter_handle(recv_sigs, output_sigs, batch_size, error_sender.clone()),
    // };
    let sig_handle = zipwriter_handle(recv_sigs, output_sigs, batch_size, error_sender.clone());

    let failures_handle = failures_handle(failed_csv, recv_failed, error_sender.clone());
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
    let dna_sig_templates = build_siginfo(&params_vec, "DNA");
    let prot_sig_templates = build_siginfo(&params_vec, "protein");

    let mut genomes_only = false;
    let mut proteomes_only = false;

    // Check if dna_sig_templates is empty and not keep_fastas
    if dna_sig_templates.is_empty() && !keep_fastas {
        eprintln!("No DNA signature templates provided, and --keep-fastas is not set.");
        proteomes_only = true;
    }
    // Check if protein_sig_templates is empty and not keep_fastas
    if prot_sig_templates.is_empty() && !keep_fastas {
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
        let semaphore_clone = Arc::clone(&semaphore);
        let client_clone = Arc::clone(&client);
        let send_sigs = send_sigs.clone();
        let send_failed = send_failed.clone();
        let download_path_clone = download_path.clone(); // Clone the path for each task
        let send_errors = error_sender.clone();

        let dna_sigs = dna_sig_templates.clone();
        let prot_sigs = prot_sig_templates.clone();

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
                dna_sigs,
                prot_sigs,
                genomes_only,
                proteomes_only,
                download_only,
            )
            .await;
            match result {
                Ok((sigs, failed_downloads)) => {
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
