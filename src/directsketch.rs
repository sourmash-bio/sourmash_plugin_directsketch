use anyhow::{anyhow, bail, Context, Result};
use async_zip::base::write::ZipFileWriter;
use async_zip::Compression;
use async_zip::{ZipDateTime, ZipEntryBuilder};
use camino::Utf8PathBuf as PathBuf;
use chrono::Utc;
use needletail::parse_fastx_reader;
use regex::Regex;
use reqwest::Client;
use std::collections::HashMap;
use std::fs::{self, create_dir_all};
use std::io::Cursor;
use std::os::unix::process;
use std::path::Path;
use tokio::fs::File;
use tokio::task;
use tokio_util::compat::Compat;

use pyo3::prelude::*;

use std::sync::Arc;
use tokio::sync::Mutex;

use tokio::sync::Semaphore;
use tokio::time::{interval, Duration};

use sourmash::manifest::{Manifest, Record};
use sourmash::signature::Signature;

use crate::utils::{build_siginfo, load_accession_info, parse_params_str};

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

    fn moltype(&self) -> String {
        match self {
            GenBankFileType::Genomic => "DNA".to_string(),
            GenBankFileType::Protein => "protein".to_string(),
            _ => "".to_string(),
        }
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
    eprintln!("base url for accession {}", base_url);
    let directory_response = client.get(&base_url).send().await?;
    if !directory_response.status().is_success() {
        eprintln!(
            "Failed to open genome directory: HTTP {}",
            directory_response.status()
        );
        return Err(anyhow!(
            "Failed to open genome directory: HTTP {}",
            directory_response.status()
        ));
    }
    eprintln!("getting response text for accession {}", base_url);
    let text = directory_response.text().await?;
    eprintln!("got response text for accession {}", base_url);
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
            eprintln!("base url: {}, clean_name: {}", base_url, clean_name);
            return Ok((format!("{}/{}", base_url, clean_name), clean_name.into()));
        }
    }

    eprintln!("No matching genome found for accession {}", accession);
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

async fn sketch_data(
    name: &str,
    filename: &str,
    compressed_data: Vec<u8>,
    mut sigs: Vec<Signature>,
    moltype: &str,
) -> Result<Vec<Signature>> {
    task::block_in_place(|| {
        let cursor = Cursor::new(compressed_data);

        let mut fastx_reader =
            parse_fastx_reader(cursor).context("Failed to parse FASTA/FASTQ data")?;

        // for each sig in template list, add sequence to sketch
        let mut set_name = false;
        while let Some(record) = fastx_reader.next() {
            let record = record.context("Failed to read record")?;
            sigs.iter_mut().for_each(|sig| {
                if !set_name {
                    sig.set_name(name);
                    sig.set_filename(filename);
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

        Ok(sigs)
    })
}

struct FailedDownload {
    accession: String,
    name: String,
    url: String,
    moltype: String,
}

#[allow(clippy::too_many_arguments)]
async fn dl_sketch_accession(
    client: &Client,
    accession: String,
    name: String,
    location: &PathBuf,
    retry: Option<u32>,
    keep_fastas: bool,
    dna_sigs: Vec<Signature>,
    prot_sigs: Vec<Signature>,
    genomes_only: bool,
    proteomes_only: bool,
    download_only: bool,
) -> Result<(Vec<Signature>, Vec<FailedDownload>)> {
    let retry_count = retry.unwrap_or(3); // Default retry count
    let mut sigs = Vec::<Signature>::new();
    let mut failed = Vec::<FailedDownload>::new();

    // keep track of any accessions for which we fail to find URLs
    let (base_url, full_name) = match fetch_genbank_filename(client, accession.as_str()).await {
        Ok(result) => result,
        Err(_err) => {
            eprintln!("download error for acc {}", accession);
            // Add accession to failed downloads with each moltype
            if !proteomes_only {
                let failed_download_dna = FailedDownload {
                    accession: accession.clone(),
                    name: name.clone(),
                    url: "".to_string(),
                    moltype: "dna".to_string(),
                };
                failed.push(failed_download_dna);
            }
            if !genomes_only {
                let failed_download_protein = FailedDownload {
                    accession: accession.clone(),
                    name: name.clone(),
                    url: "".to_string(),
                    moltype: "protein".to_string(),
                };
                failed.push(failed_download_protein);
            }
            eprintln!(
                "dl+sketched {} sigs for accession {}!",
                sigs.len(),
                accession
            );

            return Ok((sigs, failed));
        }
    };

    let mut file_types = vec![
        GenBankFileType::Genomic,
        GenBankFileType::Protein,
        // GenBankFileType::AssemblyReport,
        // GenBankFileType::Checksum, // Including standalone files like checksums here
    ];
    if genomes_only {
        file_types = vec![GenBankFileType::Genomic];
    } else if proteomes_only {
        file_types = vec![GenBankFileType::Protein];
    }

    for file_type in &file_types {
        let url = file_type.url(&base_url, &full_name);
        let data = match download_with_retry(client, &url, retry_count).await {
            Ok(data) => data,
            Err(_err) => {
                // here --> keep track of accession errors + filetype
                let failed_download = FailedDownload {
                    accession: accession.clone(),
                    name: name.clone(),
                    url: url.clone(),
                    moltype: file_type.moltype(),
                };
                failed.push(failed_download);
                continue;
            }
        };
        let file_name = file_type.filename(&accession);

        if keep_fastas {
            let path = location.join(&file_name);
            fs::write(&path, &data).context("Failed to write data to file")?;
        }
        if !download_only {
            // sketch data
            match file_type {
                GenBankFileType::Genomic => sigs.extend(
                    sketch_data(
                        name.as_str(),
                        file_name.as_str(),
                        data,
                        dna_sigs.clone(),
                        "dna",
                    )
                    .await?,
                ),
                GenBankFileType::Protein => {
                    sigs.extend(
                        sketch_data(
                            name.as_str(),
                            file_name.as_str(),
                            data,
                            prot_sigs.clone(),
                            "protein",
                        )
                        .await?,
                    );
                }
                _ => {} // Do nothing for other file types
            };
        }
    }

    Ok((sigs, failed))
}

async fn write_sig(
    sig: &Signature,
    md5sum_occurrences: &mut HashMap<String, usize>,
    manifest_rows: &mut Vec<Record>,
    // zip_writer: &mut ZipFileWriter<&mut File>,
    // zip_writer: &mut ZipFileWriter<Compat<&mut tokio::fs::File>>,
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
        .last_modification_date(ZipDateTime::from_chrono(&now));
    zip_writer
        .write_entry_whole(builder, &gzipped_buffer)
        .await
        .map_err(|e| anyhow!("Error writing zip entry: {}", e))
}

#[tokio::main]
#[allow(clippy::too_many_arguments)]
pub async fn download_and_sketch(
    py: Python,
    input_csv: String,
    output_sigs: String,
    param_str: String,
    failed_csv: String,
    retry_times: u32,
    fasta_location: String,
    keep_fastas: bool,
    genomes_only: bool,
    proteomes_only: bool,
    download_only: bool,
) -> Result<(), anyhow::Error> {
    let download_path = PathBuf::from(fasta_location);
    if !download_path.exists() {
        create_dir_all(&download_path)?;
    }
    let arc_download_path = Arc::new(download_path);

    // if sig output doesn't end in zip, bail
    if Path::new(&output_sigs)
        .extension()
        .map_or(true, |ext| ext != "zip")
    {
        bail!("Output must be a zip file.");
    }
    // start zip file; set up trackers
    let outpath: PathBuf = output_sigs.into();
    // let mut file = tokio::fs::File::create(outpath).await?;
    // let zip_f = ZipFileWriter::with_tokio(&mut file);
    // let mut zip_f = ZipFileWriter::with_tokio(&mut file);
    // let mut file = File::create(outpath).await?;
    // let mut file = tokio::fs::File::create(outpath).await?;
    // let arc_zip_writer = Arc::clone(&Arc::new(Mutex::new(ZipFileWriter::with_tokio(&mut file))));
    let file = tokio::fs::File::create(outpath).await?;
    let zip_writer = ZipFileWriter::with_tokio(file);
    let arc_zip_writer = Arc::new(Mutex::new(zip_writer));

    // let arc_zip_writer = Arc::new(Mutex::new(ZipFileWriter::with_tokio(&mut file)));
    // let arc_zip_writer = Arc::new(Mutex::new(ZipFileWriter::with_tokio(file)));
    // let file = tokio::fs::File::create(outpath).await?;
    // let compat_file = Compat::new(file);
    // let mut zip_writer = ZipFileWriter::with_tokio(&mut compat_file);

    // failures
    let fail_file = std::fs::File::create(failed_csv)?;
    let mut failed_writer = csv::Writer::from_writer(fail_file);
    failed_writer.write_record(["accession", "name", "moltype", "url"])?;
    let arc_failed_writer = Arc::new(Mutex::new(failed_writer));

    // let mut md5sum_occurrences: HashMap<String, usize> = HashMap::new();
    let mut md5sum_occurrences = Arc::new(Mutex::new(HashMap::new()));
    let mut arc_manifest_rows = Arc::new(Mutex::new(Vec::new()));
    // let mut manifest_rows: Vec<Record> = Vec::new();

    // Open the file containing the accessions synchronously
    let (accession_info, n_accs) = load_accession_info(input_csv)?;
    if n_accs == 0 {
        bail!("No accessions to download and sketch.")
    }

    // parse param string into params_vec, print error if fail
    let param_result = parse_params_str(param_str);
    let params_vec = match param_result {
        Ok(params) => params,
        Err(e) => {
            eprintln!("Error parsing params string: {}", e);
            bail!("Failed to parse params string");
        }
    };
    let dna_sig_templates = build_siginfo(&params_vec, "DNA");
    let prot_sig_templates = build_siginfo(&params_vec, "protein");

    // report every percent (or ever 1, whichever is larger)
    // let reporting_threshold = std::cmp::max(n_accs / 100, 1);

    let client = Arc::new(Client::new());
    let mut wrote_sigs = false;
    let semaphore = Arc::new(Semaphore::new(3)); // Allows up to 3 concurrent tasks
    let mut interval = tokio::time::interval(Duration::from_secs(1));
    let mut accessions = accession_info.into_iter().peekable();

    while accessions.peek().is_some() {
        interval.tick().await;
        py.check_signals()?; // If interrupted, return an Err automatically

        // Start up to 3 tasks at once per second
        for _ in 0..3 {
            if let Some(accinfo) = accessions.next() {
                let semaphore_clone = Arc::clone(&semaphore);
                let client_clone = Arc::clone(&client);
                let arc_zip_writer_clone = Arc::clone(&arc_zip_writer);
                let arc_failed_writer_clone = Arc::clone(&arc_failed_writer);
                let arc_manifest_rows_clone = Arc::clone(&arc_manifest_rows);
                let md5sum_occurrences_clone = Arc::clone(&md5sum_occurrences);
                let download_path_clone = Arc::clone(&arc_download_path);

                let dna_sigs = dna_sig_templates.clone();
                let prot_sigs = prot_sig_templates.clone();

                tokio::spawn(async move {
                    let _permit = semaphore_clone
                        .acquire()
                        .await
                        .expect("Failed to acquire semaphore permit");

                    match dl_sketch_accession(
                        &client_clone,
                        accinfo.accession,
                        accinfo.name,
                        &download_path_clone,
                        Some(retry_times),
                        keep_fastas,
                        dna_sigs,
                        prot_sigs,
                        genomes_only,
                        proteomes_only,
                        download_only,
                    )
                    .await
                    {
                        Ok((mut processed_sigs, failed_downloads)) => {
                            eprintln!("HERE, n sigs here: {:?}", processed_sigs.len());
                            let mut zip_w = arc_zip_writer_clone.lock().await;
                            let mut failed_w = arc_failed_writer_clone.lock().await;
                            let mut mf_rows = arc_manifest_rows_clone.lock().await;
                            let mut md5s = md5sum_occurrences_clone.lock().await;
                            for sig in processed_sigs.iter_mut() {
                                //keep track of whether we've written any sigs
                                if !wrote_sigs {
                                    wrote_sigs = true;
                                }
                                if let Err(e) =
                                    write_sig(sig, &mut md5s, &mut mf_rows, &mut *zip_w).await
                                {
                                    eprintln!("Error processing signature: {}", e);
                                }
                            }
                            for dl in failed_downloads {
                                if let Err(e) = failed_w.write_record(&[
                                    dl.accession,
                                    dl.name,
                                    dl.moltype,
                                    dl.url,
                                ]) {
                                    eprintln!("Error writing failed download record: {}", e);
                                }
                            }
                        }
                        Err(e) => eprintln!("Error during download and sketch: {}", e),
                    }
                });
            }
        }
    }

    // if no signatures were written, bail so user knows something went wrong
    if !wrote_sigs && !download_only {
        bail!("No signatures written.")
    }

    // Finalize the ZIP file and manifest
    let zip_writer_mutex = Arc::try_unwrap(arc_zip_writer)
        .map_err(|_| anyhow::anyhow!("Lock still has multiple owners"))?;

    let mut zip_writer = Mutex::into_inner(zip_writer_mutex);

    // write the manifest
    let manifest_filename = "SOURMASH-MANIFEST.csv".to_string();
    let mf_rows = arc_manifest_rows.lock().await;
    let manifest: Manifest = mf_rows.clone().into();
    // Create a temporary buffer to hold the manifest data
    let mut manifest_buffer = Vec::new();
    // Use manifest.to_writer to write the manifest to the buffer
    manifest.to_writer(&mut manifest_buffer)?;

    // write manifest to zipfile
    let now = Utc::now();
    let builder = ZipEntryBuilder::new(manifest_filename.into(), Compression::Stored)
        .last_modification_date(ZipDateTime::from_chrono(&now));

    zip_writer
        .write_entry_whole(builder, &manifest_buffer)
        .await?;
    // close zipfile
    zip_writer.close().await?;
    Ok(())
}

// progress report at threshold
// if (i + 1) % reporting_threshold == 0 {
//     let percent_processed = (((i + 1) as f64 / n_accs as f64) * 100.0).round();
//     println!(
//         "Starting accession {}/{} ({}%)",
//         (i + 1),
//         n_accs,
//         percent_processed
//     );
// }

// if i % 100 == 0 {
//     // Check for interrupt periodically
//     py.check_signals()?; // If interrupted, return an Err automatically
// }

// Process each accession
// let result = dl_sketch_accession(
//     &client,
//     accinfo.accession.clone(),
//     accinfo.name.clone(),
//     &download_path,
//     Some(retry_times),
//     keep_fastas,
//     dna_sig_templates.clone(),
//     prot_sig_templates.clone(),
//     genomes_only,
//     proteomes_only,
//     download_only,
// )
// .await;

//     tokio::spawn(async move {
//         let _permit = semaphore
//             .acquire()
//             .await
//             .expect("Failed to acquire semaphore permit");

//         let result = dl_sketch_accession(
//             &client,
//             accinfo.accession,
//             accinfo.name,
//             &download_path,
//             Some(retry_times),
//             keep_fastas,
//             dna_sigs,
//             prot_sigs,
//             genomes_only,
//             proteomes_only,
//             download_only,
//         )
//         .await;

//         if let Ok((mut processed_sigs, failed_downloads)) = result {
//             let mut zip_writer = arc_zip_writer.lock().await;
//             for sig in &mut processed_sigs {
//                 if !wrote_sigs {
//                     wrote_sigs = true;
//                 }
//                 write_sig(
//                     sig,
//                     &mut md5sum_occurrences,
//                     &mut manifest_rows,
//                     &mut zip_writer,
//                 )
//                 .await
//                 .map_err(|e| {
//                     eprintln!("Error processing signature: {}", e);
//                     e
//                 })?;
//             }
//             processed_sigs.clear(); // do we need this?
//             for dl in failed_downloads {
//                 failed_writer.write_record(&[dl.accession, dl.name, dl.moltype, dl.url])?;
//             }
//         }

//     });
