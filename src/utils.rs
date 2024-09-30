use anyhow::{anyhow, bail, Context, Result};
use async_zip::base::write::ZipFileWriter;
use async_zip::{Compression, ZipDateTime, ZipEntryBuilder};
use camino::Utf8PathBuf;
use chrono::Utc;
use futures::stream::{self, Stream, StreamExt};
use getset::{Getters, Setters};
use needletail::parser::SequenceRecord;
use needletail::{parse_fastx_file, parse_fastx_reader};
use reqwest::Url;
use serde::Serialize;
use sourmash::cmd::ComputeParameters;
use sourmash::collection::Collection;
use sourmash::encodings::Idx;
use sourmash::manifest::Record;
use sourmash::signature::Signature;
use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::io::{Cursor, Write};
use tokio::fs::File;
use tokio_util::compat::Compat;
// use rayon::prelude::*;

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

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Params {
    pub ksize: u32,
    pub track_abundance: bool,
    pub num: u32,
    pub scaled: u64,
    pub seed: u32,
    pub is_protein: bool,
    pub is_dayhoff: bool,
    pub is_hp: bool,
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
        self.is_dayhoff.hash(state);
        self.is_hp.hash(state);
        self.is_dna.hash(state);
    }
}

impl Params {
    pub fn from_record(record: &Record) -> Self {
        let moltype = record.moltype(); // Get the moltype (HashFunctions enum)

        Params {
            ksize: record.ksize(),
            track_abundance: record.with_abundance(),
            num: record.num().clone(),
            scaled: record.scaled().clone(),
            seed: 42,
            is_protein: moltype.protein(),
            is_dayhoff: moltype.dayhoff(),
            is_hp: moltype.hp(),
            is_dna: moltype.dna(),
        }
    }
}

pub fn build_params_hashset_from_records(records: &[&Record]) -> HashSet<u64> {
    records
        .iter()
        .map(|record| {
            let params = Params::from_record(record);
            // Create a hasher to hash the Params into a u64
            let mut hasher = DefaultHasher::new();
            params.hash(&mut hasher); // Hash the Params struct
            hasher.finish() // Get the hashed value (u64)
        })
        .collect()
}

#[derive(Debug, Default, Clone, Getters, Setters, Serialize)]
pub struct BuildRecord {
    // fields are ordered the same as Record to allow serialization to manifest
    // required fields are currently immutable once set
    #[getset(get = "pub", set = "pub")]
    internal_location: Option<Utf8PathBuf>,

    #[getset(get = "pub", set = "pub")]
    md5: Option<String>,

    #[getset(get = "pub", set = "pub")]
    md5short: Option<String>,

    #[getset(get_copy = "pub", set = "pub")]
    ksize: u32,

    moltype: String,

    #[getset(get = "pub")]
    num: u32,

    #[getset(get = "pub")]
    scaled: u64,

    #[getset(get = "pub", set = "pub")]
    n_hashes: Option<usize>,

    #[getset(get_copy = "pub", set = "pub")]
    #[serde(serialize_with = "intbool")]
    with_abundance: bool,

    #[getset(get = "pub", set = "pub")]
    name: Option<String>,

    #[getset(get = "pub", set = "pub")]
    filename: Option<String>,

    #[serde(skip)]
    pub hashed_params: u64,
}

// from sourmash (intbool is currently private there)
fn intbool<S>(x: &bool, s: S) -> std::result::Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    if *x {
        s.serialize_i32(1)
    } else {
        s.serialize_i32(0)
    }
}

impl BuildRecord {
    pub fn from_params(param: &Params, input_moltype: &str) -> Self {
        // Calculate the hash of Params
        let mut hasher = DefaultHasher::new();
        param.hash(&mut hasher);
        let hashed_params = hasher.finish();

        BuildRecord {
            ksize: param.ksize,
            moltype: input_moltype.to_string(),
            num: param.num,
            scaled: param.scaled,
            with_abundance: param.track_abundance,
            hashed_params,
            ..Default::default() // automatically set optional fields to None
        }
    }
}

#[derive(Debug, Default, Clone)]
pub struct BuildManifest {
    records: Vec<BuildRecord>,
}

impl BuildManifest {
    pub fn new() -> Self {
        BuildManifest {
            records: Vec::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn size(&self) -> usize {
        self.records.len()
    }

    // clear all records
    pub fn clear(&mut self) {
        self.records.clear();
    }

    pub fn add_record(&mut self, record: BuildRecord) {
        self.records.push(record);
    }

    pub fn extend_records(&mut self, other: impl IntoIterator<Item = BuildRecord>) {
        self.records.extend(other);
    }

    pub fn extend_from_manifest(&mut self, other: &BuildManifest) {
        self.records.extend(other.records.clone()); // Clone the records from the other manifest
    }

    pub fn to_writer<W: Write>(&self, mut wtr: W) -> Result<()> {
        // Write the manifest version as a comment
        wtr.write_all(b"# SOURMASH-MANIFEST-VERSION: 1.0\n")?;

        // Use CSV writer to serialize records
        let mut csv_writer = csv::Writer::from_writer(wtr);

        for record in &self.records {
            csv_writer.serialize(record)?; // Serialize each BuildRecord
        }

        csv_writer.flush()?; // Ensure all data is written

        Ok(())
    }

    /// Asynchronously writes manifest to a zip file
    pub async fn async_write_manifest_to_zip(
        &self,
        zip_writer: &mut ZipFileWriter<Compat<File>>,
    ) -> Result<()> {
        let manifest_filename = "SOURMASH-MANIFEST.csv".to_string();
        let mut manifest_buffer = Vec::new();

        // Serialize the manifest using `to_writer`
        self.to_writer(&mut manifest_buffer)
            .map_err(|e| anyhow!("Error serializing manifest: {}", e))?;

        // create the ZipEntryBuilder with the current time and permissions
        let now = Utc::now();
        let builder = ZipEntryBuilder::new(manifest_filename.into(), Compression::Stored)
            .last_modification_date(ZipDateTime::from_chrono(&now))
            .unix_permissions(0o644);

        // Write the manifest buffer to the zip file asynchronously
        zip_writer
            .write_entry_whole(builder, &manifest_buffer)
            .await
            .map_err(|e| anyhow!("Error writing manifest to zip: {}", e))?;

        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct BuildCollection {
    pub manifest: BuildManifest,
    pub sigs: Vec<Signature>,
    pub name: Option<String>,
    pub filename: Option<String>,
}

impl BuildCollection {
    pub fn new() -> Self {
        BuildCollection {
            manifest: BuildManifest::new(),
            sigs: Vec::new(),
            name: None,
            filename: None,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.manifest.is_empty() && self.sigs.is_empty()
    }

    pub fn size(&self) -> usize {
        self.manifest.size()
    }

    pub fn is_compatible(&self, other: &BuildCollection) -> bool {
        self.name == other.name && self.filename == other.filename
    }

    pub fn from_params(params: &[Params], input_moltype: &str) -> Self {
        let mut collection = BuildCollection::new();

        for param in params.iter().cloned() {
            collection.add_template_sig(param, input_moltype);
        }

        collection
    }

    pub fn add_template_sig(&mut self, param: Params, input_moltype: &str) {
        // Check the input_moltype against Params to decide if this should be added
        match input_moltype {
            "dna" | "DNA" if !param.is_dna => return, // Skip if it's not the correct moltype
            "protein" if !param.is_protein && !param.is_dayhoff && !param.is_hp => return,
            _ => (),
        }

        let adjusted_ksize = if param.is_protein || param.is_dayhoff || param.is_hp {
            param.ksize * 3
        } else {
            param.ksize
        };

        // Construct ComputeParameters
        let cp = ComputeParameters::builder()
            .ksizes(vec![adjusted_ksize])
            .scaled(param.scaled)
            .protein(param.is_protein)
            .dna(param.is_dna)
            .dayhoff(param.is_dayhoff)
            .hp(param.is_hp)
            .num_hashes(param.num)
            .track_abundance(param.track_abundance)
            .build();

        // Create a Signature from the ComputeParameters
        let sig = Signature::from_params(&cp);

        // Create the BuildRecord using from_param
        let template_record = BuildRecord::from_params(&param, input_moltype);

        // Add the record and signature to the collection
        self.manifest.records.push(template_record);
        self.sigs.push(sig);
    }

    pub fn filter(&mut self, params_set: &HashSet<u64>) {
        let mut index = 0;
        while index < self.manifest.records.len() {
            let record = &self.manifest.records[index];

            // filter records with matching Params
            if params_set.contains(&record.hashed_params) {
                self.manifest.records.remove(index);
                self.sigs.remove(index);
            } else {
                index += 1;
            }
        }
    }

    pub fn sigs_iter_mut(&mut self) -> impl Iterator<Item = &mut Signature> {
        self.sigs.iter_mut()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = (&mut BuildRecord, &mut Signature)> {
        // zip together mutable iterators over records and sigs
        self.manifest.records.iter_mut().zip(self.sigs.iter_mut())
    }

    // Parallel version of iter_mut using rayon's par_iter_mut
    // pub fn par_iter_mut(&mut self) -> impl ParallelIterator<Item = (&mut BuildRecord, &mut Signature)> {
    //     self.manifest.records.par_iter_mut().zip(self.sigs.par_iter_mut())
    // }

    pub fn build_sigs_from_data(
        &mut self,
        data: Vec<u8>,
        input_moltype: &str, // (protein/dna); todo - use hashfns?
        name: String,
        filename: String,
    ) -> Result<()> {
        let cursor = Cursor::new(data);
        let mut fastx_reader =
            parse_fastx_reader(cursor).context("Failed to parse FASTA/FASTQ data")?;

        // Iterate over FASTA records and add sequences/proteins to sigs
        while let Some(record) = fastx_reader.next() {
            let record = record.context("Failed to read record")?;
            self.sigs_iter_mut().for_each(|sig| {
                if input_moltype == "protein" {
                    sig.add_protein(&record.seq())
                        .expect("Failed to add protein");
                } else {
                    sig.add_sequence(&record.seq(), true)
                        .expect("Failed to add sequence");
                    // if not force, panics with 'N' in dna sequence
                }
            });
        }

        // After processing sequences, update sig, record information
        self.update_info(name, filename);

        Ok(())
    }

    pub fn build_sigs_from_file(
        &mut self,
        input_moltype: &str, // (protein/dna); todo - use hashfns?
        name: String,
        filename: String,
    ) -> Result<()> {
        let mut fastx_reader = parse_fastx_file(&filename)?;
        // Iterate over FASTA records and add sequences/proteins to sigs
        while let Some(record) = fastx_reader.next() {
            let record = record.context("Failed to read record")?;
            self.sigs_iter_mut().for_each(|sig| {
                if input_moltype == "protein" {
                    sig.add_protein(&record.seq())
                        .expect("Failed to add protein");
                } else {
                    sig.add_sequence(&record.seq(), true)
                        .expect("Failed to add sequence");
                    // if not force, panics with 'N' in dna sequence
                }
            });
        }

        // After processing sequences, update sig, record information
        self.update_info(name, filename);

        Ok(())
    }

    pub fn build_singleton_sigs(
        &mut self,
        record: SequenceRecord,
        input_moltype: &str, // (protein/dna); todo - use hashfns?
        filename: String,
    ) -> Result<()> {
        self.sigs_iter_mut().for_each(|sig| {
            if input_moltype == "protein" {
                sig.add_protein(&record.seq())
                    .expect("Failed to add protein");
            } else {
                sig.add_sequence(&record.seq(), true)
                    .expect("Failed to add sequence");
                // if not force, panics with 'N' in dna sequence
            }
        });
        let record_name = std::str::from_utf8(record.id())
            .expect("could not get record id")
            .to_string();
        // After processing sequences, update sig, record information
        self.update_info(record_name, filename);

        Ok(())
    }

    pub fn update_info(&mut self, name: String, filename: String) {
        // update the records to reflect information the signature;
        self.name = Some(name.clone());
        self.filename = Some(filename.clone());
        for (record, sig) in self.iter_mut() {
            // update signature name, filename
            sig.set_name(name.as_str());
            sig.set_filename(&filename.as_str());

            // update record: set name, filename, md5sum, n_hashes
            record.set_name(Some(name.clone()));
            record.set_filename(Some(filename.clone()));
            record.set_md5(Some(sig.md5sum().into()));
            record.set_md5short(Some(sig.md5sum()[0..8].into()));
            record.set_n_hashes(Some(sig.size()));

            // note, this needs to be set when writing sigs
            // record.set_internal_location("")
        }
    }

    pub async fn async_write_sigs_to_zip(
        &mut self, // need mutable to update records
        zip_writer: &mut ZipFileWriter<Compat<File>>,
        md5sum_occurrences: &mut HashMap<String, usize>,
    ) -> Result<()> {
        // iterate over both records and signatures
        for (record, sig) in self.iter_mut() {
            let md5sum_str = sig.md5sum();
            let count = md5sum_occurrences.entry(md5sum_str.clone()).or_insert(0);
            *count += 1;

            // Generate the signature filename
            let sig_filename = if *count > 1 {
                format!("signatures/{}_{}.sig.gz", md5sum_str, count)
            } else {
                format!("signatures/{}.sig.gz", md5sum_str)
            };

            // update record's internal_location with the signature filename
            record.internal_location = Some(sig_filename.clone().into());

            // serialize signature to JSON
            let wrapped_sig = vec![sig.clone()];
            let json_bytes = serde_json::to_vec(&wrapped_sig)
                .map_err(|e| anyhow!("Error serializing signature: {}", e))?;

            // gzip
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

            // write to zip file
            let now = Utc::now();
            let builder = ZipEntryBuilder::new(sig_filename.into(), Compression::Stored)
                .last_modification_date(ZipDateTime::from_chrono(&now))
                .unix_permissions(0o644);

            zip_writer
                .write_entry_whole(builder, &gzipped_buffer)
                .await
                .map_err(|e| anyhow!("Error writing zip entry for signature: {}", e))?;
        }

        Ok(())
    }
}

impl<'a> IntoIterator for &'a mut BuildCollection {
    type Item = (&'a mut BuildRecord, &'a mut Signature);
    type IntoIter =
        std::iter::Zip<std::slice::IterMut<'a, BuildRecord>, std::slice::IterMut<'a, Signature>>;

    fn into_iter(self) -> Self::IntoIter {
        self.manifest.records.iter_mut().zip(self.sigs.iter_mut())
    }
}

#[derive(Debug, Clone)]
pub struct MultiBuildCollection {
    pub collections: Vec<BuildCollection>,
}

impl MultiBuildCollection {
    pub fn new() -> Self {
        MultiBuildCollection {
            collections: Vec::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.collections.is_empty()
    }

    pub fn add_collection(&mut self, collection: &mut BuildCollection) {
        self.collections.push(collection.clone())
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
        let mut is_dayhoff = false;
        let mut is_hp = false;
        let mut is_dna = false;

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
                }
                "dna" => {
                    is_dna = true;
                }
                "dayhoff" => {
                    is_dayhoff = true;
                }
                "hp" => {
                    is_hp = true;
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
                is_dayhoff,
                is_hp,
            };
            unique_params.insert(param);
        }
    }

    Ok(unique_params.into_iter().collect())
}

// this should be replaced with branchwater's MultiCollection when it's ready
#[derive(Clone)]
pub struct MultiCollection {
    collections: Vec<Collection>,
}

impl MultiCollection {
    fn new(collections: Vec<Collection>, contains_revindex: bool) -> Self {
        Self { collections }
    }

    pub fn from_zipfile(sigpath: &Utf8PathBuf) -> Result<Self> {
        match Collection::from_zipfile(sigpath) {
            Ok(collection) => Ok(MultiCollection::new(vec![collection], false)),
            Err(_) => bail!("failed to load zipfile: '{}'", sigpath),
        }
    }

    // Load from multiple zip files
    pub fn from_zipfiles(sigpaths: &[Utf8PathBuf]) -> Result<Self> {
        let mut collections = Vec::new();

        for sigpath in sigpaths {
            match Collection::from_zipfile(sigpath) {
                Ok(collection) => collections.push(collection),
                Err(_) => bail!("failed to load zipfile: '{}'", sigpath),
            }
        }

        Ok(MultiCollection::new(collections, false))
    }

    pub fn stream_iter(&self) -> impl Stream<Item = (&Collection, Idx, &Record)> + '_ {
        stream::iter(&self.collections).flat_map(|collection| {
            stream::iter(
                collection
                    .iter()
                    .map(move |(idx, record)| (collection, idx, record)),
            )
        })
    }

    pub fn iter(&self) -> impl Iterator<Item = (&Collection, Idx, &Record)> + '_ {
        self.collections.iter().flat_map(|collection| {
            collection
                .iter()
                .map(move |(idx, record)| (collection, idx, record))
        })
    }
}
