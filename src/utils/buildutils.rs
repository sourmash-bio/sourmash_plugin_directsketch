//! sketching utilities

use anyhow::{anyhow, Context, Result};
use async_zip::base::write::ZipFileWriter;
use async_zip::{Compression, ZipDateTime, ZipEntryBuilder};
use camino::Utf8PathBuf;
use chrono::Utc;
use getset::{Getters, Setters};
use needletail::parser::SequenceRecord;
use needletail::{parse_fastx_file, parse_fastx_reader};
use serde::Serialize;
use sourmash::cmd::ComputeParameters;
use sourmash::encodings::HashFunctions;
use sourmash::manifest::Record;
use sourmash::signature::Signature;
use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt::Display;
use std::hash::{Hash, Hasher};
use std::io::{Cursor, Write};
use std::num::ParseIntError;
use std::str::FromStr;
use tokio::fs::File;
use tokio_util::compat::Compat;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BuildParams {
    pub ksize: u32,
    pub track_abundance: bool,
    pub num: u32,
    pub scaled: u64,
    pub seed: u32,
    pub moltype: HashFunctions,
}

impl Default for BuildParams {
    fn default() -> Self {
        BuildParams {
            ksize: 31,
            track_abundance: false,
            num: 0,
            scaled: 1000,
            seed: 42,
            moltype: HashFunctions::Murmur64Dna,
        }
    }
}

impl Hash for BuildParams {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.ksize.hash(state);
        self.track_abundance.hash(state);
        self.num.hash(state);
        self.scaled.hash(state);
        self.seed.hash(state);
        self.moltype.hash(state);
    }
}

impl BuildParams {
    pub fn calculate_hash(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.hash(&mut hasher);
        hasher.finish()
    }

    pub fn from_record(record: &Record) -> Self {
        let moltype = record.moltype(); // Get the moltype (HashFunctions enum)

        BuildParams {
            ksize: record.ksize(),
            track_abundance: record.with_abundance(),
            num: *record.num(),
            scaled: *record.scaled(),
            seed: 42,
            moltype,
        }
    }

    pub fn parse_ksize(value: &str) -> Result<u32, String> {
        value
            .parse::<u32>()
            .map_err(|_| format!("cannot parse k='{}' as a valid integer", value))
    }

    // disallow repeated values for scaled, num, seed
    pub fn parse_int_once<T>(
        value: &str,
        field: &str,
        current: &mut Option<T>,
    ) -> Result<(), String>
    where
        T: FromStr<Err = ParseIntError> + Display + Copy,
    {
        let parsed_value = value
            .parse::<T>()
            .map_err(|_| format!("cannot parse {}='{}' as a valid integer", field, value))?;

        // Check for conflicts; we don't allow multiple values for the same field.
        if let Some(old_value) = *current {
            return Err(format!(
                "Conflicting values for '{}': {} and {}",
                field, old_value, parsed_value
            ));
        }

        // Set the new value.
        *current = Some(parsed_value);

        Ok(())
    }

    pub fn parse_moltype(
        item: &str,
        current: &mut Option<HashFunctions>,
    ) -> Result<HashFunctions, String> {
        let new_moltype = match item {
            "protein" => HashFunctions::Murmur64Protein,
            "dna" => HashFunctions::Murmur64Dna,
            "dayhoff" => HashFunctions::Murmur64Dayhoff,
            "hp" => HashFunctions::Murmur64Hp,
            _ => return Err(format!("unknown moltype '{}'", item)),
        };

        // Check for conflicts and update the moltype.
        if let Some(existing) = current {
            if *existing != new_moltype {
                return Err(format!(
                    "Conflicting moltype settings in param string: '{}' and '{}'",
                    existing, new_moltype
                ));
            }
        }

        // Update the current value.
        *current = Some(new_moltype.clone());

        Ok(new_moltype)
    }

    pub fn parse_abundance(item: &str, current: &mut Option<bool>) -> Result<(), String> {
        let new_abundance = item == "abund";

        if let Some(existing) = *current {
            if existing != new_abundance {
                return Err(format!(
                    "Conflicting abundance settings in param string: '{}'",
                    item
                ));
            }
        }

        *current = Some(new_abundance);
        Ok(())
    }

    pub fn from_param_string(p_str: &str) -> Result<(Self, Vec<u32>), String> {
        let mut base_param = BuildParams::default();
        let mut ksizes = Vec::new();
        let mut moltype: Option<HashFunctions> = None;
        let mut track_abundance: Option<bool> = None;
        let mut num: Option<u32> = None;
        let mut scaled: Option<u64> = None;
        let mut seed: Option<u32> = None;

        for item in p_str.split(',') {
            match item {
                _ if item.starts_with("k=") => {
                    ksizes.push(Self::parse_ksize(&item[2..])?);
                }
                "abund" | "noabund" => {
                    Self::parse_abundance(item, &mut track_abundance)?;
                }
                "protein" | "dna" | "dayhoff" | "hp" => {
                    Self::parse_moltype(item, &mut moltype)?;
                }
                _ if item.starts_with("num=") => {
                    Self::parse_int_once(&item[4..], "num", &mut num)?;
                }
                _ if item.starts_with("scaled=") => {
                    Self::parse_int_once(&item[7..], "scaled", &mut scaled)?;
                }
                _ if item.starts_with("seed=") => {
                    Self::parse_int_once(&item[5..], "seed", &mut seed)?;
                }
                _ => return Err(format!("unknown component '{}' in params string", item)),
            }
        }

        // Ensure that num and scaled are mutually exclusive unless num is 0.
        if let (Some(n), Some(_)) = (num, scaled) {
            if n != 0 {
                return Err("Cannot specify both 'num' (non-zero) and 'scaled' in the same parameter string".to_string());
            }
        }

        // Apply parsed values to the base_param.
        if let Some(moltype) = moltype {
            base_param.moltype = moltype;
        }
        if let Some(track_abund) = track_abundance {
            base_param.track_abundance = track_abund;
        }
        if let Some(n) = num {
            base_param.num = n;
        }
        if let Some(s) = scaled {
            base_param.scaled = s;
        }
        if let Some(s) = seed {
            base_param.seed = s;
        }

        if ksizes.is_empty() {
            ksizes.push(base_param.ksize); // Use the default ksize if none were specified.
        }

        Ok((base_param, ksizes))
    }
}

#[derive(Debug)]
pub struct BuildParamsSet {
    params: HashSet<BuildParams>,
}

impl Default for BuildParamsSet {
    fn default() -> Self {
        let mut set = HashSet::new();
        set.insert(BuildParams::default());
        BuildParamsSet { params: set }
    }
}

impl BuildParamsSet {
    pub fn new() -> Self {
        Self {
            params: HashSet::new(),
        }
    }

    pub fn size(&self) -> usize {
        self.params.len()
    }

    pub fn insert(&mut self, param: BuildParams) {
        self.params.insert(param);
    }

    pub fn iter(&self) -> impl Iterator<Item = &BuildParams> {
        self.params.iter()
    }

    pub fn from_params_str(params_str: String) -> Result<Self, String> {
        if params_str.trim().is_empty() {
            return Err("Parameter string cannot be empty.".to_string());
        }

        let mut set = BuildParamsSet::new();

        for p_str in params_str.split('_') {
            let (base_param, ksizes) = BuildParams::from_param_string(p_str)?;

            for k in ksizes {
                let mut param = base_param.clone();
                param.ksize = k;
                set.insert(param);
            }
        }

        Ok(set)
    }

    pub fn get_params(&self) -> &HashSet<BuildParams> {
        &self.params
    }

    pub fn into_vec(self) -> Vec<BuildParams> {
        self.params.into_iter().collect()
    }
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
    pub fn from_buildparams(param: &BuildParams) -> Self {
        // Calculate the hash of Params
        let hashed_params = param.calculate_hash();

        BuildRecord {
            ksize: param.ksize,
            moltype: param.moltype.to_string(),
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

impl<'a> IntoIterator for &'a BuildManifest {
    type Item = &'a BuildRecord;
    type IntoIter = std::slice::Iter<'a, BuildRecord>;

    fn into_iter(self) -> Self::IntoIter {
        self.records.iter()
    }
}

impl<'a> IntoIterator for &'a mut BuildManifest {
    type Item = &'a mut BuildRecord;
    type IntoIter = std::slice::IterMut<'a, BuildRecord>;

    fn into_iter(self) -> Self::IntoIter {
        self.records.iter_mut()
    }
}

#[derive(Debug, Default, Clone)]
pub struct BuildCollection {
    pub manifest: BuildManifest,
    pub sigs: Vec<Signature>,
}

impl BuildCollection {
    pub fn new() -> Self {
        BuildCollection {
            manifest: BuildManifest::new(),
            sigs: Vec::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.manifest.is_empty() && self.sigs.is_empty()
    }

    pub fn size(&self) -> usize {
        self.manifest.size()
    }

    pub fn from_buildparams(params: &[BuildParams], input_moltype: &str) -> Self {
        let mut collection = BuildCollection::new();

        for param in params.iter().cloned() {
            collection.add_template_sig(param, input_moltype);
        }

        collection
    }

    pub fn from_buildparams_set(params_set: &BuildParamsSet, input_moltype: &str) -> Self {
        let mut collection = BuildCollection::new();

        for param in params_set.iter().cloned() {
            collection.add_template_sig(param, input_moltype);
        }

        collection
    }

    pub fn add_template_sig(&mut self, param: BuildParams, input_moltype: &str) {
        // Check the input_moltype against Params to decide if this should be added
        match input_moltype.to_lowercase().as_str() {
            "dna" if param.moltype != HashFunctions::Murmur64Dna => return, // Skip if it's not DNA
            "protein"
                if param.moltype != HashFunctions::Murmur64Protein
                    && param.moltype != HashFunctions::Murmur64Dayhoff
                    && param.moltype != HashFunctions::Murmur64Hp =>
            {
                return
            } // Skip if not a protein type
            _ => (),
        }

        // Adjust ksize for protein, dayhoff, or hp, which typically require tripling the k-mer size
        let adjusted_ksize = match param.moltype {
            HashFunctions::Murmur64Protein
            | HashFunctions::Murmur64Dayhoff
            | HashFunctions::Murmur64Hp => param.ksize * 3,
            _ => param.ksize,
        };

        // Construct ComputeParameters
        let cp = ComputeParameters::builder()
            .ksizes(vec![adjusted_ksize])
            .scaled(param.scaled)
            .protein(param.moltype == HashFunctions::Murmur64Protein)
            .dna(param.moltype == HashFunctions::Murmur64Dna)
            .dayhoff(param.moltype == HashFunctions::Murmur64Dayhoff)
            .hp(param.moltype == HashFunctions::Murmur64Hp)
            .num_hashes(param.num)
            .track_abundance(param.track_abundance)
            .build();

        // Create a Signature from the ComputeParameters
        let sig = Signature::from_params(&cp);

        // Create the BuildRecord using from_param
        let template_record = BuildRecord::from_buildparams(&param);

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
        for (record, sig) in self.iter_mut() {
            // update signature name, filename
            sig.set_name(name.as_str());
            sig.set_filename(filename.as_str());

            // update record: set name, filename, md5sum, n_hashes
            record.set_name(Some(name.clone()));
            record.set_filename(Some(filename.clone()));
            record.set_md5(Some(sig.md5sum()));
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_buildparams_consistent_hashing() {
        let params1 = BuildParams {
            ksize: 31,
            track_abundance: true,
            ..Default::default()
        };

        let params2 = BuildParams {
            ksize: 31,
            track_abundance: true,
            ..Default::default()
        };

        let hash1 = params1.calculate_hash();
        let hash2 = params2.calculate_hash();
        let hash3 = params2.calculate_hash();

        // Check that the hash for two identical Params is the same
        assert_eq!(hash1, hash2, "Hashes for identical Params should be equal");

        assert_eq!(
            hash2, hash3,
            "Hashes for the same Params should be consistent across multiple calls"
        );
    }

    #[test]
    fn test_buildparams_hashing_different() {
        let params1 = BuildParams {
            ksize: 31,
            ..Default::default()
        };

        let params2 = BuildParams {
            ksize: 21, // Changed ksize
            ..Default::default()
        };

        let params3 = BuildParams {
            ksize: 31,
            moltype: HashFunctions::Murmur64Protein,
            ..Default::default()
        };

        let hash1 = params1.calculate_hash();
        let hash2 = params2.calculate_hash();
        let hash3 = params3.calculate_hash();

        // Check that the hash for different Params is different
        assert_ne!(
            hash1, hash2,
            "Hashes for different Params should not be equal"
        );
        assert_ne!(
            hash1, hash3,
            "Hashes for different Params should not be equal"
        );
    }

    #[test]
    fn test_buildparams_generated_from_record() {
        // load signature + build record
        let mut filename = Utf8PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("tests/test-data/GCA_000175535.1.sig.gz");
        let path = filename.clone();

        let file = std::fs::File::open(filename).unwrap();
        let mut reader = std::io::BufReader::new(file);
        let sigs = Signature::load_signatures(
            &mut reader,
            Some(31),
            Some("DNA".try_into().unwrap()),
            None,
        )
        .unwrap();

        assert_eq!(sigs.len(), 1);

        let sig = sigs.get(0).unwrap();
        let record = Record::from_sig(sig, path.as_str());

        // create the expected Params based on the Record data
        let expected_params = BuildParams {
            ksize: 31,
            track_abundance: true,
            ..Default::default()
        };

        // // Generate the Params from the Record using the from_record method
        let generated_params = BuildParams::from_record(&record[0]);

        // // Assert that the generated Params match the expected Params
        assert_eq!(
            generated_params, expected_params,
            "Generated Params did not match the expected Params"
        );

        // // Calculate the hash for the expected Params
        let expected_hash = expected_params.calculate_hash();

        // // Calculate the hash for the generated Params
        let generated_hash = generated_params.calculate_hash();

        // // Assert that the hash for the generated Params matches the expected Params hash
        assert_eq!(
            generated_hash, expected_hash,
            "Hash of generated Params did not match the hash of expected Params"
        );
    }

    #[test]
    fn test_filter_removes_matching_buildparams() {
        let params1 = BuildParams {
            ksize: 31,
            track_abundance: true,
            ..Default::default()
        };

        let params2 = BuildParams {
            ksize: 21,
            track_abundance: true,
            ..Default::default()
        };

        let params3 = BuildParams {
            ksize: 31,
            scaled: 2000,
            track_abundance: true,
            ..Default::default()
        };

        let params_list = [params1.clone(), params2.clone(), params3.clone()];
        let mut build_collection = BuildCollection::from_buildparams(&params_list, "DNA");

        let mut params_set = HashSet::new();
        params_set.insert(params1.calculate_hash());
        params_set.insert(params3.calculate_hash());

        // Call the filter method
        build_collection.filter(&params_set);

        // Check that the records and signatures with matching params are removed
        assert_eq!(
            build_collection.manifest.records.len(),
            1,
            "Only one record should remain after filtering"
        );
        assert_eq!(
            build_collection.sigs.len(),
            1,
            "Only one signature should remain after filtering"
        );

        // Check that the remaining record is the one with hashed_params = 456
        let h2 = params2.calculate_hash();
        assert_eq!(
            build_collection.manifest.records[0].hashed_params, h2,
            "The remaining record should have hashed_params {}",
            h2
        );
    }

    #[test]
    fn test_build_params_set_default() {
        // Create a default BuildParamsSet.
        let default_set = BuildParamsSet::default();

        // Check that the set contains exactly one item.
        assert_eq!(
            default_set.size(),
            1,
            "Expected default BuildParamsSet to contain one default BuildParams."
        );

        // Get the default parameter from the set.
        let default_param = default_set.iter().next().unwrap();

        // Verify that the default parameter has expected default values.
        assert_eq!(default_param.ksize, 31, "Expected default ksize to be 31.");
        assert_eq!(
            default_param.track_abundance, false,
            "Expected default track_abundance to be false."
        );
        assert_eq!(
            default_param.moltype,
            HashFunctions::Murmur64Dna,
            "Expected default moltype to be DNA."
        );
        assert_eq!(default_param.num, 0, "Expected default num to be 0.");
        assert_eq!(
            default_param.scaled, 1000,
            "Expected default scaled to be 1000."
        );
        assert_eq!(default_param.seed, 42, "Expected default seed to be 42.");
    }

    #[test]
    fn test_valid_params_str() {
        let params_str = "k=31,abund,dna";
        let result = BuildParamsSet::from_params_str(params_str.to_string());

        assert!(
            result.is_ok(),
            "Expected 'k=31,abund,dna' to be valid, but got an error: {:?}",
            result
        );
        let params_set = result.unwrap();

        // Ensure the BuildParamsSet contains the expected number of parameters.
        assert_eq!(
            params_set.iter().count(),
            1,
            "Expected 1 BuildParams in the set"
        );

        // Verify that the BuildParams has the correct settings.
        let param = params_set.iter().next().unwrap();
        assert_eq!(param.ksize, 31);
        assert_eq!(param.track_abundance, true);
        assert_eq!(param.moltype, HashFunctions::Murmur64Dna);
    }

    #[test]
    fn test_multiple_valid_params_str() {
        let params_str = "k=31,abund,dna_k=21,k=31,k=51,abund_k=10,protein";
        let result = BuildParamsSet::from_params_str(params_str.to_string());

        assert!(
            result.is_ok(),
            "Param str {} is valid, but got an error: {:?}",
            params_str,
            result
        );
        let params_set: BuildParamsSet = result.unwrap();

        // Ensure the BuildParamsSet contains the expected number of parameters.
        // note that k=31,dna,abund is in two diff param strs, should only show up once.
        assert_eq!(params_set.size(), 4, "Expected 4 BuildParams in the set");
        // Define the expected BuildParams for comparison.
        let expected_params = vec![
            {
                let mut param = BuildParams::default();
                param.ksize = 31;
                param.track_abundance = true;
                param
            },
            {
                let mut param = BuildParams::default();
                param.ksize = 21;
                param.track_abundance = true;
                param
            },
            {
                let mut param = BuildParams::default();
                param.ksize = 51;
                param.track_abundance = true;
                param
            },
            {
                let mut param = BuildParams::default();
                param.ksize = 10;
                param.moltype = HashFunctions::Murmur64Protein;
                param
            },
        ];

        // Check that each expected BuildParams is present in the params_set.
        for expected_param in expected_params {
            assert!(
            params_set.get_params().contains(&expected_param),
            "Expected BuildParams with ksize: {}, track_abundance: {}, moltype: {:?} not found in the set",
            expected_param.ksize,
            expected_param.track_abundance,
            expected_param.moltype
        );
        }
    }

    #[test]
    fn test_invalid_params_str_conflicting_moltypes() {
        let params_str = "k=31,abund,dna,protein";
        let result = BuildParamsSet::from_params_str(params_str.to_string());

        assert!(
            result.is_err(),
            "Expected 'k=31,abund,dna,protein' to be invalid due to conflicting moltypes, but got a successful result"
        );

        // Check if the error message contains the expected conflict text.
        if let Err(e) = result {
            assert!(
                e.contains("Conflicting moltype settings"),
                "Expected error to contain 'Conflicting moltype settings', but got: {}",
                e
            );
        }
    }

    #[test]
    fn test_unknown_component_error() {
        // Test for an unknown component that should trigger an error.
        let result = BuildParamsSet::from_params_str("k=31,notaparam".to_string());
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "unknown component 'notaparam' in params string"
        );
    }

    #[test]
    fn test_unknown_component_error2() {
        // Test a common param string error (k=31,51 compared with valid k=31,k=51)
        let result = BuildParamsSet::from_params_str("k=31,51,abund".to_string());
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "unknown component '51' in params string"
        );
    }

    #[test]
    fn test_conflicting_num_and_scaled() {
        // Test for specifying both num and scaled, which should result in an error.
        let result = BuildParamsSet::from_params_str("k=31,num=10,scaled=1000".to_string());
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "Cannot specify both 'num' (non-zero) and 'scaled' in the same parameter string"
        );
    }

    #[test]
    fn test_conflicting_abundance() {
        // Test for providing conflicting abundance settings, which should result in an error.
        let result = BuildParamsSet::from_params_str("k=31,abund,noabund".to_string());
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "Conflicting abundance settings in param string: 'noabund'"
        );
    }

    #[test]
    fn test_invalid_ksize_format() {
        // Test for an invalid ksize format that should trigger an error.
        let result = BuildParamsSet::from_params_str("k=abc".to_string());
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "cannot parse k='abc' as a valid integer"
        );
    }

    #[test]
    fn test_invalid_num_format() {
        // Test for an invalid number format that should trigger an error.
        let result = BuildParamsSet::from_params_str("k=31,num=abc".to_string());
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "cannot parse num='abc' as a valid integer"
        );
    }

    #[test]
    fn test_invalid_scaled_format() {
        // Test for an invalid scaled format that should trigger an error.
        let result = BuildParamsSet::from_params_str("k=31,scaled=abc".to_string());
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "cannot parse scaled='abc' as a valid integer"
        );
    }

    #[test]
    fn test_invalid_seed_format() {
        // Test for an invalid seed format that should trigger an error.
        let result = BuildParamsSet::from_params_str("k=31,seed=abc".to_string());
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "cannot parse seed='abc' as a valid integer"
        );
    }

    #[test]
    fn test_repeated_values() {
        // repeated scaled
        let result = BuildParamsSet::from_params_str("k=31,scaled=1,scaled=1000".to_string());
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "Conflicting values for 'scaled': 1 and 1000"
        );

        // repeated num
        let result = BuildParamsSet::from_params_str("k=31,num=1,num=1000".to_string());
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "Conflicting values for 'num': 1 and 1000"
        );

        // repeated seed
        let result = BuildParamsSet::from_params_str("k=31,seed=1,seed=42".to_string());
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(
            result.unwrap_err(),
            "Conflicting values for 'seed': 1 and 42"
        );
    }

    #[test]
    fn test_missing_ksize() {
        // Test for a missing ksize, using default should not result in an error.
        let result = BuildParamsSet::from_params_str("abund".to_string());
        assert!(result.is_ok(), "Expected Ok but got an error.");
    }

    #[test]
    fn test_repeated_ksize() {
        // Repeated ksize settings should not trigger an error since it is valid to have multiple ksizes.
        let result = BuildParamsSet::from_params_str("k=31,k=21".to_string());
        assert!(result.is_ok(), "Expected Ok but got an error.");
    }

    #[test]
    fn test_empty_string() {
        // Test for an empty parameter string, which should now result in an error.
        let result = BuildParamsSet::from_params_str("".to_string());
        assert!(result.is_err(), "Expected an error but got Ok.");
        assert_eq!(result.unwrap_err(), "Parameter string cannot be empty.");
    }

    #[test]
    fn test_from_buildparams_abundance() {
        let mut params = BuildParams::default();
        params.track_abundance = true;

        // Create a BuildRecord using from_buildparams.
        let record = BuildRecord::from_buildparams(&params);

        // Check thqat all fields are set correctly.
        assert_eq!(record.ksize, 31, "Expected ksize to be 31.");
        assert_eq!(record.moltype, "DNA", "Expected moltype to be 'DNA'.");
        assert_eq!(record.scaled, 1000, "Expected scaled to be 1000.");
        assert_eq!(record.num, 0, "Expected num to be 0.");
        assert!(record.with_abundance, "Expected with_abundance to be true.");
        assert_eq!(
            record.hashed_params,
            params.calculate_hash(),
            "Expected the hashed_params to match the calculated hash."
        );
    }

    #[test]
    fn test_from_buildparams_protein() {
        let mut params = BuildParams::default();
        params.ksize = 10;
        params.scaled = 200;
        params.moltype = HashFunctions::Murmur64Protein;

        let record = BuildRecord::from_buildparams(&params);

        // Check that all fields are set correctly.
        assert_eq!(record.ksize, 10, "Expected ksize to be 10.");
        assert_eq!(record.moltype, "protein", "Expected moltype to be protein.");
        assert_eq!(record.scaled, 200, "Expected scaled to be 200.");
        assert_eq!(
            record.hashed_params,
            params.calculate_hash(),
            "Expected the hashed_params to match the calculated hash."
        );
    }

    #[test]
    fn test_from_buildparams_dayhoff() {
        let mut params = BuildParams::default();
        params.ksize = 10;
        params.moltype = HashFunctions::Murmur64Dayhoff;

        let record = BuildRecord::from_buildparams(&params);

        assert_eq!(record.ksize, 10, "Expected ksize to be 10.");
        assert_eq!(record.moltype, "dayhoff", "Expected moltype to be dayhoff.");
        // didn't change default scaled here, so should still be 1000
        assert_eq!(record.scaled, 1000, "Expected scaled to be 1000.");
        assert_eq!(
            record.hashed_params,
            params.calculate_hash(),
            "Expected the hashed_params to match the calculated hash."
        );
    }

    #[test]
    fn test_from_buildparams_hp() {
        let mut params = BuildParams::default();
        params.ksize = 10;
        params.moltype = HashFunctions::Murmur64Hp;

        let record = BuildRecord::from_buildparams(&params);

        assert_eq!(record.ksize, 10, "Expected ksize to be 10.");
        assert_eq!(record.moltype, "hp", "Expected moltype to be hp.");
        assert_eq!(
            record.hashed_params,
            params.calculate_hash(),
            "Expected the hashed_params to match the calculated hash."
        );
    }
}
