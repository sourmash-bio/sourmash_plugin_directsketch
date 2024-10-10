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
use sourmash::encodings::{HashFunctions, Idx};
use sourmash::errors::SourmashError;
use sourmash::manifest::Record;
use sourmash::selection::{Select, Selection};
use sourmash::signature::Signature;
use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt::Display;
use std::hash::{Hash, Hasher};
use std::io::{Cursor, Write};
use std::num::ParseIntError;
use std::ops::Index;
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

#[derive(Debug, Clone, Getters, Setters, Serialize)]
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

    #[getset(get_copy = "pub")]
    #[serde(skip)]
    pub seed: u32,

    #[serde(skip)]
    pub hashed_params: u64,

    #[serde(skip)]
    pub sequence_added: bool,
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

impl Default for BuildRecord {
    fn default() -> Self {
        // Default BuildRecord is DNA default
        BuildRecord {
            internal_location: None,
            md5: None,
            md5short: None,
            ksize: 31,
            moltype: "DNA".to_string(),
            num: 0,
            scaled: 1000,
            n_hashes: None,
            with_abundance: false,
            name: None,
            filename: None,
            seed: 42,
            hashed_params: 0,
            sequence_added: false,
        }
    }
}

impl BuildRecord {
    pub fn default_dna() -> Self {
        Self {
            ..Default::default()
        }
    }

    pub fn default_protein() -> Self {
        Self {
            moltype: "protein".to_string(),
            ksize: 10,
            scaled: 200,
            ..Default::default()
        }
    }

    pub fn default_dayhoff() -> Self {
        Self {
            moltype: "dayhoff".to_string(),
            ksize: 10,
            scaled: 200,
            ..Default::default()
        }
    }

    pub fn default_hp() -> Self {
        Self {
            moltype: "hp".to_string(),
            ksize: 10,
            scaled: 200,
            ..Default::default()
        }
    }

    pub fn moltype(&self) -> HashFunctions {
        self.moltype.as_str().try_into().unwrap()
    }

    pub fn from_record(record: &Record) -> Self {
        Self {
            ksize: record.ksize(),
            moltype: record.moltype().to_string(),
            num: *record.num(),
            scaled: *record.scaled(),
            with_abundance: record.with_abundance(),
            ..Default::default() // ignore remaining fields
        }
    }

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

impl PartialEq for BuildRecord {
    fn eq(&self, other: &Self) -> bool {
        self.ksize == other.ksize
            && self.moltype == other.moltype
            && self.with_abundance == other.with_abundance
            && self.num == other.num
            && self.scaled == other.scaled
    }
}

impl Eq for BuildRecord {}

impl Hash for BuildRecord {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.ksize.hash(state);
        self.moltype.hash(state);
        self.scaled.hash(state);
        self.num.hash(state);
        self.with_abundance.hash(state);
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

    pub fn iter(&self) -> impl Iterator<Item = &BuildRecord> {
        self.records.iter()
    }

    // clear all records
    pub fn clear(&mut self) {
        self.records.clear();
    }

    pub fn filter_manifest(&self, other: &BuildManifest) -> Self {
        // Create a HashSet of references to the `BuildRecord`s in `other`
        let pairs: HashSet<_> = other.records.iter().collect();

        // Filter `self.records` to retain only those `BuildRecord`s that are NOT in `pairs`
        let records = self
            .records
            .iter()
            .filter(|&build_record| !pairs.contains(build_record))
            .cloned()
            .collect();

        Self { records }
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

impl Select for BuildManifest {
    fn select(self, selection: &Selection) -> Result<Self, SourmashError> {
        let rows = self.records.iter().filter(|row| {
            let mut valid = true;
            valid = if let Some(ksize) = selection.ksize() {
                row.ksize == ksize
            } else {
                valid
            };
            valid = if let Some(abund) = selection.abund() {
                valid && row.with_abundance == abund
            } else {
                valid
            };
            valid = if let Some(moltype) = selection.moltype() {
                valid && row.moltype() == moltype
            } else {
                valid
            };
            valid = if let Some(scaled) = selection.scaled() {
                // num sigs have row.scaled = 0, don't include them
                valid && row.scaled != 0 && row.scaled <= scaled as u64
            } else {
                valid
            };
            valid = if let Some(num) = selection.num() {
                valid && row.num == num
            } else {
                valid
            };
            valid
        });

        Ok(BuildManifest {
            records: rows.cloned().collect(),
        })
    }
}

impl From<Vec<BuildRecord>> for BuildManifest {
    fn from(records: Vec<BuildRecord>) -> Self {
        BuildManifest { records }
    }
}

impl Index<usize> for BuildManifest {
    type Output = BuildRecord;

    fn index(&self, index: usize) -> &Self::Output {
        &self.records[index]
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
        self.manifest.is_empty()
    }

    pub fn size(&self) -> usize {
        self.manifest.size()
    }

    pub fn dna_size(&self) -> Result<usize, SourmashError> {
        let selection = Selection::builder()
            .moltype(HashFunctions::Murmur64Dna)
            .build();

        let selected_manifest = self.manifest.clone().select(&selection)?;

        Ok(selected_manifest.records.len())
    }

    pub fn protein_size(&self) -> Result<usize, SourmashError> {
        let selection = Selection::builder()
            .moltype(HashFunctions::Murmur64Protein)
            .build();

        let selected_manifest = self.manifest.clone().select(&selection)?;

        Ok(selected_manifest.records.len())
    }

    pub fn anyprotein_size(&self) -> Result<usize, SourmashError> {
        // Create selections for each protein type.
        let protein_selection = Selection::builder()
            .moltype(HashFunctions::Murmur64Protein)
            .build();
        let dayhoff_selection = Selection::builder()
            .moltype(HashFunctions::Murmur64Dayhoff)
            .build();
        let hp_selection = Selection::builder()
            .moltype(HashFunctions::Murmur64Hp)
            .build();

        // Apply each selection and sum the sizes of the selected manifests.
        let protein_count = self
            .manifest
            .clone()
            .select(&protein_selection)?
            .records
            .len();
        let dayhoff_count = self
            .manifest
            .clone()
            .select(&dayhoff_selection)?
            .records
            .len();
        let hp_count = self.manifest.clone().select(&hp_selection)?.records.len();

        Ok(protein_count + dayhoff_count + hp_count)
    }

    pub fn parse_ksize(value: &str) -> Result<u32, String> {
        value
            .parse::<u32>()
            .map_err(|_| format!("cannot parse k='{}' as a valid integer", value))
    }

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

        *current = Some(parsed_value);
        Ok(())
    }

    pub fn parse_moltype(item: &str, current: &mut Option<String>) -> Result<String, String> {
        let new_moltype = match item {
            "protein" | "dna" | "dayhoff" | "hp" => item.to_string(),
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

    pub fn parse_params(p_str: &str) -> Result<(BuildRecord, Vec<u32>), String> {
        let mut ksizes = Vec::new();
        let mut moltype: Option<String> = None;
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
                "protein" | "dna" | "DNA" | "dayhoff" | "hp" => {
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

        // Create a moltype-specific default BuildRecord.
        let mut base_record = match moltype.as_deref() {
            Some("dna") => BuildRecord::default_dna(),
            Some("DNA") => BuildRecord::default_dna(),
            Some("protein") => BuildRecord::default_protein(),
            Some("dayhoff") => BuildRecord::default_dayhoff(),
            Some("hp") => BuildRecord::default_hp(),
            _ => BuildRecord::default_dna(), // no moltype --> assume DNA
        };

        // Apply parsed values
        if let Some(track_abund) = track_abundance {
            base_record.with_abundance = track_abund;
        }
        if let Some(n) = num {
            base_record.num = n;
        }
        if let Some(s) = scaled {
            base_record.scaled = s;
        }
        if let Some(s) = seed {
            base_record.seed = s as u32;
        }

        // Use the default ksize if none were specified.
        if ksizes.is_empty() {
            ksizes.push(base_record.ksize);
        }

        // Ensure that num and scaled are mutually exclusive unless num is 0.
        if let (Some(n), Some(_)) = (num, scaled) {
            if n != 0 {
                return Err("Cannot specify both 'num' (non-zero) and 'scaled' in the same parameter string".to_string());
            }
        }

        Ok((base_record, ksizes))
    }

    pub fn from_param_str(params_str: &str) -> Result<Self, String> {
        if params_str.trim().is_empty() {
            return Err("Parameter string cannot be empty.".to_string());
        }

        let mut coll = BuildCollection::new();
        let mut seen_records = HashSet::new();

        for p_str in params_str.split('_') {
            // Use `parse_params` to get the base record and ksizes.
            let (base_record, ksizes) = Self::parse_params(p_str)?;

            // Iterate over each ksize and add a signature to the collection.
            for k in ksizes {
                let mut record = base_record.clone();
                record.ksize = k;

                // Check if the record is already in the set.
                if seen_records.insert(record.clone()) {
                    // Add the record and its associated signature to the collection.
                    // coll.add_template_sig_from_record(&record, &record.moltype);
                    coll.add_template_sig_from_record(&record);
                }
            }
        }
        Ok(coll)
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

    // pub fn from_manifest(manifest: &BuildManifest, input_moltype: &str) -> Self {
    pub fn from_manifest(manifest: &BuildManifest) -> Self {
        let mut collection = BuildCollection::new();

        // Iterate over each `BuildRecord` in the provided `BuildManifest`.
        for record in &manifest.records {
            // Add a signature to the collection using the `BuildRecord` and `input_moltype`.
            // collection.add_template_sig_from_record(record, input_moltype);
            collection.add_template_sig_from_record(record);
        }

        collection
    }

    pub fn add_template_sig_from_record(&mut self, record: &BuildRecord) {
        //, input_moltype: &str) {
        // Check the input_moltype against the `record`'s moltype to decide if this should be added.
        // this is because we don't currently allow translation --> modify to allow translate.
        // match input_moltype.to_lowercase().as_str() {
        //     "dna" if record.moltype != "DNA" => return, // Skip if it's not DNA.
        //     "protein"
        //         if record.moltype != "protein"
        //             && record.moltype != "dayhoff"
        //             && record.moltype != "hp" =>
        //     {
        //         return;
        //     } // Skip if not a protein type.
        //     _ => (),
        // }

        // Adjust ksize for protein, dayhoff, or hp, which typically require tripling the k-mer size.
        let adjusted_ksize = match record.moltype.as_str() {
            "protein" | "dayhoff" | "hp" => record.ksize * 3,
            _ => record.ksize,
        };

        // Construct ComputeParameters.
        let cp = ComputeParameters::builder()
            .ksizes(vec![adjusted_ksize])
            .scaled(record.scaled)
            .protein(record.moltype == "protein")
            .dna(record.moltype == "DNA")
            .dayhoff(record.moltype == "dayhoff")
            .hp(record.moltype == "hp")
            .num_hashes(record.num)
            .track_abundance(record.with_abundance)
            .build();

        // Create a Signature from the ComputeParameters.
        let sig = Signature::from_params(&cp);

        // Clone the `BuildRecord` and use it directly.
        let template_record = record.clone();

        // Add the record and signature to the collection.
        self.manifest.records.push(template_record);
        self.sigs.push(sig);
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

    // pub fn filter_manifest(&mut self, other: &BuildManifest) {
    //     self.manifest = self.manifest.filter_manifest(other)
    // }

    pub fn filter_by_manifest(&mut self, other: &BuildManifest) {
        // Create a HashSet for efficient filtering based on the `BuildRecord`s in `other`.
        let other_records: HashSet<_> = other.records.iter().collect();

        // Retain only the records that are not in `other_records`, filtering in place.
        let mut sig_index = 0;
        self.manifest.records.retain(|record| {
            let keep = !other_records.contains(record);
            if !keep {
                // Remove the corresponding signature at the same index.
                self.sigs.remove(sig_index);
            } else {
                sig_index += 1; // Only increment if we keep the record and signature.
            }
            keep
        });
    }

    // filter template signatures that had no sequence added
    // suggested use right before writing signatures
    pub fn filter_empty(&mut self) {
        let mut sig_index = 0;

        self.manifest.records.retain(|record| {
            // Keep only records where `sequence_added` is `true`.
            let keep = record.sequence_added;

            if !keep {
                // Remove the corresponding signature at the same index if the record is not kept.
                self.sigs.remove(sig_index);
            } else {
                sig_index += 1; // Only increment if we keep the record and signature.
            }

            keep
        });
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

    pub fn iter(&self) -> impl Iterator<Item = (Idx, &BuildRecord)> {
        self.manifest.iter().enumerate().map(|(i, r)| (i as Idx, r))
    }

    pub fn record_for_dataset(&self, dataset_id: Idx) -> Result<&BuildRecord> {
        Ok(&self.manifest[dataset_id as usize])
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
            self.iter_mut().for_each(|(rec, sig)| {
                if input_moltype == "protein"
                    && (rec.moltype() == HashFunctions::Murmur64Protein
                        || rec.moltype() == HashFunctions::Murmur64Dayhoff
                        || rec.moltype() == HashFunctions::Murmur64Hp)
                {
                    sig.add_protein(&record.seq())
                        .expect("Failed to add protein");
                    if !rec.sequence_added {
                        rec.sequence_added = true
                    }
                } else {
                    sig.add_sequence(&record.seq(), true)
                        .expect("Failed to add sequence");
                    // if not force, panics with 'N' in dna sequence
                    if !rec.sequence_added {
                        rec.sequence_added = true
                    }
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
            self.iter_mut().for_each(|(rec, sig)| {
                if input_moltype == "protein"
                    && (rec.moltype() == HashFunctions::Murmur64Protein
                        || rec.moltype() == HashFunctions::Murmur64Dayhoff
                        || rec.moltype() == HashFunctions::Murmur64Hp)
                {
                    sig.add_protein(&record.seq())
                        .expect("Failed to add protein");
                    if !rec.sequence_added {
                        rec.sequence_added = true
                    }
                } else {
                    sig.add_sequence(&record.seq(), true)
                        .expect("Failed to add sequence");
                    // if not force, panics with 'N' in dna sequence
                    if !rec.sequence_added {
                        rec.sequence_added = true
                    }
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
        self.iter_mut().for_each(|(rec, sig)| {
            if input_moltype == "protein"
                && (rec.moltype() == HashFunctions::Murmur64Protein
                    || rec.moltype() == HashFunctions::Murmur64Dayhoff
                    || rec.moltype() == HashFunctions::Murmur64Hp)
            {
                sig.add_protein(&record.seq())
                    .expect("Failed to add protein");
                if !rec.sequence_added {
                    rec.sequence_added = true
                }
            } else {
                sig.add_sequence(&record.seq(), true)
                    .expect("Failed to add sequence");
                // if not force, panics with 'N' in dna sequence
                if !rec.sequence_added {
                    rec.sequence_added = true
                }
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

    // to do -- figure out how to skip sigs or warn user if no sequence was added to the template
    // could do the following, BUT then any BuildManifest we're adding may get out of sync.
    // would need to modify the collection, e.g. filter_empty or similar
    // if !rec.sequence_added{
    //     continue
    // }
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

// impl Select for BuildCollection {
//     fn select(mut self, selection: &Selection) -> Result<Self, SourmashError> {
//         self.manifest = self.manifest.select(selection)?;
//         Ok(self)
//     }
// }

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
    fn test_valid_params_str() {
        let params_str = "k=31,abund,dna";
        let result = BuildCollection::parse_params(params_str);

        assert!(
            result.is_ok(),
            "Expected 'k=31,abund,dna' to be valid, but got an error: {:?}",
            result
        );

        let (record, ksizes) = result.unwrap();

        // Verify that the Record, ksizes have the correct settings.
        assert_eq!(record.moltype, "DNA");
        assert_eq!(record.with_abundance, true);
        assert_eq!(ksizes, vec![31]);
        assert_eq!(record.scaled, 1000, "Expected default scaled value of 1000");
        assert_eq!(record.num, 0, "Expected default num value of 0");
    }

    #[test]
    fn test_from_param_str() {
        let params_str = "k=31,abund,dna_k=21,k=31,k=51,abund_k=10,protein";
        let coll_result = BuildCollection::from_param_str(params_str);

        assert!(
            coll_result.is_ok(),
            "Param str '{}' is valid, but got an error: {:?}",
            params_str,
            coll_result
        );

        let coll = coll_result.unwrap();

        // Ensure the BuildCollection contains the expected number of records.
        // Note that "k=31,abund,dna" appears in two different parameter strings, so it should only appear once.
        assert_eq!(
            coll.manifest.records.len(),
            4,
            "Expected 4 unique BuildRecords in the collection, but found {}",
            coll.manifest.records.len()
        );

        // Define the expected BuildRecords for comparison.
        let expected_records = vec![
            BuildRecord {
                ksize: 31,
                moltype: "DNA".to_string(),
                with_abundance: true,
                ..Default::default()
            },
            BuildRecord {
                ksize: 21,
                moltype: "DNA".to_string(),
                with_abundance: true,
                ..Default::default()
            },
            BuildRecord {
                ksize: 51,
                moltype: "DNA".to_string(),
                with_abundance: true,
                ..Default::default()
            },
            BuildRecord::default_protein(),
        ];

        // Verify that each expected BuildRecord is present in the collection.
        for expected_record in expected_records {
            assert!(
                coll.manifest.records.contains(&expected_record),
                "Expected BuildRecord with ksize: {}, moltype: {}, with_abundance: {} not found in the collection",
                expected_record.ksize,
                expected_record.moltype,
                expected_record.with_abundance
            );
        }

        // Optionally, check that the corresponding signatures are present.
        assert_eq!(
            coll.sigs.len(),
            4,
            "Expected 4 Signatures in the collection, but found {}",
            coll.sigs.len()
        );
    }

    #[test]
    fn test_invalid_params_str_conflicting_moltypes() {
        let params_str = "k=31,abund,dna,protein";
        let result = BuildCollection::from_param_str(params_str);

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

    #[test]
    fn test_filter_by_manifest_with_matching_records() {
        // Create a BuildCollection with some records and signatures.

        let rec1 = BuildRecord::default_dna();
        let rec2 = BuildRecord {
            ksize: 21,
            moltype: "DNA".to_string(),
            scaled: 1000,
            ..Default::default()
        };
        let rec3 = BuildRecord {
            ksize: 31,
            moltype: "DNA".to_string(),
            scaled: 1000,
            with_abundance: true,
            ..Default::default()
        };

        let bmanifest = BuildManifest {
            records: vec![rec1.clone(), rec2.clone(), rec3.clone()],
        };
        // let mut dna_build_collection = BuildCollection::from_manifest(&bmanifest, "DNA");
        let mut dna_build_collection = BuildCollection::from_manifest(&bmanifest);

        // Create a BuildManifest with records to filter out.
        let filter_manifest = BuildManifest {
            records: vec![rec1],
        };

        // Apply the filter.
        dna_build_collection.filter_by_manifest(&filter_manifest);

        // check that the default DNA sig remains
        assert_eq!(dna_build_collection.manifest.size(), 2);

        let remaining_records = &dna_build_collection.manifest.records;

        assert!(remaining_records.contains(&rec2));
        assert!(remaining_records.contains(&rec3));
    }

    #[test]
    fn test_add_template_sig_from_record() {
        // Create a BuildCollection.
        let mut build_collection = BuildCollection::new();

        // Create a DNA BuildRecord.
        let dna_record = BuildRecord {
            ksize: 31,
            moltype: "DNA".to_string(),
            scaled: 1000,
            with_abundance: true,
            ..Default::default()
        };

        // Add the DNA record to the collection with a matching moltype.
        // build_collection.add_template_sig_from_record(&dna_record, "DNA");
        build_collection.add_template_sig_from_record(&dna_record);

        // Verify that the record was added.
        assert_eq!(build_collection.manifest.records.len(), 1);
        assert_eq!(build_collection.sigs.len(), 1);

        let added_record = &build_collection.manifest.records[0];
        assert_eq!(added_record.moltype, "DNA");
        assert_eq!(added_record.ksize, 31);
        assert_eq!(added_record.with_abundance, true);

        // Create a protein BuildRecord.
        let protein_record = BuildRecord {
            ksize: 10,
            moltype: "protein".to_string(),
            scaled: 200,
            with_abundance: false,
            ..Default::default()
        };

        // Add the protein record to the collection with a matching moltype.
        // build_collection.add_template_sig_from_record(&protein_record, "protein");
        build_collection.add_template_sig_from_record(&protein_record);

        // Verify that the protein record was added and ksize adjusted.
        assert_eq!(build_collection.manifest.records.len(), 2);
        assert_eq!(build_collection.sigs.len(), 2);

        let added_protein_record = &build_collection.manifest.records[1];
        assert_eq!(added_protein_record.moltype, "protein");
        assert_eq!(added_protein_record.ksize, 10);
        assert_eq!(added_protein_record.with_abundance, false);

        // Create a BuildRecord with a non-matching moltype.
        let non_matching_record = BuildRecord {
            ksize: 10,
            moltype: "dayhoff".to_string(),
            scaled: 200,
            with_abundance: true,
            ..Default::default()
        };

        // Attempt to add the non-matching record with "DNA" as input moltype.
        // this is because we currently don't allow translation
        // build_collection.add_template_sig_from_record(&non_matching_record, "DNA");

        // Verify that the non-matching record was not added.
        // assert_eq!(build_collection.manifest.records.len(), 2);
        // assert_eq!(build_collection.sigs.len(), 2);

        // Add the same non-matching record with a matching input moltype.
        build_collection.add_template_sig_from_record(&non_matching_record);

        // Verify that the record was added.
        assert_eq!(build_collection.manifest.records.len(), 3);
        assert_eq!(build_collection.sigs.len(), 3);

        let added_dayhoff_record = &build_collection.manifest.records[2];
        assert_eq!(added_dayhoff_record.moltype, "dayhoff");
        assert_eq!(added_dayhoff_record.ksize, 10);
        assert_eq!(added_dayhoff_record.with_abundance, true);
    }

    #[test]
    fn test_filter_empty() {
        // Create a parameter string that generates BuildRecords with different `sequence_added` values.
        let params_str = "k=31,abund,dna_k=21,protein_k=10,abund";

        // Use `from_param_str` to build a `BuildCollection`.
        let mut build_collection = BuildCollection::from_param_str(params_str)
            .expect("Failed to build BuildCollection from params_str");

        // Manually set `sequence_added` for each record to simulate different conditions.
        build_collection.manifest.records[0].sequence_added = true; // Keep this record.
        build_collection.manifest.records[1].sequence_added = false; // This record should be removed.
        build_collection.manifest.records[2].sequence_added = true; // Keep this record.

        // Check initial sizes before filtering.
        assert_eq!(
            build_collection.manifest.records.len(),
            3,
            "Expected 3 records before filtering, but found {}",
            build_collection.manifest.records.len()
        );
        assert_eq!(
            build_collection.sigs.len(),
            3,
            "Expected 3 signatures before filtering, but found {}",
            build_collection.sigs.len()
        );

        // Apply the `filter_empty` method.
        build_collection.filter_empty();

        // After filtering, only the records with `sequence_added == true` should remain.
        assert_eq!(
            build_collection.manifest.records.len(),
            2,
            "Expected 2 records after filtering, but found {}",
            build_collection.manifest.records.len()
        );

        // Check that the signatures also match the remaining records.
        assert_eq!(
            build_collection.sigs.len(),
            2,
            "Expected 2 signatures after filtering, but found {}",
            build_collection.sigs.len()
        );

        // Verify that the remaining records have `sequence_added == true`.
        assert!(
            build_collection
                .manifest
                .records
                .iter()
                .all(|rec| rec.sequence_added),
            "All remaining records should have `sequence_added == true`"
        );
    }
}
