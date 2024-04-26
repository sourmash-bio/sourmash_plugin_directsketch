/// Python interface Rust code for sourmash_plugin_directsketch.
use pyo3::prelude::*;

#[macro_use]
extern crate simple_error;

mod download;

use camino::Utf8PathBuf as PathBuf;

#[pyfunction]
fn do_gbsketch(
    input_csv: String,
    param_str: String,
    fasta_location: String,
    failed_csv: String,
    retry_times: usize,
    keep_fastas: bool,
) -> anyhow::Result<u8> {
    
    let runtime = tokio::runtime::Runtime::new().unwrap();
    
    runtime.block_on(async {
        download_accessions(input_csv, failed_csv, retry_times, location).await
    }).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Failed to process downloads: {}", e)))

    // match download::download_accessions(
    //     input_csv,
    //     failed_csv,
    //     retry_times,
    //     fasta_location,
    // ) {
    //     Ok(_) => Ok(0),
    //     Err(e) => {
    //         eprintln!("Error: {e}");
    //         Ok(1)
    //     }
    // }
}

#[pymodule]
fn download_accessions(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(do_downloads, m)?)?;
    Ok(())
}