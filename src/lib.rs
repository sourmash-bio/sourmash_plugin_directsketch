/// Python interface Rust code for sourmash_plugin_directsketch.
use pyo3::prelude::*;

// #[macro_use]
extern crate simple_error;

mod directsketch;
mod utils;

#[pyfunction]
fn set_global_thread_pool(num_threads: usize) -> PyResult<usize> {
    if std::panic::catch_unwind(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
    })
    .is_ok()
    {
        Ok(rayon::current_num_threads())
    } else {
        Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
            "Could not set the number of threads. Global thread pool might already be initialized.",
        ))
    }
}

#[pyfunction]
fn do_gbsketch(
    input_csv: String,
    param_str: String,
    failed_csv: String,
    output_sigs: String,
    retry_times: u32,
    fasta_location: String,
    keep_fastas: bool,
    genomes_only: bool,
    proteomes_only: bool,
) -> anyhow::Result<u8> {
    // let runtime = tokio::runtime::Runtime::new().unwrap();

    // runtime.block_on(async {
    //     download::download_accessions(input_csv, failed_csv, retry_times, fasta_location).await
    // }).map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Failed to process downloads: {}", e)))

    match directsketch::download_and_sketch(
        input_csv,
        output_sigs,
        param_str,
        failed_csv,
        retry_times,
        fasta_location,
        keep_fastas,
        genomes_only,
        proteomes_only,
    ) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pymodule]
fn sourmash_plugin_directsketch(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(do_gbsketch, m)?)?;
    m.add_function(wrap_pyfunction!(set_global_thread_pool, m)?)?;
    Ok(())
}
