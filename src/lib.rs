use lazy_static::lazy_static;
use pyo3::exceptions::PyRuntimeError;
/// Python interface Rust code for sourmash_plugin_directsketch.
use pyo3::prelude::*;
use std::sync::Mutex;
use tokio::runtime::{Builder, Runtime};

// #[macro_use]
extern crate simple_error;

mod directsketch;
mod utils;

lazy_static! {
    static ref GLOBAL_RUNTIME: Mutex<Option<Runtime>> = Mutex::new(None);
}

#[pyfunction]
fn set_tokio_thread_pool(num_threads: usize) -> PyResult<usize> {
    let mut rt_lock = GLOBAL_RUNTIME.lock().unwrap();
    if rt_lock.is_none() {
        // Only initialize the runtime if it has not been initialized already
        let runtime = Builder::new_multi_thread()
            .worker_threads(num_threads)
            .enable_all()
            .build()
            .map_err(|e| PyErr::new::<PyRuntimeError, _>(e))?;

        *rt_lock = Some(runtime);
    }

    // Return the number of threads, which is now guaranteed to be set
    Ok(num_threads)
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn do_gbsketch(
    py: Python,
    input_csv: String,
    param_str: String,
    failed_csv: String,
    output_sigs: String,
    retry_times: u32,
    fasta_location: String,
    keep_fastas: bool,
    genomes_only: bool,
    proteomes_only: bool,
    download_only: bool,
) -> anyhow::Result<u8> {
    match directsketch::download_and_sketch(
        py,
        input_csv,
        output_sigs,
        param_str,
        failed_csv,
        retry_times,
        fasta_location,
        keep_fastas,
        genomes_only,
        proteomes_only,
        download_only,
    ) {
        Ok(_) => Ok(0),
        Err(e) => {
            eprintln!("Error: {e}");
            Ok(1)
        }
    }
}

#[pymodule]
fn sourmash_plugin_directsketch(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(do_gbsketch, m)?)?;
    m.add_function(wrap_pyfunction!(set_tokio_thread_pool, m)?)?;
    Ok(())
}
