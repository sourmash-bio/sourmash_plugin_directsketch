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

    // Check if pytest is running
    let pytest_running = std::env::var("PYTEST_RUNNING").is_ok();

    // Check if runtime is already initialized
    if rt_lock.is_some() {
        if pytest_running {
            // If pytest is running, simply return the number of threads without error
            return Ok(num_threads);
        } else {
            // If not under pytest, return an error on reinitialization attempts
            return Err(PyErr::new::<PyRuntimeError, _>(
                "Tokio runtime is already initialized.",
            ));
        }
    }

    // Initialize the runtime if not already initialized
    let runtime = Builder::new_multi_thread()
        .worker_threads(num_threads)
        .enable_all()
        .build()
        .map_err(PyErr::new::<PyRuntimeError, _>)?;

    *rt_lock = Some(runtime);

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
