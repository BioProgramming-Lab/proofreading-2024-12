use itertools::Itertools;
use pyo3::{
    exceptions::{PyIOError, PyNotImplementedError, PyRuntimeError},
    prelude::*,
};

#[derive(thiserror::Error, Debug)]
enum Error {
    #[error("Failed to create csv reader: {0}")]
    CSVOpen(csv::Error),
    #[error("Failed to read csv record: {0}")]
    CSVRecord(csv::Error),
    #[error("Failed to find last data in csv record: {0:?}")]
    CSVRecordLastData(csv::StringRecord),
    #[error("Failed to parse data: {0}")]
    Parse(std::num::ParseFloatError),
    #[allow(unused)]
    #[error("todo")]
    Todo,
}

impl From<Error> for PyErr {
    fn from(e: Error) -> Self {
        match e {
            Error::CSVOpen(_) => PyIOError::new_err(e.to_string()),
            Error::CSVRecord(_) => PyIOError::new_err(e.to_string()),
            Error::CSVRecordLastData(_) => PyIOError::new_err(e.to_string()),
            Error::Parse(_) => PyRuntimeError::new_err(e.to_string()),
            Error::Todo => PyNotImplementedError::new_err("todo"),
        }
    }
}

/// Read a csv file of parameters for each (set of) simulation.
/// Each row is for a simulation.
/// return those data in a Vec<Vec<f32>> (list[list[float]]), same shape, n_sim * n_parameter
/// need to be converted to np.ndarray in python
#[pyfunction]
fn read_parameters(filename: &str) -> PyResult<Vec<Vec<f32>>> {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(filename)
        .map_err(Error::CSVOpen)?;
    // map and collect entire csv data, `Iterator<Result<StringRecord, Error>>` into `Result<Vec<Vec<f32>>>`
    Ok(rdr
        .records()
        // map and collect each `Result<StringRecord, Error>` into `Result<Vec<f32>>`
        .map(|res| match res {
            Ok(record) => record
                .iter()
                // convert each `&str` into `f32`
                .map(|v| v.parse().map_err(Error::Parse))
                .collect::<Result<Vec<f32>, _>>(),
            Err(e) => Err(Error::CSVRecord(e)),
        })
        .collect::<Result<Vec<_>, _>>()?)
}

/// Read a csv file of results for each simulation.
/// Each 8 rows is for a simulation.
///     c_A: Free A
///     c_B: Free B
///     c_C: Free A'
///     c_R: Receptor
///     c_AC: AA'
///     c_BC: BA'
///     c_AR: A-Receptor complex
///     c_BR: B-Receptor complex
/// Ditch all data but all species at the end of the grid.
/// return those data in a Vec<Vec<f32>> (list[float]), n_sim * n_species
/// need to be converted to np.ndarray in python
/// this function does not check if the n_species matches the data, so please check in python
#[pyfunction]
fn read_solve_at_end(filename: &str, n_species: usize) -> PyResult<Vec<Vec<f32>>> {
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(filename)
        .map_err(Error::CSVOpen)?;

    // map and collect entire csv data, `Iterator<Result<StringRecord, Error>>` into `Result<Vec<Vec<f32>>>`
    Ok(rdr
        .records()
        // map and collect each `Result<StringRecord, Error>` into `Result<Vec<f32>>`
        .map(|res| match res {
            Ok(record) => record
                // only retain the last data in each simulation
                .get(record.len() - 1)
                .map_or(Err(Error::CSVRecordLastData(record.clone())), |v| {
                    v.parse().map_err(Error::Parse)
                }),
            Err(e) => Err(Error::CSVRecord(e)),
        })
        .chunks(n_species)
        .into_iter()
        .map(|chunk| chunk.collect::<Result<Vec<f32>, _>>())
        .collect::<Result<Vec<_>, _>>()?)
}

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn rust_lib_proofreading(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_parameters, m)?)?;
    m.add_function(wrap_pyfunction!(read_solve_at_end, m)?)?;
    Ok(())
}
