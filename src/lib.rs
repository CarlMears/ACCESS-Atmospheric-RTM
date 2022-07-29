//! RTM computation
//!
//! NOTE: this module is intended for the interface between Rust and Python. The
//! real work happens in the other modules, and they do not use `pyo3`, it's
//! only used here.

use pyo3::prelude::*;

use ndarray::Array2;
use numpy::{PyArray2, PyReadonlyArray2, PyReadonlyArray1, ToPyArray};

/// Atmospheric parameters.
///
/// This is just a container of multiple numpy arrays, each dimensioned as
/// (`num_points`, `num_freq`).
#[pyclass]
struct AtmoParameters {
    tran: Array2<f32>,
    tb_up: Array2<f32>,
    tb_down: Array2<f32>,
}

/// Implement all the "getters" for the Python properties
#[pymethods]
impl AtmoParameters {
    #[getter(tran)]
    fn convert_tran<'py>(&self, py: Python<'py>) -> &'py PyArray2<f32> {
        self.tran.to_pyarray(py)
    }

    #[getter(tb_up)]
    fn convert_tb_up<'py>(&self, py: Python<'py>) -> &'py PyArray2<f32> {
        self.tb_up.to_pyarray(py)
    }

    #[getter(tb_down)]
    fn convert_tb_down<'py>(&self, py: Python<'py>) -> &'py PyArray2<f32> {
        self.tb_down.to_pyarray(py)
    }
}

impl AtmoParameters {
    fn new(num_points: usize, num_freq: usize) -> Self {
        Self {
            tran: Array2::zeros([num_points, num_freq]),
            tb_up: Array2::zeros([num_points, num_freq]),
            tb_down: Array2::zeros([num_points, num_freq]),
        }
    }
}

/// Compute the RTM.
#[pyfunction]
fn compute_rtm(
    pressure: PyReadonlyArray1<f32>,
    temperature: PyReadonlyArray2<f32>,
    height: PyReadonlyArray2<f32>,
    specific_humidity: PyReadonlyArray2<f32>,
    liquid_content: PyReadonlyArray2<f32>,
    surface_temperature: PyReadonlyArray1<f32>,
    surface_height: PyReadonlyArray1<f32>,
    surface_dewpoint: PyReadonlyArray1<f32>,
    surface_pressure: PyReadonlyArray1<f32>,
    incidence_angle: PyReadonlyArray1<f32>,
    frequency: PyReadonlyArray1<f32>,
    num_levels: usize,
    num_points: usize,
    num_freq: usize,
) -> PyResult<AtmoParameters> {
    // TODO: implement
    let output = AtmoParameters::new(num_points, num_freq);
    Ok(output)
}

/// A Python module implemented in Rust.
#[pymodule]
fn access_atmosphere(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_rtm, m)?)?;
    m.add_class::<AtmoParameters>()?;
    Ok(())
}
