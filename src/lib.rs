//! RTM computation
//!
//! NOTE: this module is intended for the interface between Rust and Python. The
//! real work happens in the other modules, and they do not use `pyo3`, it's
//! only used here.

pub(crate) mod error;
pub(crate) mod rtm;

use error::RtmError;
use ndarray::{Array2, Axis};
use numpy::{PyArray2, PyReadonlyArray1, PyReadonlyArray2, ToPyArray};
use pyo3::prelude::*;
use rayon::prelude::*;
use rtm::{RtmInputs, RtmParameters};

impl From<RtmError> for PyErr {
    fn from(e: RtmError) -> Self {
        use pyo3::exceptions::PyValueError;

        match e {
            RtmError::InconsistentInputs => PyValueError::new_err(e.to_string()),
            RtmError::NoSurface => PyValueError::new_err(e.to_string()),
            RtmError::NotContiguous => PyValueError::new_err(e.to_string()),
        }
    }
}

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
    let parameters = RtmParameters::new(frequency.as_slice()?, incidence_angle.as_slice()?)?;
    if parameters.len() != num_freq || pressure.len() != num_levels {
        return Err(RtmError::InconsistentInputs.into());
    }

    // Ensure everything is converted and contiguous
    let pressure = pressure.as_slice()?;
    let temperature = temperature.as_array();
    let height = height.as_array();
    let specific_humidity = specific_humidity.as_array();
    let liquid_content = liquid_content.as_array();
    let surface_temperature = surface_temperature.as_slice()?;
    let surface_height = surface_height.as_slice()?;
    let surface_dewpoint = surface_dewpoint.as_slice()?;
    let surface_pressure = surface_pressure.as_slice()?;

    let output = AtmoParameters::new(num_points, num_freq);

    (0..num_points)
        .into_par_iter()
        .try_for_each(|point| -> Result<(), RtmError> {
            let rtm_input = RtmInputs::new(
                &pressure,
                surface_temperature[point],
                temperature
                    .index_axis(Axis(0), point)
                    .as_slice()
                    .ok_or(RtmError::NotContiguous)?,
                surface_height[point],
                height
                    .index_axis(Axis(0), point)
                    .as_slice()
                    .ok_or(RtmError::NotContiguous)?,
                surface_dewpoint[point],
                specific_humidity
                    .index_axis(Axis(0), point)
                    .as_slice()
                    .ok_or(RtmError::NotContiguous)?,
                liquid_content
                    .index_axis(Axis(0), point)
                    .as_slice()
                    .ok_or(RtmError::NotContiguous)?,
                surface_pressure[point],
            )?;

            // TODO: remove later
            if point == 0 {
                dbg!(&rtm_input);
            }

            let output_point = rtm_input.run(&parameters);
            // output.tran;

            // TODO: update output
            // output

            Ok(())
        })?;

    // TODO: remove later
    dbg!(parameters);

    Ok(output)
}

/// A Python module implemented in Rust.
#[pymodule]
fn access_atmosphere(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_rtm, m)?)?;
    m.add_class::<AtmoParameters>()?;
    Ok(())
}
