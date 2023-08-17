//! Atmospheric radiative transfer model for the ACCESS project

mod core;

use crate::error::RtmError;
use crate::rtm::core::{abh2o_rk_modified, atm_tran, fdabsoxy_1992_modified, fdcldabs};
use smallvec::SmallVec;
use std::num::NonZeroUsize;

/// Input parameters for the RTM that are constant.
#[derive(Debug)]
pub struct RtmParameters {
    /// Microwave frequencies in GHz, with a length of `num_freqs`.
    frequency: SmallVec<[f32; 8]>,
    /// Earth incidence angle in degrees, with a length of `num_freqs`.
    incidence: SmallVec<[f32; 8]>,
}

/// Inputs for the RTM for a single point. Unlike [`RtmParameters`], these
/// values may vary over location/time.
#[derive(Debug)]
pub struct RtmInputs {
    /// Number of atmosphere profile levels.
    num_levels: NonZeroUsize,
    /// Starting index for the surface, aka `ibegin`.
    surface_index: usize,
    /// Pressure profile in hPa. This has length `num_levels+1` since the first
    /// element is for the surface.
    pressure: Vec<f32>,
    /// Temperature profile in K. This has length `num_levels+1` since the first
    /// element is for the surface.
    temperature: Vec<f32>,
    /// Water vapor pressure profile in hPa. This has length `num_levels+1` since the first
    /// element is for the surface.
    vapor_pressure: Vec<f32>,
    /// Liquid water density in g/m^3. This has length `num_levels+1` since the first
    /// element is for the surface.
    rho_l: Vec<f32>,
    /// Geometric height in m. This has length `num_levels+1` since the first
    /// element is for the surface.
    height: Vec<f32>,
}

/// Outputs from the RTM for a single point.
#[derive(Debug)]
pub struct RtmOutputs {
    /// Atmospheric transmissivity as a function of frequency index.
    pub tran: SmallVec<[f32; 8]>,
    /// Atmospheric upwelling in K as a function of frequency index.
    pub tb_up: SmallVec<[f32; 8]>,
    /// Atmospheric downwelling in K as a function of frequency index.
    pub tb_down: SmallVec<[f32; 8]>,
}

impl RtmParameters {
    pub fn new(freqs: &[f32], eia: &[f32]) -> Result<Self, RtmError> {
        if freqs.len() != eia.len() || freqs.is_empty() {
            return Err(RtmError::InconsistentInputs);
        }
        Ok(Self {
            frequency: SmallVec::from_slice(freqs),
            incidence: SmallVec::from_slice(eia),
        })
    }
}

impl RtmInputs {
    /// Prepare and convert values.
    ///
    /// The slices (`levels`, `temperature`, etc) must all be the same length.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        levels: &[f32],
        surface_temperature: f32,
        temperature: &[f32],
        surface_height: f32,
        height: &[f32],
        surface_dewpoint: f32,
        specific_humidity: &[f32],
        liquid_content: &[f32],
        surface_pressure: f32,
    ) -> Result<Self, RtmError> {
        #![allow(clippy::excessive_precision)]
        /// Mean radius of the Earth in meters
        const R_EARTH: f32 = 6371e3;

        /// Ideal gas constant (J/mol/K)
        const R: f32 = 8.3144598;
        /// Mean molar mass of dry air (g/mol)
        const M_DRY: f32 = 28.9644;
        /// Mean molar mass of water (g/mol)
        const M_H2O: f32 = 18.01528;
        /// Specific gas constant for dry air (J/g/K)
        const R_DRY: f32 = R / M_DRY;
        /// Specific gas constant for water vapor (J/g/K)
        const R_VAPOR: f32 = R / M_H2O;

        /// Coefficient for ratio between molar masses
        const EPSILON: f32 = M_H2O / M_DRY;
        /// Scaling factor using EPSILON
        const EPS_SCALE: f32 = (1. - EPSILON) / EPSILON;

        let num_levels: NonZeroUsize = levels
            .len()
            .try_into()
            .or(Err(RtmError::InconsistentInputs))?;

        // Find the starting index for the surface (aka `ibegin`). Note this
        // assumes that the levels are sorted in descending order (from high to
        // low pressure).
        let surface_index = (0..num_levels.get())
            .find(|&i| levels[i] <= surface_pressure)
            .ok_or(RtmError::NoSurface)?;

        // Prepend the surface value to these vectors
        let prepend_with = |level_data: &[f32], zero_value: f32, surface_value: f32| -> Vec<f32> {
            let mut prepended = Vec::with_capacity(num_levels.get() + 1);
            prepended.push(zero_value);
            prepended.extend_from_slice(level_data);

            prepended[surface_index] = surface_value;
            prepended
        };
        let pressure = prepend_with(levels, 0., surface_pressure);
        let temperature = prepend_with(temperature, surface_temperature, surface_temperature);
        let mut height = prepend_with(height, surface_height, surface_height);

        // Convert geopotential height to geometric height
        for z in &mut height {
            *z *= R_EARTH / (R_EARTH - *z);
        }
        if height[surface_index] >= height[surface_index + 1] {
            height[surface_index] = height[surface_index + 1] - 0.1;
        }

        // Convert specific humidity q to water vapor pressure P_v. The mass mixing
        // ratio w is:
        //
        // w = q / (1 - q)
        //
        // The vapor pressure is:
        //
        // P_v = (w P) / (R_dry/R_vapor + w)
        //
        // For the surface value, convert dewpoint to vapor pressure using the Buck equation.
        let pv = {
            let w = specific_humidity.iter().map(|q| q / (1. - q));

            let mut prepended = Vec::with_capacity(num_levels.get() + 1);
            prepended.push(buck_vap(surface_dewpoint));
            prepended.extend(
                levels
                    .iter()
                    .zip(w)
                    .map(|(p, w)| (w * p) / (R_DRY / R_VAPOR + w)),
            );

            prepended[surface_index] = prepended[0];
            prepended
        };

        // Specific liquid cloud mixing content
        let q_l = {
            let mut prepended = Vec::with_capacity(num_levels.get() + 1);
            prepended.push(0.);
            prepended.extend_from_slice(liquid_content);

            prepended[surface_index] = prepended[surface_index + 1];
            prepended
        };

        // Convert water mass mixing ratio to specific humidity
        // (https://earthscience.stackexchange.com/a/5077)
        //
        // This is a bit redundant since q_h2o (specific humidity) is an input.
        // However, the surface specific humidity is not given (it was converted
        // directly from dewpoint to vapor pressure) and also the `ibegin` above
        // modified the profile.
        let q_h2o = pressure.iter().zip(&pv).map(|(&p, &pv)| {
            if p > 0. {
                let w = (pv * R_DRY) / (R_VAPOR * (p - pv));
                w / (w + 1.)
            } else {
                0.
            }
        });

        // Convert specific cloud liquid water content (kg/kg) to liquid water
        // density (g/m^3).
        //
        // See here, section 4:
        // https://www.nwpsaf.eu/site/download/documentation/rtm/docs_rttov12/rttov_gas_cloud_aerosol_units.pdf
        // gas constant for humid air (J/gK)
        let r_moist = q_h2o.map(|q_h2o| R_DRY * (1. + EPS_SCALE * q_h2o));
        let rho_l: Vec<_> = q_l
            .iter()
            .zip(&pressure)
            .zip(&temperature)
            .zip(r_moist)
            .map(|(((q_l, p), t), r_moist)| q_l * (1e2 * p) / (r_moist * t))
            .collect();

        Ok(Self {
            num_levels,
            surface_index,
            pressure,
            temperature,
            height,
            vapor_pressure: pv,
            rho_l,
        })
    }

    /// Apply the RTM on the inputs for the given parameters.
    pub fn run(&self, parameters: &RtmParameters) -> RtmOutputs {
        // Scaling factor to convert from dB/km to Np/km
        let nep_scale = 0.1 * f32::ln(10.0);

        let mut tran = SmallVec::new();
        let mut tb_up = SmallVec::new();
        let mut tb_down = SmallVec::new();

        for (&freq, &inc) in parameters.frequency.iter().zip(&parameters.incidence) {
            // Build up total absorption coefficient profile
            let absorption_profile: SmallVec<[f32; 64]> = (self.surface_index
                ..self.num_levels.get() + 1)
                .map(|level_index| {
                    // Water vapor and oxygen absorption coefficients at this level converted to Np/km
                    let oxygen = fdabsoxy_1992_modified(
                        self.pressure[level_index],
                        self.temperature[level_index],
                        self.vapor_pressure[level_index],
                        freq,
                    ) * nep_scale;
                    let water = abh2o_rk_modified(
                        self.pressure[level_index],
                        self.temperature[level_index],
                        self.vapor_pressure[level_index],
                        freq,
                    ) * nep_scale;

                    // Cloud absorption coefficient in Np/km
                    let cloud = if self.rho_l[level_index] > 1.0e-7 {
                        fdcldabs(freq, self.temperature[level_index], self.rho_l[level_index])
                    } else {
                        0.0
                    };

                    // Total absorption coefficient at this level, converting from Np/km to Np/m
                    (water + oxygen + cloud) * 1.0e-3
                })
                .collect();

            let results = atm_tran(
                inc,
                &self.temperature[self.surface_index..],
                &self.height[self.surface_index..],
                &absorption_profile,
            );

            tran.push(results.0);
            tb_up.push(results.1);
            tb_down.push(results.2);
        }

        RtmOutputs {
            tran,
            tb_up,
            tb_down,
        }
    }
}

/// The Buck equation.
///
/// Convert `temp`, the temperature in K, into water vapor saturation pressure
/// in hPa. The equation is from [1], which cites Buck 1996.
///
/// To convert to water vapor partial pressure, multiply the result by the
/// relative humidity.
///
/// [1] https://en.wikipedia.org/wiki/Arden_Buck_equation
fn buck_vap(temp: f32) -> f32 {
    // Temperature in degrees Celsius
    let temp_c = temp - 273.15;
    6.1121 * f32::exp((18.678 - temp_c / 234.5) * (temp_c / (257.14 + temp_c)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// Check some values for the Buck equation. These are compared against the
    /// Fortran version.
    #[test]
    fn buck_vapor() {
        let input_temperatures = [273.15, 290., 300., 310.];
        let expected_outputs = [6.1121, 19.1925564, 35.3524513, 62.2872276];

        for (&temperature, &expected_output) in input_temperatures.iter().zip(&expected_outputs) {
            assert_relative_eq!(buck_vap(temperature), expected_output);
        }
    }
}
