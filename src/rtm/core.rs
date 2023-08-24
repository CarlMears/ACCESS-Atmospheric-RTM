//! Core atmosphere RTM functions.

use smallvec::SmallVec;

use super::{
    liquid_cloud::fdcldabs, oxygen::fdabsoxy_1992_modified, water_vapor::abh2o_rk_modified,
};

/// Compute the absorption coefficient for an atmospheric layer.
///
/// For a pressure (hPa), temperature (K), water vapor partial pressure (hPa),
/// liquid water density (g/m³), compute the layer absorption coefficient in
/// Np/m.
///
/// This is a wrapper to the lower-level absorption coefficient functions.
pub(crate) fn layer_absorption(
    pressure: f32,
    temperature: f32,
    vapor_pressure: f32,
    liquid_water_density: f32,
    frequency: f32,
) -> f32 {
    /// Scaling factor to convert from dB/km to Np/km: `0.1 * ln(10)`
    const NEP_SCALE: f32 = 0.1 * std::f32::consts::LN_10;

    // Water vapor and oxygen absorption coefficients at this level converted to Np/km
    let oxygen =
        fdabsoxy_1992_modified(pressure, temperature, vapor_pressure, frequency) * NEP_SCALE;
    let water = abh2o_rk_modified(pressure, temperature, vapor_pressure, frequency) * NEP_SCALE;

    // Cloud absorption coefficient in Np/km
    let cloud = if liquid_water_density > 1.0e-7 {
        fdcldabs(frequency, temperature, liquid_water_density)
    } else {
        0.0
    };

    // Total absorption coefficient at this level, converting from Np/km to Np/m
    (water + oxygen + cloud) * 1.0e-3
}

/// Compute total atmospheric parameters from level data.
///
/// For an Earth incidence angle `inc` in degrees, and profile data where `t` is
/// the temperature in K, `z` is the elevation in m, and `tabs` is the
/// atmospheric absorption coefficient in Np/m, compute the output tuple
/// (`tran`, `tb_up`, `tb_down`) for the atmospheric transmissivity, atmospheric
/// upwelling brightness temperature in K, and atmospheric downwelling
/// brightness temperature in K.
///
/// The three profile inputs (`t`, `z`, and `tabs`) all have the same length,
/// `num_levels + 1`, where the first index `0` is the value at the surface and
/// indices from `1` to `num_levels` are profile data above the surface.
pub(crate) fn atm_tran(inc: f32, t: &[f32], z: &[f32], tabs: &[f32]) -> (f32, f32, f32) {
    const DELTA: f32 = 0.00035;

    // Differential slant height
    let dsdh = (1.0 + DELTA) / f32::sqrt(inc.to_radians().cos().powi(2) + DELTA * (2.0 + DELTA));

    // Number of levels *not* including the surface
    let num_levels = t.len() - 1;

    let opacity: SmallVec<[f32; 64]> = (1..=num_levels)
        .map(|i| -dsdh * 0.5 * (tabs[i - 1] + tabs[i]) * (z[i] - z[i - 1]))
        .collect();
    let t_avg: SmallVec<[f32; 64]> = (1..=num_levels).map(|i| 0.5 * (t[i - 1] + t[i])).collect();
    let ems: SmallVec<[f32; 64]> = opacity.iter().map(|opacity| 1.0 - opacity.exp()).collect();

    let (sum_down, _sum_op) = (1..=num_levels).fold((0., 0.), |(sum_down, sum_op), i| {
        (
            sum_down + (t_avg[i - 1] - t[1]) * ems[i - 1] * f32::exp(sum_op),
            sum_op + opacity[i - 1],
        )
    });

    let (sum_up, sum_op) = (1..=num_levels)
        .rev()
        .fold((0., 0.), |(sum_up, sum_op), i| {
            (
                sum_up + (t_avg[i - 1] - t[1]) * ems[i - 1] * f32::exp(sum_op),
                sum_op + opacity[i - 1],
            )
        });

    let tran = sum_op.exp();
    let tb_avg = (1. - tran) * t[1];
    let tb_down = tb_avg + sum_down;
    let tb_up = tb_avg + sum_up;

    (tran, tb_up, tb_down)
}
