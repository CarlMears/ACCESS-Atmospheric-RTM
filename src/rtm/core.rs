//! Low-level RTM functions
//!
//! These are pretty directly re-written from the original Fortran source.

use num_complex::Complex32;
use once_cell::sync::OnceCell;
use smallvec::SmallVec;

const NLINES_O2: usize = 44;

struct OxygenCoefficients {
    f0: [f32; NLINES_O2],
    a1: [f32; NLINES_O2],
    a2: [f32; NLINES_O2],
    a3: [f32; NLINES_O2],
    a4: [f32; NLINES_O2],
    a5: [f32; NLINES_O2],
    a6: [f32; NLINES_O2],
}

/// Oxygen absorption coefficients
static O2_COEF: OnceCell<OxygenCoefficients> = OnceCell::new();

const NLINES_H2O: usize = 15;

struct WaterVaporCoefficients {
    f0: [f32; NLINES_H2O],
    b1: [f32; NLINES_H2O],
    b2: [f32; NLINES_H2O],
    b3: [f32; NLINES_H2O],
    b4: [f32; NLINES_H2O],
    b5: [f32; NLINES_H2O],
    b6: [f32; NLINES_H2O],
}

/// Water vapor absorption coefficients
static H2O_COEF: OnceCell<WaterVaporCoefficients> = OnceCell::new();

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
/// indices from `1..=num_levels` are profile data above the surface.
pub(crate) fn atm_tran(inc: f32, t: &[f32], z: &[f32], tabs: &[f32]) -> (f32, f32, f32) {
    const DELTA: f32 = 0.00035;

    // Differential slant height
    let dsdh = (1.0 + DELTA) / f32::sqrt(inc.to_radians().cos().powi(2) + DELTA * (2.0 + DELTA));

    // Number of levels *not* including the surface
    let num_levels = t.len() - 1;

    let opacity: SmallVec<[f32; 64]> = (1..=num_levels)
        .into_iter()
        .map(|i| -dsdh * 0.5 * (tabs[i - 1] + tabs[i]) * (z[i] - z[i - 1]))
        .collect();
    let t_avg: SmallVec<[f32; 64]> = (1..=num_levels)
        .into_iter()
        .map(|i| 0.5 * (t[i - 1] + t[i]))
        .collect();
    let ems: SmallVec<[f32; 64]> = opacity.iter().map(|opacity| 1.0 - opacity.exp()).collect();

    let (sum_down, _sum_op) =
        (1..=num_levels)
            .into_iter()
            .fold((0., 0.), |(sum_down, sum_op), i| {
                (
                    sum_down + (t_avg[i - 1] - t[1]) * ems[i - 1] * f32::exp(sum_op),
                    sum_op + opacity[i - 1],
                )
            });

    let (sum_up, sum_op) =
        (1..=num_levels)
            .into_iter()
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

/// Modified version of Liebe 1992 oxygen model.
///
/// For a total pressure `p` in hPa, temperature `t` in K, water vapor pressure
/// `pv` in hPa, and frequency `freq` in GHz, compute the oxygen absorption
/// coefficient in dB/km.
///
/// From: Atmospheric 60-GHz Oxygen Spectrum:.. Liebe, Rosenkranz, Hufford,
/// 1992. Modified over the years by Frank Wentz and converted from Fortran to
/// Rust by Richard Lindsley.
pub(crate) fn fdabsoxy_1992_modified(p: f32, t: f32, pv: f32, freq: f32) -> f32 {
    #![allow(clippy::excessive_precision)]
    // Many of the variables are retained from the original Fortran
    let OxygenCoefficients {
        f0,
        a1,
        a2,
        a3,
        a4,
        a5,
        a6,
    } = O2_COEF.get_or_init(|| {
        // All but the last six entries in a4 are 0
        let mut a4 = [0.; NLINES_O2];
        for a4 in a4.iter_mut().skip(NLINES_O2 - 6) {
            *a4 = 0.6;
        }

        let h1 = [
            50.474238, 50.987749, 51.503350, 52.021410, 52.542394, 53.066907, 53.595749, 54.130000,
            54.671159, 55.221367, 55.783802, 56.264775, 56.363389, 56.968206, 57.612484, 58.323877,
            58.446590, 59.164207, 59.590983, 60.306061, 60.434776, 61.150560, 61.800154, 62.411215,
            62.486260, 62.997977, 63.568518, 64.127767, 64.678903, 65.224071, 65.764772, 66.302091,
            66.836830, 67.369598, 67.900867, 68.431005, 68.960311, 118.750343, 368.498350,
            424.763124, 487.249370, 715.393150, 773.839675, 834.145330,
        ];
        let h2 = [
            0.94e-6, 2.46e-6, 6.08e-6, 14.14e-6, 31.02e-6, 64.10e-6, 124.70e-6, 228.00e-6,
            391.80e-6, 631.60e-6, 953.50e-6, 548.90e-6, 1344.00e-6, 1763.00e-6, 2141.00e-6,
            2386.00e-6, 1457.00e-6, 2404.00e-6, 2112.00e-6, 2124.00e-6, 2461.00e-6, 2504.00e-6,
            2298.00e-6, 1933.00e-6, 1517.00e-6, 1503.00e-6, 1087.00e-6, 733.50e-6, 463.50e-6,
            274.80e-6, 153.00e-6, 80.09e-6, 39.46e-6, 18.32e-6, 8.01e-6, 3.30e-6, 1.28e-6,
            945.00e-6, 67.90e-6, 638.00e-6, 235.00e-6, 99.60e-6, 671.00e-6, 180.00e-6,
        ];
        let h3 = [
            9.694, 8.694, 7.744, 6.844, 6.004, 5.224, 4.484, 3.814, 3.194, 2.624, 2.119, 0.015,
            1.660, 1.260, 0.915, 0.626, 0.084, 0.391, 0.212, 0.212, 0.391, 0.626, 0.915, 1.260,
            0.083, 1.665, 2.115, 2.620, 3.195, 3.815, 4.485, 5.225, 6.005, 6.845, 7.745, 8.695,
            9.695, 0.009, 0.049, 0.044, 0.049, 0.145, 0.130, 0.147,
        ];
        let h4 = [
            8.60e-3, 8.70e-3, 8.90e-3, 9.20e-3, 9.40e-3, 9.70e-3, 10.00e-3, 10.20e-3, 10.50e-3,
            10.79e-3, 11.10e-3, 16.46e-3, 11.44e-3, 11.81e-3, 12.21e-3, 12.66e-3, 14.49e-3,
            13.19e-3, 13.60e-3, 13.82e-3, 12.97e-3, 12.48e-3, 12.07e-3, 11.71e-3, 14.68e-3,
            11.39e-3, 11.08e-3, 10.78e-3, 10.50e-3, 10.20e-3, 10.00e-3, 9.70e-3, 9.40e-3, 9.20e-3,
            8.90e-3, 8.70e-3, 8.60e-3, 16.30e-3, 19.20e-3, 19.16e-3, 19.20e-3, 18.10e-3, 18.10e-3,
            18.10e-3,
        ];
        let h5 = [
            0.210, 0.190, 0.171, 0.144, 0.118, 0.114, 0.200, 0.291, 0.325, 0.224, -0.144, 0.339,
            -0.258, -0.362, -0.533, -0.178, 0.650, -0.628, 0.665, -0.613, 0.606, 0.090, 0.496,
            0.313, -0.433, 0.208, 0.094, -0.270, -0.366, -0.326, -0.232, -0.146, -0.147, -0.174,
            -0.198, -0.210, -0.220, -0.031, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];
        let h6 = [
            0.685, 0.680, 0.673, 0.664, 0.653, 0.621, 0.508, 0.375, 0.265, 0.295, 0.613, -0.098,
            0.655, 0.645, 0.606, 0.044, -0.127, 0.231, -0.078, 0.070, -0.282, -0.058, -0.662,
            -0.676, 0.084, -0.668, -0.614, -0.289, -0.259, -0.368, -0.500, -0.609, -0.639, -0.647,
            -0.655, -0.660, -0.665, 0.008, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ];

        let mut a1 = [0.; NLINES_O2];
        a1.iter_mut().zip(&h2).zip(&h1).for_each(|((a1, h2), h1)| {
            *a1 = h2 / h1;
        });
        let f0 = h1;
        let a2 = h3;
        let a3 = h4;
        let a5 = h5.map(|h5| 0.001 * h5);
        let a6 = h6.map(|h6| 0.001 * h6);

        OxygenCoefficients {
            f0,
            a1,
            a2,
            a3,
            a4,
            a5,
            a6,
        }
    });

    let tht = 300.0 / t;
    let pwet = 0.1 * pv;
    let pdry = 0.1 * p - pwet;
    let xterm = 1.0 - tht;

    // Rather than doing one loop over the oxygen lines (as in the original
    // Fortran), it works out better to build some lazy iterators, collect an
    // intermediate result into a stack-local array, and then finally transform
    // and sum that. This must be due to cache locality effects?
    let sum: f64 = {
        let ga = a3
            .iter()
            .zip(a4)
            .map(|(a3, a4)| a3 * (pdry * tht.powf(0.8 - a4) + 1.1 * tht * pwet));

        let delta = a5
            .iter()
            .zip(a6)
            .map(|(a5, a6)| (a5 + a6 * tht) * p * tht.powf(0.8));

        let mut ff = [0.; NLINES_O2];
        for (((ff, f0), ga), delta) in ff.iter_mut().zip(f0).zip(ga).zip(delta) {
            let rnuneg = f0 - freq;
            let rnupos = f0 + freq;
            let ga_sq = ga.powi(2);

            *ff = (ga - rnuneg * delta) / (ga_sq + rnuneg.powi(2))
                + (ga - rnupos * delta) / (ga_sq + rnupos.powi(2));
        }

        ff.iter()
            .zip(a1)
            .zip(a2)
            .map(|((ff, a1), a2)| f64::from(ff * a1 * f32::exp(a2 * xterm)))
            .sum()
    };
    let sum = sum.max(0.0);

    // add nonresonant contribution ("modification 1")
    let ga = 5.6e-3 * (pdry + 1.1 * pwet) * tht.powf(1.5);

    let zterm = ga * (1. + (freq / ga).powi(2));
    let apterm = 1.4e-10 * (1.0 - 1.2e-5 * freq.powf(1.5)) * pdry * tht.powf(1.5);
    let apterm = apterm.max(0.);
    let sftot = (f64::from(pdry * freq * tht.powi(2))
        * (f64::from(tht) * sum + 6.14e-4 / f64::from(zterm) + f64::from(apterm)))
        as f32;

    let gamoxy = 0.1820 * freq * sftot;
    if freq > 37. {
        gamoxy + 0.1820 * 26.0e-10 * pdry.powi(2) * tht.powi(3) * (freq - 37.).powf(1.8)
    } else {
        gamoxy
    }
}

/// Modified version of Rosenkranz water vapor model.
///
/// For a total pressure `p` in hPa, temperature `t` in K, water vapor pressure
/// `pv` in hPa, and frequency `freq` in GHz, compute the water vapor absorption
/// coefficient in dB/km.
///
/// From: P.W. Rosenkranz, Radio Science v.33, pp.919-928 (1998). Modified by
/// Frank Wentz over the years and converted from Fortran to Rust by Richard
/// Lindsley.
pub(crate) fn abh2o_rk_modified(p: f32, t: f32, pv: f32, freq: f32) -> f32 {
    // Many of the variables are retained from the original Fortran
    let WaterVaporCoefficients {
        f0,
        b1,
        b2,
        b3,
        b4,
        b5,
        b6,
    } = H2O_COEF.get_or_init(|| {
        #[allow(clippy::excessive_precision)]
        // line frequencies
        let f0 = [
            22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508, 443.0183, 448.0011,
            470.8890, 474.6891, 488.4911, 556.9360, 620.7008, 752.0332, 916.1712,
        ];

        // line intensities at 300 K
        let b1: [f32; NLINES_H2O] = [
            0.1310e-13, 0.2273e-11, 0.8036e-13, 0.2694e-11, 0.2438e-10, 0.2179e-11, 0.4624e-12,
            0.2562e-10, 0.8369e-12, 0.3263e-11, 0.6659e-12, 0.1531e-08, 0.1707e-10, 0.1011e-08,
            0.4227e-10,
        ];

        // t coeff. of intensities
        let b2 = [
            2.144, 0.668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405, 3.597, 2.379, 2.852, 0.159,
            2.391, 0.396, 1.441,
        ];
        // air-broadened width parameters at 300 K
        let mut b3 = [
            0.0281, 0.0281, 0.023, 0.0278, 0.0287, 0.021, 0.0186, 0.0263, 0.0215, 0.0236, 0.026,
            0.0321, 0.0244, 0.0306, 0.0267,
        ];
        // self-broadened width parameters at 300 K
        let b5 = [
            0.1349, 0.1491, 0.108, 0.135, 0.1541, 0.090, 0.0788, 0.1275, 0.0983, 0.1095, 0.1313,
            0.1320, 0.1140, 0.1253, 0.1275,
        ];
        // t-exponent of air-broadening
        let b4 = [
            0.69, 0.64, 0.67, 0.68, 0.54, 0.63, 0.60, 0.66, 0.66, 0.65, 0.69, 0.69, 0.71, 0.68,
            0.70,
        ];
        // t-exponent of self-broadening
        let b6 = [
            0.61, 0.85, 0.54, 0.74, 0.89, 0.52, 0.50, 0.67, 0.65, 0.64, 0.72, 1.0, 0.68, 0.84, 0.78,
        ];

        let mut b1_modified = [0.; NLINES_H2O];
        #[allow(clippy::excessive_precision)]
        for ((b1_new, b1), &f0) in b1_modified.iter_mut().zip(&b1).zip(&f0) {
            *b1_new = 1.8281089E+14 * b1 / f32::powi(f0, 2);
        }

        // convert b5 to Leibe notation
        let mut b5_modified = [0.; NLINES_H2O];
        for ((b5_new, b5), b3) in b5_modified.iter_mut().zip(&b5).zip(&b3) {
            *b5_new = b5 / b3;
        }

        // Modification 1
        b3[0] /= 1.040;

        WaterVaporCoefficients {
            f0,
            b1: b1_modified,
            b2,
            b3,
            b4,
            b5: b5_modified,
            b6,
        }
    });

    if pv <= 0. {
        return 0.;
    }

    let pwet = 0.1 * pv;
    let pdry = 0.1 * p - pwet;
    let tht = 300. / t;
    let xterm = 1. - tht;
    let freq_sq = freq.powi(2);

    let sum: f64 = (0..NLINES_H2O)
        .into_iter()
        .map(|i| {
            let f0sq = f0[i].powi(2);
            let ga = b3[i] * (pdry * tht.powf(b4[i]) + b5[i] * pwet * tht.powf(b6[i]));
            let ga_sq = ga.powi(2);
            let s = b1[i] * f32::exp(b2[i] * xterm);
            let rnuneg = f0[i] - freq;
            let rnupos = f0[i] + freq;

            // use clough's definition of local line contribution
            let base = ga / (562_500. + ga_sq);

            if i != 0 {
                let mut sum = 0.;
                if rnuneg.abs() < 750. {
                    sum += f64::from(s * (ga / (ga_sq + rnuneg.powi(2)) - base));
                }
                if rnupos.abs() <= 750. {
                    sum += f64::from(s * (ga / (ga_sq + rnupos.powi(2)) - base));
                }
                sum
            } else {
                // modification 2
                let chi = if freq < 19. {
                    let u = f32::clamp((freq - 19.).abs() / 16.5, 0., 1.);
                    0.07 * ga + 0.93 * ga * u.powi(2) * (3. - 2. * u)
                } else {
                    0.07 * ga
                };

                let chi_sq = chi.powi(2);
                f64::from(
                    s * 2. * ((ga - chi) * freq_sq + (ga + chi) * (f0sq + ga_sq - chi_sq))
                        / ((freq_sq - f0sq - ga_sq + chi_sq).powi(2) + 4. * freq_sq * ga_sq),
                )
            }
        })
        .sum();
    let sum = sum.max(0.);

    let ffac = if freq < 90. {
        1. + 0.1 * ((90. - freq) / 90.).powf(1.4)
    } else {
        1.
    };

    // modification 3
    let sftot = pwet
        * freq
        * tht.powf(3.5)
        * (sum
            + f64::from(ffac * 1.1 * 1.2957246e-6 * pdry / tht.sqrt())
            + f64::from(0.348 * (freq.powf(0.15)) * 4.2952193e-5 * pwet * tht.powi(4)))
            as f32;

    0.1820 * freq * sftot
}

/// Liquid cloud water absorption coefficient.
///
/// For a frequency `freq` in GHz, a temperature `t` in K, and a liquid cloud
/// water density `rhol` in g/m^3, compute the cloud water absorption
/// coefficient in Np/km.
pub(crate) fn fdcldabs(freq: f32, t: f32, rhol: f32) -> f32 {
    const C: f32 = 29.979;
    use std::f32::consts::PI;

    // Convert g/m^3 to g/cm^3
    let rhol0 = 1.0e-6 * rhol;

    let permit = meissner(freq, t, 0.0);
    let wavlen = C / freq;
    // Np/cm
    let al = (6.0 * PI * rhol0 / wavlen) * ((1.0 - permit) / (2.0 + permit)).im;

    // Convert to Np/km
    al * 1.0e5
}

/// Compute the complex dielectric constant of water.
///
/// For a frequency `freq` in GHz, SST `t` in K, salinity `s` in parts per
/// thousand, compute the complex dielectric constant of water.
///
/// The ranges for the inputs are:
///
/// - `freq`: 1 to 400 GHz
/// - `t`: -25 °C to 40 °C (or 248.16 K to 313.16 K) for pure water; -2 °C to 34
///   °C (or 271.16 K to 307.16 K) for saline water
/// - `s`: 0 to 40 ppt
///
/// From Thomas Meissner, February 2002 and October 2004.
///
/// The imaginary part is negative to be consistent with "wentz1" convention.
fn meissner(freq: f32, t: f32, s: f32) -> Complex32 {
    #![allow(clippy::excessive_precision)]
    const F0: f32 = 17.97510;

    // Convert from K to °C
    let sst = t - 273.15;
    let (e0s, e1s, e2s, n1s, n2s, sig) = dielectric_meissner_wentz(sst, s);

    // Debye law (2 relaxation wavelengths)
    let eps = (e0s - e1s) / Complex32::new(1.0, -(freq / n1s))
        + (e1s - e2s) / Complex32::new(1.0, -(freq / n2s))
        + e2s
        + Complex32::new(0., sig * F0 / freq);

    eps.conj()
}

/// Complex dielectric constant.
///
/// For an input SST, `sst` in °C and salinity `s` in parts-per-thousand,
/// compute the complex dielectric constant. Returns the tuple `(e0s, e1s, e2s,
/// n1s, n2s, sig)`.
///
/// The allowed values of the inputs are:
///
/// - `sst`: from -25°C to 40°C for pure water; -2°C to 34°C for saline water
/// - `s`: from 0 to 40 ppt
///
///
/// References: T. Meissner  and F. Wentz, IEEE TGARS, 42(9), 2004, 1836-1849.
fn dielectric_meissner_wentz(sst: f32, s: f32) -> (f32, f32, f32, f32, f32, f32) {
    #![allow(clippy::excessive_precision)]
    const X: [f32; 11] = [
        5.7230e+00,
        2.2379e-02,
        -7.1237e-04,
        5.0478e+00,
        -7.0315e-02,
        6.0059e-04,
        3.6143e+00,
        2.8841e-02,
        1.3652e-01,
        1.4825e-03,
        2.4166e-04,
    ];

    const Z: [f32; 13] = [
        -3.56417e-03,
        4.74868e-06,
        1.15574e-05,
        2.39357e-03,
        -3.13530e-05,
        2.52477e-07,
        -6.28908e-03,
        1.76032e-04,
        -9.22144e-05,
        -1.99723e-02,
        1.81176e-04,
        -2.04265e-03,
        1.57883e-04,
    ];

    const A0_COEF: [f32; 3] = [-0.33330E-02, 4.74868e-06, 0.0e0];
    const B1_COEF: [f32; 5] = [
        0.23232E-02,
        -0.79208E-04,
        0.36764E-05,
        -0.35594E-06,
        0.89795E-08,
    ];

    // protects against n1 and n2 going zero for very cold water
    let sst = sst.max(-30.16);
    let sst2 = sst.powi(2);
    let sst3 = sst.powi(3);
    let sst4 = sst.powi(4);

    let s2 = s.powi(2);

    // Pure water. e0 is from Stogryn et al.
    let e0 = (3.70886e4 - 8.2168e1 * sst) / (4.21854e2 + sst);
    let e1 = X[0] + X[1] * sst + X[2] * sst2;
    let n1 = (45.0 + sst) / (X[3] + X[4] * sst + X[5] * sst2);
    let e2 = X[6] + X[7] * sst;
    let n2 = (45.0 + sst) / (X[8] + X[9] * sst + X[10] * sst2);

    // Saline water. Conductivity [s/m] taken from Stogryn et al.
    let sig35 =
        2.903602 + 8.60700e-2 * sst + 4.738817e-4 * sst2 - 2.9910e-6 * sst3 + 4.3047e-9 * sst4;
    let r15 = s * (37.5109 + 5.45216 * s + 1.4409e-2 * s2) / (1004.75 + 182.283 * s + s2);

    let alpha0 = (6.9431 + 3.2841 * s - 9.9486e-2 * s2) / (84.850 + 69.024 * s + s2);
    let alpha1 = 49.843 - 0.2276 * s + 0.198e-2 * s2;
    let rtr15 = 1.0 + (sst - 15.0) * alpha0 / (alpha1 + sst);

    let sig = sig35 * r15 * rtr15;

    // permittivity
    let a0 = f32::exp(A0_COEF[0] * s + A0_COEF[1] * s2 + A0_COEF[2] * s * sst);
    let e0s = a0 * e0;

    let b1 = if sst <= 30. {
        1.0 + s
            * (B1_COEF[0]
                + B1_COEF[1] * sst
                + B1_COEF[2] * sst2
                + B1_COEF[3] * sst3
                + B1_COEF[4] * sst4)
    } else {
        1.0 + s * (9.1873715e-04 + 1.5012396e-04 * (sst - 30.))
    };
    let n1s = n1 * b1;

    let a1 = f32::exp(Z[6] * s + Z[7] * s2 + Z[8] * s * sst);
    let e1s = e1 * a1;

    let b2 = 1.0 + s * (Z[9] + 0.5 * Z[10] * (sst + 30.));
    let n2s = n2 * b2;

    let a2 = 1.0 + s * (Z[11] + Z[12] * sst);
    let e2s = e2 * a2;

    (e0s, e1s, e2s, n1s, n2s, sig)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// Check some values for the oxygen absorption coefficient. These values
    /// are from the Fortran version.
    #[test]
    fn oxygen() {
        let inputs_and_outputs = [
            [
                1000.00000,
                273.149994,
                1.00000000,
                2.00000000,
                7.90619664E-03,
            ],
            [
                1000.00000,
                273.149994,
                1.00000000,
                10.0000000,
                9.75634996E-03,
            ],
            [
                1000.00000,
                273.149994,
                1.00000000,
                30.0000000,
                2.46610157E-02,
            ],
            [
                1000.00000,
                273.149994,
                0.100000001,
                2.00000000,
                7.91274104E-03,
            ],
            [
                1000.00000,
                273.149994,
                0.100000001,
                10.0000000,
                9.76428948E-03,
            ],
            [
                1000.00000,
                273.149994,
                0.100000001,
                30.0000000,
                2.46806331E-02,
            ],
            [
                1000.00000,
                273.149994,
                1.00000005E-03,
                2.00000000,
                7.91346002E-03,
            ],
            [
                1000.00000,
                273.149994,
                1.00000005E-03,
                10.0000000,
                9.76515841E-03,
            ],
            [
                1000.00000,
                273.149994,
                1.00000005E-03,
                30.0000000,
                2.46827919E-02,
            ],
            [
                1000.00000,
                290.000000,
                1.00000000,
                2.00000000,
                6.51325937E-03,
            ],
            [
                1000.00000,
                290.000000,
                1.00000000,
                10.0000000,
                7.95020070E-03,
            ],
            [
                1000.00000,
                290.000000,
                1.00000000,
                30.0000000,
                2.04622969E-02,
            ],
            [
                1000.00000,
                290.000000,
                0.100000001,
                2.00000000,
                6.51863590E-03,
            ],
            [
                1000.00000,
                290.000000,
                0.100000001,
                10.0000000,
                7.95667619E-03,
            ],
            [
                1000.00000,
                290.000000,
                0.100000001,
                30.0000000,
                2.04787478E-02,
            ],
            [
                1000.00000,
                290.000000,
                1.00000005E-03,
                2.00000000,
                6.51922543E-03,
            ],
            [
                1000.00000,
                290.000000,
                1.00000005E-03,
                10.0000000,
                7.95738958E-03,
            ],
            [
                1000.00000,
                290.000000,
                1.00000005E-03,
                30.0000000,
                2.04805583E-02,
            ],
            [
                1000.00000,
                300.000000,
                1.00000000,
                2.00000000,
                5.82994567E-03,
            ],
            [
                1000.00000,
                300.000000,
                1.00000000,
                10.0000000,
                7.07963528E-03,
            ],
            [
                1000.00000,
                300.000000,
                1.00000000,
                30.0000000,
                1.84053257E-02,
            ],
            [
                1000.00000,
                300.000000,
                0.100000001,
                2.00000000,
                5.83474990E-03,
            ],
            [
                1000.00000,
                300.000000,
                0.100000001,
                10.0000000,
                7.08540762E-03,
            ],
            [
                1000.00000,
                300.000000,
                0.100000001,
                30.0000000,
                1.84202138E-02,
            ],
            [
                1000.00000,
                300.000000,
                1.00000005E-03,
                2.00000000,
                5.83527796E-03,
            ],
            [
                1000.00000,
                300.000000,
                1.00000005E-03,
                10.0000000,
                7.08604185E-03,
            ],
            [
                1000.00000,
                300.000000,
                1.00000005E-03,
                30.0000000,
                1.84218492E-02,
            ],
            [
                500.000000,
                273.149994,
                1.00000000,
                2.00000000,
                2.12392258E-03,
            ],
            [
                500.000000,
                273.149994,
                1.00000000,
                10.0000000,
                2.44359416E-03,
            ],
            [
                500.000000,
                273.149994,
                1.00000000,
                30.0000000,
                6.15879428E-03,
            ],
            [
                500.000000,
                273.149994,
                0.100000001,
                2.00000000,
                2.12739035E-03,
            ],
            [
                500.000000,
                273.149994,
                0.100000001,
                10.0000000,
                2.44757137E-03,
            ],
            [
                500.000000,
                273.149994,
                0.100000001,
                30.0000000,
                6.16860390E-03,
            ],
            [
                500.000000,
                273.149994,
                1.00000005E-03,
                2.00000000,
                2.12777173E-03,
            ],
            [
                500.000000,
                273.149994,
                1.00000005E-03,
                10.0000000,
                2.44800886E-03,
            ],
            [
                500.000000,
                273.149994,
                1.00000005E-03,
                30.0000000,
                6.16968237E-03,
            ],
            [
                500.000000,
                290.000000,
                1.00000000,
                2.00000000,
                1.72994297E-03,
            ],
            [
                500.000000,
                290.000000,
                1.00000000,
                10.0000000,
                1.99030270E-03,
            ],
            [
                500.000000,
                290.000000,
                1.00000000,
                30.0000000,
                5.11018187E-03,
            ],
            [
                500.000000,
                290.000000,
                0.100000001,
                2.00000000,
                1.73276500E-03,
            ],
            [
                500.000000,
                290.000000,
                0.100000001,
                10.0000000,
                1.99354719E-03,
            ],
            [
                500.000000,
                290.000000,
                0.100000001,
                30.0000000,
                5.11840824E-03,
            ],
            [
                500.000000,
                290.000000,
                1.00000005E-03,
                2.00000000,
                1.73307571E-03,
            ],
            [
                500.000000,
                290.000000,
                1.00000005E-03,
                10.0000000,
                1.99390342E-03,
            ],
            [
                500.000000,
                290.000000,
                1.00000005E-03,
                30.0000000,
                5.11931209E-03,
            ],
            [
                500.000000,
                300.000000,
                1.00000000,
                2.00000000,
                1.53969473E-03,
            ],
            [
                500.000000,
                300.000000,
                1.00000000,
                10.0000000,
                1.77195808E-03,
            ],
            [
                500.000000,
                300.000000,
                1.00000000,
                30.0000000,
                4.59646620E-03,
            ],
            [
                500.000000,
                300.000000,
                0.100000001,
                2.00000000,
                1.54220534E-03,
            ],
            [
                500.000000,
                300.000000,
                0.100000001,
                10.0000000,
                1.77484925E-03,
            ],
            [
                500.000000,
                300.000000,
                0.100000001,
                30.0000000,
                4.60390979E-03,
            ],
            [
                500.000000,
                300.000000,
                1.00000005E-03,
                2.00000000,
                1.54248148E-03,
            ],
            [
                500.000000,
                300.000000,
                1.00000005E-03,
                10.0000000,
                1.77516695E-03,
            ],
            [
                500.000000,
                300.000000,
                1.00000005E-03,
                30.0000000,
                4.60472750E-03,
            ],
            [
                100.000000,
                273.149994,
                1.00000000,
                2.00000000,
                8.64337417E-05,
            ],
            [
                100.000000,
                273.149994,
                1.00000000,
                10.0000000,
                9.71217014E-05,
            ],
            [
                100.000000,
                273.149994,
                1.00000000,
                30.0000000,
                2.44592986E-04,
            ],
            [
                100.000000,
                273.149994,
                0.100000001,
                2.00000000,
                8.71413431E-05,
            ],
            [
                100.000000,
                273.149994,
                0.100000001,
                10.0000000,
                9.79184406E-05,
            ],
            [
                100.000000,
                273.149994,
                0.100000001,
                30.0000000,
                2.46556592E-04,
            ],
            [
                100.000000,
                273.149994,
                1.00000005E-03,
                2.00000000,
                8.72190940E-05,
            ],
            [
                100.000000,
                273.149994,
                1.00000005E-03,
                10.0000000,
                9.80059995E-05,
            ],
            [
                100.000000,
                273.149994,
                1.00000005E-03,
                30.0000000,
                2.46772397E-04,
            ],
            [
                100.000000,
                290.000000,
                1.00000000,
                2.00000000,
                7.01206882E-05,
            ],
            [
                100.000000,
                290.000000,
                1.00000000,
                10.0000000,
                7.90928752E-05,
            ],
            [
                100.000000,
                290.000000,
                1.00000000,
                30.0000000,
                2.02932788E-04,
            ],
            [
                100.000000,
                290.000000,
                0.100000001,
                2.00000000,
                7.06947685E-05,
            ],
            [
                100.000000,
                290.000000,
                0.100000001,
                10.0000000,
                7.97426983E-05,
            ],
            [
                100.000000,
                290.000000,
                0.100000001,
                30.0000000,
                2.04579439E-04,
            ],
            [
                100.000000,
                290.000000,
                1.00000005E-03,
                2.00000000,
                7.07578365E-05,
            ],
            [
                100.000000,
                290.000000,
                1.00000005E-03,
                10.0000000,
                7.98140973E-05,
            ],
            [
                100.000000,
                290.000000,
                1.00000005E-03,
                30.0000000,
                2.04760363E-04,
            ],
            [
                100.000000,
                300.000000,
                1.00000000,
                2.00000000,
                6.22867592E-05,
            ],
            [
                100.000000,
                300.000000,
                1.00000000,
                10.0000000,
                7.04105114E-05,
            ],
            [
                100.000000,
                300.000000,
                1.00000000,
                30.0000000,
                1.82524440E-04,
            ],
            [
                100.000000,
                300.000000,
                0.100000001,
                2.00000000,
                6.27967020E-05,
            ],
            [
                100.000000,
                300.000000,
                0.100000001,
                10.0000000,
                7.09894957E-05,
            ],
            [
                100.000000,
                300.000000,
                0.100000001,
                30.0000000,
                1.84014571E-04,
            ],
            [
                100.000000,
                300.000000,
                1.00000005E-03,
                2.00000000,
                6.28527414E-05,
            ],
            [
                100.000000,
                300.000000,
                1.00000005E-03,
                10.0000000,
                7.10531240E-05,
            ],
            [
                100.000000,
                300.000000,
                1.00000005E-03,
                30.0000000,
                1.84178294E-04,
            ],
        ];

        for [p, t, pv, freq, expected_output] in inputs_and_outputs {
            assert_relative_eq!(fdabsoxy_1992_modified(p, t, pv, freq), expected_output);
        }
    }

    /// Check some values for the water vapor absorption coefficient. These values
    /// are from the Fortran version.
    #[test]
    fn water_vapor() {
        let inputs_and_outputs = [
            [
                1000.00000,
                273.149994,
                1.00000000,
                2.00000000,
                2.59790759E-05,
            ],
            [
                1000.00000,
                273.149994,
                1.00000000,
                10.0000000,
                7.06406950E-04,
            ],
            [
                1000.00000,
                273.149994,
                1.00000000,
                30.0000000,
                7.41568720E-03,
            ],
            [
                1000.00000,
                273.149994,
                0.100000001,
                2.00000000,
                2.57373767E-06,
            ],
            [
                1000.00000,
                273.149994,
                0.100000001,
                10.0000000,
                6.98733420E-05,
            ],
            [
                1000.00000,
                273.149994,
                0.100000001,
                30.0000000,
                7.33381894E-04,
            ],
            [
                1000.00000,
                273.149994,
                1.00000005E-03,
                2.00000000,
                2.57107917E-08,
            ],
            [
                1000.00000,
                273.149994,
                1.00000005E-03,
                10.0000000,
                6.97889163E-07,
            ],
            [
                1000.00000,
                273.149994,
                1.00000005E-03,
                30.0000000,
                7.32481112E-06,
            ],
            [
                1000.00000,
                290.000000,
                1.00000000,
                2.00000000,
                2.21492974E-05,
            ],
            [
                1000.00000,
                290.000000,
                1.00000000,
                10.0000000,
                6.05081732E-04,
            ],
            [
                1000.00000,
                290.000000,
                1.00000000,
                30.0000000,
                6.41573453E-03,
            ],
            [
                1000.00000,
                290.000000,
                0.100000001,
                2.00000000,
                2.19888739E-06,
            ],
            [
                1000.00000,
                290.000000,
                0.100000001,
                10.0000000,
                5.99988671E-05,
            ],
            [
                1000.00000,
                290.000000,
                0.100000001,
                30.0000000,
                6.36112352E-04,
            ],
            [
                1000.00000,
                290.000000,
                1.00000005E-03,
                2.00000000,
                2.19712266E-08,
            ],
            [
                1000.00000,
                290.000000,
                1.00000005E-03,
                10.0000000,
                5.99428290E-07,
            ],
            [
                1000.00000,
                290.000000,
                1.00000005E-03,
                30.0000000,
                6.35511469E-06,
            ],
            [
                1000.00000,
                300.000000,
                1.00000000,
                2.00000000,
                2.02244992E-05,
            ],
            [
                1000.00000,
                300.000000,
                1.00000000,
                10.0000000,
                5.53888734E-04,
            ],
            [
                1000.00000,
                300.000000,
                1.00000000,
                30.0000000,
                5.90435695E-03,
            ],
            [
                1000.00000,
                300.000000,
                0.100000001,
                2.00000000,
                2.00966224E-06,
            ],
            [
                1000.00000,
                300.000000,
                0.100000001,
                10.0000000,
                5.49828328E-05,
            ],
            [
                1000.00000,
                300.000000,
                0.100000001,
                30.0000000,
                5.86065173E-04,
            ],
            [
                1000.00000,
                300.000000,
                1.00000005E-03,
                2.00000000,
                2.00825561E-08,
            ],
            [
                1000.00000,
                300.000000,
                1.00000005E-03,
                10.0000000,
                5.49381696E-07,
            ],
            [
                1000.00000,
                300.000000,
                1.00000005E-03,
                30.0000000,
                5.85584257E-06,
            ],
            [
                500.000000,
                273.149994,
                1.00000000,
                2.00000000,
                1.31258521E-05,
            ],
            [
                500.000000,
                273.149994,
                1.00000000,
                10.0000000,
                3.61695827E-04,
            ],
            [
                500.000000,
                273.149994,
                1.00000000,
                30.0000000,
                3.94016085E-03,
            ],
            [
                500.000000,
                273.149994,
                0.100000001,
                2.00000000,
                1.28841123E-06,
            ],
            [
                500.000000,
                273.149994,
                0.100000001,
                10.0000000,
                3.53940814E-05,
            ],
            [
                500.000000,
                273.149994,
                0.100000001,
                30.0000000,
                3.85494815E-04,
            ],
            [
                500.000000,
                273.149994,
                1.00000005E-03,
                2.00000000,
                1.28575204E-08,
            ],
            [
                500.000000,
                273.149994,
                1.00000005E-03,
                10.0000000,
                3.53087728E-07,
            ],
            [
                500.000000,
                273.149994,
                1.00000005E-03,
                30.0000000,
                3.84557325E-06,
            ],
            [
                500.000000,
                290.000000,
                1.00000000,
                2.00000000,
                1.11654444E-05,
            ],
            [
                500.000000,
                290.000000,
                1.00000000,
                10.0000000,
                3.08859831E-04,
            ],
            [
                500.000000,
                290.000000,
                1.00000000,
                30.0000000,
                3.39363888E-03,
            ],
            [
                500.000000,
                290.000000,
                0.100000001,
                2.00000000,
                1.10049859E-06,
            ],
            [
                500.000000,
                290.000000,
                0.100000001,
                10.0000000,
                3.03698980E-05,
            ],
            [
                500.000000,
                290.000000,
                0.100000001,
                30.0000000,
                3.33620585E-04,
            ],
            [
                500.000000,
                290.000000,
                1.00000005E-03,
                2.00000000,
                1.09873346E-08,
            ],
            [
                500.000000,
                290.000000,
                1.00000005E-03,
                10.0000000,
                3.03131259E-07,
            ],
            [
                500.000000,
                290.000000,
                1.00000005E-03,
                30.0000000,
                3.32988748E-06,
            ],
            [
                500.000000,
                300.000000,
                1.00000000,
                2.00000000,
                1.01847809E-05,
            ],
            [
                500.000000,
                300.000000,
                1.00000000,
                10.0000000,
                2.82317982E-04,
            ],
            [
                500.000000,
                300.000000,
                1.00000000,
                30.0000000,
                3.11595527E-03,
            ],
            [
                500.000000,
                300.000000,
                0.100000001,
                2.00000000,
                1.00568730E-06,
            ],
            [
                500.000000,
                300.000000,
                0.100000001,
                10.0000000,
                2.78196840E-05,
            ],
            [
                500.000000,
                300.000000,
                0.100000001,
                30.0000000,
                3.06969858E-04,
            ],
            [
                500.000000,
                300.000000,
                1.00000005E-03,
                2.00000000,
                1.00428030E-08,
            ],
            [
                500.000000,
                300.000000,
                1.00000005E-03,
                10.0000000,
                2.77743482E-07,
            ],
            [
                500.000000,
                300.000000,
                1.00000005E-03,
                30.0000000,
                3.06460925E-06,
            ],
            [
                100.000000,
                273.149994,
                1.00000000,
                2.00000000,
                2.84019279E-06,
            ],
            [
                100.000000,
                273.149994,
                1.00000000,
                10.0000000,
                7.95400047E-05,
            ],
            [
                100.000000,
                273.149994,
                1.00000000,
                30.0000000,
                8.78573628E-04,
            ],
            [
                100.000000,
                273.149994,
                0.100000001,
                2.00000000,
                2.59843944E-07,
            ],
            [
                100.000000,
                273.149994,
                0.100000001,
                10.0000000,
                7.17569219E-06,
            ],
            [
                100.000000,
                273.149994,
                0.100000001,
                30.0000000,
                7.92029823E-05,
            ],
            [
                100.000000,
                273.149994,
                1.00000005E-03,
                2.00000000,
                2.57184651E-09,
            ],
            [
                100.000000,
                273.149994,
                1.00000005E-03,
                10.0000000,
                7.09007608E-08,
            ],
            [
                100.000000,
                273.149994,
                1.00000005E-03,
                30.0000000,
                7.82509801E-07,
            ],
            [
                100.000000,
                290.000000,
                1.00000000,
                2.00000000,
                2.37583617E-06,
            ],
            [
                100.000000,
                290.000000,
                1.00000000,
                10.0000000,
                6.66124397E-05,
            ],
            [
                100.000000,
                290.000000,
                1.00000000,
                30.0000000,
                7.41978583E-04,
            ],
            [
                100.000000,
                290.000000,
                0.100000001,
                2.00000000,
                2.21536595E-07,
            ],
            [
                100.000000,
                290.000000,
                0.100000001,
                10.0000000,
                6.14283135E-06,
            ],
            [
                100.000000,
                290.000000,
                0.100000001,
                30.0000000,
                6.83440303E-05,
            ],
            [
                100.000000,
                290.000000,
                1.00000005E-03,
                2.00000000,
                2.19771423E-09,
            ],
            [
                100.000000,
                290.000000,
                1.00000005E-03,
                10.0000000,
                6.08580493E-08,
            ],
            [
                100.000000,
                290.000000,
                1.00000005E-03,
                30.0000000,
                6.77000855E-07,
            ],
            [
                100.000000,
                300.000000,
                1.00000000,
                2.00000000,
                2.15075602E-06,
            ],
            [
                100.000000,
                300.000000,
                1.00000000,
                10.0000000,
                6.03525659E-05,
            ],
            [
                100.000000,
                300.000000,
                1.00000000,
                30.0000000,
                6.75204385E-04,
            ],
            [
                100.000000,
                300.000000,
                0.100000001,
                2.00000000,
                2.02283843E-07,
            ],
            [
                100.000000,
                300.000000,
                0.100000001,
                10.0000000,
                5.62106197E-06,
            ],
            [
                100.000000,
                300.000000,
                0.100000001,
                30.0000000,
                6.27955960E-05,
            ],
            [
                100.000000,
                300.000000,
                1.00000005E-03,
                2.00000000,
                2.00876737E-09,
            ],
            [
                100.000000,
                300.000000,
                1.00000005E-03,
                10.0000000,
                5.57549953E-08,
            ],
            [
                100.000000,
                300.000000,
                1.00000005E-03,
                30.0000000,
                6.22758421E-07,
            ],
        ];
        for [p, t, pv, freq, expected_output] in inputs_and_outputs {
            assert_relative_eq!(abh2o_rk_modified(p, t, pv, freq), expected_output);
        }
    }

    /// Check some values for the liquid cloud absorption coefficient. These
    /// values are from the Fortran version.
    #[test]
    fn liquid_cloud() {
        let inputs_and_outputs = [
            [2.00000000, 273.149994, 1.00000000, 8.60669592E-04],
            [2.00000000, 273.149994, 0.100000001, 8.60669606E-05],
            [2.00000000, 273.149994, 1.00000005E-03, 8.60669616E-07],
            [2.00000000, 290.000000, 1.00000000, 5.34837891E-04],
            [2.00000000, 290.000000, 0.100000001, 5.34837891E-05],
            [2.00000000, 290.000000, 1.00000005E-03, 5.34837909E-07],
            [2.00000000, 300.000000, 1.00000000, 4.27581428E-04],
            [2.00000000, 300.000000, 0.100000001, 4.27581472E-05],
            [2.00000000, 300.000000, 1.00000005E-03, 4.27581455E-07],
            [10.0000000, 273.149994, 1.00000000, 2.13160031E-02],
            [10.0000000, 273.149994, 0.100000001, 2.13160040E-03],
            [10.0000000, 273.149994, 1.00000005E-03, 2.13160038E-05],
            [10.0000000, 290.000000, 1.00000000, 1.33182406E-02],
            [10.0000000, 290.000000, 0.100000001, 1.33182411E-03],
            [10.0000000, 290.000000, 1.00000005E-03, 1.33182421E-05],
            [10.0000000, 300.000000, 1.00000000, 1.06627224E-02],
            [10.0000000, 300.000000, 0.100000001, 1.06627226E-03],
            [10.0000000, 300.000000, 1.00000005E-03, 1.06627231E-05],
            [30.0000000, 273.149994, 1.00000000, 0.178156421],
            [30.0000000, 273.149994, 0.100000001, 1.78156421E-02],
            [30.0000000, 273.149994, 1.00000005E-03, 1.78156421E-04],
            [30.0000000, 290.000000, 1.00000000, 0.116090909],
            [30.0000000, 290.000000, 0.100000001, 1.16090896E-02],
            [30.0000000, 290.000000, 1.00000005E-03, 1.16090901E-04],
            [30.0000000, 300.000000, 1.00000000, 9.40238163E-02],
            [30.0000000, 300.000000, 0.100000001, 9.40238219E-03],
            [30.0000000, 300.000000, 1.00000005E-03, 9.40238242E-05],
        ];

        for [freq, t, rhol, expected_output] in inputs_and_outputs {
            assert_relative_eq!(fdcldabs(freq, t, rhol), expected_output);
        }
    }

    /// Check some values for the water dielectric value. These
    /// values are from the Fortran version.
    #[test]
    fn water_dielectric() {
        let inputs_and_outputs = [
            (
                2.00000000,
                273.149994,
                0.00000000,
                Complex32::new(83.9792633, -17.5693436),
            ),
            (
                2.00000000,
                273.149994,
                15.0000000,
                Complex32::new(80.2050705, -28.2318840),
            ),
            (
                2.00000000,
                273.149994,
                30.0000000,
                Complex32::new(76.7605515, -37.6490784),
            ),
            (
                2.00000000,
                290.000000,
                0.00000000,
                Complex32::new(80.1225433, -9.69443417),
            ),
            (
                2.00000000,
                290.000000,
                15.0000000,
                Complex32::new(76.3306427, -27.7735748),
            ),
            (
                2.00000000,
                290.000000,
                30.0000000,
                Complex32::new(72.8786850, -43.5699654),
            ),
            (
                2.00000000,
                300.000000,
                0.00000000,
                Complex32::new(77.0278702, -7.13627577),
            ),
            (
                2.00000000,
                300.000000,
                15.0000000,
                Complex32::new(73.3600235, -29.7685146),
            ),
            (
                2.00000000,
                300.000000,
                30.0000000,
                Complex32::new(70.0197601, -49.4610291),
            ),
            (
                10.0000000,
                273.149994,
                0.00000000,
                Complex32::new(42.1181870, -40.8917809),
            ),
            (
                10.0000000,
                273.149994,
                15.0000000,
                Complex32::new(41.4164009, -41.5028992),
            ),
            (
                10.0000000,
                273.149994,
                30.0000000,
                Complex32::new(40.9500046, -41.8011780),
            ),
            (
                10.0000000,
                290.000000,
                0.00000000,
                Complex32::new(58.8581696, -34.6061859),
            ),
            (
                10.0000000,
                290.000000,
                15.0000000,
                Complex32::new(56.4688148, -36.5488663),
            ),
            (
                10.0000000,
                290.000000,
                30.0000000,
                Complex32::new(54.3929825, -38.0325890),
            ),
            (
                10.0000000,
                300.000000,
                0.00000000,
                Complex32::new(63.3488808, -28.8427391),
            ),
            (
                10.0000000,
                300.000000,
                15.0000000,
                Complex32::new(60.4814758, -31.9897060),
            ),
            (
                10.0000000,
                300.000000,
                30.0000000,
                Complex32::new(57.9399567, -34.5400085),
            ),
            (
                30.0000000,
                273.149994,
                0.00000000,
                Complex32::new(12.3748817, -22.6335545),
            ),
            (
                30.0000000,
                273.149994,
                15.0000000,
                Complex32::new(12.1546707, -23.0304813),
            ),
            (
                30.0000000,
                273.149994,
                30.0000000,
                Complex32::new(12.3004131, -23.3680649),
            ),
            (
                30.0000000,
                290.000000,
                0.00000000,
                Complex32::new(21.5122662, -30.7900391),
            ),
            (
                30.0000000,
                290.000000,
                15.0000000,
                Complex32::new(20.7252617, -30.8957348),
            ),
            (
                30.0000000,
                290.000000,
                30.0000000,
                Complex32::new(20.3074341, -30.8547306),
            ),
            (
                30.0000000,
                300.000000,
                0.00000000,
                Complex32::new(27.9100685, -33.4010315),
            ),
            (
                30.0000000,
                300.000000,
                15.0000000,
                Complex32::new(26.6922417, -33.5059967),
            ),
            (
                30.0000000,
                300.000000,
                30.0000000,
                Complex32::new(25.8418388, -33.4098701),
            ),
        ];

        for (freq, temp, s, expected_output) in inputs_and_outputs {
            let actual_output = meissner(freq, temp, s);
            // I'm not sure why, but the real part just needs a little higher
            // epsilon than the default, which is `f32::EPSILON`.
            assert_relative_eq!(actual_output.re, expected_output.re, epsilon = 5e-6);
            assert_relative_eq!(actual_output.im, expected_output.im);
        }
    }
}
