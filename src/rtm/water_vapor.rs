//! Water vapor absorption coefficients.
//!
//! These are pretty directly re-written from the original Fortran source.

use std::sync::OnceLock;

const NLINES: usize = 15;

struct WaterVaporCoefficients {
    f0: [f32; NLINES],
    b1: [f32; NLINES],
    b2: [f32; NLINES],
    b3: [f32; NLINES],
    b4: [f32; NLINES],
    b5: [f32; NLINES],
    b6: [f32; NLINES],
}

impl WaterVaporCoefficients {
    /// Initialize the water vapor coefficients.
    #[allow(clippy::excessive_precision)]
    fn new() -> Self {
        // line frequencies
        let f0 = [
            22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508, 443.0183, 448.0011,
            470.8890, 474.6891, 488.4911, 556.9360, 620.7008, 752.0332, 916.1712,
        ];

        // line intensities at 300 K
        let b1: [f32; NLINES] = [
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

        let mut b1_modified = [0.; NLINES];
        #[allow(clippy::excessive_precision)]
        for ((b1_new, b1), &f0) in b1_modified.iter_mut().zip(&b1).zip(&f0) {
            *b1_new = 1.8281089E+14 * b1 / f32::powi(f0, 2);
        }

        // convert b5 to Leibe notation
        let mut b5_modified = [0.; NLINES];
        for ((b5_new, b5), b3) in b5_modified.iter_mut().zip(&b5).zip(&b3) {
            *b5_new = b5 / b3;
        }

        // Modification 1
        b3[0] /= 1.040;

        Self {
            f0,
            b1: b1_modified,
            b2,
            b3,
            b4,
            b5: b5_modified,
            b6,
        }
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

    /// Ensure the coefficients are only initialized once.
    static COEF: OnceLock<WaterVaporCoefficients> = OnceLock::new();
    let WaterVaporCoefficients {
        f0,
        b1,
        b2,
        b3,
        b4,
        b5,
        b6,
    } = COEF.get_or_init(WaterVaporCoefficients::new);

    if pv <= 0. {
        return 0.;
    }

    let pwet = 0.1 * pv;
    let pdry = 0.1 * p - pwet;
    let tht = 300. / t;
    let xterm = 1. - tht;
    let freq_sq = freq.powi(2);

    let sum: f64 = (0..NLINES)
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
