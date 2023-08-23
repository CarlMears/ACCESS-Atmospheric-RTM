//! Oxygen absorption coefficients.
//!
//! These are pretty directly re-written from the original Fortran source.

use std::sync::OnceLock;

const NLINES: usize = 44;

/// Oxygen absorption coefficients
struct OxygenCoefficients {
    f0: [f32; NLINES],
    a1: [f32; NLINES],
    a2: [f32; NLINES],
    a3: [f32; NLINES],
    a4: [f32; NLINES],
    a5: [f32; NLINES],
    a6: [f32; NLINES],
}

impl OxygenCoefficients {
    /// Initialize the oxygen coefficients.
    #[allow(clippy::excessive_precision)]
    fn new() -> Self {
        // All but the last six entries in a4 are 0
        let mut a4 = [0.; NLINES];
        for a4 in a4.iter_mut().skip(NLINES - 6) {
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

        let mut a1 = [0.; NLINES];
        a1.iter_mut().zip(&h2).zip(&h1).for_each(|((a1, h2), h1)| {
            *a1 = h2 / h1;
        });
        let f0 = h1;
        let a2 = h3;
        let a3 = h4;
        let a5 = h5.map(|h5| 0.001 * h5);
        let a6 = h6.map(|h6| 0.001 * h6);

        Self {
            f0,
            a1,
            a2,
            a3,
            a4,
            a5,
            a6,
        }
    }
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
    // Many of the variables are retained from the original Fortran

    /// Ensure the coefficients are only initialized once.
    static COEF: OnceLock<OxygenCoefficients> = OnceLock::new();
    let OxygenCoefficients {
        f0,
        a1,
        a2,
        a3,
        a4,
        a5,
        a6,
    } = COEF.get_or_init(OxygenCoefficients::new);

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

        let mut ff = [0.; NLINES];
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
