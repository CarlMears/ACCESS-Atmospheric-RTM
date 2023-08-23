//! Liquid cloud water absorption coefficients.
//!
//! These are pretty directly re-written from the original Fortran source.

use num_complex::Complex32;

/// Liquid cloud water absorption coefficient.
///
/// For a frequency `freq` in GHz, a temperature `t` in K, and a liquid cloud
/// water density `rhol` in g/m³, compute the cloud water absorption
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
pub(super) fn meissner(freq: f32, t: f32, s: f32) -> Complex32 {
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
/// # References
///
/// T. Meissner and F. J. Wentz, "The complex dielectric constant of pure and
/// sea water from microwave satellite observations", in IEEE Transactions on
/// Geoscience and Remote Sensing, vol. 42, no. 9, pp. 1836-1849, Sept. 2004,
/// <https://doi.org/10.1109/TGRS.2004.831888>.
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
