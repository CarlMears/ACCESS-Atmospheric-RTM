//! Low-level RTM functions
//!
//! These are pretty directly re-written from the original Fortran source.

use once_cell::sync::OnceCell;
use smallvec::SmallVec;

const NLINES: usize = 44;

struct OxygenCoefficients {
    f0: [f32; NLINES],
    a1: [f32; NLINES],
    a2: [f32; NLINES],
    a3: [f32; NLINES],
    a4: [f32; NLINES],
    a5: [f32; NLINES],
    a6: [f32; NLINES],
}

/// Oxygen absorption coefficients
static OXY_COEF: OnceCell<OxygenCoefficients> = OnceCell::new();

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
    let OxygenCoefficients {
        f0,
        a1,
        a2,
        a3,
        a4,
        a5,
        a6,
    } = OXY_COEF.get_or_init(|| {
        let mut a4 = [0.; NLINES];
        for i in NLINES - 6..NLINES {
            a4[i] = 0.6;
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

    let sum: f64 = (0..NLINES)
        .into_iter()
        .map(|i| {
            let ga = a3[i] * (pdry * tht.powf(0.8 - a4[i]) + 1.1 * tht * pwet);
            let ga_sq = ga.powi(2);
            let delta = (a5[i] + a6[i] * tht) * p * tht.powf(0.8);

            let rnuneg = f0[i] - freq;
            let rnupos = f0[i] + freq;
            let ff = (ga - rnuneg * delta) / (ga_sq + rnuneg.powi(2))
                + (ga - rnupos * delta) / (ga_sq + rnupos.powi(2));

            f64::from(ff * a1[i] * f32::exp(a2[i] * xterm))
        })
        .sum();
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
pub(crate) fn abh2o_rk_modified(_p: f32, _t: f32, _pv: f32, _freq: f32) -> f32 {
    //   integer(int32), parameter :: nlines=15

    //   integer(int32) :: i
    //   logical, save :: first = .true.
    //   real(real32) :: s,base
    //   real(real32) :: tht,pwet,pdry,ga,gasq,sftot,xterm,rnuneg,rnupos
    //   real(real64) :: sum

    //   real(real32) :: chi,chisq,freqsq,f0sq,u,ffac

    //   !     line frequencies:
    //   real(real32), dimension(nlines), parameter :: f0 = [ &
    //        22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508, &
    //        443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360, &
    //        620.7008, 752.0332, 916.1712]
    //   !     line intensities at 300k:
    //   real(real32), dimension(nlines), save :: b1 = [ &
    //        .1310e-13, .2273e-11, .8036e-13, .2694e-11, .2438e-10, &
    //        .2179e-11, .4624e-12, .2562e-10, .8369e-12, .3263e-11, .6659e-12, &
    //        .1531e-08, .1707e-10, .1011e-08, .4227e-10]
    //   !     t coeff. of intensities:
    //   real(real32), dimension(nlines), parameter :: b2 = [ &
    //        2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405, 3.597, 2.379, 2.852, .159, 2.391, .396, 1.441]
    //   !     air-broadened width parameters at 300k:
    //   real(real32), dimension(nlines), save :: b3 = [ &
    //        .0281, .0281, .023, .0278, .0287, .021, .0186, .0263, .0215, .0236, .026, .0321, .0244, .0306, .0267]
    //   !     self-broadened width parameters at 300k:
    //   real(real32), dimension(nlines), save :: b5 = [ &
    //        .1349, .1491, .108, .135, .1541, .090, .0788, .1275, .0983, .1095, .1313, .1320, .1140, .1253, .1275]
    //   !     t-exponent of air-broadening:
    //   real(real32), dimension(nlines), parameter :: b4 = [ &
    //        .69, .64, .67, .68, .54, .63, .60, .66, .66, .65, .69, .69, .71, .68, .70]
    //   !     t-exponent of self-broadening:
    //   real(real32), dimension(nlines), parameter :: b6 = [ &
    //        .61, .85, .54, .74, .89, .52, .50, .67, .65, .64, .72, 1.0, .68, .84, .78]

    //   !$omp critical
    //   if (first) then
    //      first = .false.
    //      b1=1.8281089E+14*b1/f0**2
    //      b5=b5/b3  !convert b5 to Leibe notation
    //      !1    b1(1)=1.01*b1(1)  !this modification is no longer done
    //      b3(1)=b3(1)/1.040 !modification 1
    //   end if
    //   !$omp end critical

    //   if (pv <= 0.) then
    //      gamh2o=0
    //      return
    //   end if

    //   pwet=0.1*pv
    //   pdry=0.1*p-pwet
    //   tht = 300./t
    //   xterm=1-tht
    //   freqsq=freq*freq

    //   sum = 0.
    //   do i = 1, nlines
    //      f0sq=f0(i)*f0(i)
    //      ga=b3(i)*(pdry*tht**b4(i) + b5(i)*pwet*tht**b6(i))
    //      gasq = ga*ga
    //      s = b1(i)*exp(b2(i)*xterm)
    //      rnuneg = f0(i)-freq
    //      rnupos = f0(i)+freq
    //      base = ga/(562500. + gasq)  !use clough's definition of local line contribution

    //      if (i /= 1) then
    //         if (abs(rnuneg) < 750) sum = sum + s*(ga/(gasq + rnuneg**2) - base)
    //         if (abs(rnupos) < 750) sum = sum + s*(ga/(gasq + rnupos**2) - base)
    //      else
    //         chi=0.07*ga   !modification 2

    //         if (freq < 19) then
    //            u=abs(freq-19.)/16.5
    //            if (u < 0) u=0
    //            if (u > 1) u=1
    //            chi=0.07*ga + 0.93*ga*u*u*(3-2*u)  !modification 2
    //         end if

    //         chisq=chi*chi
    //         sum=sum + s*2*((ga-chi)*freqsq + (ga+chi)*(f0sq+gasq-chisq))/((freqsq-f0sq-gasq+chisq)**2 + 4*freqsq*gasq)
    //      end if
    //   end do
    //   if (sum < 0) sum=0

    //   ffac=1
    //   if (freq < 90) ffac=ffac + 0.1*((90. - freq)/90.)**1.4
    //   !1    sftot=pwet*freq*tht**3.5*(sum +      1.1*1.2957246e-6*pdry/tht**0.5 + 0.425*(freq**0.10)*4.2952193e-5*pwet*tht**4) !previous version
    //   sftot=real(pwet*freq*tht**3.5*(sum + ffac*1.1*1.2957246e-6*pdry/tht**0.5 + 0.348*(freq**0.15)*4.2952193e-5*pwet*tht**4), real32) !modification 3

    //   gamh2o=0.1820*freq*sftot

    // TODO
    0.0
}

/// Liquid cloud water absorption coefficient.
///
/// For a frequency `freq` in GHz, a temperature `t` in K, and a liquid cloud
/// water density `rhol` in g/m^3, compute the cloud water absorption
/// coefficient in Np/km.
pub(crate) fn fdcldabs(_freq: f32, _t: f32, _rhol: f32) -> f32 {
    //   real(real32) :: rhol0, wavlen
    //   real(real32), parameter :: c=29.979, pi=3.14159
    //   complex(real32) :: permit

    //   rhol0 = 1.0e-6*rhol ![g/cm**3]
    //   call meissner(freq,t,0.0,   permit)
    //   wavlen = c/freq
    //   al = (6.0*pi*rhol0/wavlen)*aimag((1.0-permit)/(2.0+permit))  ! Np/cm
    //   al = 1.0e5*al ! Np/km

    // TODO
    0.0
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
