//! Atmospheric radiative transfer model for the ACCESS project

use std::num::NonZeroUsize;

use crate::error::RtmError;

/// Input parameters for the RTM that are constant.
#[derive(Debug)]
pub struct RtmParameters {
    /// Number of frequencies to compute.
    num_freqs: NonZeroUsize,
    /// Microwave frequencies in GHz, with a length of `num_freqs`.
    frequency: Vec<f32>,
    /// Earth incidence angle in degrees, with a length of `num_freqs`.
    incidence: Vec<f32>,
}

/// Inputs for the RTM for a single point. Unlike [`RtmParameters`], these
/// values may vary over location/time.
#[derive(Debug)]
pub struct RtmInputs {
    /// Number of atmosphere profile levels.
    num_levels: NonZeroUsize,
    /// Starting index for the surface, aka `ibegin`.
    surface_index: usize,
    /// Pressure profile in hPa. This has length `num_levels`+1 since the first
    /// element is for the surface.
    pressure: Vec<f32>,
    /// Temperature profile in K. This has length `num_levels`+1 since the first
    /// element is for the surface.
    temperature: Vec<f32>,
    /// Water vapor pressure profile in hPa. This has length `num_levels`+1 since the first
    /// element is for the surface.
    vapor_pressure: Vec<f32>,
    /// Liquid water density in g/m^3. This has length `num_levels`+1 since the first
    /// element is for the surface.
    rho_l: Vec<f32>,
    /// Geometric height in m. This has length `num_levels`+1 since the first
    /// element is for the surface.
    height: Vec<f32>,
}

/// Outputs from the RTM for a single point.
#[derive(Debug)]
pub struct RtmOutputs {
    /// Atmospheric transmissivity as a function of frequency index.
    tran: Vec<f32>,
    /// Atmospheric upwelling in K as a function of frequency index.
    tb_up: Vec<f32>,
    /// Atmospheric downwelling in K as a function of frequency index.
    tb_down: Vec<f32>,
}

impl RtmParameters {
    pub fn new(freqs: &[f32], eia: &[f32]) -> Result<Self, RtmError> {
        if freqs.len() != eia.len() || freqs.is_empty() {
            return Err(RtmError::InconsistentInputs);
        }
        // The unwrap() here is okay since the length of `freqs` is already
        // checked
        Ok(Self {
            num_freqs: freqs.len().try_into().unwrap(),
            frequency: freqs.to_vec(),
            incidence: eia.to_vec(),
        })
    }

    /// Return the number of frequencies used.
    pub fn len(&self) -> usize {
        self.num_freqs.into()
    }
}

impl RtmInputs {
    /// Prepare and convert values.
    ///
    /// The slices (`levels`, `temperature`, etc) must all be the same length.
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
            .into_iter()
            .find(|&i| levels[i] <= surface_pressure)
            .ok_or(RtmError::NoSurface)?;

        // Prepend the surface value to these vectors
        let prepend_with = |level_data: &[f32], zero_value: f32, surface_value: f32| -> Vec<f32> {
            let mut prepended = Vec::with_capacity(num_levels.get() + 1);
            prepended.push(zero_value);
            prepended.extend_from_slice(&level_data);

            prepended[surface_index] = surface_value;
            prepended
        };
        let pressure = prepend_with(levels, 0., surface_pressure);
        let temperature = prepend_with(temperature, surface_temperature, surface_temperature);
        let mut height = prepend_with(height, surface_height, surface_height);

        // Convert geopotential height to geometric height
        for z in height.iter_mut() {
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

        let mut tran = Vec::with_capacity(parameters.num_freqs.get());
        let mut tb_up = Vec::with_capacity(parameters.num_freqs.get());
        let mut tb_down = Vec::with_capacity(parameters.num_freqs.get());

        for (&freq, &inc) in parameters.frequency.iter().zip(&parameters.incidence) {
            // Build up total absorption coefficient profile
            let absorption_profile: Vec<_> = (self.surface_index..self.num_levels.get())
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

/// Modified version of Liebe 1992 oxygen model.
///
/// For a total pressure `p` in hPa, temperature `t` in K, water vapor pressure
/// `pv` in hPa, and frequency `freq` in GHz, compute the oxygen absorption
/// coefficient in dB/km.
///
/// From: Atmospheric 60-GHz Oxygen Spectrum:.. Liebe, Rosenkranz, Hufford,
/// 1992. Modified over the years by Frank Wentz and converted from Fortran to
/// Rust by Richard Lindsley.
fn fdabsoxy_1992_modified(_p: f32, _t: f32, _pv: f32, _freq: f32) -> f32 {
    //   integer(int32), parameter :: nlines=44
    //   integer(int32) :: i
    //   real(real32) :: tht,pwet,pdry,ga,gasq,delta,rnuneg,rnupos,ff,zterm,apterm,sftot,xterm

    //   logical, save :: first = .true.
    //   real(real32), dimension(6, nlines), save :: h
    //   real(real32), dimension(nlines), save :: f0, a1, a2, a3, a4, a5, a6

    //   real(real64) :: sum

    //   data a4/38*0., 6*0.6/

    //   !          freq          a1      a2       a3       a5          a6
    //   data h/ &
    //        50.474238,    0.94e-6,  9.694,  8.60e-3,  0.210,  0.685, &
    //        50.987749,    2.46e-6,  8.694,  8.70e-3,  0.190,  0.680, &
    //        51.503350,    6.08e-6,  7.744,  8.90e-3,  0.171,  0.673, &
    //        52.021410,   14.14e-6,  6.844,  9.20e-3,  0.144,  0.664, &
    //        52.542394,   31.02e-6,  6.004,  9.40e-3,  0.118,  0.653, &
    //        53.066907,   64.10e-6,  5.224,  9.70e-3,  0.114,  0.621, &
    //        53.595749,  124.70e-6,  4.484, 10.00e-3,  0.200,  0.508, &
    //        54.130000,  228.00e-6,  3.814, 10.20e-3,  0.291,  0.375, &
    //        54.671159,  391.80e-6,  3.194, 10.50e-3,  0.325,  0.265, &
    //        55.221367,  631.60e-6,  2.624, 10.79e-3,  0.224,  0.295, &
    //        55.783802,  953.50e-6,  2.119, 11.10e-3, -0.144,  0.613, &
    //        56.264775,  548.90e-6,  0.015, 16.46e-3,  0.339, -0.098, &
    //        56.363389, 1344.00e-6,  1.660, 11.44e-3, -0.258,  0.655, &
    //        56.968206, 1763.00e-6,  1.260, 11.81e-3, -0.362,  0.645, &
    //        57.612484, 2141.00e-6,  0.915, 12.21e-3, -0.533,  0.606, &
    //        58.323877, 2386.00e-6,  0.626, 12.66e-3, -0.178,  0.044, &
    //        58.446590, 1457.00e-6,  0.084, 14.49e-3,  0.650, -0.127, &
    //        59.164207, 2404.00e-6,  0.391, 13.19e-3, -0.628,  0.231, &
    //        59.590983, 2112.00e-6,  0.212, 13.60e-3,  0.665, -0.078, &
    //        60.306061, 2124.00e-6,  0.212, 13.82e-3, -0.613,  0.070, &
    //        60.434776, 2461.00e-6,  0.391, 12.97e-3,  0.606, -0.282, &
    //        61.150560, 2504.00e-6,  0.626, 12.48e-3,  0.090, -0.058, &
    //        61.800154, 2298.00e-6,  0.915, 12.07e-3,  0.496, -0.662, &
    //        62.411215, 1933.00e-6,  1.260, 11.71e-3,  0.313, -0.676, &
    //        62.486260, 1517.00e-6,  0.083, 14.68e-3, -0.433,  0.084, &
    //        62.997977, 1503.00e-6,  1.665, 11.39e-3,  0.208, -0.668, &
    //        63.568518, 1087.00e-6,  2.115, 11.08e-3,  0.094, -0.614, &
    //        64.127767,  733.50e-6,  2.620, 10.78e-3, -0.270, -0.289, &
    //        64.678903,  463.50e-6,  3.195, 10.50e-3, -0.366, -0.259, &
    //        65.224071,  274.80e-6,  3.815, 10.20e-3, -0.326, -0.368, &
    //        65.764772,  153.00e-6,  4.485, 10.00e-3, -0.232, -0.500, &
    //        66.302091,   80.09e-6,  5.225,  9.70e-3, -0.146, -0.609, &
    //        66.836830,   39.46e-6,  6.005,  9.40e-3, -0.147, -0.639, &
    //        67.369598,   18.32e-6,  6.845,  9.20e-3, -0.174, -0.647, &
    //        67.900867,    8.01e-6,  7.745,  8.90e-3, -0.198, -0.655, &
    //        68.431005,    3.30e-6,  8.695,  8.70e-3, -0.210, -0.660, &
    //        68.960311,    1.28e-6,  9.695,  8.60e-3, -0.220, -0.665, &
    //        118.750343,  945.00e-6,  0.009, 16.30e-3, -0.031,  0.008, &
    //        368.498350,   67.90e-6,  0.049, 19.20e-3,  0.0,    0.0,   &
    //        424.763124,  638.00e-6,  0.044, 19.16e-3,  0.0,    0.0,   &
    //        487.249370,  235.00e-6,  0.049, 19.20e-3,  0.0,    0.0,   &
    //        715.393150,   99.60e-6,  0.145, 18.10e-3,  0.0,    0.0,   &
    //        773.839675,  671.00e-6,  0.130, 18.10e-3,  0.0,    0.0,   &
    //        834.145330,  180.00e-6,  0.147, 18.10e-3,  0.0,    0.0/

    //   !$omp critical
    //   if (first) then
    //      first=.false.
    //      f0(:)=h(1,:)
    //      a1(:)=h(2,:)/h(1,:)
    //      a2(:)=h(3,:)
    //      a3(:)=h(4,:)
    //      a5(:)=0.001*h(5,:)
    //      a6(:)=0.001*h(6,:)
    //   end if
    //   !$omp end critical

    //   tht = 300/t
    //   pwet=0.1*pv
    //   pdry=0.1*p-pwet
    //   xterm=1-tht

    //   sum = 0.
    //   !$omp simd reduction(+:sum) private(ga, gasq, delta, rnuneg, rnupos, ff)
    //   do i=1, nlines
    //      ga = a3(i)*(pdry*tht**(0.8-a4(i)) + 1.1*tht*pwet)
    //      gasq=ga*ga
    //      delta=(a5(i) + a6(i)*tht)*p*tht**0.8
    //      rnuneg = f0(i)-freq
    //      rnupos = f0(i)+freq
    //      ff = (ga-rnuneg*delta)/(gasq+rnuneg**2) +  (ga-rnupos*delta)/(gasq+rnupos**2)
    //      sum = sum + ff*a1(i)*exp(a2(i)*xterm)
    //   end do
    //   if (sum < 0) sum=0

    //   !     add nonresonant contribution

    //   !     ga=5.6e-3*(pdry+1.1*pwet)*tht**0.8
    //   ga=5.6e-3*(pdry+1.1*pwet)*tht**1.5  !modification 1

    //   zterm=ga*(1.+(freq/ga)**2)
    //   apterm=1.4e-10*(1-1.2e-5*freq**1.5)*pdry*tht**1.5
    //   if (apterm < 0) apterm=0
    //   sftot=real(pdry*freq*tht**2 * (tht*sum + 6.14e-4/zterm + apterm), real32)

    //   gamoxy=0.1820*freq*sftot
    //   !x    if(freq.gt.37) gamoxy=gamoxy + 0.1820*43.e-10 *pdry**2*tht**3*(freq-37.)**1.7  !prior to 7/17/2015
    //   if (freq > 37) gamoxy=gamoxy + 0.1820*26.e-10 *pdry**2*tht**3*(freq-37.)**1.8  !implemented 7/17/2015.

    // TODO
    0.0
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
fn abh2o_rk_modified(_p: f32, _t: f32, _pv: f32, _freq: f32) -> f32 {
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
fn fdcldabs(_freq: f32, _t: f32, _rhol: f32) -> f32 {
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
/// `num_levels`, where the first index `0` is the value at the surface and
/// indices from `1..num_levels` are profile data above the surface.
fn atm_tran(_inc: f32, _t: &[f32], _z: &[f32], _tabs: &[f32]) -> (f32, f32, f32) {
    //   integer(int32), intent(in) :: nlev
    //   real(real32), intent(in) :: tht
    //   real(real32), dimension(0:nlev), intent(in) :: t, z, tabs
    //   real(real32), intent(out) :: tran, tbdw, tbup

    //   ! real(real32), parameter :: re=6378.135
    //   real(real32), parameter :: delta=0.00035

    //   integer(int32) :: i
    //   real(real32) :: opacty(nlev),tavg(nlev),ems(nlev)
    //   real(real32) :: sumop, sumdw, sumup, tbavg, dsdh

    //   dsdh = (1.0+delta)/sqrt(cosd(tht)**2 + delta*(2+delta))

    //   do i=1,nlev
    //      opacty(i)=-dsdh*0.5*(tabs(i-1)+tabs(i))*(z(i)-z(i-1))
    //      tavg(i)  =0.5*(t(i-1)+t(i))
    //      ems(i)   =1.-exp(opacty(i))
    //   end do

    //   sumop=0
    //   sumdw=0
    //   do i=1,nlev
    //      sumdw=sumdw+(tavg(i)-t(1))*ems(i)*exp(sumop)
    //      sumop=sumop+opacty(i)
    //   end do

    //   sumop=0
    //   sumup=0.
    //   do i=nlev,1,-1
    //      sumup=sumup+(tavg(i)-t(1))*ems(i)*exp(sumop)
    //      sumop=sumop+opacty(i)
    //   end do

    //   tran=exp(sumop)
    //   tbavg=(1.-tran)*t(1)
    //   tbdw=tbavg+sumdw
    //   tbup=tbavg+sumup

    // TODO
    (0.0, 0.0, 0.0)
}
