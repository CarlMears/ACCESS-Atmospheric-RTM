! This is the public interface to the RTM
!
! Sadly, due to the limitations of f2py, derived types are not supported.

module access_rtm
  use, intrinsic :: iso_fortran_env, only: int32, real32, real64, ERROR_UNIT
  use, intrinsic :: iso_c_binding, only: c_int, c_int32_t, c_ptr, c_f_pointer
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use atms_abs_routines, only: atm_tran, fdcldabs, fdabscoeff
  use wvap_convert, only: goff_gratch_vap
  implicit none
  private
  public :: NMAX, compute_rtm

  ! ! Number of hours used per day (TODO: probably increase to 24 later)
  ! integer, parameter :: NUM_HR = 2
  ! integer(int32), dimension(NUM_HR), parameter :: HOURS = [0, 12]

  ! ! All data is on a 0.25-degree grid
  ! integer, parameter :: NUM_LAT = 721, NUM_LON = 1440
  ! real(real64), parameter :: DLAT = 0.25, DLON = 0.25
  ! real(real64), parameter :: LAT0 = -90, LON0 = -180

  ! Maximum number of pressure levels to use
  integer, parameter :: NMAX = 26

contains

  ! ----------------------------------------------------------------------
  ! Compute the RTM
  !
  ! Returns 0 if all okay
  subroutine compute_rtm(num_points, num_freq, &
    levels, temperature, height, relative_humidity, liquid_content, &
    surface_temperature, surface_height, surface_relative_humidity, surface_pressure, &
    eia, freq, tran, tb_up, tb_down)
    integer, intent(in) :: num_points, num_freq
    real(real32), dimension(NMAX), intent(in) :: levels
    real(real32), dimension(num_points), intent(in) :: surface_temperature, surface_height, &
      surface_relative_humidity, surface_pressure
    real(real32), dimension(NMAX, num_points), intent(in) :: temperature, height, &
      relative_humidity, liquid_content
    real(real32), dimension(num_freq), intent(in) :: eia, freq
    real(real32), dimension(num_freq, num_points), intent(out) :: tran, tb_up, tb_down

    integer :: ibegin, ifreq, i
    real(real32), dimension(0:NMAX) :: p, t, pv, rhol, z

    !$omp parallel do private(ibegin, p, t, pv, rhol, z)
    do i = 1, num_points
      call prepare_parameters(levels(:), surface_temperature(i), temperature(:, i), &
        surface_height(i), height(:, i), &
        surface_relative_humidity(i), relative_humidity(:, i), &
        liquid_content(:, i), surface_pressure(i), &
        ibegin, p, t, pv, rhol, z)
      do ifreq = 1, num_freq
        call atmo_params(p, t, pv, rhol, z, ibegin, eia(ifreq), freq(ifreq), &
          tran(ifreq, i), tb_up(ifreq, i), tb_down(ifreq, i))
      end do
    end do
  end subroutine compute_rtm

  ! ----------------------------------------------------------------------
  ! From the ERA5 data, prepare these profile/surface parameters for the RTM
  pure subroutine prepare_parameters(levels, surface_temperature, temperature, &
    surface_height, height, surface_relative_humidity, relative_humidity, &
    liquid_content, surface_pressure, &
    ibegin, p, t, pv, rhol, z)
    real(real32), intent(in) :: surface_temperature, surface_height, surface_relative_humidity, surface_pressure
    real(real32), dimension(NMAX), intent(in) :: levels, temperature, height, relative_humidity, liquid_content
    integer, intent(out) :: ibegin
    real(real32), dimension(0:NMAX), intent(out) :: p, t, pv, rhol, z

    integer :: ipr
    real(real32), dimension(0:NMAX) :: hgt, rh, rhov, q_l
    real(real32), dimension(0:NMAX) :: R_moist, q_h2o, w

    ! Ideal gas constant (J/mol/K)
    real(real32), parameter :: R = 8.3144598
    ! Mean molar mass of dry air (g/mol)
    real(real32), parameter :: M_dry = 28.9644
    ! Mean molar mass of water (g/mol)
    real(real32), parameter :: M_h2o = 18.01528
    ! Specific gas constant for dry air (J/g/K)
    real(real32), parameter :: R_dry = R / M_dry
    ! Specific gas constant for water vapor (J/g/K)
    real(real32), parameter :: R_vapor = R / M_h2o
    ! Mean radius of the Earth in meters
    real(real32), parameter :: R_e = 6371e3
    ! Coefficient for ratio between molar masses
    real(real32), parameter :: epsilon = M_h2o / M_dry
    ! Scaling factor using epsilon
    real(real32), parameter :: eps_scale = (1 - epsilon) / epsilon

    p(0) = 0.
    p(1:NMAX) = levels(:)

    t(0) = surface_temperature
    t(1:NMAX) = temperature(1:NMAX)

    hgt(0) = surface_height
    hgt(1:NMAX) = height(1:NMAX)

    rh(0) = surface_relative_humidity
    rh(1:NMAX) = relative_humidity(1:NMAX)

    q_l(0) = 0.
    q_l(1:NMAX) = liquid_content(1:NMAX)

    ! Find "ibegin", or the starting index for the surface
    ibegin = -1
    do ipr = 1, NMAX
      if (p(ipr) <= surface_pressure) then
        ibegin = ipr - 1
        exit
      end if
    end do
    if (ibegin < 0) error stop "Couldn't find ibegin"

    p(ibegin) = surface_pressure
    t(ibegin) = t(0)
    hgt(ibegin) = hgt(0)
    rh(ibegin) = rh(0)
    q_l(ibegin) = q_l(ibegin+1)

    ! Convert geopotential height to geometric height
    z = hgt * R_e / (R_e - hgt)
    if (z(ibegin) >= z(ibegin+1)) z(ibegin) = z(ibegin+1) - 0.1

    ! Find the vapor pressure and water vapor density
    call goff_gratch_vap(t, rh, p, pv, rhov)

    ! Convert relative humidity to specific humidity
    ! (https://earthscience.stackexchange.com/a/5077)
    !
    ! w is the mass mixing ratio of the water vapor to dry air
    where (p > 0)
      w = (pv * R_dry) / (R_vapor * (p - pv))
      q_h2o = w / (w + 1)
    elsewhere
      q_h2o = 0
    end where

    ! Convert specific cloud liquid water content (kg/kg) to liquid water
    ! density (g/m^3).
    !
    ! See here, section 4:
    ! https://www.nwpsaf.eu/site/download/documentation/rtm/docs_rttov12/rttov_gas_cloud_aerosol_units.pdf
    ! gas constant for humid air (J/gK)
    R_moist(:) = R_dry * (1 + eps_scale * q_h2o)
    rhol(:) = q_l * (1e2 * p) / (R_moist * T)
  end subroutine prepare_parameters

  ! ----------------------------------------------------------------------
  ! Apply the RTM to obtain atmospheric terms
  subroutine atmo_params(p, t, pv, rhol, z, ibegin, &
       eia, freq, tran, tb_up, tb_down)
    real(real32), dimension(0:NMAX), intent(in) :: p, t, pv, rhol, z
    integer, intent(in) :: ibegin
    real(real32), intent(in) :: eia, freq
    real(real32), intent(out) :: tran, tb_up, tb_down

    real(real32) :: abh2o,abo2,abcld
    real(real32), dimension(0:NMAX) :: tabs
    integer :: ipr

    do ipr = ibegin,NMAX
      ! Water/oxygen absorption coefficients in Np/km
      call fdabscoeff(freq,p(ipr),t(ipr),pv(ipr), abh2o,abo2)

      ! Cloud absorption coefficients in Np/km
      if (rhol(ipr) > 1.0e-7) then
          call fdcldabs(freq,t(ipr),rhol(ipr), abcld)
       else
          abcld = 0.0
       endif

       ! Total absorption coefficient, converting from Np/km to Np/m
       tabs(ipr) = (abh2o + abo2 + abcld) * 1.e-3
    end do

    call atm_tran(NMAX-ibegin,eia,t(ibegin:NMAX),z(ibegin:NMAX),tabs(ibegin:NMAX), &
         tran,tb_down,tb_up)
  end subroutine atmo_params

end module access_rtm
