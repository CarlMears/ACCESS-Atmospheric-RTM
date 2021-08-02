module find_abs_coeff
  use, intrinsic :: iso_fortran_env, only: real32, ERROR_UNIT
  use atms_abs_routines, only: fdabsoxy_1992_modified, abh2o_rk_modified
  implicit none
  private
  public :: fdabscoeff

contains

  !   input:
  !     oxygen absorption from rosenkranz
  !     freq  frequency [in ghz]
  !     p      pressure [in h pa]
  !     t      temperature [in k]
  !     pv     water vapor pressure  [in hpa]
  !
  !     output:
  !     av          water vapor absorption coefficients [neper/km]
  !     ao          oxygen absortption coefficient        [neper/km]
  subroutine fdabscoeff(freq,p,t,pv, av,ao)
    real(real32), intent(in) :: freq,p,t,pv
    real(real32), intent(out) :: av,ao

    real(real32), parameter :: nep_scale = 0.1 * log(10.0)
    real(real32) gamoxy,gamh2o

    call fdabsoxy_1992_modified(p,t,pv,freq, gamoxy)
    call abh2o_rk_modified(p,t,pv,freq, gamh2o)

    ! The result values above are in dB/km, so convert to Np/km
    ao=nep_scale * gamoxy
    av=nep_scale * gamh2o
  end subroutine fdabscoeff

end module find_abs_coeff
