module dielectric_meissner
  use, intrinsic :: iso_fortran_env, only: real32
  use geomod, only: dielectric_meissner_wentz
  implicit none
  private
  public :: meissner

contains

  !    complex dielectric constant: eps
  !    t. meissner, february 2002
  !    updated oct 2004

  !   input:
  !   name   parameter  unit  range
  !   freq   frequency  [ghz] 1 to 400
  !   t      sst        [k]   248.16 k (-25 c) to 313.16 k (40 c) for pure water
  !                           271.16 k (-2  c) to 307.16 k (34 c) for saline water
  !   s      salinity   [ppt]  0 to 40
  !
  !   output:
  !   eps    complex dielectric constant
  !          negative imaginary part to be consistent with wentz1 conventionc
  !
  pure subroutine  meissner(freq,t,s,   eps)
    real(real32), intent(in) :: freq,t,s
    complex(real32), intent(out) :: eps

    real(real32), parameter :: f0=17.97510
    real(real32) sst
    real(real32) e0s,e1s,e2s,n1s,n2s,sig
    complex(real32), parameter :: j = (0.0,1.0)

    sst = t - 273.15 ! [celsius]
    call dielectric_meissner_wentz(sst,s,  e0s,e1s,e2s,n1s,n2s,sig)

    !     debye law (2 relaxation wavelengths)
    eps = (e0s - e1s)/(1.0 - j*(freq/n1s)) + &
         (e1s - e2s)/(1.0 - j*(freq/n2s)) + e2s + &
         j*sig*f0/freq

    eps = conjg(eps)
  end subroutine meissner
end module dielectric_meissner
