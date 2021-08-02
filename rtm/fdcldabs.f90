module cloud_abs
  use, intrinsic :: iso_fortran_env, only: real32
  use dielectric_meissner, only: meissner
  implicit none
  private
  public :: fdcldabs

contains

  !     liquid cloud water absorption
  !     rayleigh
  !     freq:      frequency [ghz]
  !     t:         temperature [k]
  !     rhol:      liquid cloud water density [g/m**3]
  !     output:
  !     al:        cloud water absorption coefficient [neper/km]
  pure subroutine fdcldabs(freq,t,rhol,   al)
    real(real32), intent(in) :: freq, t, rhol
    real(real32), intent(out) :: al
    real(real32) :: rhol0, wavlen
    real(real32), parameter :: c=29.979, pi=3.14159
    complex(real32) :: permit

    rhol0 = 1.0e-6*rhol ![g/cm**3]
    call meissner(freq,t,0.0,   permit)
    wavlen = c/freq
    al = (6.0*pi*rhol0/wavlen)*aimag((1.0-permit)/(2.0+permit))  ! nepers/cm
    al = 1.0e5*al ! nepers/km
  end subroutine fdcldabs
end module cloud_abs
