! ---------------------------------------------------------------------------
!
! Wrappers to call trigonometric functions in terms of degrees instead
! of radians
module trig_degrees
  use, intrinsic :: iso_fortran_env, only: real32, real64
  implicit none
  private
  public :: sind, cosd, tand, dsind, dcosd, dtand
  public :: asind, acosd, atand, dasind, dacosd, datand
  public :: atan2d, datan2d
  public :: RAD2DEG_F32, RAD2DEG_F64
  public :: DEG2RAD_F32, DEG2RAD_F64

  real(real32), parameter :: PI_F32 = 4.0_real32 * atan(1.0_real32)
  real(real64), parameter :: PI_F64 = 4.0_real64 * atan(1.0_real64)

  real(real32), parameter :: DEG2RAD_F32 = PI_F32 / 180.0_real32
  real(real64), parameter :: DEG2RAD_F64 = PI_F64 / 180.0_real64

  real(real32), parameter :: RAD2DEG_F32 = 180.0_real32 / PI_F32
  real(real64), parameter :: RAD2DEG_F64 = 180.0_real64 / PI_F64

contains

  ! ------------------------------------------

  pure elemental function sind(x)
    real(real32), intent(in) :: x
    real(real32) :: sind
    sind = sin(x * DEG2RAD_F32)
  end function sind

  pure elemental function cosd(x)
    real(real32), intent(in) :: x
    real(real32) :: cosd
    cosd = cos(x * DEG2RAD_F32)
  end function cosd

  pure elemental function tand(x)
    real(real32), intent(in) :: x
    real(real32) :: tand
    tand = tan(x * DEG2RAD_F32)
  end function tand

  ! ------------------------------------------

  pure elemental function dsind(x)
    real(real64), intent(in) :: x
    real(real64) :: dsind
    dsind = sin(x * DEG2RAD_F64)
  end function dsind

  pure elemental function dcosd(x)
    real(real64), intent(in) :: x
    real(real64) :: dcosd
    dcosd = cos(x * DEG2RAD_F64)
  end function dcosd

  pure elemental function dtand(x)
    real(real64), intent(in) :: x
    real(real64) :: dtand
    dtand = tan(x * DEG2RAD_F64)
  end function dtand

  ! ------------------------------------------

  pure elemental function asind(x)
    real(real32), intent(in) :: x
    real(real32) :: asind
    asind = RAD2DEG_F32 * asin(x)
  end function asind

  pure elemental function acosd(x)
    real(real32), intent(in) :: x
    real(real32) :: acosd
    acosd = RAD2DEG_F32 * acos(x)
  end function acosd

  pure elemental function atand(x)
    real(real32), intent(in) :: x
    real(real32) :: atand
    atand = RAD2DEG_F32 * atan(x)
  end function atand

  ! ------------------------------------------

  pure elemental function dasind(x)
    real(real64), intent(in) :: x
    real(real64) :: dasind
    dasind = RAD2DEG_F64 * asin(x)
  end function dasind

  pure elemental function dacosd(x)
    real(real64), intent(in) :: x
    real(real64) :: dacosd
    dacosd = RAD2DEG_F64 * acos(x)
  end function dacosd

  pure elemental function datand(x)
    real(real64), intent(in) :: x
    real(real64) :: datand
    datand = RAD2DEG_F64 * atan(x)
  end function datand

  ! ------------------------------------------

  pure elemental function atan2d(x, y)
    real(real32), intent(in) :: x, y
    real(real32) :: atan2d
    atan2d = RAD2DEG_F32 * atan2(x, y)
  end function atan2d

  pure elemental function datan2d(x, y)
    real(real64), intent(in) :: x, y
    real(real64) :: datan2d
    datan2d = RAD2DEG_F64 * atan2(x, y)
  end function datan2d

  ! ------------------------------------------

end module trig_degrees
