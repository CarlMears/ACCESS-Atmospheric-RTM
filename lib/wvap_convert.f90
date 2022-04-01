module wvap_convert
  use, intrinsic :: iso_fortran_env, only: real32
  implicit none
  private
  public :: buck_vap

contains

  ! Use the Buck equation to convert temperature into water vapor saturation
  ! pressure. The equation is from [1], which cites Buck 1996.
  !
  ! To convert to water vapor partial pressure, multiply the result by the
  ! relative humidity.
  !
  ! [1] https://en.wikipedia.org/wiki/Arden_Buck_equation
  elemental function buck_vap(temp)
    real(real32), intent(in) :: temp ! temperature [K]
    real(real32) :: buck_vap ! saturation vapor pressure [hPa]

    real(real32) :: temp_c ! temperature in degrees Celsius

    temp_c = temp - 273.15
    buck_vap = 6.1121 * exp((18.678 - temp_c / 234.5) * (temp_c / (257.14 + temp_c)))
  end function buck_vap

end module wvap_convert
