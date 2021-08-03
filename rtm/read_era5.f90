! Read the ERA5 surface/profile netCDF files
module read_era5
  use, intrinsic :: iso_fortran_env, only: real32, int32
  use netcdf
  implicit none
  private
  public :: Era5DailyData

  ! Daily ERA5 surface/profile data
  type Era5DailyData
    ! The day this data is valid for
    integer(int32) :: year, month, day
    
    ! Lengths of the dimensions
    integer(int32) :: num_hours, num_lats, num_lons, num_levels

    ! Pressure level in hPa, dimensioned as num_levels
    integer(int32), dimension(:), allocatable :: levels

    ! Latitude in degrees North, dimensioned as num_lats
    real(real32), dimension(:), allocatable :: lats

    ! Longitude in degrees East, dimensioned as num_lons
    real(real32), dimension(:), allocatable :: lons

    ! Hours since 1900-01-01, dimensioned as num_hours
    integer(int32), dimension(:), allocatable :: hours

    ! Profile air temperature in Kelvin, dimensioned as (lons, lats, levels, hours)
    real(real32), dimension(:, :, :, :), allocatable :: temperature

    ! Profile relative humidity in percentage (0 to 100), dimensioned as (lons,
    ! lats, levels, hours)
    real(real32), dimension(:, :, :, :), allocatable :: relative_humidity

    ! Profile height in meters, dimensioned as (lons, lats, levels, hours)
    real(real32), dimension(:, :, :, :), allocatable :: height

    ! TODO: cloud water mixing ratio (profile)? Looks like this is only used to
    ! convert to "rhocwat" (cloud water density), and then with temperature, to
    ! "rhol" and "rhoi" (liquid and ice components of the density)

    ! Surface pressure in hPa, dimensioned as (lons, lats, hours)
    real(real32), dimension(:, :, :), allocatable :: surface_pressure

    ! TODO: precipitable water (surface) and columnar cloud water (surface)?

    contains
      procedure :: load => read_era5_data
      final :: free_era5_data
  end type Era5DailyData

contains

  subroutine read_era5_data(self, filename)
    class(Era5DailyData), intent(inout) :: self
    character(len=*), intent(in) :: filename

    ! TODO
  end subroutine read_era5_data

  subroutine free_era5_data(self)
    type(Era5DailyData), intent(inout) :: self

    deallocate(self%levels, self%lats, self%lons, self%hours, &
      self%temperature, self%relative_humidity, self%height, &
      self%surface_pressure)
  end subroutine free_era5_data

end module read_era5