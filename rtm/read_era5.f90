! Read the ERA5 surface/profile netCDF files
module read_era5
  use, intrinsic :: iso_fortran_env, only: real32, int16, int32, ERROR_UNIT
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use wvap_convert, only: buck_vap
  use netcdf
  implicit none
  private
  public :: Era5DailyData

  ! Daily ERA5 surface/profile data.
  !
  ! The data represented here is the combination of the surface and profile data
  ! for a day with a little bit of post-processing.
  type Era5DailyData
    ! Lengths of the dimensions
    integer(int32) :: num_time, num_lats, num_lons, num_levels

    ! Pressure level in hPa, dimensioned as num_levels. They should be in
    ! descending order (e.g., 1000 to 10).
    integer(int32), dimension(:), allocatable :: levels

    ! Latitude in degrees North, dimensioned as num_lats. They should be in
    ! ascending order (e.g., -90 to 90).
    real(real32), dimension(:), allocatable :: lats

    ! Longitude in degrees East, dimensioned as num_lons. They should be in
    ! ascending order (e.g., -180 to 180).
    real(real32), dimension(:), allocatable :: lons

    ! Hours since 1900-01-01, dimensioned as num_time
    integer(int32), dimension(:), allocatable :: time

    ! Profile air temperature in Kelvin, dimensioned as (lons, lats, levels, time)
    real(real32), dimension(:, :, :, :), allocatable :: temperature

    ! Profile relative humidity in percentage (0 to 100), dimensioned as (lons,
    ! lats, levels, time)
    real(real32), dimension(:, :, :, :), allocatable :: relative_humidity

    ! Geopotential height profile in meters, dimensioned as (lons, lats, levels, time)
    real(real32), dimension(:, :, :, :), allocatable :: height

    ! Profile specific liquid water content (from clouds) in kg/kg, dimensioned
    ! as (lons, lats, levels, time)
    real(real32), dimension(:, :, :, :), allocatable :: liquid_content

    ! Surface pressure in hPa, dimensioned as (lons, lats, time)
    real(real32), dimension(:, :, :), allocatable :: surface_pressure

    ! 2-meter air temperature in kelvin, dimensioned as (lons, lats, time)
    real(real32), dimension(:, :, :), allocatable :: surface_temperature

    ! 2-meter relative humidity in percentage, dimensioned as (lons, lats, time)
    real(real32), dimension(:, :, :), allocatable :: surface_relative_humidity

    ! Geopotential height at the surface in meters, dimensioned as (lons, lats, time)
    real(real32), dimension(:, :, :), allocatable :: surface_height

    ! Total column water vapor in kg/m^2, dimensioned as (lons, lats, time)
    real(real32), dimension(:, :, :), allocatable :: columnar_water_vapor

    ! Total column cloud liquid water in kg/m^2, dimensioned as (lons, lats, time)
    real(real32), dimension(:, :, :), allocatable :: columnar_cloud_liquid

    contains
      procedure :: load => read_era5_data
      final :: free_era5_data
  end type Era5DailyData

contains

! -------------------------------------------------------------------------------
  subroutine read_era5_data(self, filename_profiles, filename_surface)
    class(Era5DailyData), intent(inout) :: self
    character(len=*), intent(in) :: filename_profiles, filename_surface

    integer :: ncid, varid
    integer :: lat_len, lon_len, level_len, time_len
    integer(int16), dimension(:, :, :, :), allocatable :: packed_short_4d
    integer(int16), dimension(:, :, :), allocatable :: packed_short_3d

    ! The reciprocal of the standard gravity, in units of s^2 / m
    ! https://en.wikipedia.org/wiki/Standard_gravity
    real(real32), parameter :: INV_STANDARD_GRAVITY = 1 / 9.80665

    call check_nc(nf90_open(filename_profiles, NF90_NOWRITE, ncid))
    lat_len = dim_len(ncid, "latitude")
    lon_len = dim_len(ncid, "longitude")
    level_len = dim_len(ncid, "level")
    time_len = dim_len(ncid, "time")

    ! If the arrays are already allocated with the right size then reuse them,
    ! otherwise allocate
    if (allocated(self%temperature)) then
      if (self%num_time /= time_len &
        .or. self%num_levels /= level_len &
        .or. self%num_lats /= lat_len &
        .or. self%num_lons /= lon_len) then
          error stop "Allocated size is different than expected"
        end if
    else
      self%num_time = time_len
      self%num_levels = level_len
      self%num_lats = lat_len
      self%num_lons = lon_len
      allocate(self%levels(level_len), &
        self%lats(lat_len), &
        self%lons(lon_len), &
        self%time(time_len), &
        self%temperature(lon_len, lat_len, level_len, time_len), &
        self%relative_humidity(lon_len, lat_len, level_len, time_len), &
        self%height(lon_len, lat_len, level_len, time_len), &
        self%liquid_content(lon_len, lat_len, level_len, time_len), &
        self%surface_pressure(lon_len, lat_len, time_len), &
        self%surface_temperature(lon_len, lat_len, time_len), &
        self%surface_relative_humidity(lon_len, lat_len, time_len), &
        self%surface_height(lon_len, lat_len, time_len), &
        self%columnar_water_vapor(lon_len, lat_len, time_len), &
        self%columnar_cloud_liquid(lon_len, lat_len, time_len))
    end if

    ! Read the coordinate variables
    call check_nc(nf90_inq_varid(ncid, "longitude", varid))
    call check_nc(nf90_get_var(ncid, varid, self%lons))

    call check_nc(nf90_inq_varid(ncid, "latitude", varid))
    call check_nc(nf90_get_var(ncid, varid, self%lats))

    call check_nc(nf90_inq_varid(ncid, "level", varid))
    call check_nc(nf90_get_var(ncid, varid, self%levels))

    call check_nc(nf90_inq_varid(ncid, "time", varid))
    call check_nc(nf90_get_var(ncid, varid, self%time))

    ! The remaining data is packed, so allocate a temporary working array and
    ! unpack the data
    allocate(packed_short_4d(lon_len, lat_len, level_len, time_len))
    call unpack_variable_4d(ncid, "t", packed_short_4d, self%temperature)
    call unpack_variable_4d(ncid, "r", packed_short_4d, self%relative_humidity)
    call unpack_variable_4d(ncid, "z", packed_short_4d, self%height)
    call unpack_variable_4d(ncid, "clwc", packed_short_4d, self%liquid_content)
    deallocate(packed_short_4d)
    call check_nc(nf90_close(ncid))

    call check_nc(nf90_open(filename_surface, NF90_NOWRITE, ncid))
    ! Perhaps the lat/lon/time should be checked to ensure it matches...but for
    ! now we'll just be really trusting

    allocate(packed_short_3d(lon_len, lat_len, time_len))
    call unpack_variable_3d(ncid, "sp", packed_short_3d, self%surface_pressure)
    call unpack_variable_3d(ncid, "t2m", packed_short_3d, self%surface_temperature)
    call unpack_variable_3d(ncid, "d2m", packed_short_3d, self%surface_relative_humidity)
    call unpack_variable_3d(ncid, "z", packed_short_3d, self%surface_height)
    call unpack_variable_3d(ncid, "tcwv", packed_short_3d, self%columnar_water_vapor)
    call unpack_variable_3d(ncid, "tclw", packed_short_3d, self%columnar_cloud_liquid)
    deallocate(packed_short_3d)
    call check_nc(nf90_close(ncid))

    ! Convert geopotential to geopotential height
    ! (https://apps.ecmwf.int/codes/grib/param-db?id=129)
    self%height = self%height * INV_STANDARD_GRAVITY
    self%surface_height = self%surface_height * INV_STANDARD_GRAVITY

    ! Convert surface pressure from Pa to hPa
    self%surface_pressure = self%surface_pressure * 1e-2

    ! Convert dewpoint at 2 m to relative humidity at 2 m
    ! http://bmcnoldy.rsmas.miami.edu/Humidity.html
    ! However, rather than using the "Magnus approximation" I use the Buck equation.
    self%surface_relative_humidity = 100. &
      * buck_vap(self%surface_relative_humidity) / buck_vap(self%surface_temperature)

    ! The latitudes/longitudes need to be adjusted. In the ERA5 files, the
    ! latitudes are in *descending* order from 90 to -90, and while the longitudes
    ! are in ascending order, they go from 0 to 360. The desired output is that
    ! the latitudes go from -90 to 90 and longitudes from -180 to 180.
    self%lats = -self%lats
    self%lons = self%lons - 180.
    ! Flip latitudes and shift longitudes
    self%temperature = cshift(self%temperature(:, lat_len:1:-1, :, :), lon_len / 2, 1)
    self%relative_humidity = cshift(self%relative_humidity(:, lat_len:1:-1, :, :), lon_len / 2, 1)
    self%height = cshift(self%height(:, lat_len:1:-1, :, :), lon_len / 2, 1)
    self%liquid_content = cshift(self%liquid_content(:, lat_len:1:-1, :, :), lon_len / 2, 1)
    self%surface_pressure = cshift(self%surface_pressure(:, lat_len:1:-1, :), lon_len / 2, 1)
    self%surface_temperature = cshift(self%surface_temperature(:, lat_len:1:-1, :), lon_len / 2, 1)
    self%surface_relative_humidity = cshift(self%surface_relative_humidity(:, lat_len:1:-1, :), lon_len / 2, 1)
    self%surface_height = cshift(self%surface_height(:, lat_len:1:-1, :), lon_len / 2, 1)
    self%columnar_water_vapor = cshift(self%columnar_water_vapor(:, lat_len:1:-1, :), lon_len / 2, 1)
    self%columnar_cloud_liquid = cshift(self%columnar_cloud_liquid(:, lat_len:1:-1, :), lon_len / 2, 1)
  end subroutine read_era5_data

  ! -------------------------------------------------------------------------------
  subroutine free_era5_data(self)
    type(Era5DailyData), intent(inout) :: self

    deallocate(self%levels, self%lats, self%lons, self%time, &
      self%temperature, self%relative_humidity, self%height, &
      self%liquid_content, self%surface_pressure, self%surface_temperature, &
      self%surface_relative_humidity, self%surface_height, &
      self%columnar_water_vapor, self%columnar_cloud_liquid)
  end subroutine free_era5_data

  ! -------------------------------------------------------------------------------
  subroutine unpack_variable_3d(ncid, varname, packed_short, unpacked_float)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varname
    integer(int16), dimension(:, :, :), intent(inout) :: packed_short
    real(real32), dimension(:, :, :), intent(out) :: unpacked_float

    integer :: varid, no_fill
    real(real32) :: scale_factor, add_offset
    integer(int16) :: fill_value

    ! Read the packed data
    call check_nc(nf90_inq_varid(ncid, varname, varid))
    call check_nc(nf90_get_var(ncid, varid, packed_short))

    ! Get the needed attributes
    call check_nc(nf90_inq_var_fill(ncid, varid, no_fill, fill_value))
    call check_nc(nf90_get_att(ncid, varid, "scale_factor", scale_factor))
    call check_nc(nf90_get_att(ncid, varid, "add_offset", add_offset))

    ! And unpack the data
    where (packed_short /= fill_value)
      unpacked_float = real(packed_short, real32) * scale_factor + add_offset
    elsewhere
      unpacked_float = ieee_value(0._real32, ieee_quiet_nan)
    end where
  end subroutine unpack_variable_3d

  ! -------------------------------------------------------------------------------
  subroutine unpack_variable_4d(ncid, varname, packed_short, unpacked_float)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: varname
    integer(int16), dimension(:, :, :, :), intent(inout) :: packed_short
    real(real32), dimension(:, :, :, :), intent(out) :: unpacked_float

    integer :: varid, no_fill
    real(real32) :: scale_factor, add_offset
    integer(int16) :: fill_value

    ! Read the packed data
    call check_nc(nf90_inq_varid(ncid, varname, varid))
    call check_nc(nf90_get_var(ncid, varid, packed_short))

    ! Get the needed attributes
    call check_nc(nf90_inq_var_fill(ncid, varid, no_fill, fill_value))
    call check_nc(nf90_get_att(ncid, varid, "scale_factor", scale_factor))
    call check_nc(nf90_get_att(ncid, varid, "add_offset", add_offset))

    ! And unpack the data
    where (packed_short /= fill_value)
      unpacked_float = real(packed_short, real32) * scale_factor + add_offset
    elsewhere
      unpacked_float = ieee_value(0._real32, ieee_quiet_nan)
    end where
  end subroutine unpack_variable_4d

  ! -------------------------------------------------------------------------------
  ! Query the dimension length
  function dim_len(ncid, dimname)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: dimname
    integer :: dim_len

    integer :: dimid
    call check_nc(nf90_inq_dimid(ncid, dimname, dimid))
    call check_nc(nf90_inquire_dimension(ncid, dimid, len=dim_len))
  end function dim_len

  ! -------------------------------------------------------------------------------
  ! Check the netCDF return code and halt if there's an error
  subroutine check_nc(status)
    integer, intent(in) :: status

    if (status /= NF90_NOERR) then
       write (ERROR_UNIT, *) "NetCDF ERROR: " // nf90_strerror(status)
       error stop 1
    end if
  end subroutine check_nc

end module read_era5