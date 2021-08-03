! Read the ERA5 surface/profile netCDF files
module read_era5
  use, intrinsic :: iso_fortran_env, only: real32, int32, ERROR_UNIT
  use netcdf
  implicit none
  private
  public :: Era5DailyData

  ! Daily ERA5 surface/profile data
  type Era5DailyData
    ! The day this data is valid for
    integer(int32) :: year, month, day

    ! Lengths of the dimensions
    integer(int32) :: num_time, num_lats, num_lons, num_levels

    ! Pressure level in hPa, dimensioned as num_levels
    integer(int32), dimension(:), allocatable :: levels

    ! Latitude in degrees North, dimensioned as num_lats
    real(real32), dimension(:), allocatable :: lats

    ! Longitude in degrees East, dimensioned as num_lons
    real(real32), dimension(:), allocatable :: lons

    ! Hours since 1900-01-01, dimensioned as num_time
    integer(int32), dimension(:), allocatable :: time

    ! Profile air temperature in Kelvin, dimensioned as (lons, lats, levels, time)
    real(real32), dimension(:, :, :, :), allocatable :: temperature

    ! Profile relative humidity in percentage (0 to 100), dimensioned as (lons,
    ! lats, levels, time)
    real(real32), dimension(:, :, :, :), allocatable :: relative_humidity

    ! Profile height in meters, dimensioned as (lons, lats, levels, time)
    real(real32), dimension(:, :, :, :), allocatable :: height

    ! TODO: cloud water mixing ratio (profile)? Looks like this is only used to
    ! convert to "rhocwat" (cloud water density), and then with temperature, to
    ! "rhol" and "rhoi" (liquid and ice components of the density)

    ! Surface pressure in hPa, dimensioned as (lons, lats, time)
    real(real32), dimension(:, :, :), allocatable :: surface_pressure

    ! TODO: precipitable water (surface) and columnar cloud water (surface)?

    contains
      procedure :: load => read_era5_data
      final :: free_era5_data
  end type Era5DailyData

contains

  subroutine read_era5_data(self, filename_profiles, filename_surface)
    class(Era5DailyData), intent(inout) :: self
    character(len=*), intent(in) :: filename_profiles, filename_surface

    integer :: ncid, varid
    integer :: lat_len, lon_len, level_len, time_len

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
        self%surface_pressure(lon_len, lat_len, time_len))
    end if

    call check_nc(nf90_inq_varid(ncid, "longitude", varid))
    call check_nc(nf90_get_var(ncid, varid, self%lons))

    call check_nc(nf90_inq_varid(ncid, "latitude", varid))
    call check_nc(nf90_get_var(ncid, varid, self%lats))

    call check_nc(nf90_inq_varid(ncid, "level", varid))
    call check_nc(nf90_get_var(ncid, varid, self%levels))

    call check_nc(nf90_inq_varid(ncid, "time", varid))
    call check_nc(nf90_get_var(ncid, varid, self%time))

    ! TODO: read remaining data

    call check_nc(nf90_close(ncid))

    call check_nc(nf90_open(filename_surface, NF90_NOWRITE, ncid))
    call check_nc(nf90_close(ncid))
  end subroutine read_era5_data

  subroutine free_era5_data(self)
    type(Era5DailyData), intent(inout) :: self

    deallocate(self%levels, self%lats, self%lons, self%time, &
      self%temperature, self%relative_humidity, self%height, &
      self%surface_pressure)
  end subroutine free_era5_data

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