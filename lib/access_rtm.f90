module access_rtm
  use, intrinsic :: iso_fortran_env, only: int32, real32, real64, ERROR_UNIT
  use, intrinsic :: iso_c_binding, only: c_int, c_int32_t, c_ptr, c_f_pointer
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use atms_abs_routines, only: atm_tran, fdcldabs, fdabscoeff
  use time_conversions, only: yo_to_ymd
  use wvap_convert, only: goff_gratch_vap
  implicit none
  private
  public :: Era5DailyData, RtmDailyData, compute_rtm

  ! ! Number of hours used per day (TODO: probably increase to 24 later)
  ! integer, parameter :: NUM_HR = 2
  ! integer(int32), dimension(NUM_HR), parameter :: HOURS = [0, 12]

  ! ! All data is on a 0.25-degree grid
  ! integer, parameter :: NUM_LAT = 721, NUM_LON = 1440
  ! real(real64), parameter :: DLAT = 0.25, DLON = 0.25
  ! real(real64), parameter :: LAT0 = -90, LON0 = -180

  ! ! Reference frequencies (in GHz) to use
  ! integer, parameter :: NUM_FREQ = 7
  ! real(real32), dimension(NUM_FREQ), parameter :: REF_FREQ = [1.41, 6.8, 10.7, 18.7, 23.8, 37.0, 89.0]

  ! ! Reference Earth incidence angle to use for each reference frequency (in degrees)
  ! real(real32), dimension(NUM_FREQ), parameter :: REF_EIA = [40., 53., 53., 53., 53., 53., 53.]

  ! Maximum number of pressure levels to use
  integer, parameter :: NMAX = 26

  ! This is the C-compatible derived type for the ERA5 daily data
  type, bind(C) :: Era5DailyData
    ! Lengths of the dimensions
    integer(c_int32_t) :: num_time, num_lats, num_lons, num_levels

    ! Pressure level in hPa, dimensioned as num_levels. They should be in
    ! descending order (e.g., 1000 to 10).
    type(c_ptr) :: levels

    ! Latitude in degrees North, dimensioned as num_lats. They should be in
    ! ascending order (e.g., -90 to 90).
    type(c_ptr) :: lats

    ! Longitude in degrees East, dimensioned as num_lons. They should be in
    ! ascending order (e.g., -180 to 180).
    type(c_ptr) :: lons

    ! Hours since 1900-01-01, dimensioned as num_time
    type(c_ptr) :: time

    ! Profile air temperature in Kelvin, dimensioned as (lons, lats, levels, time)
    type(c_ptr) :: temperature

    ! Profile relative humidity in percentage (0 to 100), dimensioned as (lons,
    ! lats, levels, time)
    type(c_ptr) :: relative_humidity

    ! Geopotential height profile in meters, dimensioned as (lons, lats, levels, time)
    type(c_ptr) :: height

    ! Profile specific liquid water content (from clouds) in kg/kg, dimensioned
    ! as (lons, lats, levels, time)
    type(c_ptr) :: liquid_content

    ! Surface pressure in hPa, dimensioned as (lons, lats, time)
    type(c_ptr) :: surface_pressure

    ! 2-meter air temperature in kelvin, dimensioned as (lons, lats, time)
    type(c_ptr) :: surface_temperature

    ! 2-meter relative humidity in percentage, dimensioned as (lons, lats, time)
    type(c_ptr) :: surface_relative_humidity

    ! Geopotential height at the surface in meters, dimensioned as (lons, lats, time)
    type(c_ptr) :: surface_height

    ! Total column water vapor in kg/m^2, dimensioned as (lons, lats, time)
    type(c_ptr) :: columnar_water_vapor

    ! Total column cloud liquid water in kg/m^2, dimensioned as (lons, lats, time)
    type(c_ptr) :: columnar_cloud_liquid
  end type Era5DailyData

  type :: Era5DailyData_f
    ! Lengths of the dimensions
    integer(int32) :: num_time, num_lats, num_lons, num_levels

    ! Pressure level in hPa, dimensioned as num_levels. They should be in
    ! descending order (e.g., 1000 to 10).
    integer(int32), dimension(:), pointer :: levels

    ! Latitude in degrees North, dimensioned as num_lats. They should be in
    ! ascending order (e.g., -90 to 90).
    real(real32), dimension(:), pointer :: lats

    ! Longitude in degrees East, dimensioned as num_lons. They should be in
    ! ascending order (e.g., -180 to 180).
    real(real32), dimension(:), pointer :: lons

    ! Hours since 1900-01-01, dimensioned as num_time
    integer(int32), dimension(:), pointer :: time

    ! Profile air temperature in Kelvin, dimensioned as (lons, lats, levels, time)
    real(real32), dimension(:, :, :, :), pointer :: temperature

    ! Profile relative humidity in percentage (0 to 100), dimensioned as (lons,
    ! lats, levels, time)
    real(real32), dimension(:, :, :, :), pointer :: relative_humidity

    ! Geopotential height profile in meters, dimensioned as (lons, lats, levels, time)
    real(real32), dimension(:, :, :, :), pointer :: height

    ! Profile specific liquid water content (from clouds) in kg/kg, dimensioned
    ! as (lons, lats, levels, time)
    real(real32), dimension(:, :, :, :), pointer :: liquid_content

    ! Surface pressure in hPa, dimensioned as (lons, lats, time)
    real(real32), dimension(:, :, :), pointer :: surface_pressure

    ! 2-meter air temperature in kelvin, dimensioned as (lons, lats, time)
    real(real32), dimension(:, :, :), pointer :: surface_temperature

    ! 2-meter relative humidity in percentage, dimensioned as (lons, lats, time)
    real(real32), dimension(:, :, :), pointer :: surface_relative_humidity

    ! Geopotential height at the surface in meters, dimensioned as (lons, lats, time)
    real(real32), dimension(:, :, :), pointer :: surface_height

    ! Total column water vapor in kg/m^2, dimensioned as (lons, lats, time)
    real(real32), dimension(:, :, :), pointer :: columnar_water_vapor

    ! Total column cloud liquid water in kg/m^2, dimensioned as (lons, lats, time)
    real(real32), dimension(:, :, :), pointer :: columnar_cloud_liquid
  end type Era5DailyData_f

  ! Output values after computing the atmospheric RTM.
  type, bind(C) :: RtmDailyData
    ! Lengths of the dimensions
    integer(c_int32_t) :: num_hour, num_lat, num_lon, num_freq

    ! Coordinate variables. Both freq/eia have the same length.
    type(c_ptr) :: hour
    type(c_ptr) :: lat
    type(c_ptr) :: lon
    type(c_ptr) :: freq
    type(c_ptr) :: eia

    ! These are dimensioned as (lon, lat, hour)
    type(c_ptr) :: col_vapor, col_water

    ! These are dimensioned as (lon, lat, freq, hour)
    type(c_ptr) :: tran, tb_up, tb_down
  end type RtmDailyData

  ! Output values after computing the atmospheric RTM.
  type RtmDailyData_f
    ! Lengths of the dimensions
    integer(int32) :: num_hour, num_lat, num_lon, num_freq

    ! Coordinate variables. Both freq/eia have the same length.
    integer(int32), dimension(:), pointer :: hour
    real(real32), dimension(:), pointer :: lat
    real(real32), dimension(:), pointer :: lon
    real(real32), dimension(:), pointer :: freq
    real(real32), dimension(:), pointer :: eia

    ! These are dimensioned as (lon, lat, hour)
    real(real32), dimension(:, :, :), pointer :: col_vapor, col_water

    ! These are dimensioned as (lon, lat, freq, hour)
    real(real32), dimension(:, :, :, :), pointer :: tran, tb_up, tb_down
  end type RtmDailyData_f

contains

  ! ----------------------------------------------------------------------
  ! Compute the RTM
  !
  ! Returns 0 if all okay
  function compute_rtm(atmo_data_c, rtm_data_c) bind(C)
    type(Era5DailyData), intent(in) :: atmo_data_c
    type(RtmDailyData), intent(inout) :: rtm_data_c
    integer(c_int) :: compute_rtm

    integer :: ihour, ilat, ilon, ifreq
    integer :: ibegin

    type(Era5DailyData_f) :: atmo_data
    type(RtmDailyData_f) :: rtm_data

    real(real32), dimension(0:NMAX) :: p, t, pv, rhol, z
    real(real32) :: tran, tb_up, tb_down

    integer :: num_freq

    ! Convert between C/Fortran inputs/outputs
    compute_rtm = convert_atmo_data(atmo_data_c, atmo_data)
    if (compute_rtm /= 0) then
      write(ERROR_UNIT, *) "Couldn't convert input atmosphere data"
      return
    end if
    compute_rtm = convert_rtm_data(rtm_data_c, rtm_data)
    if (compute_rtm /= 0) then
      write(ERROR_UNIT, *) "Couldn't convert output RTM data"
      return
    end if

    ! Check that the output is properly allocated
    if (size(rtm_data%hour) /= atmo_data%num_time &
      .or. size(rtm_data%lat) /= atmo_data%num_lats &
      .or. size(rtm_data%lon) /= atmo_data%num_lons) then
        write(ERROR_UNIT, *) "Output sizes don't match input sizes"
        compute_rtm = 1
        return
    else if (size(rtm_data%freq) /= size(rtm_data%eia)) then
      write(ERROR_UNIT, *) "Output freq/eia sizes don't match"
      compute_rtm = 1
      return
    ! else if (.not. all([allocated(rtm_data%col_vapor), ...]))
    !   write(ERROR_UNIT, *) "Output values aren't allocated"
    !   compute_rtm = 1
    !   return
    end if
    
    ! These coordinate variables are copied from the input, but the "freq" and
    ! "eia" coordinate variables are presumed to be already set
    rtm_data%hour = atmo_data%time
    rtm_data%lat = atmo_data%lats
    rtm_data%lon = atmo_data%lons

    num_freq = size(rtm_data%freq)
    
    ! allocate(self%col_vapor(NUM_LON, NUM_LAT, NUM_HR), &
    !      self%col_water(NUM_LON, NUM_LAT, NUM_HR), &
    !      self%tran(NUM_LON, NUM_LAT, NUM_FREQ, NUM_HR), &
    !      self%tb_up(NUM_LON, NUM_LAT, NUM_FREQ, NUM_HR), &
    !      self%tb_down(NUM_LON, NUM_LAT, NUM_FREQ, NUM_HR))

    do ihour = 1, atmo_data%num_time
      !$omp parallel do collapse(2) private(ibegin, p, t, pv, rhol, z, tran, tb_up, tb_down)
      do ilat = 1, atmo_data%num_lats
        do ilon = 1, atmo_data%num_lons
          rtm_data%col_vapor(ilon, ilat, ihour) = atmo_data%columnar_water_vapor(ilon, ilat, ihour)
          rtm_data%col_water(ilon, ilat, ihour) = atmo_data%columnar_cloud_liquid(ilon, ilat, ihour)

          call prepare_parameters(atmo_data, ilat, ilon, ihour, ibegin, p, t, pv, rhol, z)
          do ifreq = 1, num_freq
            call atmo_params(p, t, pv, rhol, z, ibegin, &
              rtm_data%eia(ifreq), rtm_data%freq(ifreq), &
              tran, tb_up, tb_down)

              rtm_data%tran(ilon, ilat, ifreq, ihour) = tran
              rtm_data%tb_up(ilon, ilat, ifreq, ihour) = tb_up
              rtm_data%tb_down(ilon, ilat, ifreq, ihour) = tb_down
           end do
        end do
      end do
    end do
  end function compute_rtm

  ! Convert between the C and Fortran data
  function convert_atmo_data(atmo_data_c, atmo_data)
    type(Era5DailyData), intent(in) :: atmo_data_c
    type(Era5DailyData_f), intent(out) :: atmo_data
    integer :: convert_atmo_data

    convert_atmo_data = 0
    atmo_data%num_time = atmo_data_c%num_time
    atmo_data%num_lats = atmo_data_c%num_lats
    atmo_data%num_lons = atmo_data_c%num_lons
    atmo_data%num_levels = atmo_data_c%num_levels

    associate(num_levels => atmo_data_c%num_levels, &
      num_time => atmo_data_c%num_time, &
      num_lats => atmo_data_c%num_lats, &
      num_lons => atmo_data_c%num_lons)
      call c_f_pointer(atmo_data_c%levels, atmo_data%levels, [num_levels])
      call c_f_pointer(atmo_data_c%lats, atmo_data%lats, [num_lats])
      call c_f_pointer(atmo_data_c%lons, atmo_data%lons, [num_lons])
      call c_f_pointer(atmo_data_c%time, atmo_data%time, [num_time])

      call c_f_pointer(atmo_data_c%temperature, atmo_data%temperature, [num_lons, num_lats, num_levels, num_time])
      call c_f_pointer(atmo_data_c%relative_humidity, atmo_data%relative_humidity, [num_lons, num_lats, num_levels, num_time])
      call c_f_pointer(atmo_data_c%height, atmo_data%height, [num_lons, num_lats, num_levels, num_time])
      call c_f_pointer(atmo_data_c%liquid_content, atmo_data%liquid_content, [num_lons, num_lats, num_levels, num_time])

      call c_f_pointer(atmo_data_c%surface_pressure, atmo_data%surface_pressure, [num_lons, num_lats, num_time])
      call c_f_pointer(atmo_data_c%surface_temperature, atmo_data%surface_temperature, [num_lons, num_lats, num_time])
      call c_f_pointer(atmo_data_c%surface_relative_humidity, atmo_data%surface_relative_humidity, [num_lons, num_lats, num_time])
      call c_f_pointer(atmo_data_c%surface_height, atmo_data%surface_height, [num_lons, num_lats, num_time])
      call c_f_pointer(atmo_data_c%columnar_water_vapor, atmo_data%columnar_water_vapor, [num_lons, num_lats, num_time])
      call c_f_pointer(atmo_data_c%columnar_cloud_liquid, atmo_data%columnar_cloud_liquid, [num_lons, num_lats, num_time])
    end associate
  end function convert_atmo_data

  ! Convert between the C and Fortran data
  function convert_rtm_data(rtm_data_c, rtm_data)
    type(RtmDailyData), intent(in) :: rtm_data_c
    type(RtmDailyData_f), intent(out) :: rtm_data
    integer :: convert_rtm_data

    convert_rtm_data = 0
    rtm_data%num_hour = rtm_data_c%num_hour
    rtm_data%num_lat = rtm_data_c%num_lat
    rtm_data%num_lon = rtm_data_c%num_lon
    rtm_data%num_freq = rtm_data_c%num_freq

    associate(num_freq => rtm_data_c%num_freq, &
      num_hour => rtm_data_c%num_hour, &
      num_lat => rtm_data_c%num_lat, &
      num_lon => rtm_data_c%num_lon)
      call c_f_pointer(rtm_data_c%freq, rtm_data%freq, [num_freq])
      call c_f_pointer(rtm_data_c%eia, rtm_data%eia, [num_freq])
      call c_f_pointer(rtm_data_c%lat, rtm_data%lat, [num_lat])
      call c_f_pointer(rtm_data_c%lon, rtm_data%lon, [num_lon])
      call c_f_pointer(rtm_data_c%hour, rtm_data%hour, [num_hour])

      call c_f_pointer(rtm_data_c%col_vapor, rtm_data%col_vapor, [num_lon, num_lat, num_hour])
      call c_f_pointer(rtm_data_c%col_water, rtm_data%col_water, [num_lon, num_lat, num_hour])

      call c_f_pointer(rtm_data_c%tran, rtm_data%tran, [num_lon, num_lat, num_freq, num_hour])
      call c_f_pointer(rtm_data_c%tb_up, rtm_data%tb_up, [num_lon, num_lat, num_freq, num_hour])
      call c_f_pointer(rtm_data_c%tb_down, rtm_data%tb_down, [num_lon, num_lat, num_freq, num_hour])
    end associate
  end function convert_rtm_data

  ! ! ----------------------------------------------------------------------
  ! ! Write a day's standard scene data to a netCDF compatible file
  ! subroutine write_scene(scene, filename_out)
  !   class(DailySceneDataEra5), intent(in) :: scene
  !   character(len=*), intent(in) :: filename_out

  !   integer :: ncid, varid
  !   integer :: dim_hour, dim_lat, dim_lon, dim_freq
  !   character(len=80) :: nc_version_full
  !   character(len=10) :: nc_version
  !   character(len=5) :: time_zone
  !   integer, dimension(8) :: time_vals
  !   character(len=24) :: timestamp
  !   character(len=32) :: epoch
  !   character(len=20) :: time_start, time_end
  !   character(len=:), allocatable :: cmdline, history
  !   integer :: cmdline_len
  !   integer :: month, day

  !   nc_version_full = trim(nf90_inq_libvers())
  !   ! Only use the first whitespace-delimited word, which is the version number
  !   nc_version = nc_version_full(1:index(nc_version_full, " "))

  !   call date_and_time(zone=time_zone, values=time_vals)
  !   write(timestamp, '(I4, "-", I2.2, "-", I2.2, " ", I2.2, ":", I2.2, ":", I2.2, A5)') &
  !        time_vals(1), time_vals(2), time_vals(3), time_vals(5), time_vals(6), time_vals(7), time_zone

  !   call yo_to_ymd(scene%year, scene%doy, month, day)
  !   write(time_start, '(I4, "-", I2.2, "-", I2.2, " 00:00:00Z")') scene%year, month, day
  !   write(time_end, '(I4, "-", I2.2, "-", I2.2, " 23:59:59Z")') scene%year, month, day
  !   write(epoch, '(I4, "-", I2.2, "-", I2.2, " 00:00:00 +0000")') scene%year, month, day

  !   call get_command(length=cmdline_len)
  !   allocate(character(len=cmdline_len) :: cmdline)
  !   allocate(character(len=cmdline_len + 180) :: history)
  !   call get_command(cmdline)
  !   write(history, '(A24, " created: ", A)') timestamp, trim(cmdline)

  !   call handle_nc_err(nf90_create(filename_out, ior(NF90_CLOBBER, NF90_NETCDF4), ncid))

  !   ! Define global attributes
  !   call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "CF-1.9,ACDD-1.3"))
  !   call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "title", "ACCESS RTM output"))
  !   call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "institution", "REMSS"))
  !   call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "history", trim(history)))
  !   call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "netcdf_version_id", trim(nc_version)))
  !  !  call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "git_revision", trim(MWI_SIM_VERSION)))
  !  !  call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "git_date", trim(MWI_SIM_DATE)))
  !   call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "date_created", timestamp))
  !   call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_name", "Remote Sensing Systems"))
  !   call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_email", "support@remss.com"))
  !   call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_url", "http://www.remss.com"))
  !   call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lat_min", -90.0_real64))
  !   call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lat_max", 90.0_real64))
  !   call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lon_min", -180.0_real64))
  !   call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lon_max", 180.0_real64))
  !   call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "time_coverage_start", time_start))
  !   call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "time_coverage_end", time_end))
  !   call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "standard_name_vocabulary", "CF Standard Name Table v78"))

  !   ! Define dimensions
  !   call handle_nc_err(nf90_def_dim(ncid, "hour", NUM_HR, dim_hour))
  !   call handle_nc_err(nf90_def_dim(ncid, "lat", NUM_LAT, dim_lat))
  !   call handle_nc_err(nf90_def_dim(ncid, "lon", NUM_LON, dim_lon))
  !   call handle_nc_err(nf90_def_dim(ncid, "freq", NUM_FREQ, dim_freq))

  !   ! Define and write coordinate variables and their attributes
  !   call handle_nc_err(nf90_def_var(ncid, "hour", NF90_INT, dim_hour, varid))
  !   call handle_nc_err(nf90_put_var(ncid, varid, scene%hour))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "time"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "synoptic hour"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "axis", "T"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "units", "hours since " // trim(epoch)))

  !   call handle_nc_err(nf90_def_var(ncid, "lat", NF90_DOUBLE, dim_lat, varid))
  !   call handle_nc_err(nf90_put_var(ncid, varid, scene%lat))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "latitude"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "units", "degrees_north"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "axis", "Y"))

  !   call handle_nc_err(nf90_def_var(ncid, "lon", NF90_DOUBLE, dim_lon, varid))
  !   call handle_nc_err(nf90_put_var(ncid, varid, scene%lon))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "longitude"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "units", "degrees_east"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "axis", "X"))

  !   call handle_nc_err(nf90_def_var(ncid, "freq", NF90_DOUBLE, dim_freq, varid))
  !   call handle_nc_err(nf90_put_var(ncid, varid, scene%freq))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "sensor_band_central_radiation_frequency"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "frequency"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "units", "GHz"))

  !   call handle_nc_err(nf90_def_var(ncid, "eia", NF90_DOUBLE, dim_freq, varid))
  !   call handle_nc_err(nf90_put_var(ncid, varid, scene%eia))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "sensor_zenith_angle"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "incidence angle"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "units", "degree"))

  !   ! Define and write the datasets and their attributes
  !   call handle_nc_err(nf90_def_var(ncid, "col_vapor", NF90_FLOAT, [dim_lon, dim_lat, dim_hour], varid, &
  !        deflate_level=2, shuffle=.true.))
  !   call handle_nc_err(nf90_put_var(ncid, varid, scene%col_vapor))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "atmosphere_mass_content_of_water_vapor"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "columnar water vapor"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "units", "kg m-2"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

  !   call handle_nc_err(nf90_def_var(ncid, "col_water", NF90_FLOAT, [dim_lon, dim_lat, dim_hour], varid, &
  !        deflate_level=2, shuffle=.true.))
  !   call handle_nc_err(nf90_put_var(ncid, varid, scene%col_water))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "atmosphere_mass_content_of_cloud_liquid_water"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "columnar liquid cloud content"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "units", "kg m-2"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

  !   call handle_nc_err(nf90_def_var(ncid, "tran", NF90_FLOAT, [dim_lon, dim_lat, dim_freq, dim_hour], varid, &
  !        deflate_level=2, shuffle=.true.))
  !   call handle_nc_err(nf90_put_var(ncid, varid, scene%tran))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "atmospheric transmissivity"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

  !   call handle_nc_err(nf90_def_var(ncid, "tb_up", NF90_FLOAT, [dim_lon, dim_lat, dim_freq, dim_hour], varid, &
  !        deflate_level=2, shuffle=.true.))
  !   call handle_nc_err(nf90_put_var(ncid, varid, scene%tb_up))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "upwelling brightness temperature"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "units", "kelvin"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

  !   call handle_nc_err(nf90_def_var(ncid, "tb_down", NF90_FLOAT, [dim_lon, dim_lat, dim_freq, dim_hour], varid, &
  !        deflate_level=2, shuffle=.true.))
  !   call handle_nc_err(nf90_put_var(ncid, varid, scene%tb_down))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "downwelling brightness temperature"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "units", "kelvin"))
  !   call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

  !   call handle_nc_err(nf90_close(ncid))
  ! end subroutine write_scene

  ! ! ----------------------------------------------------------------------
  ! ! This is a very simple netCDF error handler: if there's an error,
  ! ! then display it and exit with a nonzero code.
  ! subroutine handle_nc_err(status)
  !   integer, intent(in) :: status
  !   if (status /= NF90_NOERR) then
  !      write (ERROR_UNIT, *) "ERROR: ", nf90_strerror(status)
  !      error stop "NetCDF error"
  !   end if
  ! end subroutine handle_nc_err

  ! ----------------------------------------------------------------------
  ! From the ERA5 data, prepare these profile/surface parameters for the RTM
  pure subroutine prepare_parameters(era5_data, ilat, ilon, itime, ibegin, p, t, pv, rhol, z)
    type(Era5DailyData_f), intent(in) :: era5_data
    integer, intent(in) :: ilat, ilon, itime
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
    p(1:NMAX) = era5_data%levels(:)

    t(0) = era5_data%surface_temperature(ilon, ilat, itime)
    t(1:NMAX) = era5_data%temperature(ilon, ilat, 1:NMAX, itime)

    hgt(0) = era5_data%surface_height(ilon, ilat, itime)
    hgt(1:NMAX) = era5_data%height(ilon, ilat, 1:NMAX, itime)

    rh(0) = era5_data%surface_relative_humidity(ilon, ilat, itime)
    rh(1:NMAX) = era5_data%relative_humidity(ilon, ilat, 1:NMAX, itime)

    q_l(0) = 0.
    q_l(1:NMAX) = era5_data%liquid_content(ilon, ilat, 1:NMAX, itime)

    ! Find "ibegin", or the starting index for the surface
    ibegin = -1
    do ipr = 1, NMAX
      if (p(ipr) <= era5_data%surface_pressure(ilon, ilat, itime)) then
        ibegin = ipr - 1
        exit
      end if
    end do
    if (ibegin < 0) error stop "Couldn't find ibegin"

    p(ibegin) = era5_data%surface_pressure(ilon, ilat, itime)
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
    ! w is the mass mixing ratio of the water vapor to dry air; note that
    ! pressure is converted from hPa to Pa
    where (p > 0)
      w = (pv * R_dry) / (R_vapor * (p * 1e2 - pv))
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
      ! Water/oxygen absorption coefficients in neper/km
      call fdabscoeff(freq,p(ipr),t(ipr),pv(ipr), abh2o,abo2)

      ! Cloud absorption coefficients in neper/km
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
