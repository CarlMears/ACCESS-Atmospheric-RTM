module daily_scene_ncep
  use, intrinsic :: iso_fortran_env, only: int32, real32, real64, ERROR_UNIT
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use atms_abs_routines, only: atm_tran, fdcldabs, fdabscoeff
  use read_ncep, only: findncep
  use time_conversions, only: yo_to_ymd
  use netcdf
  implicit none
  private
  public :: DailySceneDataNcep

  ! Number of hours used per day
  integer, parameter :: NUM_HR = 4
  integer(int32), dimension(NUM_HR), parameter :: HOURS = [0, 6, 12, 18]

  ! All data is on a 1-degree grid
  integer, parameter :: NUM_LAT = 181, NUM_LON = 360
  real(real64), parameter :: DLAT = 1.0, DLON = 1.0
  real(real64), parameter :: LAT0 = -90, LON0 = -180

  ! Reference frequencies (in GHz) to use
  integer, parameter :: NUM_FREQ = 7
  real(real32), dimension(NUM_FREQ), parameter :: REF_FREQ = [1.41, 6.8, 10.7, 18.7, 23.8, 37.0, 89.0]

  ! Reference Earth incidence angle to use for each reference frequency (in degrees)
  real(real32), dimension(NUM_FREQ), parameter :: REF_EIA = [40., 53., 53., 53., 53., 53., 53.]

  ! Daily scene data
  type DailySceneDataNcep
     integer(int32) :: year
     integer(int32) :: doy

     ! Coordinate variables
     integer(int32), dimension(NUM_HR) :: hour
     real(real64), dimension(NUM_LAT) :: lat
     real(real64), dimension(NUM_LON) :: lon
     real(real64), dimension(NUM_FREQ) :: freq
     real(real64), dimension(NUM_FREQ) :: eia

     ! These are dimensioned as (lon, lat, hour)
     real(real32), dimension(:, :, :), allocatable :: col_vapor, col_water

     ! These are dimensioned as (lon, lat, freq, hour)
     real(real32), dimension(:, :, :, :), allocatable :: tran, tb_up, tb_down

   contains
     procedure :: load => read_daily_scene_data
     procedure :: save => write_scene
     final :: free_daily_scene_data
  end type DailySceneDataNcep

contains

  ! ----------------------------------------------------------------------  
  subroutine read_daily_scene_data(self, year, doy)
    class(DailySceneDataNcep), intent(inout) :: self
    integer(int32), intent(in) :: year, doy

    integer :: ihour, ilat, ilon, ifreq
    integer :: month, day
    integer :: ibegin
    integer, parameter :: NMAX = 26

    real(real32), dimension(0:NMAX) :: p, t, pv, rhov, rhol, z
    real(real32) :: colvap, colwat, pwat, cwat
    real(real32) :: lat, lon, tran, tb_up, tb_down

    self%year = year
    self%doy = doy
    self%hour = HOURS
    self%freq = REF_FREQ
    self%eia = REF_EIA
    self%lat = [(LAT0 + DLAT * ilat, ilat = 0, NUM_LAT - 1)]
    self%lon = [(LON0 + DLON * ilon, ilon = 0, NUM_LON - 1)]

    allocate(self%col_vapor(NUM_LON, NUM_LAT, NUM_HR), &
         self%col_water(NUM_LON, NUM_LAT, NUM_HR), &
         self%tran(NUM_LON, NUM_LAT, NUM_FREQ, NUM_HR), &
         self%tb_up(NUM_LON, NUM_LAT, NUM_FREQ, NUM_HR), &
         self%tb_down(NUM_LON, NUM_LAT, NUM_FREQ, NUM_HR))

    call yo_to_ymd(year, doy, month, day)

    ! Read 1-degree NCEP surface/profile data and apply the RTM
    do ihour = 1, NUM_HR
       write (*, '("   surface/profile and RTM, hour ", i0, "/", i0)') ihour, NUM_HR
       do ilat = 1, NUM_LAT
          do ilon = 1, NUM_LON
             lat = real(self%lat(ilat), real32)
             lon = real(self%lon(ilon), real32)
             if (lon < 0) lon = lon + 360

             call findncep(year, month, day, (ihour - 1) * 6, &
                  lat, lon, &
                  colvap, colwat, pwat, cwat, p, t, pv, rhov, rhol, z, ibegin)

             self%col_water(ilon, ilat, ihour) = colwat
             self%col_vapor(ilon, ilat, ihour) = colvap

             do ifreq = 1, NUM_FREQ
                call atmo_params(colvap, pwat, p, t, pv, rhol, z, ibegin, &
                     REF_EIA(ifreq), REF_FREQ(ifreq), tran, tb_up, tb_down)

                self%tran(ilon, ilat, ifreq, ihour) = tran
                self%tb_up(ilon, ilat, ifreq, ihour) = tb_up
                self%tb_down(ilon, ilat, ifreq, ihour) = tb_down
             end do
          end do
       end do
    end do
  end subroutine read_daily_scene_data

  subroutine free_daily_scene_data(self)
    type(DailySceneDataNcep), intent(inout) :: self

    deallocate(self%col_vapor, self%col_water, &
         self%tran, self%tb_up, self%tb_down)

    self%hour(:) = 0
    self%freq(:) = 0
    self%lat(:) = 0
    self%lon(:) = 0
    self%year = -1
    self%doy = -1
  end subroutine free_daily_scene_data

  ! ----------------------------------------------------------------------  
  ! Write a day's standard scene data to a netCDF compatible file
  subroutine write_scene(scene, filename_out)
    class(DailySceneDataNcep), intent(in) :: scene
    character(len=*), intent(in) :: filename_out

    integer :: ncid, varid
    integer :: dim_hour, dim_lat, dim_lon, dim_freq
    character(len=80) :: nc_version_full
    character(len=10) :: nc_version
    character(len=5) :: time_zone
    integer, dimension(8) :: time_vals
    character(len=24) :: timestamp
    character(len=32) :: epoch
    character(len=20) :: time_start, time_end
    character(len=:), allocatable :: cmdline, history
    integer :: cmdline_len
    integer :: month, day

    nc_version_full = trim(nf90_inq_libvers())
    ! Only use the first whitespace-delimited word, which is the version number
    nc_version = nc_version_full(1:index(nc_version_full, " "))

    call date_and_time(zone=time_zone, values=time_vals)
    write(timestamp, '(I4, "-", I2.2, "-", I2.2, " ", I2.2, ":", I2.2, ":", I2.2, A5)') &
         time_vals(1), time_vals(2), time_vals(3), time_vals(5), time_vals(6), time_vals(7), time_zone

    call yo_to_ymd(scene%year, scene%doy, month, day)
    write(time_start, '(I4, "-", I2.2, "-", I2.2, " 00:00:00Z")') scene%year, month, day
    write(time_end, '(I4, "-", I2.2, "-", I2.2, " 23:59:59Z")') scene%year, month, day
    write(epoch, '(I4, "-", I2.2, "-", I2.2, " 00:00:00 +0000")') scene%year, month, day

    call get_command(length=cmdline_len)
    allocate(character(len=cmdline_len) :: cmdline)
    allocate(character(len=cmdline_len + 180) :: history)
    call get_command(cmdline)
    write(history, '(A24, " created: ", A)') timestamp, trim(cmdline)

    call handle_nc_err(nf90_create(filename_out, ior(NF90_CLOBBER, NF90_NETCDF4), ncid))

    ! Define global attributes
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "CF-1.8,ACDD-1.3"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "title", "ACCESS RTM output"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "institution", "REMSS"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "history", trim(history)))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "netcdf_version_id", trim(nc_version)))
   !  call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "git_revision", trim(MWI_SIM_VERSION)))
   !  call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "git_date", trim(MWI_SIM_DATE)))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "date_created", timestamp))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_name", "Remote Sensing Systems"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_email", "support@remss.com"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "creator_url", "http://www.remss.com"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lat_min", -90.0_real64))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lat_max", 90.0_real64))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lon_min", -180.0_real64))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "geospatial_lon_max", 180.0_real64))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "time_coverage_start", time_start))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "time_coverage_end", time_end))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "standard_name_vocabulary", "CF Standard Name Table v64"))

    ! Define dimensions
    call handle_nc_err(nf90_def_dim(ncid, "hour", NUM_HR, dim_hour))
    call handle_nc_err(nf90_def_dim(ncid, "lat", NUM_LAT, dim_lat))
    call handle_nc_err(nf90_def_dim(ncid, "lon", NUM_LON, dim_lon))
    call handle_nc_err(nf90_def_dim(ncid, "freq", NUM_FREQ, dim_freq))

    ! Define and write coordinate variables and their attributes
    call handle_nc_err(nf90_def_var(ncid, "hour", NF90_INT, dim_hour, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%hour))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "time"))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "synoptic hour"))
    call handle_nc_err(nf90_put_att(ncid, varid, "axis", "T"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "hours since " // trim(epoch)))

    call handle_nc_err(nf90_def_var(ncid, "lat", NF90_DOUBLE, dim_lat, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%lat))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "latitude"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "degrees_north"))
    call handle_nc_err(nf90_put_att(ncid, varid, "axis", "Y"))

    call handle_nc_err(nf90_def_var(ncid, "lon", NF90_DOUBLE, dim_lon, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%lon))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "longitude"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "degrees_east"))
    call handle_nc_err(nf90_put_att(ncid, varid, "axis", "X"))

    call handle_nc_err(nf90_def_var(ncid, "freq", NF90_DOUBLE, dim_freq, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%freq))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "sensor_band_central_radiation_frequency"))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "frequency"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "GHz"))

    call handle_nc_err(nf90_def_var(ncid, "eia", NF90_DOUBLE, dim_freq, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%eia))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "sensor_zenith_angle"))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "incidence angle"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "degree"))

    ! Define and write the datasets and their attributes
    call handle_nc_err(nf90_def_var(ncid, "col_vapor", NF90_FLOAT, [dim_lon, dim_lat, dim_hour], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%col_vapor))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "atmosphere_mass_content_of_water_vapor"))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "columnar water vapor"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "kg m-2"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

    call handle_nc_err(nf90_def_var(ncid, "col_water", NF90_FLOAT, [dim_lon, dim_lat, dim_hour], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%col_water))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "atmosphere_mass_content_of_cloud_liquid_water"))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "columnar liquid cloud content"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "kg m-2"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

    call handle_nc_err(nf90_def_var(ncid, "tran", NF90_FLOAT, [dim_lon, dim_lat, dim_freq, dim_hour], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%tran))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "atmospheric transmissivity"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

    call handle_nc_err(nf90_def_var(ncid, "tb_up", NF90_FLOAT, [dim_lon, dim_lat, dim_freq, dim_hour], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%tb_up))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "upwelling brightness temperature"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "kelvin"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

    call handle_nc_err(nf90_def_var(ncid, "tb_down", NF90_FLOAT, [dim_lon, dim_lat, dim_freq, dim_hour], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%tb_down))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "downwelling brightness temperature"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "kelvin"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

    call handle_nc_err(nf90_close(ncid))
  end subroutine write_scene

  ! ----------------------------------------------------------------------  
  ! This is a very simple netCDF error handler: if there's an error,
  ! then display it and exit with a nonzero code.
  subroutine handle_nc_err(status)
    integer, intent(in) :: status
    if (status /= NF90_NOERR) then
       write (ERROR_UNIT, *) "ERROR: ", nf90_strerror(status)
       error stop "NetCDF error"
    end if
  end subroutine handle_nc_err

  ! ----------------------------------------------------------------------  
  ! Apply the RTM to obtain atmospheric terms
  subroutine atmo_params(colvap, pwat, p, t, pv, rhol, z, ibegin, &
       eia, freq, tran, tb_up, tb_down)
    real(real32), intent(in) :: colvap, pwat
    integer, parameter :: NMAX = 26
    real(real32), dimension(0:NMAX), intent(in) :: p, t, pv, rhol, z
    integer, intent(in) :: ibegin
    real(real32), intent(in) :: eia, freq
    real(real32), intent(out) :: tran, tb_up, tb_down

    real(real32) :: abh2o,abo2,abcld
    real(real32), dimension(0:NMAX) :: tabs
    real(real32) :: pv_fix
    integer :: ipr

    do ipr = ibegin,NMAX
       !     for 2002-2007, there is a big dif between vap and pwat, with vap being near 36 mm for amazon and pwat being 50 mm
       !     in 2008 vap increases and is very similar to pwat.  I do not understand why vap changes.
       !     to make things consistent, here is normalized to pwat
       pv_fix=(pwat/colvap)*pv(ipr)

       call fdabscoeff(freq,p(ipr),t(ipr),pv_fix,  abh2o,abo2) ! neper/km

       if (rhol(ipr) > 1.0e-7) then
          call fdcldabs(freq,t(ipr),rhol(ipr),   abcld) ! nepers/km
       else
          abcld = 0.0
       endif

       ! total absorption, convert from Np/km to Np/m
       tabs(ipr) = (abh2o + abo2 + abcld) * 1.e-3
    enddo

    call atm_tran(NMAX-ibegin,eia,t(ibegin:NMAX),z(ibegin:NMAX),tabs(ibegin:NMAX), &
         tran,tb_down,tb_up)

  end subroutine atmo_params

end module daily_scene_ncep
