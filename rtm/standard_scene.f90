module standard_scene
  use, intrinsic :: iso_fortran_env, only: int32, real32, real64, ERROR_UNIT
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use atmos_tables, only: atm_tran
  use cloud_abs, only: fdcldabs
  use find_abs_coeff, only: fdabscoeff
  use month_day, only: find_month_day
  use ncep, only: findncep
  use netcdf
  use trig_degrees, only: atan2d
  implicit none
  private
  public :: write_scene, daily_scene_data

  integer, parameter :: NUM_HR = 4, NUM_FREQ = 6
  integer, parameter :: NUM_LAT_LO = 181, NUM_LON_LO = 360
  integer, parameter :: NUM_LAT_HI = 721, NUM_LON_HI = 1440

  real(real64), parameter :: DLAT_HI = 0.25, DLON_HI = 0.25
  real(real64), parameter :: DLAT_LO = 1, DLON_LO = 1
  real(real64), parameter :: LAT0 = -90, LON0 = -180

  ! Nominal MWI Earth incidence angle
  real(real32), parameter :: EIA_NOMINAL = 53.

  ! Daily scene data
  !
  ! Two spatial resolutions are used: most of the data is on a
  ! 1-degree grid. However, the wind data is on a 0.25-degree grid.
  ! When writing to the output file, each low-resolution (1-degree
  ! grid) variable is bilinearly interpolated to the high-resolution
  ! grid.
  type daily_scene_data
     integer(int32) :: year
     integer(int32) :: doy

     ! Coordinate variables
     integer(int32), dimension(NUM_HR) :: hour
     real(real64), dimension(NUM_LAT_LO) :: lat_lo
     real(real64), dimension(NUM_LAT_HI) :: lat_hi
     real(real64), dimension(NUM_LON_LO) :: lon_lo
     real(real64), dimension(NUM_LON_HI) :: lon_hi
     real(real64), dimension(NUM_FREQ) :: freq
     integer(int32), dimension(2) :: pol

     ! These are dimensioned as (lon_lo, lat_lo, hour)
     real(real32), dimension(:, :, :), allocatable :: sst, col_vapor, col_water, rain

     ! These are dimensioned as (lon_lo, lat_lo, freq, hour)
     real(real32), dimension(:, :, :, :), allocatable :: tran, tb_up, tb_down

     ! These are dimensioned as (lon_hi, lat_hi, hour)
     real(real32), dimension(:, :, :), allocatable :: wspd, wdir, tec

     ! These are dimensioned as (lon_hi, lat_hi)
     real(real32), dimension(:, :), allocatable :: ice_conc, ice_fraction_new, ice_fraction_first, ice_fraction_multi

     ! These are dimensioned as (lon_hi, lat_hi, freq, pol)
     real(real32), dimension(:, :, :, :), allocatable :: ice_emissivity

   contains
     procedure :: load => read_daily_scene_data
     procedure :: free => free_daily_scene_data
  end type daily_scene_data

contains

  subroutine read_daily_scene_data(self, year, doy)
    class(daily_scene_data), intent(inout) :: self
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
    self%hour = [0, 6, 12, 18]
    self%freq = [10.85, 18.85, 23.8, 36.75, 37.3, 89.]
    self%pol = [1, 2]
    self%lat_lo = [(LAT0 + DLAT_LO * ilat, ilat = 0, NUM_LAT_LO - 1)]
    self%lon_lo = [(LON0 + DLON_LO * ilon, ilon = 0, NUM_LON_LO - 1)]
    self%lat_hi = [(LAT0 + DLAT_HI * ilat, ilat = 0, NUM_LAT_HI - 1)]
    self%lon_hi = [(LON0 + DLON_HI * ilon, ilon = 0, NUM_LON_HI - 1)]

    allocate(self%wspd(NUM_LON_HI, NUM_LAT_HI, NUM_HR), &
         self%wdir(NUM_LON_HI, NUM_LAT_HI, NUM_HR), &
         self%sst(NUM_LON_LO, NUM_LAT_LO, NUM_HR), &
         self%col_vapor(NUM_LON_LO, NUM_LAT_LO, NUM_HR), &
         self%col_water(NUM_LON_LO, NUM_LAT_LO, NUM_HR), &
         self%rain(NUM_LON_LO, NUM_LAT_LO, NUM_HR), &
         self%tec(NUM_LON_HI, NUM_LAT_HI, NUM_HR), &
         self%tran(NUM_LON_LO, NUM_LAT_LO, NUM_FREQ, NUM_HR), &
         self%tb_up(NUM_LON_LO, NUM_LAT_LO, NUM_FREQ, NUM_HR), &
         self%tb_down(NUM_LON_LO, NUM_LAT_LO, NUM_FREQ, NUM_HR), &
         self%ice_conc(NUM_LON_HI, NUM_LAT_HI), &
         self%ice_fraction_new(NUM_LON_HI, NUM_LAT_HI), &
         self%ice_fraction_first(NUM_LON_HI, NUM_LAT_HI), &
         self%ice_fraction_multi(NUM_LON_HI, NUM_LAT_HI), &
         self%ice_emissivity(NUM_LON_HI, NUM_LAT_HI, NUM_FREQ, 2))

    call find_month_day(year, doy, month, day)

    ! Read 1-degree NCEP surface/profile data and apply the RTM
    do ihour = 1, 4
       write (*, '("   surface/profile and RTM, hour ", i0, "/4")') ihour
       do ilat = 1, NUM_LAT_LO
          do ilon = 1, NUM_LON_LO
             lat = real(self%lat_lo(ilat), real32)
             lon = real(self%lon_lo(ilon), real32)
             if (lon < 0) lon = lon + 360

             call findncep(year, month, day, (ihour - 1) * 6, &
                  lat, lon, &
                  colvap, colwat, pwat, cwat, p, t, pv, rhov, rhol, z, ibegin)

             self%col_water(ilon, ilat, ihour) = colwat
             self%col_vapor(ilon, ilat, ihour) = colvap
             self%sst(ilon, ilat, ihour) = t(0)

             do ifreq = 1, NUM_FREQ
                call atmo_params(colvap, pwat, p, t, pv, rhol, z, ibegin, &
                     EIA_NOMINAL, real(self%freq(ifreq), real32), tran, tb_up, tb_down)

                self%tran(ilon, ilat, ifreq, ihour) = tran
                self%tb_up(ilon, ilat, ifreq, ihour) = tb_up
                self%tb_down(ilon, ilat, ifreq, ihour) = tb_down
             end do
          end do
       end do
    end do

    ! Transform from integrated liquid cloud to equivalent rain rate
    write (*, *) "  rain"
    self%rain(:, :, :) = cld_to_rain(self%sst, self%col_water)
  end subroutine read_daily_scene_data

  subroutine free_daily_scene_data(self)
    class(daily_scene_data), intent(inout) :: self

    deallocate(self%wspd, self%wdir, self%sst, self%col_vapor, self%col_water, self%rain, self%tec, &
         self%tran, self%tb_up, self%tb_down, self%ice_conc, &
         self%ice_fraction_new, self%ice_fraction_first, &
         self%ice_fraction_multi, self%ice_emissivity)

    self%hour(:) = 0
    self%freq(:) = 0
    self%pol(:) = 0
    self%lat_lo(:) = 0
    self%lon_lo(:) = 0
    self%lat_hi(:) = 0
    self%lon_hi(:) = 0
    self%year = -1
    self%doy = -1
  end subroutine free_daily_scene_data

  ! Write a day's standard scene data to a netCDF compatible file
  subroutine write_scene(scene, filename_out)
    type(daily_scene_data), intent(in) :: scene
    character(len=*), intent(in) :: filename_out

    integer :: ncid, varid
    integer :: dim_hour, dim_lat, dim_lon, dim_freq, dim_pol
    ! integer :: pol_enum_id
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

    real(real32), dimension(:,:,:), allocatable :: upsampled_data_3d
    real(real32), dimension(:,:,:,:), allocatable :: upsampled_data_4d

    nc_version_full = trim(nf90_inq_libvers())
    ! Only use the first whitespace-delimited word, which is the version number
    nc_version = nc_version_full(1:index(nc_version_full, " "))

    call date_and_time(zone=time_zone, values=time_vals)
    write(timestamp, '(I4, "-", I2.2, "-", I2.2, " ", I2.2, ":", I2.2, ":", I2.2, A5)') &
         time_vals(1), time_vals(2), time_vals(3), time_vals(5), time_vals(6), time_vals(7), time_zone

    call find_month_day(scene%year, scene%doy, month, day)
    write(time_start, '(I4, "-", I2.2, "-", I2.2, " 00:00:00Z")') scene%year, month, day
    write(time_end, '(I4, "-", I2.2, "-", I2.2, " 23:59:59Z")') scene%year, month, day
    write(epoch, '(I4, "-", I2.2, "-", I2.2, " 00:00:00 +0000")') scene%year, month, day

    call get_command(length=cmdline_len)
    allocate(character(len=cmdline_len) :: cmdline)
    allocate(character(len=cmdline_len + 180) :: history)
    call get_command(cmdline)
    write(history, '(A24, " created: ", A)') timestamp, trim(cmdline)

    call handle_nc_err(nf90_create(filename_out, ior(NF90_NOCLOBBER, NF90_NETCDF4), ncid))

    ! Define global attributes
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "CF-1.8,ACDD-1.3"))
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "title", "WSF-M MWI standard scene"))
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
    call handle_nc_err(nf90_put_att(ncid, NF90_GLOBAL, "nominal_mwi_incidence_angle", EIA_NOMINAL))

    ! Define dimensions
    call handle_nc_err(nf90_def_dim(ncid, "hour", NUM_HR, dim_hour))
    call handle_nc_err(nf90_def_dim(ncid, "lat", NUM_LAT_HI, dim_lat))
    call handle_nc_err(nf90_def_dim(ncid, "lon", NUM_LON_HI, dim_lon))
    call handle_nc_err(nf90_def_dim(ncid, "freq", NUM_FREQ, dim_freq))
    call handle_nc_err(nf90_def_dim(ncid, "pol", 2, dim_pol))

    ! ! Define enum
    ! call handle_nc_err(nf90_def_enum(ncid, NF90_INT, "pol", pol_enum_id))
    ! call handle_nc_err(nf90_insert_enum(ncid, pol_enum_id, "v", 1))
    ! call handle_nc_err(nf90_insert_enum(ncid, pol_enum_id, "h", 2))

    ! Define and write coordinate variables and their attributes
    call handle_nc_err(nf90_def_var(ncid, "hour", NF90_INT, dim_hour, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%hour))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "time"))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "synoptic hour"))
    call handle_nc_err(nf90_put_att(ncid, varid, "axis", "T"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "hours since " // trim(epoch)))

    call handle_nc_err(nf90_def_var(ncid, "lat", NF90_DOUBLE, dim_lat, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%lat_hi))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "latitude"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "degrees_north"))
    call handle_nc_err(nf90_put_att(ncid, varid, "axis", "Y"))

    call handle_nc_err(nf90_def_var(ncid, "lon", NF90_DOUBLE, dim_lon, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%lon_hi))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "longitude"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "degrees_east"))
    call handle_nc_err(nf90_put_att(ncid, varid, "axis", "X"))

    call handle_nc_err(nf90_def_var(ncid, "freq", NF90_DOUBLE, dim_freq, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%freq))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "sensor_band_central_radiation_frequency"))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "frequency"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "GHz"))

    ! call handle_nc_err(nf90_def_var(ncid, "pol", pol_enum_id, dim_pol, varid))
    call handle_nc_err(nf90_def_var(ncid, "pol", NF90_INT, dim_pol, varid))
    call handle_nc_err(nf90_put_var(ncid, varid, [1, 2]))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "polarization"))

    ! Define and write the high-resolution datasets and their attributes
    call handle_nc_err(nf90_def_var(ncid, "wspd", NF90_FLOAT, [dim_lon, dim_lat, dim_hour], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%wspd))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "wind_speed"))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "wind speed"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "m s-1"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

    call handle_nc_err(nf90_def_var(ncid, "wdir", NF90_FLOAT, [dim_lon, dim_lat, dim_hour], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%wdir))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "wind_to_direction"))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "wind direction"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "degree"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))
    call handle_nc_err(nf90_put_att(ncid, varid, "comment", &
         "Oceanographic convention: degrees clockwise from North, pointing along mass flow"))

    call handle_nc_err(nf90_def_var(ncid, "tec", NF90_FLOAT, [dim_lon, dim_lat, dim_hour], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%tec))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "total electron count"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "tecu"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

    call handle_nc_err(nf90_def_var(ncid, "ice_conc", NF90_FLOAT, [dim_lon, dim_lat], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%ice_conc))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "sea_ice_area_fraction"))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "sea ice concentration"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

    call handle_nc_err(nf90_def_var(ncid, "ice_fraction_new", NF90_FLOAT, [dim_lon, dim_lat], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%ice_fraction_new))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "fraction new/thin ice"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

    call handle_nc_err(nf90_def_var(ncid, "ice_fraction_first", NF90_FLOAT, [dim_lon, dim_lat], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%ice_fraction_first))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "fraction first-year ice"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

    call handle_nc_err(nf90_def_var(ncid, "ice_fraction_multi", NF90_FLOAT, [dim_lon, dim_lat], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%ice_fraction_multi))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "fraction multi-year ice"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

    call handle_nc_err(nf90_def_var(ncid, "ice_emissivity", NF90_FLOAT, [dim_lon, dim_lat, dim_freq, dim_pol], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, scene%ice_emissivity))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "sea ice emissivity"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

    ! Define and write low-resolution datasets, interpolating to the high-resolution grid
    allocate(upsampled_data_3d(NUM_LON_HI, NUM_LAT_HI, NUM_HR))

    call bilin_scene_surface(scene%sst, upsampled_data_3d)
    call handle_nc_err(nf90_def_var(ncid, "sst", NF90_FLOAT, [dim_lon, dim_lat, dim_hour], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, upsampled_data_3d))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "sea surface temperature"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "kelvin"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

    call bilin_scene_surface(scene%col_vapor, upsampled_data_3d)
    call handle_nc_err(nf90_def_var(ncid, "col_vapor", NF90_FLOAT, [dim_lon, dim_lat, dim_hour], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, upsampled_data_3d))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "atmosphere_mass_content_of_water_vapor"))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "columnar water vapor"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "kg m-2"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

    call bilin_scene_surface(scene%col_water, upsampled_data_3d)
    call handle_nc_err(nf90_def_var(ncid, "col_water", NF90_FLOAT, [dim_lon, dim_lat, dim_hour], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, upsampled_data_3d))
    call handle_nc_err(nf90_put_att(ncid, varid, "standard_name", "atmosphere_mass_content_of_cloud_liquid_water"))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "columnar liquid cloud content"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "kg m-2"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

    call bilin_scene_surface(scene%rain, upsampled_data_3d)
    call handle_nc_err(nf90_def_var(ncid, "rain", NF90_FLOAT, [dim_lon, dim_lat, dim_hour], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, upsampled_data_3d))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "rain rate"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "mm hr-1"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))

    deallocate(upsampled_data_3d)
    allocate(upsampled_data_4d(NUM_LON_HI, NUM_LAT_HI, NUM_FREQ, NUM_HR))

    call bilin_scene_profile(scene%tran, upsampled_data_4d)
    call handle_nc_err(nf90_def_var(ncid, "tran", NF90_FLOAT, [dim_lon, dim_lat, dim_freq, dim_hour], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, upsampled_data_4d))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "atmospheric transmissivity"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))
    call handle_nc_err(nf90_put_att(ncid, varid, "nominal_mwi_incidence_angle", EIA_NOMINAL))

    call bilin_scene_profile(scene%tb_up, upsampled_data_4d)
    call handle_nc_err(nf90_def_var(ncid, "tb_up", NF90_FLOAT, [dim_lon, dim_lat, dim_freq, dim_hour], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, upsampled_data_4d))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "upwelling brightness temperature"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "kelvin"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))
    call handle_nc_err(nf90_put_att(ncid, varid, "nominal_mwi_incidence_angle", EIA_NOMINAL))

    call bilin_scene_profile(scene%tb_down, upsampled_data_4d)
    call handle_nc_err(nf90_def_var(ncid, "tb_down", NF90_FLOAT, [dim_lon, dim_lat, dim_freq, dim_hour], varid, &
         deflate_level=2, shuffle=.true.))
    call handle_nc_err(nf90_put_var(ncid, varid, upsampled_data_4d))
    call handle_nc_err(nf90_put_att(ncid, varid, "long_name", "downwelling brightness temperature"))
    call handle_nc_err(nf90_put_att(ncid, varid, "units", "kelvin"))
    call handle_nc_err(nf90_put_att(ncid, varid, "coordinates", "lat lon"))
    call handle_nc_err(nf90_put_att(ncid, varid, "nominal_mwi_incidence_angle", EIA_NOMINAL))

    deallocate(upsampled_data_4d)
    call handle_nc_err(nf90_close(ncid))
  end subroutine write_scene

  ! This is a very simple netCDF error handler: if there's an error,
  ! then display it and exit with a nonzero code.
  subroutine handle_nc_err(status)
    integer, intent(in) :: status
    if (status /= NF90_NOERR) then
       write (ERROR_UNIT, *) "ERROR: ", nf90_strerror(status)
       error stop "NetCDF error"
    end if
  end subroutine handle_nc_err

  ! Bilinearly interpolate the surface data
  pure subroutine bilin_scene_surface(data_lo, data_hi)
    real(real32), dimension(NUM_LON_LO, NUM_LAT_LO, NUM_HR), intent(in) :: data_lo
    real(real32), dimension(NUM_LON_HI, NUM_LAT_HI, NUM_HR), intent(out) :: data_hi

    integer :: ilat, ilon, ihour, i
    real(real32), dimension(4) :: weights, points
    integer, dimension(4) :: ilons_lo, ilats_lo

    ! Iterate over the output pixels. The "hour" dimension is
    ! independent, so it's bilinear interpolation over lat/lon. For
    ! boundary conditions, the longitudes wrap, but the latitudes do
    ! not (it's treated as a reflection).
    do ilat = 1, NUM_LAT_HI
       do ilon = 1, NUM_LON_HI
          call lo_to_hi(ilon, ilat, ilons_lo, ilats_lo, weights)
          do ihour=1, NUM_HR
             do i = 1, 4
                points(i) = data_lo(ilons_lo(i), ilats_lo(i), ihour)
             end do
             data_hi(ilon, ilat, ihour) = dot_product(weights, points)
          end do
       end do
    end do

  end subroutine bilin_scene_surface

  ! Bilinearly interpolate the profile data
  pure subroutine bilin_scene_profile(data_lo, data_hi)
    real(real32), dimension(NUM_LON_LO, NUM_LAT_LO, NUM_FREQ, NUM_HR), intent(in) :: data_lo
    real(real32), dimension(NUM_LON_HI, NUM_LAT_HI, NUM_FREQ, NUM_HR), intent(out) :: data_hi

    integer :: ilat, ilon, ifreq, ihour, i
    real(real32), dimension(4) :: weights, points
    integer, dimension(4) :: ilons_lo, ilats_lo

    ! Iterate over the output pixels. The "hour" and "freq" dimensions
    ! are independent, so it's bilinear interpolation over lat/lon.
    ! For boundary conditions, the longitudes wrap, but the latitudes
    ! do not (it's treated as a reflection).
    do ilat = 1, NUM_LAT_HI
       do ilon = 1, NUM_LON_HI
          call lo_to_hi(ilon, ilat, ilons_lo, ilats_lo, weights)
          do ihour=1, NUM_HR
             do ifreq=1, NUM_FREQ
                do i = 1, 4
                   points(i) = data_lo(ilons_lo(i), ilats_lo(i), ifreq, ihour)
                end do
                data_hi(ilon, ilat, ifreq, ihour) = dot_product(weights, points)
             end do
          end do
       end do
    end do

  end subroutine bilin_scene_profile

  ! Convert indices of high-resolution data to low-resolution data for
  ! bilinear interpolation
  !
  ! ilon_hi, ilat_hi: the indices on the high-resolution grid
  !
  ! ilats_lo, ilons_hi: the four neighboring indices on the
  ! low-resolution grid (note that longitudes wrap but latitudes
  ! reflect)
  !
  ! weights: the corresponding weight for each point. The weights sum
  ! to 1.
  !
  pure subroutine lo_to_hi(ilon_hi, ilat_hi, ilons_lo, ilats_lo, weights)
    integer, intent(in) :: ilon_hi, ilat_hi
    integer, dimension(4), intent(out) :: ilons_lo, ilats_lo
    real(real32), dimension(4), intent(out) :: weights

    real(real64) :: x_point, y_point
    ! real(real64), dimension(NUM_LAT_LO), parameter :: lat_lo = [(LAT0 + DLAT_LO * i, i = 0, NUM_LAT_LO - 1)]
    ! real(real64), dimension(NUM_LON_LO), parameter :: lon_lo = [(LON0 + DLON_LO * i, i = 0, NUM_LON_LO - 1)]
    ! real(real64), dimension(NUM_LAT_HI), parameter :: lat_hi = [(LAT0 + DLAT_HI * i, i = 0, NUM_LAT_HI - 1)]
    ! real(real64), dimension(NUM_LON_HI), parameter :: lon_hi = [(LON0 + DLON_HI * i, i = 0, NUM_LON_HI - 1)]

    ! Convert from high-resolution lat/lon to floating-point
    ! low-resolution bin "coordinates". As a consequence, the "units"
    ! in this system are 1. So the weights below do not need to be
    ! divided by the bin spacing (since it's 1). Also, since the two
    ! grids are nested (same LAT0/LON0), there is a special case to
    ! consider: if the low-resolution lat/lon is exactly the same as
    ! the high resolution lat/lon, then the weight array only has 1
    ! non-zero value. (This happens when floor(x) = ceiling(x), or in
    ! other words there is no fractional part to x.)

    x_point = (ilon_hi - 1) * DLON_HI / DLON_LO
    y_point = (ilat_hi - 1) * DLAT_HI / DLAT_LO

    ! Extract the neighboring bin indices (0-based indexing!) and
    ! weights
    ilons_lo(1) = floor(x_point)
    ilats_lo(1) = floor(y_point)

    ilons_lo(2) = floor(x_point) + 1
    ilats_lo(2) = floor(y_point)

    ilons_lo(3) = floor(x_point)
    ilats_lo(3) = floor(y_point) + 1

    ilons_lo(4) = floor(x_point) + 1
    ilats_lo(4) = floor(y_point) + 1

    weights(4) = real((x_point - ilons_lo(1)) * (y_point - ilats_lo(1)), real32)
    weights(3) = real((ilons_lo(2) - x_point) * (y_point - ilats_lo(2)), real32)
    weights(2) = real((x_point - ilons_lo(3)) * (ilats_lo(3) - y_point), real32)
    weights(1) = real((ilons_lo(4) - x_point) * (ilats_lo(4) - y_point), real32)

    ! Convert indexing to 1-based
    ilons_lo(:) = ilons_lo(:) + 1
    ilats_lo(:) = ilats_lo(:) + 1

    ! Wrap longitudes
    where (ilons_lo > NUM_LON_LO)
       ilons_lo = ilons_lo - NUM_LON_LO
    end where

    ! Reflect latitudes
    where (ilats_lo == NUM_LAT_LO + 1)
       ilats_lo = NUM_LAT_LO
    end where
  end subroutine lo_to_hi

  ! Convert liquid cloud content to equivalent rain rate
  pure elemental function cld_to_rain(sst, cld)
    real(real32), intent(in) :: sst ! kelvin
    real(real32), intent(in) :: cld ! mm
    real(real32) :: cld_to_rain ! mm/h

    real(real32) :: ht, rr
    real(real32) :: sst_c

    ! Internally, SST is in degrees Celsius
    sst_c = sst - 273.15

    if (sst_c < 0.) then
       ht=0.46
    else if (sst_c > 30.) then
       ht=5.26
    else
       ht=0.46 + 0.16*sst_c
    end if

    if (cld <= 0.18) then
       rr=0.
    else
       rr=(cld-0.18)/0.18
       rr=rr*rr
    endif
    cld_to_rain=rr/ht
  end function cld_to_rain

  ! Apply the RTM to obtain atmospheric terms
  !
  ! This is basically the same thing as atmos_tables::get_atm but
  ! without the lat/lon dependence.
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

end module standard_scene
