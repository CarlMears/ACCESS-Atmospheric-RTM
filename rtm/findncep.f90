!     12/10/2016 changed 8/22/2016.  this subroutine was modified so that it can do two years at once
!     without excessive reads of the ncep data.  for the first call use standard lyear.
!     for the second call append a minus sign to lyear.
!     this change has not effect on previous or future usage of the routine

!     same as 'o:\aquarius\mk_atmos_tables\findncep.f' dated march 4 2013 except this version does not call
!     call sleepqq(10000) !wait 10 seconds for possible file download to complete


!     03/04/2013
!   changed path from Nasserver5 to ops1p
!
!     01/05/2006
!     was: if t(ipr) < 100 k : ierr=-10, bail out
!     changed: if t(ipr) < 100 k, set rhocwat=0.0
!
!     01/28/2003
!     rhol = liquid cloud water density
!
!     ordered ncep profiles, surtep, wind, pwat ,colvap
!
!     input:
!     nmax: maximum number of levels
!     year:   long year (e.g. 1999)
!     mon:    month
!     day:    day of the month
!     hour: uct hour
!     lat: latitude  [between -90 and 90]
!     lon: longitude   [between 0 and 360]


!    output:
!     colvap: ncep water vapor integrated [mm]
!     colwat: ncep columnar liquid cloud water integrated [mm]
!     pwat:   ncep value for precipitable water [mm]
!     cwat:   ncep value for columnar cloud water   (liquid + ice) [mm]
!     p: air pressure profile      [mb] (0:nmax)   ordered
!     t: air temperature profile [k]  (0:nmax)
!     z: elevation profile [m]
!     pv: water vapor pressure profile [mb]  (0:nmax)
!     rhov: water vapor density [g /m**3]
!     rhol: liquid water density [g /m**3]
!     ibegin: index of surface level

module ncep
  use, intrinsic :: iso_fortran_env, only: real32, real64, int32, ERROR_UNIT
  use wvap_convert, only: goff_gratch_vap
  use columnar_int, only: column
  implicit none
  private
  public :: findncep

  integer(int32), parameter :: NMAX = 26

  ! Hourly NCEP data
  type NcepHourly
     integer(int32) :: year, month, day, hour
     ! These  are dimensioned as (360, 181, 0:NMAX)
     real(real32), dimension(:,:,:), allocatable :: t_map, hgt_map, rh_map,clwmr_map
     ! These are dimensioned as (360, 181)
     real(real32), dimension(:,:), allocatable :: p_sfc_map,pwat_map,cwat_map
  end type NcepHourly

contains

  subroutine init_cache(cache)
    type(NcepHourly), intent(inout) :: cache

    cache%year = -999
    cache%month = -999
    cache%day = -999
    cache%hour = -999
    allocate(cache%t_map(360,181,0:NMAX), &
         cache%hgt_map(360,181,0:NMAX), &
         cache%rh_map(360,181,0:NMAX), &
         cache%clwmr_map(360,181,0:NMAX), &
         cache%p_sfc_map(360, 181), &
         cache%pwat_map(360, 181), &
         cache%cwat_map(360, 181))
  end subroutine init_cache

  subroutine findncep(year,month,day,hour,lat,lon,  colvap,colwat,pwat,cwat,p,t,pv,rhov,rhol,z,ibegin)
    integer(int32), intent(in) :: year, month, day, hour
    real(real32), intent(in) :: lat, lon
    integer(int32), intent(out) :: ibegin
    real(real32), dimension(0:NMAX), intent(out) :: t, p
    real(real32), dimension(0:NMAX), intent(out) :: z,pv,rhov,rhol
    real(real32), intent(out) :: pwat,cwat
    real(real32), intent(out) :: colvap,colwat

    real(real32), dimension(0:NMAX) :: hgt,rh,clwmr
    real(real32), dimension(0:NMAX) :: rhocwat,rhoi

    real(real32) :: p_sfc
    real(real32) :: p0,rhoair
    integer(int32) :: ipr
    real(real32), parameter :: r_e = 6371000, rd=287.05, epsilon = 0.622

    logical, save :: first = .true.
    ! Two cache types are used for year-wrappping cases
    type(NcepHourly), save :: ncep_positive, ncep_negative

    p(0:NMAX) = [ 0., &
         1000.,975.,950.,925.,900.,850.,800.,750.,700.,650.,600., &
         550.,500.,450.,400.,350.,300.,250.,200.,150.,100., &
         70., 50., 30., 20. ,10. ]  ! [hPa]

    if (first) then
       call init_cache(ncep_positive)
       call init_cache(ncep_negative)
       first = .false.
    end if

    if (year > 0) then
       call fnd_ncep(ncep_positive, year,month,day,hour,lat,lon,  t,hgt,rh,clwmr,p_sfc,pwat,cwat)
    else
       call fnd_ncep(ncep_negative, -year,month,day,hour,lat,lon,  t,hgt,rh,clwmr,p_sfc,pwat,cwat)
    end if

    p0 = 0.01*p_sfc    ! hpa

    ibegin=-1
    do ipr=1,NMAX
       if (p(ipr).le.p0) then  !le is used rather than lt to be consistent with old sorting method
          ibegin=ipr-1
          exit
       end if
    end do
    if (ibegin.eq.-1) error stop 'ibegin.eq.-1 in findncep'

    p(ibegin) =p0
    t(ibegin) =t(0)
    hgt(ibegin) =hgt(0)
    rh(ibegin) =rh(0)
    clwmr(ibegin) =clwmr(0)

    !     write(*,*) t(0),hgt(0),rh(0),clwmr(0)
    !     write(*,*) ibegin,p(ibegin),p0
    !   do ipr=0,nmax
    !   write(*,'(i5,12f10.2)') ipr,p(ipr),t(ipr),hgt(ipr),rh(ipr),clwmr(ipr)
    !     end do
    !   stop 'aaa'

    z = hgt * r_e / (r_e - hgt)
    if (z(ibegin) .ge. z(ibegin+1)) z(ibegin) = z(ibegin+1) - 0.1

    !   transform rh -> water vapor pressure and density
    call goff_gratch_vap(t,rh,p,  pv,rhov)

    !     convert ncep clwmr to rhocwat
    clwmr(ibegin) = clwmr(ibegin+1) ! clwmr at surface

    do ipr = 0,NMAX
       if (t(ipr) >= 100 .and. p(ipr) > 0) then
          rhoair = p(ipr)/(rd*t(ipr)) ! density of dry air
          rhoair = rhoair*(1.0 - (1.0-epsilon)*pv(ipr)/p(ipr))
          ! density of moist air (without cloud), unusual ncep definition for mixing ratio

          rhocwat(ipr) = 1.0e5*clwmr(ipr)*rhoair    ! [g/m**3]
          call cldwatice(t(ipr),rhocwat(ipr),  rhol(ipr),rhoi(ipr))
       else ! unphysical temperature
          rhocwat(ipr)=0.0
          rhol(ipr) = 0.0
          rhoi(ipr) = 0.0
       end if

    end do

    where (rhov < 0) rhov=0.0
    where (rhol < 0) rhol=0.0

    call column(NMAX-ibegin,z(ibegin:NMAX),rhov(ibegin:NMAX),2, colvap)
    call column(NMAX-ibegin,z(ibegin:NMAX),rhol(ibegin:NMAX),1, colwat)

    colvap = colvap *1.e-3   ! g/m**2 -> mm = kg/m**2
    colwat = colwat *1.e-3   ! g/m**2 -> mm = kg/m**2
  end subroutine findncep

  pure subroutine cldwatice(t,rhoc,  rhol,rhoi)
    !     t: temperature [k]
    !     rhoc:   cloud water density [arbitrary unit]
    !
    !     rhol : liquid part [same unit as rhoc]
    !     rhoi : ice    part [same unit as rhoc]
    real(real32), intent(in) :: t,rhoc
    real(real32), intent(out) :: rhol,rhoi

    real(real32), parameter :: twat = 273.15, tice = 253.15

    !     water or ice
    if (t >= twat) then
       rhol = rhoc ! all liquid water
       rhoi = 0.0
    else if (t <= tice) then
       rhol = 0.0       ! all ice
       rhoi = rhoc
    else
       rhol = rhoc*(t-tice)/(twat-tice)
       rhoi = rhoc - rhol
       ! linear temperature interpolation in between
    end if
  end subroutine cldwatice

  !    01/23/2006
  !    kyle changed subroutine yread
  !    read fileheader between each profile level
  !
  !
  !    ncep variables
  !    1 deg resolution
  !    4 times daily: 00z 06z 12z 18z
  !
  !    thomas meissner : apr 1999
  !
  !    needs to be linked with time_routines.f
  !
  !
  !         varaible    level (subfolder)
  !           sst               sfc    (0)
  !           t                 sfc    (0) 2 m above ground
  !           rh                sfc    (0)
  !           hgt               sfc    (0)
  !           p                 sfc    (0)
  !           pwat              col
  !           cwat              col
  !           t                 1:nmax  (prf)
  !           rh                1:nmax   (prf)
  !           hgt               1:nmax  (prf)
  !           clwmr             1:nmax  (prf)
  !
  !    named as year_month_day_xxz.dat :
  !    binary file:
  !    year month day hour (integer(4))
  !    real*4 array
  !    grid: lat from +90 to -90 by -1 deg
  !          lon from 0 to 360 by 1 deg
  !    total size 4*4 + 360*181*4 = 260656
  !
  !    input: lyear  (year, long form, e.g. 2005)   integer
  !           imon   (month, 1-12)                  integer
  !           iday (day of month, 1-31)             integer
  !           hour   (hour of day)                  integer
  !           xlat   (latitude)                     real(4)
  !           xlon   (longitude)                    real(4)
  !
  !    output: var (interpolated variable) real(4)
  !

  subroutine fnd_ncep(cache, year,month,day,hour,lat,lon,  t,hgt,rh,clwmr,p_sfc,pwat,cwat)
    type(NcepHourly), intent(inout) :: cache
    integer(int32), intent(in) :: year, month, day, hour
    real(real32), intent(in) :: lat, lon

    real(real32), dimension(0:NMAX), intent(out) :: t,hgt,rh,clwmr
    real(real32), intent(out) :: p_sfc,pwat,cwat

    integer(int32) :: x0, x1, y0, y1
    real(real32) :: xd, yd
    integer :: ilevel

    if (year /= cache%year .or. month /= cache%month .or. day /= cache%day .or. hour /= cache%hour) then
       cache%year = year
       cache%month = month
       cache%day = day
       cache%hour = hour
       call yread(year,month,day,hour, &
            cache%t_map,cache%hgt_map,cache%rh_map,cache%clwmr_map,&
            cache%p_sfc_map,cache%pwat_map,cache%cwat_map)
    end if

    ! Bilinearly interpolate over lat/lon. The notation below is from:
    ! https://en.wikipedia.org/wiki/Trilinear_interpolation
    ! https://en.wikipedia.org/wiki/Linear_interpolation

    call lookup_bin_inds_ncep(lat, lon, x0, x1, xd, y0, y1, yd)

    ! For each dataset, find the 4 neighboring values or profiles and
    ! then apply the bilinear interpolation
    do concurrent (ilevel=0:NMAX)
       t(ilevel) = interp_2d(cache%t_map(:,:,ilevel), x0, x1, xd, y0, y1, yd)
       hgt(ilevel) = interp_2d(cache%hgt_map(:,:,ilevel), x0, x1, xd, y0, y1, yd)
       rh(ilevel) = interp_2d(cache%rh_map(:,:,ilevel), x0, x1, xd, y0, y1, yd)
       clwmr(ilevel) = interp_2d(cache%clwmr_map(:,:,ilevel), x0, x1, xd, y0, y1, yd)
    end do

    p_sfc = interp_2d(cache%p_sfc_map, x0, x1, xd, y0, y1, yd)
    pwat = interp_2d(cache%pwat_map, x0, x1, xd, y0, y1, yd)
    cwat = interp_2d(cache%cwat_map, x0, x1, xd, y0, y1, yd)

  end subroutine fnd_ncep

  ! Find the indices of the neighboring bins that a value belongs to.
  !
  ! This is specific to the NCEP data where the maps are indexed using
  ! (longitude, latitude, level), where each 1-degree bin center
  ! ranges from 0 to 359, for longitude, and from +90 to -90, for
  ! latitude.
  !
  ! Since each bin is 1 degree wide, this permits a constant-time
  ! lookup to return the left and right bin indices for a given value,
  ! as well as the fractional contribution between the two bins.
  !
  ! x0, x1, and xd correspond to longitude
  ! y0, y1, and yd correspond to latitutde
  !
  ! x0 and x1 vary between 1 and 360; y0 and y1 vary between 1 and
  ! 181.
  pure subroutine lookup_bin_inds_ncep(lat, lon, &
       x0, x1, xd, y0, y1, yd)
    real(real32), intent(in) :: lat, lon
    integer(int32), intent(out) :: x0, x1, y0, y1
    real(real32), intent(out) :: xd, yd

    ! The latitude bin centers vary from +90 to -90
    y0 = floor(90 - lat) + 1
    y1 = y0 + 1
    yd = (90 - lat) - floor(90 - lat)

    ! The longitude bin centers vary from 0 to 359
    x0 = floor(lon) + 1
    x1 = x0 + 1
    xd = lon - floor(lon)

    if (x0 == 361) then
       ! If the longitude is between 359.5 and 360, then x0 is 361, which
       ! is out of range, so it wraps around
       x0 = 1
       x1 = 2
    else if (x1 == 361) then
       ! If the longitude is between 358.5 and 359.5, then although x0
       ! is 360, which is in range, x1 is 361, which is not. Since
       ! it's periodic, it's set to 1.
       x1 = 1
    end if
  end subroutine lookup_bin_inds_ncep

  subroutine yread(year,month,day,hour,  t,hgt,rh,clwmr,p_sfc,pwat,cwat)
    integer(int32), intent(in) :: year,month,day,hour
    real(real32), dimension(360,181,0:NMAX), intent(out) :: t,hgt,rh,clwmr
    real(real32), dimension(360,181), intent(out) :: p_sfc,pwat,cwat

    integer(int32), parameter :: NRH = 21
    character(len=100) :: filename

    integer(int32) :: kyear,kmon,kday,khour
    integer(int32) :: ilevel,icase

    integer :: ioerr, lu
    character(len=80) :: iomsg

9001 format('/mnt/ops1p/n/data/model/NCEP/PROFILES/TMP/Y',i4.4,'/',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
9002 format('/mnt/ops1p/n/data/model/NCEP/SURFACE/TMP_2m/Y',i4.4,'/',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
9003 format('/mnt/ops1p/n/data/model/NCEP/PROFILES/RH/Y',i4.4,'/',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
9004 format('/mnt/ops1p/n/data/model/NCEP/SURFACE/RH/Y',i4.4,'/',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
9005 format('/mnt/ops1p/n/data/model/NCEP/PROFILES/HGT/Y',i4.4,'/',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
9006 format('/mnt/ops1p/n/data/model/NCEP/SURFACE/HGT/Y',i4.4,'/',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
9007 format('/mnt/ops1p/n/data/model/NCEP/SURFACE/PRES/Y',i4.4,'/',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
9008 format('/mnt/ops1p/n/data/model/NCEP/SURFACE/PWAT/Y',i4.4,'/',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
9009 format('/mnt/ops1p/n/data/model/NCEP/PROFILES/CLWMR/Y',i4.4,'/',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')
9010 format('/mnt/ops1p/n/data/model/NCEP/SURFACE/CWAT/Y',i4.4,'/',i4.4,'_',i2.2,'_',i2.2,'_',i2.2,'z.dat')

    cases: do icase=1,10
       select case (icase)
       case (1)
          write(filename,9001) year,year,month,day,hour
       case (2)
          write(filename,9002) year,year,month,day,hour
       case (3)
          write(filename,9003) year,year,month,day,hour
       case (4)
          write(filename,9004) year,year,month,day,hour
       case (5)
          write(filename,9005) year,year,month,day,hour
       case (6)
          write(filename,9006) year,year,month,day,hour
       case (7)
          write(filename,9007) year,year,month,day,hour
       case (8)
          write(filename,9008) year,year,month,day,hour
       case (9)
          write(filename,9009) year,year,month,day,hour
       case (10)
          write(filename,9010) year,year,month,day,hour
       end select

       open(newunit=lu, file=filename, status='old', action='read', &
            form='unformatted', access='stream', iostat=ioerr, iomsg=iomsg)
       if (ioerr /= 0) then
          write(ERROR_UNIT, *) iomsg
          error stop 1
       end if

       read(lu) kyear, kmon, kday, khour
       if (kyear /= year .or. kmon /= month .or. kday /= day .or. khour /= hour) then
          write(*,*) filename
          write(*,*) kyear,kmon,kday,khour
          write(*,*) year,month,day,hour
          error stop 'header error in yread'
       end if

       select case (icase)
       case (1)
          do ilevel=1,NMAX
             read(lu) t(:,:,ilevel)
          end do
       case (2)
          read(lu) t(:,:,0)
       case (3)
          do ilevel=1,NRH
             read(lu) rh(:,:,ilevel)
          end do
          rh(:,:,NRH+1:NMAX) = 0.0
       case (4)
          read(lu) rh(:,:,0)
       case (5)
          do ilevel=1,NMAX
             read(lu) hgt(:,:,ilevel)
          end do
       case (6)
          read(lu) hgt(:,:,0)
       case (7)
          read(lu) p_sfc
       case (8)
          read(lu) pwat
       case (9)
          do ilevel=1,NRH
             read(lu) clwmr(:,:,ilevel)
          end do
          clwmr(:,:,NRH+1:NMAX) = 0.0
          clwmr(:,:,0) = 0.0
       case (10)
          read(lu) cwat
       end select

       close(lu)

    end do cases
  end subroutine yread

  ! Linear interpolation
  !
  ! t: ranges between 0 and 1
  ! v0: value when t == 0
  ! v1: value when t == 1
  pure function lerp(v0, v1, t)
    real(real32), intent(in) :: v0, v1, t
    real(real32) :: lerp
    lerp = (1 - t) * v0 + t * v1
  end function lerp

  ! 2D interpolation
  pure function interp_2d(map, x0, x1, xd, y0, y1, yd)
    real(real32), dimension(:, :), intent(in) :: map
    integer(int32), intent(in) :: x0, x1, y0, y1
    real(real32), intent(in) :: xd, yd
    real(real32) :: interp_2d

    real(real32) :: c00, c01, c10, c11, c0, c1

    c00 = map(x0,y0)
    c01 = map(x0,y1)
    c10 = map(x1,y0)
    c11 = map(x1,y1)
    c0 = lerp(c00, c10, xd)
    c1 = lerp(c01, c11, xd)
    interp_2d = lerp(c0, c1, yd)
  end function interp_2d

end module ncep
