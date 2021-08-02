!     1/5/2016 changed 2/20/2016.  water vapor normalized to pwat

!     icld = 0: no clouds
!     icld = 1: ncep cloud density, only liquid part
!

module atmos_tables
  use, intrinsic :: iso_fortran_env, only: real32, int32
  use find_abs_coeff, only: fdabscoeff
  use columnar_int, only: column
  use cloud_abs, only: fdcldabs
  use ncep, only: findncep
  use trig_degrees, only: cosd
  implicit none
  private
  public :: get_atm, atm_tran

contains
  subroutine get_atm(ncep_cloud,year,month,day,hour,lat,lon,nfreq,eia,freq, &
       surtep,vap,cld,pwat,cwat,tbup,tbdw,tran)
    logical, intent(in) :: ncep_cloud
    integer(int32), intent(in) :: year,month,day,hour
    real(real32), intent(in) :: lat, lon
    integer(int32), intent(in) :: nfreq

    real(real32), intent(in) :: eia
    real(real32), dimension(:), intent(in) :: freq
    real(real32), intent(out) :: surtep,vap,cld,pwat,cwat
    real(real32), dimension(:), intent(out) :: tbup,tbdw,tran

    integer(int32), parameter :: nmax = 26
    integer(int32) :: ipr,ibegin
    integer(int32) :: ifreq
    real(real32), dimension(0:nmax) :: t,p,pv,z,rhov,rhol,rhol0
    real(real32), dimension(0:nmax) :: abh2o,abo2,abcld,tabs
    real(real32) :: ao,av,al
    real(real32) :: pv_fix

    !     ncep parameters
    call findncep(year,month,day,hour,lat,lon, vap,cld,pwat,cwat,p,t,pv,rhov,rhol0,z,ibegin)

    surtep=t(0)

    if (ncep_cloud) then
       rhol= rhol0
    else
       rhol= 0.0
    endif

    do ifreq=1,nfreq
       do ipr  = ibegin,nmax

          !x    call fdabscoeff(freq(ifreq),p(ipr),t(ipr),pv(ipr),  abh2o(ipr),abo2(ipr)) ! neper/km
          !     for 2002-2007, there is a big dif between vap and pwat, with vap being near 36 mm for amazon and pwat being 50 mm
          !     in 2008 vap increases and is very similar to pwat.  I do not understand why vap changes.
          !     to make things consistent, here is normalized to pwat
          pv_fix=(pwat/vap)*pv(ipr)

          call fdabscoeff(freq(ifreq),p(ipr),t(ipr),pv_fix,  abh2o(ipr),abo2(ipr)) ! neper/km

          if (rhol(ipr) > 1.0e-7) then
             call fdcldabs(freq(ifreq),t(ipr),rhol(ipr),   abcld(ipr)) ! nepers/km
          else
             abcld(ipr) = 0.0
          endif
          ! Convert Np/km to Np/m
          abh2o(ipr)=abh2o(ipr) * 1e-3
          abo2(ipr)=abo2(ipr) * 1e-3
          abcld(ipr)=abcld(ipr) * 1e-3
       enddo ! ipr

       !     vertical integrals
       call column(nmax-ibegin,z(ibegin:nmax),abh2o(ibegin:nmax),2, av)
       call column(nmax-ibegin,z(ibegin:nmax), abo2(ibegin:nmax),2, ao)
       call column(nmax-ibegin,z(ibegin:nmax),abcld(ibegin:nmax),1, al)

       !     total absorption
       tabs(ibegin:nmax) = abh2o(ibegin:nmax) + abo2(ibegin:nmax) + abcld(ibegin:nmax)

       call atm_tran(nmax-ibegin,eia,t(ibegin:nmax),z(ibegin:nmax),tabs(ibegin:nmax), &
            tran(ifreq),tbdw(ifreq),tbup(ifreq))

    enddo  !ifreq
  end subroutine get_atm


  !   compute atmospheric downwelling and upwelling brightness temperatures
  !   and upward transmittance at each pressure level (altitude)

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !   input:
  !     nlev           number of atmosphere levels
  !     tht            earth incidence angle [in deg]
  !     tabs(0:nlev)   atmosphric absorptrion coefficients [nepers/m]
  !     t(0:nlev)      temperature profile[in k]
  !   z(0:nlev)      elevation (m)


  !     output:
  !     tran          total atmospheric transmission
  !     tbdw          downwelling brightness temperature t_bd [in k]
  !     tbup          upwelling   brightness temperature t_bu [in k]
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  pure subroutine atm_tran(nlev,tht,t,z,tabs,  tran,tbdw,tbup)
    integer(int32), intent(in) :: nlev
    real(real32), intent(in) :: tht
    real(real32), dimension(0:nlev), intent(in) :: t, z, tabs
    real(real32), intent(out) :: tran, tbdw, tbup

    ! real(real32), parameter :: re=6378.135
    real(real32), parameter :: delta=0.00035

    integer(int32) :: i
    real(real32) :: opacty(nlev),tavg(nlev),ems(nlev)
    real(real32) :: sumop, sumdw, sumup, tbavg, dsdh

    dsdh = (1.0+delta)/sqrt(cosd(tht)**2 + delta*(2+delta))

    do i=1,nlev
       opacty(i)=-dsdh*0.5*(tabs(i-1)+tabs(i))*(z(i)-z(i-1))
       tavg(i)  =0.5*(t(i-1)+t(i))
       ems(i)   =1.-exp(opacty(i))
    enddo

    sumop=0
    sumdw=0
    do i=1,nlev
       sumdw=sumdw+(tavg(i)-t(1))*ems(i)*exp(sumop)
       sumop=sumop+opacty(i)
    enddo

    sumop=0
    sumup=0.
    do i=nlev,1,-1
       sumup=sumup+(tavg(i)-t(1))*ems(i)*exp(sumop)
       sumop=sumop+opacty(i)
    enddo

    tran=exp(sumop)
    tbavg=(1.-tran)*t(1)
    tbdw=tbavg+sumdw
    tbup=tbavg+sumup
  end subroutine atm_tran
end module atmos_tables
