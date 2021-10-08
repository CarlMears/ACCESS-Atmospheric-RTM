!     july 17 2015  changed oct 11 2016.  changes made to vapor model. see 'memo4.txt'
!     july 17 2015 veresion save in O:\skytemp3\version_07172015

!     june 11 2015 changed july 17 2015.  oxyopc adjustment above 37 ghz changed slightly. see 'O:\gmi\abs_cal\memo20.txt'

!     this module contains the oxygen absoprtion routine, the vapor absorption routine and fdaray
!     there were used for the april-june 2009 skytemp update.  see 'memo8.txt'

module atms_abs_routines
  use, intrinsic :: iso_fortran_env, only: real32, int32, real64
  use dielectric_meissner, only: meissner
  use trig_degrees, only: cosd
  implicit none
  private
  public :: atm_tran, fdabscoeff, fdcldabs

contains

  !     ===========================================================
  !     ====== modified version of Liebe 1992 oxygen model ========
  !     ===========================================================

  !     This is from Atmospheric 60-GHz Oxygen Spectrum:.. Liebe, Rosenkranz, Hufford, 1992
  !     coded: June 2009 1992 by f.wentz
  !           inputs: t, temperature (K)
  !                   p, total pressure (mb)
  !                   pv, water vapor pressure (mb)
  !                   freq, frequency (GHz)
  !           output: gamoxy, oxygen absorption coefficient (dB/km)


  !     It is the same as fdabsoxy_1989 except for the a5 and a6 coefs for finding delta have different values.
  !     Also this 1992 versions says delta is proprotional to total pressure p rather than pdry.
  !     Also note in this routine, the 1.e-3 scaling is done in the startup block.

  !     compared to abo2_rk (the code Rosenkranz sent us), this routine gives very similar results if you set the apterm to 0.
  !     for my freqs 6-85 ghz, the largest dif was 0.0003 at the coldest vapor at 85.5 GHz.
  !     the apterm adds 0.003 at 85 ghz
  !     Apart from the apterm, you can essentially say fdabsoxy_1992 and abo2_rk are the same for my purposes.

  !     this routine has been modified in the following ways (see 'memo8.txt')
  !     1.  non-resonance continuum temperature coef changed from 0.8 to 1.5
  !     2.  a p*p continuum was added
  !     these modifications were done June 22 2009

  subroutine fdabsoxy_1992_modified(p,t,pv,freq, gamoxy)
    real(real32), intent(in) :: p,t,pv,freq
    real(real32), intent(out) :: gamoxy

    integer(int32), parameter :: nlines=44
    integer(int32) :: i
    real(real32) :: tht,pwet,pdry,ga,gasq,delta,rnuneg,rnupos,ff,zterm,apterm,sftot,xterm

    logical, save :: first = .true.
    real(real32), dimension(6, nlines), save :: h
    real(real32), dimension(nlines), save :: f0, a1, a2, a3, a4, a5, a6

    real(real64) :: sum

    data a4/38*0., 6*0.6/

    !          freq          a1      a2       a3       a5          a6
    data h/ &
         50.474238,    0.94e-6,  9.694,  8.60e-3,  0.210,  0.685, &
         50.987749,    2.46e-6,  8.694,  8.70e-3,  0.190,  0.680, &
         51.503350,    6.08e-6,  7.744,  8.90e-3,  0.171,  0.673, &
         52.021410,   14.14e-6,  6.844,  9.20e-3,  0.144,  0.664, &
         52.542394,   31.02e-6,  6.004,  9.40e-3,  0.118,  0.653, &
         53.066907,   64.10e-6,  5.224,  9.70e-3,  0.114,  0.621, &
         53.595749,  124.70e-6,  4.484, 10.00e-3,  0.200,  0.508, &
         54.130000,  228.00e-6,  3.814, 10.20e-3,  0.291,  0.375, &
         54.671159,  391.80e-6,  3.194, 10.50e-3,  0.325,  0.265, &
         55.221367,  631.60e-6,  2.624, 10.79e-3,  0.224,  0.295, &
         55.783802,  953.50e-6,  2.119, 11.10e-3, -0.144,  0.613, &
         56.264775,  548.90e-6,  0.015, 16.46e-3,  0.339, -0.098, &
         56.363389, 1344.00e-6,  1.660, 11.44e-3, -0.258,  0.655, &
         56.968206, 1763.00e-6,  1.260, 11.81e-3, -0.362,  0.645, &
         57.612484, 2141.00e-6,  0.915, 12.21e-3, -0.533,  0.606, &
         58.323877, 2386.00e-6,  0.626, 12.66e-3, -0.178,  0.044, &
         58.446590, 1457.00e-6,  0.084, 14.49e-3,  0.650, -0.127, &
         59.164207, 2404.00e-6,  0.391, 13.19e-3, -0.628,  0.231, &
         59.590983, 2112.00e-6,  0.212, 13.60e-3,  0.665, -0.078, &
         60.306061, 2124.00e-6,  0.212, 13.82e-3, -0.613,  0.070, &
         60.434776, 2461.00e-6,  0.391, 12.97e-3,  0.606, -0.282, &
         61.150560, 2504.00e-6,  0.626, 12.48e-3,  0.090, -0.058, &
         61.800154, 2298.00e-6,  0.915, 12.07e-3,  0.496, -0.662, &
         62.411215, 1933.00e-6,  1.260, 11.71e-3,  0.313, -0.676, &
         62.486260, 1517.00e-6,  0.083, 14.68e-3, -0.433,  0.084, &
         62.997977, 1503.00e-6,  1.665, 11.39e-3,  0.208, -0.668, &
         63.568518, 1087.00e-6,  2.115, 11.08e-3,  0.094, -0.614, &
         64.127767,  733.50e-6,  2.620, 10.78e-3, -0.270, -0.289, &
         64.678903,  463.50e-6,  3.195, 10.50e-3, -0.366, -0.259, &
         65.224071,  274.80e-6,  3.815, 10.20e-3, -0.326, -0.368, &
         65.764772,  153.00e-6,  4.485, 10.00e-3, -0.232, -0.500, &
         66.302091,   80.09e-6,  5.225,  9.70e-3, -0.146, -0.609, &
         66.836830,   39.46e-6,  6.005,  9.40e-3, -0.147, -0.639, &
         67.369598,   18.32e-6,  6.845,  9.20e-3, -0.174, -0.647, &
         67.900867,    8.01e-6,  7.745,  8.90e-3, -0.198, -0.655, &
         68.431005,    3.30e-6,  8.695,  8.70e-3, -0.210, -0.660, &
         68.960311,    1.28e-6,  9.695,  8.60e-3, -0.220, -0.665, &
         118.750343,  945.00e-6,  0.009, 16.30e-3, -0.031,  0.008, &
         368.498350,   67.90e-6,  0.049, 19.20e-3,  0.0,    0.0,   &
         424.763124,  638.00e-6,  0.044, 19.16e-3,  0.0,    0.0,   &
         487.249370,  235.00e-6,  0.049, 19.20e-3,  0.0,    0.0,   &
         715.393150,   99.60e-6,  0.145, 18.10e-3,  0.0,    0.0,   &
         773.839675,  671.00e-6,  0.130, 18.10e-3,  0.0,    0.0,   &
         834.145330,  180.00e-6,  0.147, 18.10e-3,  0.0,    0.0/

    if (first) then
       first=.false.
       f0(:)=h(1,:)
       a1(:)=h(2,:)/h(1,:)
       a2(:)=h(3,:)
       a3(:)=h(4,:)
       a5(:)=0.001*h(5,:)
       a6(:)=0.001*h(6,:)
    end if

    tht = 300/t
    pwet=0.1*pv
    pdry=0.1*p-pwet
    xterm=1-tht

    sum = 0.
    !$omp simd reduction(+:sum) private(ga, gasq, delta, rnuneg, rnupos, ff)
    do i=1, nlines
       ga = a3(i)*(pdry*tht**(0.8-a4(i)) + 1.1*tht*pwet)
       gasq=ga*ga
       delta=(a5(i) + a6(i)*tht)*p*tht**0.8
       rnuneg = f0(i)-freq
       rnupos = f0(i)+freq
       ff = (ga-rnuneg*delta)/(gasq+rnuneg**2) +  (ga-rnupos*delta)/(gasq+rnupos**2)
       sum = sum + ff*a1(i)*exp(a2(i)*xterm)
    end do
    if (sum < 0) sum=0

    !     add nonresonant contribution

    !     ga=5.6e-3*(pdry+1.1*pwet)*tht**0.8
    ga=5.6e-3*(pdry+1.1*pwet)*tht**1.5  !modification 1

    zterm=ga*(1.+(freq/ga)**2)
    apterm=1.4e-10*(1-1.2e-5*freq**1.5)*pdry*tht**1.5
    if (apterm < 0) apterm=0
    sftot=real(pdry*freq*tht**2 * (tht*sum + 6.14e-4/zterm + apterm), real32)

    gamoxy=0.1820*freq*sftot
    !x    if(freq.gt.37) gamoxy=gamoxy + 0.1820*43.e-10 *pdry**2*tht**3*(freq-37.)**1.7  !prior to 7/17/2015
    if (freq > 37) gamoxy=gamoxy + 0.1820*26.e-10 *pdry**2*tht**3*(freq-37.)**1.8  !implemented 7/17/2015.
  end subroutine fdabsoxy_1992_modified


  !     ====================================================================
  !     ========= modified version of Rosenkranz water vapor model =========
  !     ====================================================================


  ! purpose- compute absorption coef in atmosphere due to water vapor
  !
  !  calling sequence parameters-
  !    specifications
  !      name    units    i/o  descripton            valid range
  !      t       kelvin    i   temperature
  !      p       millibar  i   pressure              .1 to 1000
  !      f       GHz       i   frequency             0 to 800
  !      gamh2o  dB/km     o   absorption coefficient
  !
  !   references-
  !    p.w. rosenkranz, radio science v.33, pp.919-928 (1998).
  !
  !   line intensities selection threshold=
  !     half of continuum absorption at 1000 mb.
  !   widths measured at 22,183,380 ghz, others calculated.
  !     a.bauer et al.asa workshop (sept. 1989) (380ghz).
  !
  !   revision history-
  !    date- oct.6, 1988  p.w.rosenkranz - eqs as publ. in 1993.
  !          oct.4, 1995  pwr- use clough's definition of local line
  !                   contribution,  hitran intensities, add 7 lines.
  !          oct. 24, 95  pwr -add 1 line.
  !          july 7, 97   pwr -separate coeff. for self-broadening,
  !                       revised continuum.
  !          dec. 11, 98  pwr - added comments

  !     the routine is a modified version of abh2o_rk_reformat.
  !     this routine has been modified in the following three ways (see 'memo8.txt')
  !     1.  b1(1)=1.01*b1(1)  :22 ghz line strength increase slightly
  !     2.  22 ghz line shape below 22 ghz has been modified
  !     3.  foreign and self broadening continuum has been adjusted
  !     these modification were done June 22 2009

  subroutine abh2o_rk_modified(p,t,pv,freq,  gamh2o)
    real(real32), intent(in) :: p,t,pv,freq
    real(real32), intent(out) :: gamh2o

    integer(int32), parameter :: nlines=15

    integer(int32) :: i
    logical, save :: first = .true.
    real(real32) :: s,base
    real(real32) :: tht,pwet,pdry,ga,gasq,sftot,xterm,rnuneg,rnupos
    real(real64) :: sum

    real(real32) :: chi,chisq,freqsq,f0sq,u,ffac

    !     line frequencies:
    real(real32), dimension(nlines), parameter :: f0 = [ &
         22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508, &
         443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360, &
         620.7008, 752.0332, 916.1712]
    !     line intensities at 300k:
    real(real32), dimension(nlines), save :: b1 = [ &
         .1310e-13, .2273e-11, .8036e-13, .2694e-11, .2438e-10, &
         .2179e-11, .4624e-12, .2562e-10, .8369e-12, .3263e-11, .6659e-12, &
         .1531e-08, .1707e-10, .1011e-08, .4227e-10]
    !     t coeff. of intensities:
    real(real32), dimension(nlines), parameter :: b2 = [ &
         2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405, 3.597, 2.379, 2.852, .159, 2.391, .396, 1.441]
    !     air-broadened width parameters at 300k:
    real(real32), dimension(nlines), save :: b3 = [ &
         .0281, .0281, .023, .0278, .0287, .021, .0186, .0263, .0215, .0236, .026, .0321, .0244, .0306, .0267]
    !     self-broadened width parameters at 300k:
    real(real32), dimension(nlines), save :: b5 = [ &
         .1349, .1491, .108, .135, .1541, .090, .0788, .1275, .0983, .1095, .1313, .1320, .1140, .1253, .1275]
    !     t-exponent of air-broadening:
    real(real32), dimension(nlines), parameter :: b4 = [ &
         .69, .64, .67, .68, .54, .63, .60, .66, .66, .65, .69, .69, .71, .68, .70]
    !     t-exponent of self-broadening:
    real(real32), dimension(nlines), parameter :: b6 = [ &
         .61, .85, .54, .74, .89, .52, .50, .67, .65, .64, .72, 1.0, .68, .84, .78]
    !
    if (first) then
       first = .false.
       b1=1.8281089E+14*b1/f0**2
       b5=b5/b3  !convert b5 to Leibe notation
       !1    b1(1)=1.01*b1(1)  !this modification is no longer done
       b3(1)=b3(1)/1.040 !modification 1
    end if

    if (pv <= 0.) then
       gamh2o=0
       return
    end if

    pwet=0.1*pv
    pdry=0.1*p-pwet
    tht = 300./t
    xterm=1-tht
    freqsq=freq*freq

    sum = 0.
    do i = 1, nlines
       f0sq=f0(i)*f0(i)
       ga=b3(i)*(pdry*tht**b4(i) + b5(i)*pwet*tht**b6(i))
       gasq = ga*ga
       s = b1(i)*exp(b2(i)*xterm)
       rnuneg = f0(i)-freq
       rnupos = f0(i)+freq
       base = ga/(562500. + gasq)  !use clough's definition of local line contribution

       if (i /= 1) then
          if (abs(rnuneg) < 750) sum = sum + s*(ga/(gasq + rnuneg**2) - base)
          if (abs(rnupos) < 750) sum = sum + s*(ga/(gasq + rnupos**2) - base)
       else
          chi=0.07*ga   !modification 2

          if (freq < 19) then
             u=abs(freq-19.)/16.5
             if (u < 0) u=0
             if (u > 1) u=1
             chi=0.07*ga + 0.93*ga*u*u*(3-2*u)  !modification 2
          end if

          chisq=chi*chi
          sum=sum + s*2*((ga-chi)*freqsq + (ga+chi)*(f0sq+gasq-chisq))/((freqsq-f0sq-gasq+chisq)**2 + 4*freqsq*gasq)
       end if
    end do
    if (sum < 0) sum=0

    ffac=1
    if (freq < 90) ffac=ffac + 0.1*((90. - freq)/90.)**1.4
    !1    sftot=pwet*freq*tht**3.5*(sum +      1.1*1.2957246e-6*pdry/tht**0.5 + 0.425*(freq**0.10)*4.2952193e-5*pwet*tht**4) !previous version
    sftot=real(pwet*freq*tht**3.5*(sum + ffac*1.1*1.2957246e-6*pdry/tht**0.5 + 0.348*(freq**0.15)*4.2952193e-5*pwet*tht**4), real32) !modification 3

    gamh2o=0.1820*freq*sftot
  end subroutine abh2o_rk_modified

  !   compute atmospheric downwelling and upwelling brightness temperatures
  !   and upward transmittance at each pressure level (altitude)

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !   input:
  !     nlev           number of atmosphere levels
  !     tht            earth incidence angle [degree]
  !     tabs(0:nlev)   atmospheric absorption coefficients [Np/m]
  !     t(0:nlev)      temperature profile [K]
  !   z(0:nlev)      elevation [m]


  !     output:
  !     tran          total atmospheric transmission
  !     tbdw          downwelling brightness temperature t_bd [K]
  !     tbup          upwelling   brightness temperature t_bu [K]
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
    end do

    sumop=0
    sumdw=0
    do i=1,nlev
       sumdw=sumdw+(tavg(i)-t(1))*ems(i)*exp(sumop)
       sumop=sumop+opacty(i)
    end do

    sumop=0
    sumup=0.
    do i=nlev,1,-1
       sumup=sumup+(tavg(i)-t(1))*ems(i)*exp(sumop)
       sumop=sumop+opacty(i)
    end do

    tran=exp(sumop)
    tbavg=(1.-tran)*t(1)
    tbdw=tbavg+sumdw
    tbup=tbavg+sumup
  end subroutine atm_tran

  !   input:
  !     oxygen absorption from rosenkranz
  !     freq  frequency [GHz]
  !     p      pressure [hPa]
  !     t      temperature [k]
  !     pv     water vapor pressure  [hPa]
  !
  !     output:
  !     av          water vapor absorption coefficient [Np/km]
  !     ao          oxygen absorption coefficient        [Np/km]
  subroutine fdabscoeff(freq,p,t,pv, av,ao)
    real(real32), intent(in) :: freq,p,t,pv
    real(real32), intent(out) :: av,ao

    real(real32), parameter :: nep_scale = 0.1 * log(10.0)
    real(real32) gamoxy,gamh2o

    call fdabsoxy_1992_modified(p,t,pv,freq, gamoxy)
    call abh2o_rk_modified(p,t,pv,freq, gamh2o)

    ! The result values above are in dB/km, so convert to Np/km
    ao=nep_scale * gamoxy
    av=nep_scale * gamh2o
  end subroutine fdabscoeff

  !     liquid cloud water absorption
  !     rayleigh
  !     freq:      frequency [GHz]
  !     t:         temperature [K]
  !     rhol:      liquid cloud water density [g/m**3]
  !     output:
  !     al:        cloud water absorption coefficient [Np/km]
  pure subroutine fdcldabs(freq,t,rhol,   al)
    real(real32), intent(in) :: freq, t, rhol
    real(real32), intent(out) :: al
    real(real32) :: rhol0, wavlen
    real(real32), parameter :: c=29.979, pi=3.14159
    complex(real32) :: permit

    rhol0 = 1.0e-6*rhol ![g/cm**3]
    call meissner(freq,t,0.0,   permit)
    wavlen = c/freq
    al = (6.0*pi*rhol0/wavlen)*aimag((1.0-permit)/(2.0+permit))  ! Np/cm
    al = 1.0e5*al ! Np/km
  end subroutine fdcldabs

end module atms_abs_routines
