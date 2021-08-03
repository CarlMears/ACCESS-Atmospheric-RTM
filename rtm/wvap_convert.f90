module wvap_convert
  use, intrinsic :: iso_fortran_env, only: real32
  implicit none
  private
  public :: goff_gratch_vap

contains
  elemental subroutine goff_gratch_vap(T,RH,P,  P_V,rho_V)
    !
    !
    !     calculates water vapor pressure P_V and density rho_V
    !     from relative humidity RH
    !     S. Cruz Pol, C. Ruf and S. Keihm, Radio Science 33 (5),1319 (1998)
    !
    !    water vapor pressure: P_v = P_s * RH
    !    water vapor density   rho_v = F_w *(P_v * eps) / (R_d * T)
    !
    real(real32), intent(in) :: T,RH,P
    real(real32), intent(out) :: P_v,rho_v

    real(real32) :: P_s
    real(real32) :: F_w,xi1,xi2,xi3

    real(real32), parameter :: p_standard = 1013.246  ! hPa
    real(real32), parameter :: TS = 373.14     ! K  (water boiling)
    real(real32), parameter :: T0 = 273.14  ! K
    real(real32), parameter :: eps = 1./1.607795  ! M(H2O)/M(air)
    real(real32), parameter :: R_d = 287.05 ! J /(K*kg)  gas constant for dry air
    !
    F_w = 1.0 + 1.E-4 * (5.92854 + 3.740346e-2 * P + &
         1.971198E-4 * (T-T0) * (800-P)  + &
         6.045511E-6 * P * (T-T0)**2 ) ! deviation from ideal gas
    !
    xi1 = -7.90298*(Ts/T - 1) + 5.02808*Log10(TS/T)
    xi2 = -1.3816E-7 * 10**(11.344*(1.-T/TS) -1 )
    xi3 = 8.1328E-3 * (10**(-3.49149*(TS/T - 1)) - 1)
    !
    P_s = p_standard * 10**(xi1 + xi2 + xi3)
    P_v = P_s * RH* 0.01                                           !  mbar
    rho_v = ((F_w * P_v * eps) / (R_d * T)) * 1.E5  !  [g/m**3]
  end subroutine goff_gratch_vap
end module wvap_convert
