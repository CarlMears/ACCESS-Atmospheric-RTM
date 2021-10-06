module columnar_int
  implicit none
  private
  public :: column

contains
  !     columnar integral
  !     input:
  !     nlevel: number of profiles (nlevel = 0 -> sfc)
  !     z:      altitude levels [in m]
  !     rho:    density [in unit/m]     (unit arbitrary)
  !     ip:     =1 : linear varying profile -> arithmetic mean for integration
  !             =2 : approximately exponentially varying profile ->
  !                  use average of arithmetic and geometric mean for integration
  !             this has smaller error in integration than either
  !             arithmetic (a) or geometric (g) mean
  !             (if profile is varying exactly exponentially with z, then the integration error
  !              is minimal for the combination: 2/3 g + 1/3 a)
  !
  !     output:
  !     col [in unit/m]
  !
  subroutine column(nlevel,z,rho,ip,  col)
    integer(4), intent(in) :: nlevel, ip
    real(4), dimension(0:nlevel), intent(in) :: z, rho
    real(4), intent(out) :: col

    integer(4) :: i
    real(4) :: dz, avg

    col = 0.
    do i=1,nlevel
       if (z(i) <= z(i-1)) error stop 'error in column'

       dz = z(i) - z(i-1)

       if (ip == 1) then
          avg = 0.5 * (rho(i) + rho(i-1) )
       else
          avg = 0.25 * ( rho(i) + rho(i-1) + 2.*sqrt(rho(i-1)*rho(i)) )
       endif

       col = col + avg*dz
    end do
  end subroutine column
end module columnar_int
