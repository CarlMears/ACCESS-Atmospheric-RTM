module month_day
  use, intrinsic :: iso_fortran_env, only: int32
  implicit none
  private
  public :: find_month_day

contains

  pure SUBROUTINE FIND_MONTH_DAY(LYEAR,IDAYJL, IMON,IDAY)
    integer, intent(in) :: lyear, idayjl
    integer, intent(out) :: imon, iday

    integer :: ileap, jmon

    integer, dimension(12, 0:1), parameter :: idayfx = reshape( &
         [1,32,60,91,121,152,182,213,244,274,305,335, &
         1 ,32,61,92,122,153,183,214,245,275,306,336], &
         [12, 2])

    ILEAP=0
    IF(LYEAR.EQ.4*INT(LYEAR/4)) ILEAP=1

    DO JMON=2,12
       IF(IDAYFX(JMON,ILEAP).GT.IDAYJL) THEN
          IMON=JMON-1
          GO TO 20
       ENDIF
    end do
    IMON=12
20  CONTINUE

    IDAY=1+IDAYJL-IDAYFX(IMON,ILEAP)
  END SUBROUTINE FIND_MONTH_DAY
end module month_day
