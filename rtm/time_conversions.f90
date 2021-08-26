module time_conversions
  implicit none
  private
  public :: yo_to_ymd, is_leap

contains

  ! Convert an ordinal date to month and day
  pure subroutine yo_to_ymd(year, day_of_year, month, day_of_month)
    integer, intent(in) :: year, day_of_year
    integer, intent(out) :: month, day_of_month

    ! The first column is for non-leap years, the second column for leap years
    integer, dimension(12, 2), parameter :: DAYS_PER_MONTH = &
         reshape([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
         31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], [12, 2])

    integer :: ileap, imonth

    if (is_leap(year)) then
       ileap = 2
    else
       ileap = 1
    end if

    day_of_month = day_of_year
    do imonth = 1, 12
       if (day_of_month <= DAYS_PER_MONTH(imonth, ileap)) then
          month = imonth
          exit
       else
          day_of_month = day_of_month - DAYS_PER_MONTH(imonth, ileap)
       end if
    end do
  end subroutine yo_to_ymd

  ! Is this a leap year?
  pure function is_leap(year)
    integer, intent(in) :: year
    logical :: is_leap

    ! https://en.wikipedia.org/w/index.php?title=Leap_year&oldid=852830613#Algorithm
    if (modulo(year, 4) /= 0) then
       is_leap = .false.
    else if (modulo(year, 100) /= 0) then
       is_leap = .true.
    else if (modulo(year, 400) /= 0) then
       is_leap = .false.
    else
       is_leap = .true.
    end if
  end function is_leap

end module time_conversions
