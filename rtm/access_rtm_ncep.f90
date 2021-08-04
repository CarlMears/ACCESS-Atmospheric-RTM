program access_rtm_ncep
  use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
  use daily_scene_ncep, only: DailySceneDataNcep
  implicit none

  character(len=*), parameter :: PROG_NAME = "access_rtm_ncep"

  character(len=120) :: filename_out
  integer :: year, doy
  type(DailySceneDataNcep) :: data

  ! ----------------------------------------------------------------------
  call parse_args(year, doy, filename_out)

  write (*,'(" Applying RTM using NCEP data for: ", i4.4, "-", i3.3)') year, doy
  call data%load(year, doy)

  write (*, *) "Writing output: ", filename_out
  call data%save(filename_out)

contains

  ! Parse the command-line arguments, stopping if there's an error
  subroutine parse_args(year, doy, filename_out)
    integer, intent(out) :: year, doy
    character(len=*), intent(out) :: filename_out

    integer :: nargs, iarg, argstat, readstat, num_parsed
    character(len=80) :: arg_val

    nargs = command_argument_count()
    if (nargs == 0) then
      call print_help()
      stop
    end if

    num_parsed = 0
    args: do iarg = 1, nargs
      call get_command_argument(iarg, arg_val, status=argstat)
      if (argstat /= 0) error stop "Could not get an arg properly"

      if (trim(arg_val) == "-h" .or. trim(arg_val) == "--help") then
          call print_help()
          stop
      else if (trim(arg_val) == "-V" .or. trim(arg_val) == "--version") then
          write (*,*) PROG_NAME
          stop
      else
         if (num_parsed == 0) then
             read(arg_val, *, iostat=readstat) year
             if (readstat /= 0) then
                write (ERROR_UNIT, *) "Cannot parse to integer: ", arg_val
                error stop "Argument parsing failed"
             end if
             num_parsed = num_parsed + 1
         else if (num_parsed == 1) then
             read(arg_val, *, iostat=readstat) doy
             if (readstat /= 0) then
                write (ERROR_UNIT, *) "Cannot parse to integer: ", arg_val
                error stop "Argument parsing failed"
             end if
             num_parsed = num_parsed + 1
          elseif (num_parsed == 2) then
             filename_out = trim(arg_val)
             num_parsed = num_parsed + 1
          else
             call print_help()
             error stop "Too many positional arguments"
          end if
       end if
    end do args

    if (num_parsed /= 3) then
       call print_help()
       error stop "Didn't parse enough positional arguments"
    end if
  end subroutine parse_args

  subroutine print_help
    write (*, *) " usage: " // PROG_NAME // " [-h] [-V] " &
         // "year doy out_file"
    write (*, *)
    write (*, *) " Apply RTM for ACCESS for one day of NCEP data"
    write (*, *)
    write (*, *) " positional arguments:"
    write (*, *) "   year      Year to use for geophysical input data"
    write (*, *) "   doy       Ordinal day for geophysical input data"
    write (*, *) "   out_file  Output file to write (netCDF format)"
    write (*, *)
    write (*, *) " optional arguments:"
    write (*, *) "   -h, --help      show this help message and exit"
    write (*, *) "   -V, --version   show the version info and exit"
  end subroutine print_help

end program access_rtm_ncep
