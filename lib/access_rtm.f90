program access_rtm
  use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
  use daily_scene_era5, only: DailySceneDataEra5
  implicit none

  character(len=*), parameter :: PROG_NAME = "access_rtm"

  character(len=250) :: filename_out, input_dir
  integer :: year, doy
  type(DailySceneDataEra5) :: data

  ! ----------------------------------------------------------------------
  call parse_args(year, doy, filename_out, input_dir)

  write (*,'(" Applying RTM using ERA5 data for: ", i4.4, "-", i3.3)') year, doy
  write (*,*) "ERA5 data directory: ", trim(input_dir)
  call data%load(input_dir, year, doy)

  write (*, *) "Writing output: ", filename_out
  call data%save(filename_out)

contains

  ! Parse the command-line arguments, stopping if there's an error
  subroutine parse_args(year, doy, filename_out, input_dir)
    integer, intent(out) :: year, doy
    character(len=*), intent(out) :: filename_out, input_dir

    integer :: nargs, iarg, argstat, readstat, num_parsed, mode
    character(len=80) :: arg_val

    nargs = command_argument_count()
    if (nargs == 0) then
      call print_help()
      stop
    end if

    num_parsed = 0
    mode = 0
    input_dir = "."
    args: do iarg = 1, nargs
      call get_command_argument(iarg, arg_val, status=argstat)
      if (argstat /= 0) error stop "Could not get an arg properly"

      select case(trim(arg_val))
      case ("-h", "--help")
        call print_help()
        stop
      case ("-V", "--version")
        write (*,*) PROG_NAME
        stop
      case ("--in-dir")
        mode = 1
      case default
        if (mode == 0) then
          ! Mode 0: checking for top-level args
          select case (num_parsed)
          case (0)
            read(arg_val, *, iostat=readstat) year
            if (readstat /= 0) then
              write (ERROR_UNIT, *) "Cannot parse to integer: ", arg_val
              error stop "Argument parsing failed"
            end if
            num_parsed = num_parsed + 1
          case (1)
            read(arg_val, *, iostat=readstat) doy
            if (readstat /= 0) then
              write (ERROR_UNIT, *) "Cannot parse to integer: ", arg_val
              error stop "Argument parsing failed"
            end if
            num_parsed = num_parsed + 1
          case (2)
            filename_out = trim(arg_val)
            num_parsed = num_parsed + 1
          case default
            call print_help()
            error stop "Too many positional arguments"
          end select

        else if (mode == 1) then
          ! Mode 1: checking for arg for --in-dir
          input_dir = trim(arg_val)
          mode = 0
        end if
      end select
    end do args

    if (num_parsed /= 3) then
       call print_help()
       error stop "Didn't parse enough positional arguments"
    end if
    if (mode /= 0) then
      call print_help()
      error stop "Didn't finish parsing arguments"
   end if
  end subroutine parse_args

  subroutine print_help
    write (*, *) " usage: " // PROG_NAME // " [-h] [-V] " &
         // "[--in-dir DIR] year doy out_file"
    write (*, *)
    write (*, *) " Apply RTM for ACCESS for one day of ERA5 data"
    write (*, *)
    write (*, *) " positional arguments:"
    write (*, *) "   year      Year to use for geophysical input data"
    write (*, *) "   doy       Ordinal day for geophysical input data"
    write (*, *) "   out_file  Output file to write (netCDF format)"
    write (*, *)
    write (*, *) " optional arguments:"
    write (*, *) "   --in-dir DIR    read ERA5 data files from this directory"
    write (*, *) "   -h, --help      show this help message and exit"
    write (*, *) "   -V, --version   show the version info and exit"
  end subroutine print_help

end program access_rtm
