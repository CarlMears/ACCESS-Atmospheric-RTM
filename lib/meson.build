# ACCESS RTM project
#
# Note that some suggested examples include:
#
# meson setup build_debug/ rtm/
# meson setup --buildtype release build_release/ rtm/
# meson setup -D b_sanitize=address,undefined build_debug_sanitize/ rtm/
# meson setup --buildtype release -D b_lto=true build_release_lto/ rtm/

project('access_rtm',
    'fortran',
    version : '0.1',
    default_options : [
        'buildtype=debugoptimized',
        'warning_level=3',
        'b_pie=true',
        'fortran_std=f2018',
    ],
    meson_version : '>=0.49')

if meson.get_compiler('fortran').get_argument_syntax() == 'gcc'
    add_project_arguments(
        '-Wimplicit-interface',
        '-ffpe-trap=zero,overflow',
        '-fno-plt',
        # The code itself should already satisfy these next two, but this is
        # just to check that I didn't miss any
        '-fmodule-private',
        '-fimplicit-none',
    language : 'fortran')

    if meson.version().version_compare('<0.53')
        # This is instead handled by the "fortran_std" option for newer
        # versions of meson
        add_project_arguments('-std=f2018', language : 'fortran')
    endif

    add_project_link_arguments(
        '-Wl,-z,now',
        '-Wl,-z,relro',
        language : 'fortran')

    if get_option('debug')
        add_project_arguments('-fcheck=all', language : 'fortran')
    endif
    if get_option('optimization') == '3'
        add_project_arguments('-ftree-vectorize', language : 'fortran')
    endif

else
    warning('Not adding any compiler options')
endif

openmp = dependency('openmp', required : false)

# -----------------------
static_library('access_rtm',
    sources : [
        'access_rtm.f90',
        'atms_abs_routines.f90',
        'dielectric_meissner.f90',
        'trig_degrees.f90',
        'wvap_convert.f90',
    ],
    dependencies : [openmp])
