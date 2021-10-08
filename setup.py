"""Setuptools config for access-atmosphere."""

from numpy.distutils.core import Extension, setup

extensions = [
    Extension(
        name="access_atmosphere.rtm_f90",
        # Note the order here is important
        sources=[
            "lib/access_rtm.pyf",
            "lib/trig_degrees.f90",
            "lib/dielectric_meissner.f90",
            "lib/wvap_convert.f90",
            "lib/atms_abs_routines.f90",
            "lib/access_rtm.f90",
        ],
        extra_f90_compile_args=["-std=f2018", "-fopenmp"],
        extra_link_args=["-fopenmp"],
    )
]

setup(
    name="access-atmosphere",
    version="0.1.0a0",
    packages=["access_atmosphere"],
    python_requires=">=3.7",
    install_requires=[
        "cdsapi",
        "numpy",
        "netCDF4",
    ],
    author="Richard Lindsley",
    author_email="lindsley@remss.com",
    url="http://gitlab.remss.com/access/atmospheric-rtm/",
    include_package_data=True,
    ext_modules=extensions,
)
