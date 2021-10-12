# ACCESS Atmospheric RTM

For the ACCESS project, download ERA5 atmosphere geophysical parameters and apply the radiative transfer model to obtain atmospheric microwave terms.

## Overview

This is a Python package, `access-atmosphere`, that currently performs two main
tasks:

- Download relevant ERA5 data
- Compute atmospheric RTM based on the ERA5 data

The RTM is implemented in Fortran and compiled using `f2py` into a Python extension.

## Building

Assuming a Fortran compiler is available, the package can be built using
[build](https://github.com/pypa/build). Starting from a scratch in a new
[virtual environment](https://docs.python.org/3/library/venv.html):

```bash
python3 -m venv --upgrade-deps .venv
.venv/bin/pip install build
.venv/bin/python -m build
```

This creates a source distribution as `dist/access-atmosphere-$VERSION.tar.gz`
and a wheel as `dist/access-atmosphere-$VERSION-*.whl`. (The extra information
at the end of the filename contains the Python version and platform
information.) Either file (sdist or wheel) can be used to install the package,
but note that a big benefit for the wheel is that it's ready to go and no
Fortran compiler is needed.

To install the wheel in a new venv and download the dependencies from PyPI:

```bash
python3 -m venv --upgrade-deps .venv
.venv/bin/pip install access-atmosphere-$VERSION-*.whl
```

Alternately, Gitlab CI jobs are set up to build
[manylinux](https://github.com/pypa/manylinux) wheels for multiple versions of
Python. The wheels can be downloaded as CI job artifacts, and for every release,
are generated and saved in the [local package
registry](http://gitlab.remss.com/access/atmospheric-rtm/-/packages).

As another option, a [`Dockerfile`](Dockerfile) is provided. Build it with
`docker` or `podman`:

```bash
podman build -t access_atmosphere -f Dockerfile
```

## Running

### Downloading ERA5 data

To download data from the [Climate Data
Store](https://cds.climate.copernicus.eu/cdsapp), an account needs to be
registered. Note the UID and API key for later.

To download the ERA5 datasets of interest for some time range, the script can be
run with two environment variables set for CDS authentication:

```bash
env CDS_UID=xxx CDS_API_KEY=xxx python3 -m access_atmosphere.download 2020-01-01 2020-01-31 --out-dir era5
```

For each day, two files are created: one for the surface data and one for the
atmospheric profiles. The netCDF files are written to the directory specified by
`--out-dir`, or the current working directory if it's not specified. Any
existing files are not re-downloaded from CDS.

### Applying the RTM

A Python interface to the RTM is provided in the `access_atmosphere.rtm` module.

As an example usage, the `access_atmosphere.process` module reads in a pair of
ERA5 daily files and runs the RTM for every point (latitude/longitude/time) and
writes the results in a new netCDF file. It uses reference values for the
incidence angles and microwave frequencies. The outputs are on the same
0.25-degree grid that ERA5 uses.

As an example using the 2020-01-01 data downloaded above:

```bash
python -m access_atmosphere.process \
    era5_surface_2020-01-01.nc \
    era5_levels_2020-01-01.nc \
    access_era5_2020-01-01.nc
```

The Fortran code is compiled with [OpenMP](https://www.openmp.org/) in order to
process each profile in parallel using a pool of worker threads. By default,
this will be as many threads as there are logical CPUs detected on the machine.
To modify this at runtime, use the `OMP_NUM_THREADS` environment variable. For
instance, to set exactly 4 threads:

```bash
env OMP_NUM_THREADS=4 python3 -m access_atmosphere.process ...
```