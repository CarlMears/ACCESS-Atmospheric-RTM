# ACCESS Atmospheric RTM

For the ACCESS project, download ERA5 atmosphere geophysical parameters and apply the radiative transfer model to obtain atmospheric microwave terms.

## Overview

This is a Python package, `access-atmosphere`, that currently performs two main
tasks:

- Download relevant ERA5 data
- Compute atmospheric RTM based on the ERA5 data

The RTM is implemented in Rust and compiled into a Python extension using
[maturin](https://maturin.rs/) and [pyo3](https://pyo3.rs/).

Following the [NumPy policy of supported Python
versions](https://numpy.org/neps/nep-0029-deprecation_policy.html#drop-schedule),
Python 3.8 is the minimum version supported.

## Building

To build the package locally, both [Python](https://www.python.org/) and
[Rust](https://www.python.org/) are required. Rust can be installed using using
[`rustup`](https://rustup.rs/). [Maturin](https://maturin.rs/) is used to build
everything.

On Linux:

```bash
# Install maturin in a virtual environment
python3 -m venv --upgrade-deps .venv
source .venv/bin/activate
pip install maturin

# Build the wheel
maturin build --release
# The wheel is in target/wheels/ and can be installed using "pip install"

# For local development, the wheel can be built and installed in the current venv
maturin develop
```

On Windows (noting that `Set-ExecutionPolicy -ExecutionPolicy RemoteSigned
-Scope CurrentUser` is run first, as mentioned in the [venv
documentation](https://docs.python.org/3/library/venv.html)):

```powershell
py.exe -m venv --upgrade-deps venv
.\venv\Scripts\Activate.ps1
pip install maturin

# Build the wheel. The "abi3" feature is used to avoid looking for the Python
# development libraries since they're bundled in this way.
maturin build --release --features abi3

# For local development, building and installing the wheel in the venv
maturin develop --features abi3
```

Alternately, the GitLab CI automatically builds wheels. The built wheels are:

- x86_64 [manylinux](https://github.com/pypa/manylinux) wheels specifically for
  Python versions 3.8 through 3.10
- x86_64 manylinux abi3 wheel for Python 3.8 and later
- x86_64 Windows abi3 wheel for Python 3.8 and later

The wheels can be downloaded as CI job artifacts, and for every release, are
generated and saved in the [local package
registry](http://gitlab.remss.com/access/atmospheric-rtm/-/packages).

Note that with a wheel built, it's ready to install without needed a local Rust
compiler. To install the wheel in a new venv and download the dependencies from
PyPI:

```bash
python3 -m venv --upgrade-deps .venv
.venv/bin/pip install access-atmosphere-*.whl
```

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

<!-- 
# TODO: update this for however I expose the parallelism from Rust

The Fortran code is compiled with [OpenMP](https://www.openmp.org/) in order to
process each profile in parallel using a pool of worker threads. By default,
this will be as many threads as there are logical CPUs detected on the machine.
To modify this at runtime, use the `OMP_NUM_THREADS` environment variable. For
instance, to set exactly 4 threads:

```bash
env OMP_NUM_THREADS=4 python3 -m access_atmosphere.process ...
``` -->