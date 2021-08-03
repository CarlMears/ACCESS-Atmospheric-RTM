# ACCESS Atmospheric RTM

For the ACCESS project, download ERA5 atmosphere geophysical parameters and apply the radiative transfer model to obtain atmospheric microwave terms.

## Overview

Currently, this contains two main projects:

- `download_era5.py`: Python script to download ERA5 data
- `rtm/`: Fortran executable to apply RTM to ERA5 data

## Downloading ERA5 data

To download data from the [Climate Data
Store](https://cds.climate.copernicus.eu/cdsapp), an account needs to be
registered. Note the UID and API key for later.

Before running the Python download script, its dependencies can be installed into a local [virtual environment](https://docs.python.org/3/library/venv.html):

```bash
python3 -m venv --upgrade-deps .venv
source .venv/bin/activate
pip install cdsapi
```

To download the ERA5 datasets of interest for some time range, the script can be
run with two environment variables set for CDS authentication:

```bash
# (assuming the virtualenv is activated)
env CDS_UID=xxx CDS_API_KEY=xxx python3 ./download_era5.py 2020-01-01 2020-01-31 --out-dir era5
```

By default the netCDF files are written to the current working directory.

## Applying the RTM

A Fortran executable, `access_rtm`, is in the `rtm/` directory. Besides a
Fortran compiler, it also depends on the
[netCDF-Fortran](https://github.com/Unidata/netcdf-fortran) library. On a
[RHEL](https://en.wikipedia.org/wiki/Red_Hat_Enterprise_Linux)-like system this can all be set up using:

```bash
sudo dnf install -y gcc-gfortran netcdf-fortran-devel meson ninja-build
```

And then it can be built using [meson](https://mesonbuild.com/):

```bash
meson setup build rtm
meson compile -C build
```

An example run for 2020-01-01:

```bash
./build/access_rtm 1 2020 1 access_era_2020-01-01.nc
```

It reads the ERA5 data for a single day and runs the RTM for some hard-coded
incidence angle and microwave frequency configuration. The outputs are on a
0.25-degree grid.

Currently it also supports reading NCEP GDAS data (as post-processed by RSS) for
debugging. This, however, is on a 1-degree grid.

```bash
./build/access_rtm 2 2020 1 access_ncep_2020-01-01.nc
```