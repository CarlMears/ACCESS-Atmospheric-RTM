# ACCESS Atmospheric RTM

For the ACCESS project, download ERA5 atmosphere geophysical parameters and apply the radiative transfer model to obtain atmospheric microwave terms.

## Example usage

Before running the Python download script:

```bash
python3 -m venv --upgrade-deps .venv
source .venv/bin/activate
pip install cdsapi
```

To download the ERA5 datasets of interest for some time range:

```bash
# (assuming the virtualenv is activated)
python3 ./download_era5.py 2020-01-01 2020-01-31 --out-dir era5
```

By default the netCDF files are written to the current working directory.