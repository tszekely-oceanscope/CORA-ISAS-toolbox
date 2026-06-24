"""
calc_SST_CORA.py

Compute area-weighted mean SST anomalies (relative to a 2010-2020 climatology)
from CORA5.2 objective-analysis temperature fields, for the global ocean and
four sub-basins (Pacific, Atlantic, Indian, Southern).

Usage:
    python calc_SST_CORA.py --Y0 <start_year> --Y1 <end_year> --inputpath <data_dir>

Arguments:
    --Y0          First year to process (inclusive)
    --Y1          Last year to process (inclusive)
    --inputpath   Root directory containing per-year CORA5.2 netCDF files

Output:
    Mean_SST-CORA.csv  — monthly time series of area-weighted SST anomalies
"""

import argparse
import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

import matplotlib
matplotlib.use("Agg")  # Pas besoin de X11

parser = argparse.ArgumentParser(
    description="Compute area-weighted SST anomalies from CORA5.2 fields."
)
parser.add_argument("--Y0", type=int, required=True, help="Start year")
parser.add_argument("--Y1", type=int, required=True, help="End year")
parser.add_argument("--inputpath", type=str, required=True, help="Input data directory")
args = parser.parse_args()


def cell_area(lat, lon):
    """
    Compute spherical cell areas (in units of Earth-radius squared) on a regular grid.

    Parameters
    ----------
    lat : array-like, shape (nlat,)
    lon : array-like, shape (nlon,)

    Returns
    -------
    dA : ndarray, shape (nlat, nlon)
    """
    lat = np.asarray(lat).flatten()
    lon = np.asarray(lon).flatten()

    R = 6.371  # Earth radius in Mm (factors out; result is dimensionless ratio)
    dlat = np.deg2rad(np.gradient(lat))
    dlon = np.deg2rad(np.gradient(lon))
    lat_rad = np.deg2rad(lat)

    dA = (R**2) * np.outer(dlat, dlon) * np.cos(lat_rad[:, None])  # (nlat, nlon)
    return dA


data_dir = args.inputpath

records = []
isfirst = 1

for yy in range(args.Y0, args.Y1 + 1):
    for mm in range(1, 13):
        year = str(yy)
        month = '%02d' % mm

        # Build path to monthly CORA5.2 temperature field
        files = data_dir + f"{year}/OA_CORA5.2_{year}{month}15_fld_TEMP.nc"
        print(files)
        ds = xr.open_dataset(files)

        date_str = year + month
        date = pd.to_datetime(date_str, format="%Y%m")

        # Load monthly climatology for anomaly computation
        clim = f"../Steric/clim/CORA_OA_CLIM_TEMP_2010-2020_{month}.nc"
        dscl = xr.open_dataset(clim)

        lat = ds.latitude
        lon = ds.longitude.values

        # SST anomaly = observation minus climatology (surface level, first time step)
        sal = ds["TEMP"].values[0, 0, :, :] - dscl["TEMP"].values[0, :, :]

        # Build cell areas and basin masks once on the first iteration
        if isfirst == 1:
            isfirst = 0
            area = cell_area(lat, lon)  # (nlat, nlon)
            Lon, Lat = np.meshgrid(lon, lat)
            mask_pac = (
                ((Lon >= -99) & (Lon <= 120) & (Lat > 0) & (Lat < 60)) |
                ((Lon >= -70) & (Lon <= 120) & (Lat > -60) & (Lat < 0))
            )
            mask_atl = (
                ((Lon >= -99) & (Lon <= 42) & (Lat > 0) & (Lat < 60)) |
                ((Lon >= -70) & (Lon <= 22) & (Lat > -60) & (Lat < 0))
            )
            mask_ind = ((Lon > 22) & (Lon <= 120)) & (Lat > -60) & (Lat < 31)
            mask_antarctic = (Lat <= -30)

        # Global mean: restrict to 60S–60N to exclude sea-ice-affected poles
        lat_mask = (lat > -60) & (lat < 60)
        sal_masked = sal[lat_mask, :]
        area_masked = area[lat_mask, :]

        valid = ~np.isnan(sal_masked)
        global_mean = np.sum(sal_masked[valid] * area_masked[valid]) / np.sum(area_masked[valid])

        def weighted_mean_region(sal, mask):
            """Area-weighted mean of `sal` over grid points selected by `mask`."""
            valid = ~np.isnan(sal) & mask
            if np.sum(valid) == 0:
                return np.nan
            return np.sum(sal[valid] * area[valid]) / np.sum(area[valid])

        mean_pac   = weighted_mean_region(sal, mask_pac)
        mean_atl   = weighted_mean_region(sal, mask_atl)
        mean_ind   = weighted_mean_region(sal, mask_ind)
        mean_antar = weighted_mean_region(sal, mask_antarctic)

        records.append((date, global_mean, mean_pac, mean_atl, mean_ind, mean_antar))

# Assemble time series and write output
df = pd.DataFrame(
    records,
    columns=["date", "Global", "Pacific", "Atlantic", "Indian", "South"]
).sort_values("date")
df.to_csv('Mean_SST-CORA.csv', index=False)
