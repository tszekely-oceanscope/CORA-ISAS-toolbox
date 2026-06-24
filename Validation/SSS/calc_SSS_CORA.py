"""
calc_SSS_CORA.py
----------------
Compute area-weighted mean Sea Surface Salinity (SSS) from CORA OA monthly
fields and write a CSV time series per ocean basin.

Usage
-----
    python calc_SSS_CORA.py --inputpath <root_dir> [--Y0 1960] [--Y1 2024]

Arguments
---------
--inputpath : str
    Root directory containing yearly sub-folders of OA field NetCDF files,
    e.g. /data/ISAS_RESU/field/  ->  .../field/2005/OA_CORA5.2_200508_fld_PSAL.nc
--Y0 : int
    First year to process (inclusive). Default: 1960.
--Y1 : int
    Last year to process (inclusive). Default: 2024.

Output
------
Mean_SSS_CORA.csv — columns: date, Global, Pacific, Atlantic, Indian, South
"""

import argparse
import xarray as xr
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import gsw
from matplotlib import pyplot as plt
import sys
import glob
from scipy.interpolate import griddata
import pandas as pd

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Compute mean SSS from CORA OA fields")
parser.add_argument("--inputpath", required=True, help="Root directory of OA field files")
parser.add_argument("--Y0", type=int, default=1960, help="Start year (inclusive)")
parser.add_argument("--Y1", type=int, default=2024, help="End year (inclusive)")
args = parser.parse_args()


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def cell_area(lat, lon):
    """Return 2-D array of grid-cell areas (in Earth-radius units, R=6.371).

    Parameters
    ----------
    lat, lon : array-like
        1-D latitude and longitude vectors (degrees).

    Returns
    -------
    dA : ndarray, shape (len(lat), len(lon))
    """
    lat = np.asarray(lat).flatten()
    lon = np.asarray(lon).flatten()

    R = 6.371  # Earth radius (arbitrary units; consistent across script)
    dlat = np.deg2rad(np.gradient(lat))
    dlon = np.deg2rad(np.gradient(lon))
    lat_rad = np.deg2rad(lat)

    dA = (R**2) * np.outer(dlat, dlon) * np.cos(lat_rad[:, None])
    return dA


def weighted_mean_region(sal, mask, area):
    """Area-weighted mean salinity over a masked region.

    Parameters
    ----------
    sal : ndarray
        2-D salinity field (lat × lon).
    mask : ndarray of bool
        Region mask, same shape as sal.
    area : ndarray
        Grid-cell areas, same shape as sal.

    Returns
    -------
    float or np.nan if the masked region contains no valid data.
    """
    valid = ~np.isnan(sal) & mask
    if np.sum(valid) == 0:
        return np.nan
    return np.sum(sal[valid] * area[valid]) / np.sum(area[valid])


# ---------------------------------------------------------------------------
# Bathymetry — loaded once, interpolated to OA grid on first file
# ---------------------------------------------------------------------------
bathy_path = "../../bathy_isas15.nc"
ds_bathy = xr.open_dataset(bathy_path)
bathy = ds_bathy["bathymetry"].astype("float32")
bathy = bathy.where(bathy != 32767)  # mask fill value
latbath = ds_bathy["latitude"]
lonbath = ds_bathy["longitude"]

# ---------------------------------------------------------------------------
# Main loop — iterate over months, accumulate basin means
# ---------------------------------------------------------------------------
rep_OA = args.inputpath.rstrip("/") + "/"
records = []
isfirst = True

for year in range(args.Y0, args.Y1 + 1):
    for mm in range(1, 13):
        month = "%02d" % mm

        fpath = rep_OA + str(year) + "/OA_CORA5.2_" + str(year) + month + "15_fld_PSAL.nc"
        if not os.path.exists(fpath):
            continue

        print(fpath)
        ds_sss = xr.open_dataset(fpath)

        # Surface salinity (depth level = 1 m)
        sal = ds_sss.PSAL.sel(depth=1).squeeze().values
        lon = ds_sss.longitude.values
        lat = ds_sss.latitude.values

        # One-time initialisation: grid areas and basin masks
        if isfirst:
            isfirst = False

            area = cell_area(lat, lon)  # shape (lat, lon)

            # Interpolate bathymetry onto the OA grid
            points = np.array([(la, lo) for la in latbath for lo in lonbath])
            values = bathy.values.flatten()
            grid_lat, grid_lon = np.meshgrid(lat, lon, indexing="ij")
            bathy_interp = griddata(points, values, (grid_lat, grid_lon), method="linear")

            Lon, Lat = np.meshgrid(lon, lat)
            deep = bathy_interp > 1500  # open-ocean depth threshold (m)

            # Basin masks (longitude/latitude boxes + depth filter)
            mask_pac = (
                ((Lon >= -99) & (Lon <= 120) & (Lat > 0) & (Lat < 60))
                | ((Lon >= -70) & (Lon <= 120) & (Lat > -60) & (Lat < 0))
            ) & deep

            mask_atl = (
                ((Lon >= -99) & (Lon <= 42) & (Lat > 0) & (Lat < 60))
                | ((Lon >= -70) & (Lon <= 22) & (Lat > -60) & (Lat < 0))
            ) & deep

            mask_ind = ((Lon > 22) & (Lon <= 120) & (Lat > -60) & (Lat < 31)) & deep

            mask_antarctic = (Lat <= -35) & deep

            lat_mask = (Lat > -60) & (Lat < 60) & deep  # global (60S–60N)

        date = pd.to_datetime(str(year) + month, format="%Y%m")

        global_mean = weighted_mean_region(sal, lat_mask, area)
        mean_pac    = weighted_mean_region(sal, mask_pac, area)
        mean_atl    = weighted_mean_region(sal, mask_atl, area)
        mean_ind    = weighted_mean_region(sal, mask_ind, area)
        mean_antar  = weighted_mean_region(sal, mask_antarctic, area)

        records.append((date, global_mean, mean_pac, mean_atl, mean_ind, mean_antar))

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------
df = pd.DataFrame(
    records,
    columns=["date", "Global", "Pacific", "Atlantic", "Indian", "South"],
).sort_values("date")

df.to_csv("Mean_SSS_CORA.csv", index=False)
print("Saved Mean_SSS_CORA.csv")
