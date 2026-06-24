import xarray as xr
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import argparse


def cell_area(lat, lon):
    """Compute spherical surface area of each grid cell [km²].

    Parameters
    ----------
    lat : array-like, shape (lat,)
        Latitude vector [°].
    lon : array-like, shape (lon,)
        Longitude vector [°].

    Returns
    -------
    dA : ndarray, shape (lat, lon)
        Cell area [km²] using Earth radius R = 6.371 km.
    """
    lat = np.asarray(lat).flatten()
    lon = np.asarray(lon).flatten()
    R = 6.371
    dlat = np.deg2rad(np.gradient(lat))
    dlon = np.deg2rad(np.gradient(lon))
    lat_rad = np.deg2rad(lat)
    return (R ** 2) * np.outer(dlat, dlon) * np.cos(lat_rad[:, None])


def main():
    """Entry point: parse CLI arguments and build the halosteric height time series.

    For each month between Y0 and Y1, reads the corresponding halosteric
    height field (CORA_HSSL_{year}{month}.nc), computes the area-weighted
    mean between 60°S and 60°N, and appends the result to a time series.
    The full series is written to a CSV; the DataFrame is then filtered to
    [Y0-01, Y1-12] for downstream use.

    Usage
    -----
    python calc_timeseries.py --Y0 1993 --Y1 2024 \\
        --inputpath /path/to/ssh/files --bathypath /path/to/bathy_isas15.nc
    """
    parser = argparse.ArgumentParser(description='Compute area-weighted steric height time series from SSH files.')
    parser.add_argument('--Y0', type=int, required=True, help='Start year (inclusive)')
    parser.add_argument('--Y1', type=int, required=True, help='End year (inclusive)')
    parser.add_argument('--inputpath', required=True,
                        help='Directory containing CORA_HSSL_{year}{month}.nc files')
    parser.add_argument('--bathypath', default='../../bathy_isas15.nc',
                        help='Path to bathymetry file (default: ../../bathy_isas15.nc)')
    args = parser.parse_args()

    ds_bathy = xr.open_dataset(args.bathypath)
    bathy = ds_bathy["bathymetry"].astype("float32")
    bathy = bathy.where(bathy != 32767)

    records = []
    for yy in range(args.Y0, args.Y1 + 1):
        for mm in range(1, 13):
            year = str(yy)
            month = f"{mm:02d}"
            filepath = os.path.join(args.inputpath, f"CORA_HSSL_{year}{month}.nc")
            if not os.path.exists(filepath):
                continue

            ds = xr.open_dataset(filepath)
            date = pd.to_datetime(year + month, format="%Y%m")

            lat = ds.latitude
            lon = ds.longitude.values
            area = cell_area(lat, lon)

            sal = ds["steric_height_anomaly"].values

            lat_mask = (lat > -60) & (lat < 60)
            sal_masked = sal[lat_mask, :]
            area_masked = area[lat_mask, :]

            valid = ~np.isnan(sal_masked)
            weighted_mean = np.sum(sal_masked[valid] * area_masked[valid]) / np.sum(area_masked[valid])
            records.append((date, weighted_mean))

    df = pd.DataFrame(records, columns=["date", "HSH"]).sort_values("date")

    output_csv = f"Mean_HSH_CORA_{args.Y0}_{args.Y1}.csv"
    df.to_csv(output_csv, index=False)
    print(f"Written: {output_csv}")

    df = df[(df["date"] >= f"{args.Y0}-01") & (df["date"] <= f"{args.Y1}-12")]


if __name__ == '__main__':
    main()
