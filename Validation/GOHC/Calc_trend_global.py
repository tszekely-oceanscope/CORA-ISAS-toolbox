"""
Calc_OHC_trend_global_NRT.py
-----------------------------
Compute the global Ocean Heat Content (OHC) time series from monthly OHC
NetCDF files produced by Calc_OHC_CORA.py and write the result to a
single output NetCDF file.

Usage
-----
    python Calc_OHC_trend_global_NRT.py \
        --inputdir /mydir/ \
        --outputdir /my/outdir/ \
        --file_name my_file.nc

Input files expected at:
    {inputdir}/OHC_ISASgrid_yearly_yearlyclim_{year}_{month:02d}.nc

Output file written to:
    {outputdir}/{file_name}
"""

import argparse
import os
import datetime as dt

import matplotlib
matplotlib.use("Agg")  # Use a non-GUI backend
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import xarray as xr


def build_ohc_timeseries(inputdir, outputdir, file_name):
    """
    Build the global OHC time series from pre-computed monthly OHC fields.

    For each month (1-12) the function:
      1. Loads all yearly OHC files for that month and computes a long-term climatology.
      2. Subtracts the climatology from each individual year to obtain anomalies.
      3. Sums the anomalies over all (longitude, latitude) grid cells to produce
         one scalar value per (year, month) pair.

    The resulting time series is chronologically sorted and saved as a NetCDF file.

    Parameters
    ----------
    inputdir : str
        Directory containing the monthly OHC NetCDF files
        (``OHC_ISASgrid_yearly_yearlyclim_{year}_{month:02d}.nc``).
    outputdir : str
        Directory where the output NetCDF file will be written (created if absent).
    file_name : str
        Name of the output NetCDF file (e.g. ``my_file.nc``).
    """
    print2d = 0

    if print2d == 0:
        plt.figure(figsize=(10, 6))

    vectime = []
    vecOHC = []
    for month in range(1, 13):
        print(month)
        file_paths = []
        for yy in range(1960, 2026):
            file_paths.append(os.path.join(
                inputdir,
                f"OHC_ISASgrid_yearly_yearlyclim_{yy}_{month:02d}.nc"
            ))

        datasets = [xr.open_dataset(fp, decode_times=False) for fp in file_paths]
        combined = xr.concat(datasets, dim="time")
        climatology = combined["OHC"].mean(dim="time", skipna=True)

        for yy in range(1960, 2026):
            print('   ', yy)
            file_path = os.path.join(
                inputdir,
                f"OHC_ISASgrid_yearly_yearlyclim_{yy}_{month:02d}.nc"
            )
            ds = xr.open_dataset(file_path, decode_times=False)
            temperature = ds["OHC"]

            delta = temperature - climatology
            temperature_zonal_mean = delta.sum(dim=("longitude", "latitude"))
            vecOHC.append(temperature_zonal_mean)
            vectime.append(dt.datetime(yy, month, 15))

    vecOHC = np.array(vecOHC)
    vectime = np.array(vectime)
    vectime2 = []
    vecOHC2 = []
    for yy in range(1960, 2026):
        for mm in range(1, 13):
            idg = np.where(vectime == dt.datetime(yy, mm, 15))
            vectime2.append(dt.datetime(yy, mm, 15))
            vecOHC2.append(vecOHC[idg[0][0]])

    vecglo = vecOHC2.copy()

    ohc_dataset = xr.Dataset(
        {
            "OHC_Global": (["time"], np.array(vecglo).ravel()),
        },
        coords={
            "time": np.array(vectime2).ravel(),
        }
    )

    ohc_dataset["OHC_Global"].attrs = {
        "long_name": "Ocean Heat Content",
        "units": "Joules",
        "description": "Total ocean heat content, integrated over depth 0-2000m"
    }

    ohc_dataset.attrs = {
        "title": "Ocean Heat Content Calculation",
        "source": "Computed from temperature field using rho=1030 kg/m3 and cp=3980 J/kg/K",
        "clim": "WOA23 interpolated on ISAS grid"
    }

    os.makedirs(outputdir, exist_ok=True)
    out_file = os.path.join(outputdir, file_name)
    ohc_dataset.to_netcdf(out_file)
    print(f"Saved: {out_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute global OHC time series from monthly CORA/ISAS OHC fields."
    )
    parser.add_argument("--inputdir",  required=True,
                        help="Directory containing the monthly OHC NetCDF files and basin mask")
    parser.add_argument("--outputdir", required=True,
                        help="Output directory for the resulting NetCDF file")
    parser.add_argument("--file_name", required=True,
                        help="Output NetCDF file name (e.g. my_file.nc)")
    args = parser.parse_args()

    build_ohc_timeseries(args.inputdir, args.outputdir, args.file_name)
