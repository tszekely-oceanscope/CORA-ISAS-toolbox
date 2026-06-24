"""
Calc_OHC_CORA.py
----------------
Compute the Ocean Heat Content (OHC) from a CORA/ISAS temperature field
and write the result as a NetCDF file.

Usage
-----
    python Calc_OHC_CORA.py --inputdir /mydir/ --outputdir /my/outdir/ \
                             --year 2010 --month 01

Input file expected at:
    {inputdir}/{year}/OA_CORA5.2_{year}{month}15_fld_TEMP.nc

Output file written to:
    {outputdir}/OHC_ISASgrid_yearly_yearlyclim_{year}_{month}.nc
"""

import argparse
import xarray as xr
import numpy as np
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def calculate_ohc(temp_file, zmax):
    """
    Compute the column-integrated Ocean Heat Content from a CORA temperature file.

    The OHC per grid cell is defined as:
        OHC = rho * cp * integral(T(z) dz, 0, zmax) * cell_area

    Parameters
    ----------
    temp_file : str
        Path to the input NetCDF file containing the TEMP variable on the ISAS grid.
    zmax : float
        Maximum integration depth in metres (e.g. 2000).

    Returns
    -------
    xr.Dataset
        Dataset with a single variable ``OHC`` (Joules) on the (latitude, longitude) grid.
    """
    ds = xr.open_dataset(temp_file, decode_times=False)
    
    temp = ds["TEMP"]  
    depth = ds["depth"]  # Profondeur (en metres)
    lat = ds["latitude"]
    lon = ds["longitude"]
    # Physical constants
    RHO = 1030  
    CP = 3980   
    temp = temp.data[:,:,:,:]


    dz = np.diff(depth, prepend=0)  
        
    lat_values = lat.values
    lon_values = lon.values

    lat_diff = np.abs(np.diff(lat_values, append=lat_values[-1]))
    lon_diff = np.abs(np.diff(lon_values, append=lon_values[-1]))

    earth_radius = 6371000  # Rayon moyen de la Terre en metres
    cell_area = (
        (np.radians(lat_diff)[:, None]) *  
        (np.radians(lon_diff)[None, :]) *  # Difference de longitude (rad)
        (earth_radius ** 2) *             # Rayon de la Terre
        np.cos(np.radians(lat_values))[:, None]  # Correction de la courbure terrestre
    )
    
    ohc = RHO * CP * (temp * dz[:, np.newaxis, np.newaxis]) * cell_area

    ohc_total = sum(ohc[0,0:np.where(depth==zmax)[0][0],:,:])
    

    ohc_dataset = xr.Dataset(
        {
            "OHC": (["latitude", "longitude"], ohc_total),
        },
        coords={
            "latitude": lat.values,
            "longitude": lon.values,
        }
    )

    ohc_dataset["OHC"].attrs = {
        "long_name": "Ocean Heat Content anomaly",
        "units": "Joules",
        "description": "Total ocean heat content per grid cell, integrated over depth"
    }

    ohc_dataset.attrs = {
        "title": "Ocean Heat Content Calculation",
        "source": "Computed from temperature field using rho=1030 kg/m3 and cp=3980 J/kg/K",
        "clim": "WOA23 interpolated on ISAS grid"
    }

    return ohc_dataset


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate Ocean Heat Content from CORA temperature fields.")
    parser.add_argument("--inputdir",  required=True, help="Root input directory (year subdirectories expected inside)")
    parser.add_argument("--outputdir", required=True, help="Output directory for OHC NetCDF files")
    parser.add_argument("--year",  required=True, type=int, help="Year to process (e.g. 2010)")
    parser.add_argument("--month", required=True, help="Month to process, zero-padded (e.g. 01)")
    args = parser.parse_args()

    zmax = 2000
    os.makedirs(args.outputdir, exist_ok=True)

    temp_file = os.path.join(args.inputdir, str(args.year),
                             f"OA_CORA5.2_{args.year}{args.month}15_fld_TEMP.nc")
    out_file  = os.path.join(args.outputdir,
                             f"OHC_ISASgrid_yearly_yearlyclim_{args.year}_{args.month}.nc")

    print(f"Processing: {temp_file}")
    ohc = calculate_ohc(temp_file, zmax)
    ohc.to_netcdf(out_file)
    print(f"Saved: {out_file}")
