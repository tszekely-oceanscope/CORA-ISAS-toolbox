import numpy as np
import gsw
import netCDF4 as nc
import os
import argparse

CORA_VER = "OA_CORA5.2_"
OUTPUT_DIR = "/home/datawork-coriolis-cora-s/CORA_release/CORA-production/Datarmor/ISAS/Validation/Steric/SSH"


def compute_steric_height(TEMP, PSAL, depth, latitude):
    """Integrate halosteric height anomaly from surface to 2000 m.

    Uses a reference water mass of SP=35 psu, CT=0°C at each depth level.
    The density anomaly (rhoref - rho) is integrated vertically using
    layer thicknesses derived from the depth grid (central differences,
    with 2 m at the surface and 10 m at the deepest level).

    Parameters
    ----------
    TEMP : ndarray, shape (depth, lat, lon)
        Conservative temperature [°C] from climatology.
    PSAL : ndarray, shape (depth, lat, lon)
        Practical salinity [psu] from CORA monthly field.
    depth : ndarray, shape (depth,)
        Depth levels [m], positive downward.
    latitude : ndarray, shape (lat,)
        Latitude vector [°N].

    Returns
    -------
    SSH : ndarray, shape (lat, lon)
        Halosteric height anomaly [mm], integrated 0–2000 m.
    """
    dz = np.empty_like(depth)
    dz[1:-1] = (depth[2:] - depth[:-2]) / 2
    dz[0] = 2
    dz[-1] = 10

    SSH = None
    for zz in range(len(depth)):
        if depth[zz] >= 2000:
            continue
        if depth[zz] == 1980:
            dz[zz] = 10

        temp = TEMP[zz, :, :]
        psal = PSAL[zz, :, :]

        lat2d = np.repeat(latitude[:, np.newaxis], psal.shape[1], axis=1)
        p = gsw.p_from_z(-depth[zz], lat2d)

        SR = gsw.SR_from_SP(psal)
        CT = gsw.CT_from_t(SR, temp, p)

        rhoref = gsw.rho(gsw.SR_from_SP(35.0), 0.0, p)
        rho = gsw.rho(SR, CT, p)

        ts_ssh = ((rhoref - rho) * dz[zz] / rhoref) * 1000
        ts_ssh[np.isnan(ts_ssh)] = 0

        SSH = ts_ssh if SSH is None else SSH + ts_ssh

    return SSH


def write_output(output_file, SSH, latitude, longitude):
    """Write halosteric height field to a compressed NetCDF4 file.

    Parameters
    ----------
    output_file : str
        Destination path for the output NetCDF file.
    SSH : ndarray, shape (lat, lon)
        Halosteric height anomaly [mm].
    latitude : ndarray, shape (lat,)
        Latitude coordinates [°N].
    longitude : ndarray, shape (lon,)
        Longitude coordinates [°E].
    """
    with nc.Dataset(output_file, 'w', format='NETCDF4') as ds:
        ds.createDimension('latitude', SSH.shape[0])
        ds.createDimension('longitude', SSH.shape[1])

        lat_var = ds.createVariable('latitude', 'f4', ('latitude',))
        lon_var = ds.createVariable('longitude', 'f4', ('longitude',))
        lat_var[:] = latitude
        lon_var[:] = longitude
        lat_var.units = "degree_north"
        lon_var.units = "degree_east"

        ssh_var = ds.createVariable('steric_height_anomaly', 'f4',
                                    ('latitude', 'longitude'),
                                    zlib=True, complevel=4, fill_value=np.nan)
        ssh_var.units = "mm"
        ssh_var.long_name = "Hauteur stérique intégrée (0–2000 m)"
        ssh_var.coordinates = "latitude longitude"
        ssh_var[:, :] = SSH

        ds.title = 'Hauteur halostérique (CORA, 2005–2023)'
        ds.institution = 'OceanScope'
        ds.source = 'Calculs basés sur PSAL mensuel et température climatologique'


def main():
    """Entry point: parse CLI arguments and run the steric height computation.

    Loops over every month between Y0 and Y1 (inclusive). For each month,
    reads the corresponding CORA PSAL field and the climatological TEMP,
    computes the halosteric height anomaly integrated from 0 to 2000 m,
    and writes the result to OUTPUT_DIR as a NetCDF4 file.

    Usage
    -----
    python Calc_steric_height_CORA.py --Y0 2010 --Y1 2020 \\
        --inputpath /path/to/psal/fields --climref /path/to/clim_temp.nc
    """
    parser = argparse.ArgumentParser(description='Compute steric height from CORA PSAL and a climatological TEMP.')
    parser.add_argument('--Y0', type=int, required=True, help='Start year (inclusive)')
    parser.add_argument('--Y1', type=int, required=True, help='End year (inclusive)')
    parser.add_argument('--inputpath', required=True,
                        help='Root directory for PSAL files (expects {inputpath}/{year}/OA_CORA5.2_{year}{month}15_fld_PSAL.nc)')
    parser.add_argument('--climref', required=True,
                        help='Path to climatological TEMP file (.nc). '
                             'TEMP must have shape (depth, lat, lon) or (12, depth, lat, lon).')
    args = parser.parse_args()

    # Load climatology once
    with nc.Dataset(args.climref) as ds_clim:
        TEMP_CLIM = ds_clim.variables['TEMP'][...]
    has_month_dim = TEMP_CLIM.ndim == 4

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for year in range(args.Y0, args.Y1 + 1):
        for m in range(1, 13):
            month = f"{m:02d}"

            psal_file = os.path.join(args.inputpath, str(year),
                                     f"{CORA_VER}{year}{month}15_fld_PSAL.nc")
            print(psal_file)

            if not os.path.exists(psal_file):
                print(f"  -> not found, skipping.")
                continue

            TEMP = TEMP_CLIM[m - 1].copy() if has_month_dim else TEMP_CLIM.copy()
            TEMP[TEMP > 50] = np.nan

            with nc.Dataset(psal_file) as ds:
                latitude  = ds.variables['latitude'][:]
                longitude = ds.variables['longitude'][:]
                PSAL      = ds.variables['PSAL'][0, :, :, :].astype(float)
                depth     = ds.variables['depth'][:]

            PSAL[PSAL > 50] = np.nan

            SSH = compute_steric_height(TEMP, PSAL, depth, latitude)

            output_file = os.path.join(OUTPUT_DIR,
                                       f"CORA_halosteric_height_{year}{month}.nc")
            write_output(output_file, SSH, latitude, longitude)
            print(f"  -> written: {output_file}")


if __name__ == '__main__':
    main()
