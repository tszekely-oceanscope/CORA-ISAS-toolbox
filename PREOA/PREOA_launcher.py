import sys
import os
import argparse
import time
import logging
from datetime import datetime, timezone
from pathlib import Path

from types import SimpleNamespace
import numpy as np
from xml.etree import ElementTree as ET
import xarray as xr
import glob
import re
import numpy as np
from netCDF4 import Dataset
from scipy.io import loadmat

from netCDF4 import Dataset

# ---------------------------------------------------------------------------
#  Descripteur de variables (équivalent de import_format / vardesc)
#  Structure : name, long_name, standard_name, units, nc_type,
#              valid_min, valid_max, offset, scale_factor
# ---------------------------------------------------------------------------
VARDESC = {
    'TEMP_ANO': dict(
        long_name='anomaly',
        standard_name='sea_water_temperature_anomaly',
        units='degree_Celsius', nc_type='NC_SHORT',
        valid_min=-10.0, valid_max=10.0,
        offset=0.0, scale_factor=0.001),
    'PSAL_ANO': dict(
        long_name='anomaly',
        standard_name='sea_water_salinity_anomaly',
        units='PSS-78', nc_type='NC_SHORT',
        valid_min=-5.0, valid_max=5.0,
        offset=0.0, scale_factor=0.001),
    'PSAL_ANO_ERR': dict(
        long_name='Practical salinity Error',
        standard_name='',
        units='PSS-78', nc_type='NC_SHORT',
        valid_min=0.0, valid_max=15.0,
        offset=0.0, scale_factor=0.001),
    'PSAL_ANO_PCTVAR': dict(
        long_name='Error on salinity anomaly  (% variance)',
        standard_name='', units='%', nc_type='NC_BYTE',
        valid_min=0, valid_max=100,
        offset=0.0, scale_factor=1.0),
    'TEMP_ANO_ERR': dict(
        long_name='Temperature anomaly error',
        standard_name='',
        units='degree_Celsius', nc_type='NC_SHORT',
        valid_min=0.0, valid_max=10.0,
        offset=0.0, scale_factor=0.001),
    'TEMP_ANO_PCTVAR': dict(
        long_name='Error on temperature anomaly  (% variance)',
        standard_name='', units='%', nc_type='NC_BYTE',
        valid_min=0, valid_max=100,
        offset=0.0, scale_factor=1.0),
}

# Mapping type MATLAB → dtype numpy et type netCDF4
NC_TYPE_MAP = {
    'NC_FLOAT': ('f4', np.float32),
    'NC_DOUBLE': ('f8', np.float64),
    'NC_INT':   ('i4', np.int32),
    'NC_SHORT': ('i2', np.int16),
    'NC_BYTE':  ('i1', np.int8),
}

FILLVAL_MAP = {
    'NC_INT':   np.int32(2**31 - 1),
    'NC_SHORT': np.int16(2**15 - 1),
    'NC_BYTE':  np.int8(2**7 - 1),
    'NC_FLOAT': np.float32(9999999.0),
    'NC_DOUBLE':np.float32(9999999.0),
}

def NCW_data_hdr_and_var(ncfile_name, GLOB_ATT, STDHDR, PARAM, STDVAR):
    """
    Create a NetCDF file and write header + variables in a single open.

    N_LEVELS and N_PROF are the only dimensions. 2D variables are declared
    as (N_PROF, N_LEVELS) to compensate for the dimid inversion by
    ifx/netcdf-fortran; data arrays are transposed accordingly.

    Parameters
    ----------
    ncfile_name : str            - output NetCDF file path
    GLOB_ATT    : dict/namespace - global attributes
    STDHDR      : dict           - header metadata (juld, deph, latitude, longitude, etc.)
    PARAM       : str            - parameter name ('TEMP', 'PSAL', or 'DOXY')
    STDVAR      : dict           - variable arrays: proc, qc, data, clmn, clsd, erme,
                                   erur (optional), resi (optional)

    Returns
    -------
    msg : str - 'OK: NCW_data_hdr_and_var' or error string
    """
    msg_err = 'ERROR NCW_data_hdr_and_var'
    msg_ok  = 'OK: NCW_data_hdr_and_var'

    ref_date_time_iso, msg_iso = iso_date_time(
        STDHDR.get('reference_date_time', ''))
    if msg_iso.startswith('er'):
        return f'{msg_err}, \n   {msg_iso}'

    for ch in ('-', ':', 'T', 'Z'):
        ref_date_time_iso = ref_date_time_iso.replace(ch, '')

    from datetime import datetime, timezone
    str_now_iso = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + 'L'

    n_prof = len(np.asarray(STDHDR['juld']).flatten())
    n_lev  = len(np.asarray(STDHDR['deph']).flatten())

    FILLVALUE = np.float32(99999.0)

    PAR_PROC = f'{PARAM}_PROC'
    PAR_QC   = f'{PARAM}_QC'
    PAR_CLMN = f'{PARAM}_CLMN'
    PAR_CLSD = f'{PARAM}_CLSD'
    PAR_ERME = f'{PARAM}_ERME'
    PAR_ERUR = f'{PARAM}_ERUR'
    PAR_RESI = f'{PARAM}_RESI'

    PARAM_DEFS = {
        'TEMP': dict(
            units='degree_Celsius',
            long_name='Temperature (T90) (interpolated on Z_levels)',
            std_name='sea_water_temperature',
            par_valid=[-3.0, 40.0], err_valid=[0.001, 10.0],
            std_valid=[0.001, 10.0]),
        'PSAL': dict(
            units='none',
            long_name='Salinity (S78 - PSS) (interpolated on Z_levels)',
            std_name='sea_water_salinity',
            par_valid=[3.0, 60.0], err_valid=[0.001, 10.0],
            std_valid=[0.001, 10.0]),
        'DOXY': dict(
            units='micromole/kg',
            long_name='Disolved oxygen (interpolated on Z_levels)',
            std_name='Disolved oxygen',
            par_valid=[0.0, 1400.0], err_valid=[0.01, 100.0],
            std_valid=[0.01, 1400.0]),
    }

    if PARAM not in PARAM_DEFS:
        return f'W: NCW_data_hdr_and_var: Parameter {PARAM} not defined'

    p         = PARAM_DEFS[PARAM]
    UNITS     = p['units']
    LONG_NAME = p['long_name']
    STND_NAME = p['std_name']
    PAR_VALID = p['par_valid']
    ERR_VALID = p['err_valid']
    STD_VALID = p['std_valid']

    has_optional = ('erur' in STDVAR and STDVAR['erur'] is not None
                    and np.asarray(STDVAR['erur']).size > 0
                    and 'resi' in STDVAR and STDVAR['resi'] is not None
                    and np.asarray(STDVAR['resi']).size > 0)
    # ------------------------------------------------------------------
    #  Prépare les tableaux — force tous à shape (n_lev, n_prof)
    # ------------------------------------------------------------------
    def ensure_shape(arr, n_lev, n_prof):
        """Garantit shape (n_lev, n_prof) en paddant avec NaN si nécessaire."""
        a = np.asarray(arr, dtype=float)
        if a.ndim == 1:
            a = a[:, np.newaxis]
        nz, np_ = a.shape
        if nz < n_lev:
            pad = np.full((n_lev - nz, np_), np.nan)
            a = np.vstack([a, pad])
        elif nz > n_lev:
            a = a[:n_lev, :]
        return a

    tab_data = ensure_shape(STDVAR['data'], n_lev, n_prof)
    tab_erme = ensure_shape(STDVAR['erme'], n_lev, n_prof)
    tab_clmn = ensure_shape(STDVAR['clmn'], n_lev, n_prof)
    tab_clsd = ensure_shape(STDVAR['clsd'], n_lev, n_prof)
    tab_qcc  = convqc(ensure_shape(STDVAR['qc'], n_lev, n_prof))

    isok_dat = np.isfinite(tab_data)
    isok_cli = np.isfinite(tab_clmn)
    isok_csd = np.isfinite(tab_clsd)
    inok     = ~(isok_dat & isok_cli & isok_csd)

    if np.any(inok):
        tab_qcc[inok]  = 9
        tab_data[inok] = FILLVALUE
        tab_erme[inok] = FILLVALUE
        tab_clmn[inok] = FILLVALUE
        tab_clsd[inok] = FILLVALUE

    if has_optional:
        tab_erur = ensure_shape(STDVAR['erur'], n_lev, n_prof)
        tab_resi = ensure_shape(STDVAR['resi'], n_lev, n_prof)
        if np.any(inok):
            tab_erur[inok] = FILLVALUE
            tab_resi[inok] = FILLVALUE
    """
    # ------------------------------------------------------------------
    #  Prépare les tableaux
    # ------------------------------------------------------------------
    """

    proc = np.asarray(STDVAR['proc']).flatten()

    def gatt(name, default=''):
        if hasattr(GLOB_ATT, name):
            return getattr(GLOB_ATT, name)
        return GLOB_ATT.get(name, default) if isinstance(GLOB_ATT, dict) \
               else default

    try:
        with Dataset(ncfile_name, 'w', format='NETCDF3_64BIT') as nc:

            # ----------------------------------------------------------
            #  Dimensions — N_LEVELS=1, N_PROF=2
            # ----------------------------------------------------------
            nc.createDimension('N_LEVELS', n_lev)   # ID 1
            nc.createDimension('N_PROF',   n_prof)  # ID 2
# ----------------------------------------------------------
            #  Dimensions STRING — après N_LEVELS et N_PROF
            # ----------------------------------------------------------
            nc.createDimension('DATE_TIME', 16)
            nc.createDimension('STRING64',  64)
            nc.createDimension('STRING32',  32)
            nc.createDimension('STRING16',  16)
            nc.createDimension('STRING8',    8)
            nc.createDimension('STRING4',    4)
            nc.createDimension('STRING2',    2)

            # ----------------------------------------------------------
            #  Variables char header
            # ----------------------------------------------------------
            v = nc.createVariable('REFERENCE_DATE_TIME', 'S1', ('DATE_TIME',))
            v.long_name   = 'UTC-Date-time of reference for Julian days'
            v.conventions = 'ISO 8601'

            v = nc.createVariable('DATA_TYPE', 'S1', ('STRING16',))
            v.long_name = 'Data type'

            v = nc.createVariable('PLATFORM_NUMBER', 'S1', ( 'N_PROF','STRING8'))
            v.long_name   = 'Float unique identifier'
            v.conventions = 'WMO float identifier'

            v = nc.createVariable('WMO_INST_TYPE', 'S1', ( 'N_PROF','STRING4'))
            v.long_name   = 'Coded instrument type'
            v.conventions = 'WMO code table 1770 - instrument type'

            v = nc.createVariable('DC_REFERENCE', 'S1', ( 'N_PROF','STRING32'))
            v.long_name   = 'Unique identifier for profile/in data centre'
            v.conventions = 'Data centre convention'
            # ----------------------------------------------------------
            #  Variables 1D
            # ----------------------------------------------------------
            v = nc.createVariable('DEPH', 'f8', ('N_LEVELS',))
            v.standard_name = 'depth'
            v.units         = 'm'
            v.positive      = 'down'
            v.valid_min     = 0.0
            v.valid_max     = 10000.0

            v = nc.createVariable('JULD', 'f8', ('N_PROF',))
            v.standard_name = 'time'
            v.long_name     = 'Julian day (UTC) relative to REFERENCE_DATE_TIME'
            v.units         = 'days since REFERENCE_DATE_TIME'
            v.conventions   = 'Relative julian days with decimal part (as parts of day)'

            v = nc.createVariable('LATITUDE', 'f8', ('N_PROF',))
            v.standard_name = 'latitude'
            v.long_name     = 'Latitude of the station, best estimate'
            v.units         = 'degree_north'
            v.valid_min     = -90.0
            v.valid_max     =  90.0

            v = nc.createVariable('LONGITUDE', 'f8', ('N_PROF',))
            v.standard_name = 'longitude'
            v.long_name     = 'Longitude of the station, best estimate'
            v.units         = 'degree_east'
            v.valid_min     = -180.0
            v.valid_max     =  180.0

            v = nc.createVariable(PAR_PROC, 'i1', ('N_PROF',))
            v.long_name   = 'profile processing level'
            v.conventions = 'processing status: 1=raw, 2=adjusted'
            v.valid_min   = np.int8(0)
            v.valid_max   = np.int8(5)

            # ----------------------------------------------------------
            #  Variables 2D déclarées en (N_PROF, N_LEVELS)
            #  ifx inverse les dimids → dimids(1)=N_LEVELS ✓
            # ----------------------------------------------------------
            v = nc.createVariable(PAR_QC, 'i1', ('N_PROF', 'N_LEVELS'),
                                  fill_value=np.int8(0))
            v.long_name     = 'Quality flag on interpolated variable'
            v.flag_values   = np.array([1,2,3,4,5,6,7,8,9], dtype=np.int8)
            v.flag_meanings = ('1--4:good_to_acceptable;'
                               '5:bad_quality_interpolation;'
                               '6:bad-manual_flag;7--8:not_used;'
                               '9:not_interpolated;')

            v = nc.createVariable(PARAM, 'f4', ('N_PROF', 'N_LEVELS'),
                                  fill_value=FILLVALUE)
            v.long_name     = LONG_NAME
            v.units         = UNITS
            v.standard_name = STND_NAME
            v.valid_min     = np.float32(PAR_VALID[0])
            v.valid_max     = np.float32(PAR_VALID[1])

            v = nc.createVariable(PAR_CLMN, 'f4', ('N_PROF', 'N_LEVELS'),
                                  fill_value=FILLVALUE)
            v.long_name = 'Climatology mean for profile'
            v.units     = UNITS
            v.valid_min = np.float32(PAR_VALID[0])
            v.valid_max = np.float32(PAR_VALID[1])

            v = nc.createVariable(PAR_CLSD, 'f4', ('N_PROF', 'N_LEVELS'),
                                  fill_value=FILLVALUE)
            v.long_name = 'Climatology standard deviation for profile'
            v.units     = UNITS
            v.valid_min = np.float32(STD_VALID[0])
            v.valid_max = np.float32(STD_VALID[1])

            v = nc.createVariable(PAR_ERME, 'f4', ('N_PROF', 'N_LEVELS'),
                                  fill_value=FILLVALUE)
            v.long_name = 'Measurement error'
            v.units     = UNITS
            v.valid_min = np.float32(ERR_VALID[0])
            v.valid_max = np.float32(ERR_VALID[1])

            if 1:
                v = nc.createVariable(PAR_ERUR, 'f4', ('N_PROF', 'N_LEVELS'),
                                      fill_value=FILLVALUE)
                v.long_name = 'Error from unresolved scales'
                v.units     = UNITS
                v.valid_min = np.float32(ERR_VALID[0])
                v.valid_max = np.float32(ERR_VALID[1])

                v = nc.createVariable(PAR_RESI, 'f4', ('N_PROF', 'N_LEVELS'),
                                      fill_value=FILLVALUE)
                v.long_name = 'Residual'
                v.units     = UNITS
                v.valid_min = np.float32(-2 * ERR_VALID[1])
                v.valid_max = np.float32(+2 * ERR_VALID[1])

            # ----------------------------------------------------------
            #  Attributs globaux
            # ----------------------------------------------------------
            nc.Conventions            = 'CF-1.4'
            nc.title                  = gatt('TITLE')
            nc.institution            = gatt('INSTITUTION')
            nc.project_name           = gatt('PROJECT_NAME')
            nc.data_manager           = gatt('DATA_MANAGER')
            nc.software_version       = gatt('SOFTWARE_VERSION')
            nc.references             = gatt('REFERENCES')
            nc.history                = f'{str_now_iso} : Creation'
            nc.reference_date_time    = ref_date_time_iso[:16].ljust(16)
            nc.data_type              = str(STDHDR.get('data_type', ''))

            # ----------------------------------------------------------
            #  Écriture 1D
            # ----------------------------------------------------------
            nc.variables['DEPH'][:]      = np.asarray(STDHDR['deph']).flatten()
            nc.variables['JULD'][:]      = np.asarray(STDHDR['juld']).flatten()
            nc.variables['LATITUDE'][:]  = np.asarray(STDHDR['latitude']).flatten()
            nc.variables['LONGITUDE'][:] = np.asarray(STDHDR['longitude']).flatten()
            nc.variables[PAR_PROC][:]    = proc.astype(np.int8)
            # Char variables
            rdt = ref_date_time_iso[:16].ljust(16)
            nc.variables['REFERENCE_DATE_TIME'][:] = np.array(
                list(rdt), dtype='S1')
            dc_ref = STDHDR.get('dc_reference', '')
            nc.variables['PLATFORM_NUMBER'][:] = data_hdr_format(
                STDHDR.get('platform_number', ''), 8, n_prof).T   
            nc.variables['WMO_INST_TYPE'][:] = data_hdr_format(
                STDHDR.get('wmo_inst_type', ''), 4, n_prof).T   
            nc.variables['DC_REFERENCE'][:] = data_hdr_format(
                STDHDR.get('dc_reference', ''), 32, n_prof).T   
            #result = data_hdr_format(dc_ref, 32, n_prof)
            # ----------------------------------------------------------
            #  Écriture 2D — transpose (N_LEVELS, N_PROF) → (N_PROF, N_LEVELS)
            #  car les variables sont déclarées en (N_PROF, N_LEVELS)
            # ----------------------------------------------------------
            nc.variables[PAR_QC][:]   = tab_qcc.T.astype(np.int8)
            nc.variables[PARAM][:]    = tab_data.T.astype(np.float32)
            nc.variables[PAR_ERME][:] = tab_erme.T.astype(np.float32)
            nc.variables[PAR_CLMN][:] = tab_clmn.T.astype(np.float32)
            nc.variables[PAR_CLSD][:] = tab_clsd.T.astype(np.float32)

            if has_optional:
                nc.variables[PAR_ERUR][:] = tab_erur.T.astype(np.float32)
                nc.variables[PAR_RESI][:] = tab_resi.T.astype(np.float32)

    except Exception as e:
        return f'{msg_err}, on file: {ncfile_name} — {e}'

    return msg_ok


def _get_gatt(GLOB_ATT, name, default=''):
    """Unified getter for GLOB_ATT supporting both dict and SimpleNamespace."""
    if isinstance(GLOB_ATT, dict):
        return GLOB_ATT.get(name, default)
    return getattr(GLOB_ATT, name, default)


def _iso_date_time(s):
    """Parse various date string formats and return an ISO 8601 'yyyy-mm-ddTHH:MM:SSZ' string."""
    s = str(s).strip().replace('-', '').replace(':', '').replace('T', '') \
               .replace('Z', '').replace(' ', '')[:14].ljust(14, '0')
    try:
        datetime.strptime(s, '%Y%m%d%H%M%S')
        return f'{s[0:4]}-{s[4:6]}-{s[6:8]}T{s[8:10]}:{s[10:12]}:{s[12:14]}Z', 'ok'
    except ValueError as e:
        return '', f'error: {e}'


def _get_vardesc(var_name):
    """Return the variable descriptor dict from VARDESC, or None if not found."""
    return VARDESC.get(var_name)


def _make_fillvalue(nc_type):
    """Return the fill value scalar for a given NetCDF type string."""
    return FILLVAL_MAP.get(nc_type, np.float32(9999999.0))


def _apply_scaling_and_fill(data, valid_min, valid_max, fillvalue,
                             nc_type, offset, scale_factor, scaling):
    """Apply range clipping, optional linear scaling, dtype cast, and fill-value masking."""
    out = np.array(data, dtype=float)
    inok = ~np.isfinite(out) | (out > valid_max) | (out < valid_min)

    if scaling:
        out = (out - offset) / scale_factor
        nc_t, np_t = NC_TYPE_MAP.get(nc_type, ('f4', np.float32))
        out = out.astype(np_t)
    else:
        _, np_t = NC_TYPE_MAP.get(nc_type, ('f4', np.float32))
        out = out.astype(np_t)

    out[inok] = fillvalue
    return out, inok

def _write_field_var(nc, var_name, GLOB_ATT, dimids,
                     data, COMPRESS, suffix=''):
    """
    Create and write a scaled/filled variable in an open NetCDF dataset.

    Parameters
    ----------
    nc       : netCDF4.Dataset - open dataset in write mode
    var_name : str             - variable name as defined in VARDESC
    GLOB_ATT : dict/namespace  - global attributes (kept for API symmetry)
    dimids   : tuple           - NetCDF dimension names for the variable
    data     : array-like      - raw data to scale and write; None → variable created but not filled
    COMPRESS : dict            - compression kwargs passed to createVariable
    suffix   : str             - unused, reserved for future use

    Returns
    -------
    None on success, or an error string on failure.
    """
    vd = _get_vardesc(var_name)
    if vd is None:
        return f'ERROR: variable {var_name} not in VARDESC'

    nc_type      = vd['nc_type']
    nc_str, _    = NC_TYPE_MAP.get(nc_type, ('f4', np.float32))
    fillvalue    = _make_fillvalue(nc_type)
    offset       = vd['offset']
    scale_factor = vd['scale_factor']
    scaling      = not (scale_factor == 1.0 and offset == 0.0)

    valid_min = vd['valid_min']
    valid_max = vd['valid_max']
    if scaling:
        valid_min = int(np.floor((valid_min - offset) / scale_factor))
        valid_max = int(np.ceil( (valid_max - offset) / scale_factor))

    long_name = 'anomaly' if var_name.endswith('_ANO') else vd['long_name']

    v = nc.createVariable(var_name, nc_str, dimids,
                          fill_value=fillvalue, **COMPRESS)
    v.long_name     = long_name
    v.standard_name = vd.get('standard_name', '')
    v.units         = vd['units']
    v.valid_min     = valid_min
    v.valid_max     = valid_max
    if scaling:
        v.add_offset   = np.float32(offset)
        v.scale_factor = np.float32(scale_factor)

    if data is not None and np.asarray(data).size > 0:
        arr = np.asarray(data, dtype=float)

        # Transpose (lon, lat, dep) → (time, dep, lat, lon) pour stockage C
        # Le Fortran écrit en (lon, lat, dep, time) = ordre Fortran
        if arr.ndim == 3:
            arr = arr.transpose(2, 1, 0)        # (lon,lat,dep) → (dep,lat,lon)
            arr = arr[np.newaxis, :, :, :]       # → (1, dep, lat, lon)
        elif arr.ndim == 4:
            arr = arr.transpose(3, 2, 1, 0)     # (lon,lat,dep,t) → (t,dep,lat,lon)

        out, _ = _apply_scaling_and_fill(
            arr, valid_min, valid_max, fillvalue,
            nc_type, offset, scale_factor, scaling)
        nc.variables[var_name][:] = out

    return None

def NCW_OA_field(ncfile_name, PARAM, GLOB_ATT,
                 longitude, latitude, depth, DATE_EST,
                 FIELD=None, ERROR=None, PCTVAR=None):
    """
    Crée et écrit un fichier '.fld' NetCDF contenant un champ 3D grillé
    (lon, lat, depth, time).

    Parameters
    ----------
    ncfile_name : str
    PARAM       : str   - nom du paramètre (ex: 'TEMP_ANO')
    GLOB_ATT    : dict/namespace - attributs globaux
    longitude   : array 1D
    latitude    : array 1D
    depth       : array 1D
    DATE_EST    : str   - date 'yyyymmdd'
    FIELD       : array 3D optionnel (lon, lat, depth)
    ERROR       : array 3D optionnel
    PCTVAR      : array 3D optionnel

    Returns
    -------
    msg_error : str
    """

    fct_name = 'NCW_OA_field'
    msg_err  = f'ERROR:{fct_name}'
    msg_ok   = f'OK:{fct_name}'

    COMPRESS = dict(zlib=True, complevel=9, shuffle=True)

    # ------------------------------------------------------------------
    #  Date de référence
    # ------------------------------------------------------------------
    ref_str = _get_gatt(GLOB_ATT, 'reference_date_time', '19500101000000')
    ref_date_time_iso, msg_iso = _iso_date_time(ref_str)
    if msg_iso.startswith('er'):
        return f'{msg_err}\n   {msg_iso}'

    # Date courte pour datenum
    ref_short = ref_date_time_iso[:4] + ref_date_time_iso[5:7] + ref_date_time_iso[8:10]

    from datetime import date as _date
    def _datenum(yyyymmdd_str):
        d = datetime.strptime(yyyymmdd_str, '%Y%m%d').date()
        return (_date(d.year, d.month, d.day) - _date(1, 1, 1)).days + 367

    jul_rel = float(_datenum(DATE_EST) - _datenum(ref_short))

    # ------------------------------------------------------------------
    #  Vérifie le descripteur de variable
    # ------------------------------------------------------------------
    var_name = PARAM
    vd = _get_vardesc(var_name)
    if vd is None:
        return (f'{msg_err}\n  on file {ncfile_name}\n'
                f'  Variable {PARAM} not documented in VARDESC')

    # ------------------------------------------------------------------
    #  Crée le fichier
    # ------------------------------------------------------------------
    try:
        with Dataset(ncfile_name, 'w', format='NETCDF4') as nc:

            # ----------------------------------------------------------
            #  Attributs globaux
            # ----------------------------------------------------------
            str_now = datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')
            nc.Conventions      = 'CF-1.4'
            nc.title            = _get_gatt(GLOB_ATT, 'TITLE')
            nc.institution      = _get_gatt(GLOB_ATT, 'INSTITUTION')
            nc.project_name     = _get_gatt(GLOB_ATT, 'PROJECT_NAME')
            nc.data_manager     = _get_gatt(GLOB_ATT, 'DATA_MANAGER')
            nc.software_version = _get_gatt(GLOB_ATT, 'SOFTWARE_VERSION')
            nc.references       = _get_gatt(GLOB_ATT, 'REFERENCES')
            nc.reference_date_time = ref_date_time_iso
            nc.history          = f'{str_now} : Creation'
            # ----------------------------------------------------------
            #  Dimensions
            # ----------------------------------------------------------
            lon = np.asarray(longitude).flatten()
            lat = np.asarray(latitude).flatten()
            dep = np.asarray(depth).flatten()

            if len(lat) == 0 or len(lon) == 0:
                return f'{msg_err}: empty lat ({len(lat)}) or lon ({len(lon)}) for {ncfile_name}'

            # time seul en unlimited — les autres toutes fixes
            nc.createDimension('time',      None)      # unlimited en premier
            nc.createDimension('depth',     len(dep))
            nc.createDimension('latitude',  len(lat))
            nc.createDimension('longitude', len(lon))
            # ----------------------------------------------------------
            #  Coordonnée X — longitude
            # ----------------------------------------------------------
            v = nc.createVariable('longitude', 'f4', ('longitude',), **COMPRESS)
            v.standard_name = 'longitude'
            v.units         = 'degrees_east'
            v.valid_min     = np.float32(-180)
            v.valid_max     = np.float32(180)
            v.axis          = 'X'
            nc.variables['longitude'][:] = lon.astype(np.float32)

            # ----------------------------------------------------------
            #  Coordonnée Y — latitude
            # ----------------------------------------------------------
            v = nc.createVariable('latitude', 'f4', ('latitude',), **COMPRESS)
            v.standard_name = 'latitude'
            v.units         = 'degrees_north'
            v.valid_min     = np.float32(-90)
            v.valid_max     = np.float32(90)
            v.axis          = 'Y'
            nc.variables['latitude'][:] = lat.astype(np.float32)

            # ----------------------------------------------------------
            #  Coordonnée Z — depth ou sig0
            # ----------------------------------------------------------
            type_zcoo = _get_gatt(GLOB_ATT, 'ZCOO', 'depth')

            if type_zcoo == 'sig0':
                v = nc.createVariable('sig0', 'f4', ('depth',), **COMPRESS)
                v.standard_name = 'potential density'
                v.units         = 'kg*m-3'
                v.valid_min     = np.float32(0)
                v.valid_max     = np.float32(70)
                v.axis          = 'Z'
            else:
                v = nc.createVariable('depth', 'f4', ('depth',), **COMPRESS)
                v.standard_name = 'depth'
                v.units         = 'm'
                v.positive      = 'down'
                v.valid_min     = np.float32(0)
                v.valid_max     = np.float32(12000)
                v.axis          = 'Z'

            nc.variables[type_zcoo if type_zcoo == 'sig0' else 'depth'][:] = \
                dep.astype(np.float32)

            # ----------------------------------------------------------
            #  Coordonnée T — time
            # ----------------------------------------------------------
            v = nc.createVariable('time', 'f4', ('time',), **COMPRESS)
            v.standard_name = 'time'
            v.units         = f'days since {ref_date_time_iso}'
            v.axis          = 'T'
            nc.variables['time'][0] = np.float32(jul_rel)

            # ----------------------------------------------------------
            #  Variable principale PARAM
            # ----------------------------------------------------------
            dimids_4d = ('time', 'depth', 'latitude', 'longitude')
            err = _write_field_var(nc, PARAM, GLOB_ATT, dimids_4d,
                                   FIELD, COMPRESS)
            if err:
                return errVARDESC

            # ----------------------------------------------------------
            #  Variable ERROR (optionnelle)
            # ----------------------------------------------------------
            #if ERROR is not None:
            err_var_name = f'{PARAM}_ERR'
            vd_err = _get_vardesc(err_var_name)
            if vd_err is not None:
                err = _write_field_var(nc, err_var_name, GLOB_ATT,
                                        dimids_4d, ERROR, COMPRESS)
                if err:
                    return err

            # ----------------------------------------------------------
            #  Variable PCTVAR (optionnelle)
            # ----------------------------------------------------------
            pct_var_name = f'{PARAM}_PCTVAR'
            vd_pct = _get_vardesc(pct_var_name)
            if vd_pct is not None:
                err = _write_field_var(nc, pct_var_name, GLOB_ATT,
                                        dimids_4d, PCTVAR, COMPRESS)
                if err:
                    return err

    except Exception as e:
        return f'{msg_err} on file {ncfile_name} — {e}'

    return msg_ok


def normalize_tab(tab, n_prof):
    """Ensure array shape is (nlevels, nprofs), transposing if needed."""
    tab = np.ma.atleast_2d(tab)
    if tab.shape[0] == n_prof and tab.shape[1] != n_prof:
        tab = tab.T
    elif tab.shape[0] != n_prof and tab.shape[1] == n_prof:
        pass  # déjà correct
    return tab

def NCR_data(fnam_data, PARAM, iopt_nan=0, pltnum='0',
             jlim=None, tab_msk=None, lat_msk=None, lon_msk=None):
    """
    Read profiles from a STD or dat NetCDF file, filtered by platform, time, and position.

    Parameters
    ----------
    fnam_data : str          - full path to the NetCDF file
    PARAM     : str          - parameter name ('TEMP', 'PSAL', 'DOXY', ...)
    iopt_nan  : int          - if 1, replace fill_value with NaN in output arrays
    pltnum    : str          - platform number to select; '0' means no filter
    jlim      : array-like   - [jmin, jmax] temporal bounds in Julian days
    tab_msk   : array-like   - 2D geographic mask (lat × lon)
    lat_msk   : array-like   - latitude axis of the mask
    lon_msk   : array-like   - longitude axis of the mask

    Returns
    -------
    STDHDR : dict - header metadata and global attributes
    STDVAR : dict - variable arrays for selected profiles
    """

    # ------------------------------------------------------------------
    #  Initialisation
    # ------------------------------------------------------------------
    STDHDR = {'nbdat': 0}
    STDVAR = {'data': np.array([]),
              'erur': np.array([]),
              'resi': np.array([])}

    with Dataset(fnam_data, 'r') as nc:
        var_list = list(nc.variables.keys())

        PAR_DATA = PARAM
        PAR_PROC = f'{PARAM}_PROC'
        PAR_QC   = f'{PARAM}_QC'

        # ------------------------------------------------------------------
        #  Vérifie que le paramètre existe
        # ------------------------------------------------------------------
        if PARAM not in var_list:
            return STDHDR, STDVAR

        # ------------------------------------------------------------------
        #  Détermine le type de fichier (RAW ou STD)
        # ------------------------------------------------------------------
        file_type = 'RAW' if 'STATION_PARAMETERS' in var_list else 'STD'

        fldHDR_0 = ['reference_date_time', 'data_type']
        fldHDR_pos, fldHDR_dat, fldVAR_dat = flddat_init_copernicus_PREOA()
        nb_fldHDRdat = len(fldHDR_dat)

        if file_type == 'RAW':
            fldHDR_pos = fldHDR_pos[:3]   # pas de depth (c'est une variable)
            fldVAR_dat = ['qc', 'data']
            nampar     = [PAR_QC, PAR_DATA]
            fldHDR_datall = fldHDR_0 + fldHDR_pos + fldHDR_dat
            nam_varHDR    = [f.upper() for f in fldHDR_datall]
            varid_ERUR = -1
            varid_RESI = -1
            is_V6 = 0

        else:
            # Détecte version V6 (présence de PAR_PROC)
            is_V6 = 1 if PAR_PROC in var_list else 0

            if is_V6:
                PAR_CLMN = f'{PARAM}_CLMN'
                PAR_CLSD = f'{PARAM}_CLSD'
                PAR_ERME = f'{PARAM}_ERME'
                PAR_ERUR = f'{PARAM}_ERUR'
                PAR_RESI = f'{PARAM}_RESI'
            else:
                PAR_CLMN = f'{PARAM}_CLIM'
                PAR_CLSD = f'{PARAM}_CLIM_STD'
                PAR_ERME = f'{PARAM}_ERR_ME'
                PAR_ERUR = f'{PARAM}_ERR_UR'
                PAR_RESI = f'{PARAM}_RESID'

            nampar = [PAR_PROC, PAR_QC, PAR_DATA, PAR_CLMN,
                      PAR_CLSD, PAR_ERME, PAR_ERUR, PAR_RESI]

            fldHDR_datall = fldHDR_0 + fldHDR_pos + fldHDR_dat
            nam_varHDR    = [f.upper() for f in fldHDR_datall]
            if not is_V6:
                # index 5 correspond à 'deph' → 'DEPTH' pour ancienne version
                nam_varHDR[5] = 'DEPTH'

            varid_ERUR = 1 if PAR_ERUR in var_list else -1
            varid_RESI = 1 if PAR_RESI in var_list else -1

        nb_fldvar = len(fldVAR_dat)

        # ------------------------------------------------------------------
        #  Lecture du header
        # ------------------------------------------------------------------
        n_prof_read = len(nc.variables['LATITUDE'][:])
        for ifi, fld in enumerate(fldHDR_datall):
            nam_var = nam_varHDR[ifi]
            if nam_var and nam_var in nc.variables:
                if  nc.variables[nam_var][:].dtype !='S1':
                    STDHDR[fld] = nc.variables[nam_var][:]
                elif nam_var not in ['REFERENCE_DATE_TIME','DATA_TYPE'] :
                        
                        arr = nc.variables[nam_var][:]
                        
                        final_array=[]
                        for zz in range(0,n_prof_read):
                            final_array.append(arr[zz,:].tobytes().decode('utf-8'))
                        # ma      sked_array possible


                        # reconstruire une string par ligne
                        STDHDR[fld] =np.array(final_array)
                else: 
                    arr = nc.variables[nam_var][:].tobytes().decode('utf-8')
                    STDHDR[fld] =arr

        

        n_prof_read = len(STDHDR['latitude'])
        # Nettoyage de reference_date_time (retire '-', ':', 'T', 'Z')
        ref = STDHDR.get('reference_date_time', np.array([]))
        if hasattr(ref, 'tobytes'):
            ref_str = ref.tobytes().decode('utf-8', errors='ignore')
        else:
            ref_str = str(ref)
        ref_clean = ''.join(c for c in ref_str if c not in '-:TZ \x00')
        STDHDR['reference_date_time'] = ref_clean

        # Transpose les tableaux de caractères (sauf data_mode)
        for fld in fldHDR_dat:
            val = STDHDR.get(fld)
            if val is None:
                continue
            #if val.dtype.kind in ('S', 'U') and fld.lower() != 'data_mode':
            #    STDHDR[fld] = val.T if val.ndim > 1 else val

        # ------------------------------------------------------------------
        #  Sélection plateforme
        # ------------------------------------------------------------------
        if pltnum == '0':
            listok = np.ones(n_prof_read, dtype=bool)
        else:
            pn = STDHDR.get('platform_number', np.array([]))
            # Convertit en liste de strings nettoyées
            if pn.ndim == 2:
                pn_str = [''.join(row.astype(str)).strip() for row in pn]
            else:
                pn_str = [str(pn).strip()]
            listok = np.array([p == pltnum.strip() for p in pn_str])

        if not np.any(listok):
            STDHDR['nbdat'] = 0
            return STDHDR, STDVAR

        # ------------------------------------------------------------------
        #  Sélection temporelle
        # ------------------------------------------------------------------
        if jlim is not None:
            juld = STDHDR['juld']
            isok_juld = (juld >= jlim[0]) & (juld < jlim[1])
            listok = listok & isok_juld

        if not np.any(listok):
            STDHDR['nbdat'] = 0
            return STDHDR, STDVAR

        # ------------------------------------------------------------------
        #  Sélection géographique + masque
        # ------------------------------------------------------------------
        if tab_msk is not None and lat_msk is not None and lon_msk is not None:
            lat = STDHDR['latitude'].copy()
            lon = STDHDR['longitude'].copy()

            isok_lat = (lat >= lat_msk[0]) & (lat <= lat_msk[-1])

            if np.any(isok_lat):
                if np.any(lon_msk > 180):
                    ineg = lon < 0
                    lon[ineg] += 360
                elif np.any(lon_msk < -180):
                    ipos = lon > 0
                    lon[ipos] -= 360
                isok_lon  = (lon >= lon_msk[0]) & (lon <= lon_msk[-1])
                isok_area = isok_lat & isok_lon
            else:
                isok_area = isok_lat

            listok = listok & isok_area

            # Vérification sur le masque
            if np.any(listok):
                isok_idx = np.where(listok)[0]
                for i in isok_idx:
                    ii = int(np.argmin(np.abs(lat_msk - lat[i])))
                    jj = int(np.argmin(np.abs(lon_msk - lon[i])))
                    if tab_msk[ii, jj] == 0:
                        listok[i] = False

        if not np.any(listok):
            STDHDR['nbdat'] = 0
            return STDHDR, STDVAR

        STDHDR['nbdat'] = int(np.sum(listok))
        listok = np.asarray(listok).flatten().astype(bool)
        # ------------------------------------------------------------------
        #  Lecture des données
        # ------------------------------------------------------------------
        for ifi, nam_fld in enumerate(fldVAR_dat):

            if nam_fld == 'erur':
                tab = nc.variables[PAR_ERUR][:] if varid_ERUR > 0 else np.array([])

            elif nam_fld == 'resi':
                tab = nc.variables[PAR_RESI][:] if varid_RESI > 0 else np.array([])

            elif nam_fld == 'proc':
                tab = nc.variables[PAR_PROC][:] if is_V6 else np.array([])

            else:
                var_name = nampar[ifi]
                if var_name in nc.variables:
                    tab = nc.variables[var_name][:]
                else:
                    tab = np.array([])

            # Affectation selon le type de champ
            if nam_fld == 'proc':
                STDVAR[nam_fld] = tab[listok] if tab.size > 0 else np.array([])

            elif nam_fld == 'qc':
                if tab.size > 0:
                    tab = normalize_tab(tab, n_prof_read)
                    if file_type == 'RAW':
                        STDVAR[nam_fld] = tab[:, listok].astype(float) - 48
                    else:
                        STDVAR[nam_fld] = tab[:, listok]
                else:
                    STDVAR[nam_fld] = np.array([])

            else:
                if tab.size > 0:
                    tab = normalize_tab(tab, n_prof_read)
                    if iopt_nan:
                        var_name = nampar[ifi]
                        if var_name in nc.variables:
                            fv = nc.variables[var_name]._FillValue
                            tab = tab.astype(float)
                            tab[tab == fv] = np.nan
                    STDVAR[nam_fld] = tab[:, listok].astype(float)
                else:
                    STDVAR[nam_fld] = np.array([])
    # ------------------------------------------------------------------
    #  Extraction du header pour les profils sélectionnés
    # ------------------------------------------------------------------

    for fld in fldHDR_datall:
        val = STDHDR.get(fld)
        if val is None:
            continue
        # Ne sélectionne que les arrays NumPy de taille n_prof_read
        if not isinstance(val, (np.ndarray, np.ma.MaskedArray)):
            continue
        val_flat = np.asarray(val)
        if val_flat.shape[0] == n_prof_read:
            STDHDR[fld] = val[listok]

    STDHDR['num_sel'] = np.where(listok)[0]

    return STDHDR, STDVAR


def PREOA_select(OAHDR, OAVAR, liste_ok):
    """
    Subset OAHDR/OAVAR to the profiles indicated by liste_ok.

    Parameters
    ----------
    OAHDR    : dict                    - header and metadata
    OAVAR    : dict                    - variable arrays
    liste_ok : bool array or int array - profiles to keep (boolean mask or indices)

    Returns
    -------
    OAHDR_sel : dict
    OAVAR_sel : dict
    """

    liste_ok = np.asarray(liste_ok).flatten()

    # Champs à copier tels quels (scalaires / non-indexables par profil)
    fields_copy = {'nbdat', 'deph', 'reference_date_time', 'data_type'}

    OAHDR_sel = {}
    OAVAR_sel = {}

    for key, val in OAHDR.items():
        if key in fields_copy:
            # Copie directe sans sélection
            OAHDR_sel[key] = val
        else:
            if val is None:
                OAHDR_sel[key] = val
                continue
            arr = np.asarray(val)
            if arr.ndim == 0:
                # Scalaire
                OAHDR_sel[key] = val
            elif arr.shape[0] == len(liste_ok) or arr.dtype.kind not in ('S', 'U'):
                # Sélection sur l'axe 0 (nprofs)
                try:
                    OAHDR_sel[key] = arr[liste_ok]
                except (IndexError, TypeError):
                    OAHDR_sel[key] = val
            else:
                OAHDR_sel[key] = val

    # ------------------------------------------------------------------
    #  Sélection des variables
    # ------------------------------------------------------------------
    for nam_fld, val in OAVAR.items():
        if val is None or (hasattr(val, 'size') and val.size == 0):
            OAVAR_sel[nam_fld] = val
            continue

        arr = np.asarray(val)

        if nam_fld == 'proc':
            # 1D : (nprofs,)
            OAVAR_sel[nam_fld] = arr[liste_ok]
        else:
            # 2D : (nlevels, nprofs)
            if arr.ndim == 2:
                OAVAR_sel[nam_fld] = arr[:, liste_ok]
            elif arr.ndim == 1:
                OAVAR_sel[nam_fld] = arr[liste_ok]
            else:
                OAVAR_sel[nam_fld] = val

    # Recalcule nbdat
    OAHDR_sel['nbdat'] = int(np.asarray(OAHDR_sel.get('juld', [])).flatten().size)

    return OAHDR_sel, OAVAR_sel


def PREOA_append(OAHDR, OAVAR, STDHDR, STDVAR):
    """
    Append new profiles (STDHDR/STDVAR) to the area accumulator (OAHDR/OAVAR).

    Handles depth-level mismatches by padding the shorter array with NaN rows
    before concatenation along the profile axis.

    Parameters
    ----------
    OAHDR  : dict - accumulated header (modified in place and returned)
    OAVAR  : dict - accumulated variables (modified in place and returned)
    STDHDR : dict - new header to append
    STDVAR : dict - new variables to append

    Returns
    -------
    OAHDR : dict
    OAVAR : dict
    """
    # Case: empty initial structure
    if OAHDR["nbdat"] == 0:
        return STDHDR.copy(), STDVAR.copy()

    # --- Fields handling ---
    fieldsHDR = [k for k in STDHDR.keys()
                 if k not in ["nbdat", "deph", "reference_date_time"]]
    fieldsVAR = list(STDVAR.keys())

    # --- Append header fields ---
    for fld in fieldsHDR:
        try: #isinstance(STDHDR[fld], str):
            # MATLAB char concatenation → list of strings or join
            arr = np.hstack((OAHDR[fld], STDHDR[fld]))
            

            OAHDR[fld] =arr
        except:
            OAHDR[fld] = np.concatenate([OAHDR[fld], STDHDR[fld]], axis=0)
                # --- Depth consistency ---
    nb_dep_st = len(np.asarray(STDHDR['deph']).flatten())
    nb_dep_oa = len(np.asarray(OAHDR['deph']).flatten())
    n_prof_oa  = int(np.asarray(OAHDR.get('juld', [])).flatten().size)
    """
    if nb_dep_st > nb_dep_oa:
        OAHDR['deph'] = np.asarray(STDHDR['deph']).flatten()
        nb_add  = nb_dep_st - nb_dep_oa

        for fld in fieldsVAR:
            val = OAVAR.get(fld)
            if val is None or np.asarray(val).size == 0:
                continue
            arr = np.asarray(val)
            if fld == 'proc':
                pass  # 1D — pas de nivaux à ajouter
            elif arr.ndim == 2:
                # shape actuelle (N_LEVELS_old, N_PROF)
                n_prof_cur = arr.shape[1]
                tab_add = np.full((nb_add, n_prof_cur), np.nan)
                if fld == 'qc':
                    tab_add = np.zeros((nb_add, n_prof_cur), dtype=np.int8)
                OAVAR[fld] = np.vstack([arr, tab_add])

    elif nb_dep_st < nb_dep_oa:
        nb_add = nb_dep_oa - nb_dep_st

        for fld in fieldsVAR:
            val = STDVAR.get(fld)
            if val is None or np.asarray(val).size == 0:
                continue
            arr = np.asarray(val)
            if fld == 'proc':
                pass
            elif arr.ndim == 2:
                n_prof_cur = arr.shape[1]
                tab_add = np.full((nb_add, n_prof_cur), np.nan)
                if fld == 'qc':
                    tab_add = np.zeros((nb_add, n_prof_cur), dtype=np.int8)
                STDVAR[fld] = np.vstack([arr, tab_add])  """  
    # --- Depth consistency ---
    nb_dep_st = len(np.asarray(STDHDR['deph']).flatten())
    nb_dep_oa = len(np.asarray(OAHDR['deph']).flatten())
    n_prof_oa  = int(np.asarray(OAHDR.get('juld', [])).flatten().size)

    if nb_dep_st > nb_dep_oa:
        OAHDR['deph'] = np.asarray(STDHDR['deph']).flatten()
        nb_add  = nb_dep_st - nb_dep_oa

        for fld in fieldsVAR:
            val = OAVAR.get(fld)
            if val is None or np.asarray(val).size == 0:
                continue
            arr = np.asarray(val)
            if fld == 'proc':
                pass  # 1D — pas de nivaux à ajouter
            elif arr.ndim == 2:
                # shape actuelle (N_LEVELS_old, N_PROF)
                n_prof_cur = arr.shape[1]
                tab_add = np.full((nb_add, n_prof_cur), np.nan)
                if fld == 'qc':
                    tab_add = np.zeros((nb_add, n_prof_cur), dtype=np.int8)
                OAVAR[fld] = np.vstack([arr, tab_add])

    elif nb_dep_st < nb_dep_oa:
        nb_add = nb_dep_oa - nb_dep_st

        for fld in fieldsVAR:
            val = STDVAR.get(fld)
            if val is None or np.asarray(val).size == 0:
                continue
            arr = np.asarray(val)
            if fld == 'proc':
                pass
            elif arr.ndim == 2:
                n_prof_cur = arr.shape[1]
                tab_add = np.full((nb_add, n_prof_cur), np.nan)
                if fld == 'qc':
                    tab_add = np.zeros((nb_add, n_prof_cur), dtype=np.int8)
                STDVAR[fld] = np.vstack([arr, tab_add])
    # --- Append variable fields ---
    for fld in fieldsVAR:
        if STDVAR[fld] is None or np.asarray(STDVAR[fld]).size == 0:
            continue

        oa_arr  = np.asarray(OAVAR[fld])
        std_arr = np.asarray(STDVAR[fld])

        if fld == 'proc':
            # 1D — concatène sur axe 0
            OAVAR[fld] = np.concatenate(
                [oa_arr.flatten(), std_arr.flatten()])
        else:
            # 2D (N_LEVELS, N_PROF) — concatène sur axe 1
            if oa_arr.ndim == 1:
                oa_arr = oa_arr[:, np.newaxis]
            if std_arr.ndim == 1:
                std_arr = std_arr[:, np.newaxis]
            OAVAR[fld] = np.concatenate([oa_arr, std_arr], axis=1)
    # --- Update count ---
    OAHDR["nbdat"] = len(OAHDR["juld"])

    return OAHDR, OAVAR

def iso_date_time(s):
    """
    Parse a date string and return it in ARGO/ISO format 'yyyymmddTHHMMSSZ' (16 chars).
    """
    # Nettoie la chaîne
    clean = str(s).strip().replace('-', '').replace(':', '') \
                  .replace('T', '').replace('Z', '').replace(' ', '') \
                  [:14].ljust(14, '0')
    try:
        datetime.strptime(clean, '%Y%m%d%H%M%S')
        # Format ARGO : yyyymmddTHHMMSSZ
        iso = f"{clean[:8]}T{clean[8:14]}Z"
        return iso, 'ok'
    except ValueError as e:
        return '', f'error: {e}'


def data_hdr_format(val, str_dim, n_prof):
    """Format a header string field into a (str_dim, n_prof) array of S1 bytes."""
    out = np.full((str_dim, n_prof), b' ', dtype='S1')

    if val is None:
        return out

    arr = np.asarray(val)

    # 1D array de strings — un élément par profil
    if arr.ndim == 1 and arr.dtype.kind in ('U', 'S', 'O'):
        for j in range(min(n_prof, len(arr))):
            s = str(arr[j]).strip()[:str_dim]
            for k, c in enumerate(s):
                out[k, j] = c.encode('utf-8')

    # 2D array de chars
    elif arr.ndim == 2:
        # Normalise en (str_dim, n_prof)
        if arr.shape[0] == n_prof and arr.shape[1] != n_prof:
            arr = arr.T  # (n_prof, nchars) → (nchars, n_prof)
        for j in range(min(n_prof, arr.shape[1])):
            col = arr[:, j]
            if col.dtype.kind == 'S':
                s = b''.join(col).decode('utf-8', errors='ignore').strip()[:str_dim]
            else:
                s = ''.join(str(c) for c in col).strip()[:str_dim]
            for k, c in enumerate(s):
                out[k, j] = c.encode('utf-8')

    return out  # shape (str_dim, n_prof)Sonnet 4.6


def NCW_data_hdr_PREOA(ncfile_name, GLOB_ATT, STDHDR):
    """
    Create and write the header (metadata) for STD/dat data files with maximum compression.

    Parameters
    ----------
    ncfile_name : str            - full path of the file to create
    GLOB_ATT    : dict/namespace - global attributes
    STDHDR      : dict           - header and metadata

    Returns
    -------
    msg_error : str - 'OK: ...' or 'ERROR: ...'
    """

    msg_err = 'ERROR NCW_data_hdr_PREOA'
    msg_ok  = 'OK: NCW_data_hdr_PREOA'

    # ------------------------------------------------------------------
    #  Formate la date de référence
    # ------------------------------------------------------------------
    ref_date_time_iso, msg_iso = iso_date_time(STDHDR.get('reference_date_time', ''))
    if msg_iso.startswith('er'):
        return f'{msg_err}, \n   {msg_iso}'

    # Supprime tirets et ':' → format ARGO 'yyyymmddHHMMSS'
    for ch in ('-', ':'):
        ref_date_time_iso = ref_date_time_iso.replace(ch, '')

    str_now_iso = datetime.utcnow().strftime('%Y%m%dT%H%M%SZ') + 'L'

    # ------------------------------------------------------------------
    #  Paramètres de compression (zlib, niveau 9)
    # ------------------------------------------------------------------
    COMPRESS = dict(zlib=True, complevel=9, shuffle=True)

    # ------------------------------------------------------------------
    #  Crée le fichier NetCDF
    # ------------------------------------------------------------------
    try:
        with Dataset(ncfile_name, 'w', format='NETCDF4_CLASSIC') as nc:

            # ----------------------------------------------------------
            #  Dimensions fixes
            # ----------------------------------------------------------
            nc.createDimension('DATE_TIME', 16)
            nc.createDimension('STRING64',  64)
            nc.createDimension('STRING32',  32)
            nc.createDimension('STRING16',  16)
            nc.createDimension('STRING8',    8)
            nc.createDimension('STRING4',    4)
            nc.createDimension('STRING2',    2)

            # ----------------------------------------------------------
            #  Dimensions variables
            # ----------------------------------------------------------
            n_prof = len(np.asarray(STDHDR['juld']).flatten())
            n_lev  = len(np.asarray(STDHDR['deph']).flatten())

            nc.createDimension('N_LEVELS', n_lev)
            nc.createDimension('N_PROF',   n_prof)

            # ----------------------------------------------------------
            #  Attributs globaux
            # ----------------------------------------------------------
            def gatt(name):
                if hasattr(GLOB_ATT, name):
                    return getattr(GLOB_ATT, name)
                return GLOB_ATT.get(name, '') if isinstance(GLOB_ATT, dict) else ''

            nc.Conventions      = 'CF-1.4'
            nc.title            = gatt('TITLE')
            nc.institution      = gatt('INSTITUTION')
            nc.project_name     = gatt('PROJECT_NAME')
            nc.data_manager     = gatt('DATA_MANAGER')
            nc.software_version = gatt('SOFTWARE_VERSION')
            nc.references       = gatt('REFERENCES')
            nc.history          = f'{str_now_iso} : Creation'

            # ----------------------------------------------------------
            #  Définition des variables
            # ----------------------------------------------------------

            # REFERENCE_DATE_TIME — char (DATE_TIME,)
            v = nc.createVariable('REFERENCE_DATE_TIME', 'S1', ('DATE_TIME',),
                                  **COMPRESS)
            v.long_name   = 'UTC-Date-time of reference for Julian days'
            v.conventions = 'ISO 8601'

            # DATA_TYPE — char (STRING16,)
            v = nc.createVariable('DATA_TYPE', 'S1', ('STRING16',),
                                  **COMPRESS)
            v.long_name = 'Data type'

            # PLATFORM_NUMBER — char (STRING8, N_PROF)
            v = nc.createVariable('PLATFORM_NUMBER', 'S1',
                                ('N_PROF', 'STRING8'), **COMPRESS)
            v.long_name   = 'Float unique identifier'
            v.conventions = 'WMO float identifier'


            # WMO_INST_TYPE — char (STRING4, N_PROF)
            v = nc.createVariable('WMO_INST_TYPE', 'S1',
                                ('N_PROF', 'STRING4'), **COMPRESS)
            v.long_name   = 'Coded instrument type'
            v.conventions = 'WMO code table 1770 - instrument type'
            # DC_REFERENCE — char (STRING32, N_PROF)
            v = nc.createVariable('DC_REFERENCE', 'S1',
                                ('N_PROF', 'STRING32'), **COMPRESS)
            v.long_name   = 'Unique identifier for profile/in data centre'
            v.conventions = 'Data centre convention'

            # JULD — double (N_PROF,)
            v = nc.createVariable('JULD', 'f8', ('N_PROF',), **COMPRESS)
            v.standard_name = 'time'
            v.long_name     = 'Julian day (UTC) relative to REFERENCE_DATE_TIME'
            v.units         = 'days since REFERENCE_DATE_TIME'
            v.conventions   = 'Relative julian days with decimal part (as parts of day)'

            # LATITUDE — double (N_PROF,)
            v = nc.createVariable('LATITUDE', 'f8', ('N_PROF',), **COMPRESS)
            v.standard_name = 'latitude'
            v.long_name     = 'Latitude of the station, best estimate'
            v.units         = 'degree_north'
            v.valid_min     = -90.0
            v.valid_max     =  90.0

            # LONGITUDE — double (N_PROF,)
            v = nc.createVariable('LONGITUDE', 'f8', ('N_PROF',), **COMPRESS)
            v.standard_name = 'longitude'
            v.long_name     = 'Longitude of the station, best estimate'
            v.units         = 'degree_east'
            v.valid_min     = -180.0
            v.valid_max     =  180.0

            # DEPH — double (N_LEVELS,)
            v = nc.createVariable('DEPH', 'f8', ('N_LEVELS',), **COMPRESS)
            v.standard_name = 'depth'
            v.units         = 'm'
            v.positive      = 'down'
            v.valid_min     =     0.0
            v.valid_max     = 10000.0

            # ----------------------------------------------------------
            #  Écriture des variables
            # ----------------------------------------------------------

            # reference_date_time : pad à 16 chars
            rdt = ref_date_time_iso[:16].ljust(16)
            nc.variables['REFERENCE_DATE_TIME'][:] = np.array(
                list(rdt), dtype='S1')

            # Tableaux de caractères header
            nc.variables['PLATFORM_NUMBER'][:] = data_hdr_format(
                STDHDR['platform_number'], 8, n_prof)

            nc.variables['WMO_INST_TYPE'][:] = data_hdr_format(
                STDHDR['wmo_inst_type'], 4, n_prof)

            nc.variables['DC_REFERENCE'][:] = data_hdr_format(
                STDHDR['dc_reference'], 32, n_prof)
            result = nc.variables['DC_REFERENCE'][:]
            print(f'data_hdr_format result shape: {result.shape}')
            print(f'data_hdr_format result[0,0]: {result[0,0]}')
            print(f'data_hdr_format result[:,0]: {result[:,0]}')
            # Variables numériques
            nc.variables['JULD'][:]      = np.asarray(STDHDR['juld']).flatten()
            nc.variables['LATITUDE'][:]  = np.asarray(STDHDR['latitude']).flatten()
            nc.variables['LONGITUDE'][:] = np.asarray(STDHDR['longitude']).flatten()

            deph = np.asarray(STDHDR.get('deph', [])).flatten()
            if deph.size > 0:
                nc.variables['DEPH'][:] = deph

    except Exception as e:
        return f'{msg_err}, on file: {ncfile_name} — {e}'

    return msg_ok
def PREOA_newdata(TYPE, preoa_mess, msg_remove,
                  STDHDR, STDVAR, OAHDR, OAVAR,
                  DM_platform, nb_platform):
    """
    Add new profiles to the area accumulators, removing CO duplicates already present in DM.

    Parameters
    ----------
    TYPE        : str  - data type prefix ('DM', 'CO', ...)
    preoa_mess  : file - open log file
    msg_remove  : str  - printf-style format string for the removal message
    STDHDR      : dict - header of the new data batch
    STDVAR      : dict - variables of the new data batch
    OAHDR       : dict - accumulated area header
    OAVAR       : dict - accumulated area variables
    DM_platform : list - list of DM platform numbers already seen
    nb_platform : int  - number of DM platforms registered so far

    Returns
    -------
    OAHDR, OAVAR, DM_platform, nb_platform
    """

    txt = TYPE[:2]

    # ------------------------------------------------------------------
    #  Récupère les numéros de plateforme sous forme de liste de strings
    # ------------------------------------------------------------------
    def get_platform_numbers(STDHDR):
        pn = STDHDR.get('platform_number')
        if pn is None:
            return []
        pn = np.asarray(pn)
        if pn.ndim == 2:
            # tableau de caractères (nprofs, nchars)
            return [''.join(row.astype(str)).strip() for row in pn]
        elif pn.ndim == 1:
            if pn.dtype.kind in ('S', 'U'):
                return [str(pn).strip()]
            # array 1D de chars → un seul profil
            return [''.join(pn.astype(str)).strip()]
        return []

    # ------------------------------------------------------------------
    #  Cas DM : enregistre les plateformes
    # ------------------------------------------------------------------
    if TYPE[:2] == 'DM':
        plt_list = get_platform_numbers(STDHDR)

        if nb_platform == 0:
            # Premier passage DM : initialise la liste
            unique_plts = list(dict.fromkeys(plt_list))  # unique, ordre préservé
            for p in unique_plts:
                DM_platform[nb_platform] = p
                nb_platform += 1
        else:
            # Passages suivants : ajoute les nouvelles plateformes
            for p in plt_list:
                if p not in DM_platform[:nb_platform]:
                    DM_platform[nb_platform] = p
                    nb_platform += 1

    # ------------------------------------------------------------------
    #  Cas non-DM : supprime les profils déjà présents en DM
    # ------------------------------------------------------------------
    else:
        if nb_platform > 0:
            plt_list = get_platform_numbers(STDHDR)
            dm_set   = set(DM_platform[:nb_platform])

            list_ok = np.array([p not in dm_set for p in plt_list], dtype=bool)
            nb_ok      = int(list_ok.sum())
            nb_removed = STDHDR['nbdat'] - nb_ok

            print(msg_remove % (nb_removed, txt), end='')
            preoa_mess.write(msg_remove % (nb_removed, txt))

            if nb_ok > 0:
                if nb_removed > 0:
                    STDHDR, STDVAR = PREOA_select(STDHDR, STDVAR, list_ok)
            else:
                STDHDR['nbdat'] = 0

    # ------------------------------------------------------------------
    #  Agrège dans la structure de zone
    # ------------------------------------------------------------------
    OAHDR, OAVAR = PREOA_append(OAHDR, OAVAR, STDHDR, STDVAR)

    return OAHDR, OAVAR, DM_platform, nb_platform

def flddat_init_copernicus_PREOA():
    """
    Return the standard field-name lists for Copernicus PREOA STD/dat files.

    Returns
    -------
    fldHDR_pos : list - position header fields: juld, latitude, longitude, deph
    fldHDR_dat : list - data header fields: dc_reference, wmo_inst_type, platform_number
    fldVAR_dat : list - variable fields: proc, qc, data, clmn, clsd, erme, erur, resi
    """

    fldHDR_pos = ['juld', 'latitude', 'longitude', 'deph']

    fldHDR_dat = ['dc_reference', 'wmo_inst_type', 'platform_number']

    fldVAR_dat = ['proc', 'qc', 'data', 'clmn', 'clsd', 'erme', 'erur', 'resi']

    return fldHDR_pos, fldHDR_dat, fldVAR_dat

def expand_ranges(range_str):
    """
    Expand a MATLAB-style range string into a flat list of integers.

    Parameters
    ----------
    range_str : str - e.g. "[1:3,5,7:9]" → [1, 2, 3, 5, 7, 8, 9]

    Returns
    -------
    result : list of int
    """
    # Remove brackets
    range_str = range_str.strip('[]')

    result = []
    
    # Split by comma
    parts = range_str.split(',')

    for part in parts:
        if ':' in part:
            start, end = map(int, part.split(':'))
            result.extend(range(start, end + 1))  # inclusive like MATLAB
        else:
            result.append(int(part))

    return result

def list_ZSTFiles_coriolis(DirZst, VarName, date_files=None):
    """
    List and index Coriolis ZST NetCDF files for a given variable.

    Scans DirZst for files matching ZST_CO_DMQCGL01_ISAS*<VarName>.nc and
    extracts spatial/temporal bounds from each file.

    Parameters
    ----------
    DirZst     : str           - directory containing ZST files
    VarName    : str           - parameter name ('TEMP', 'PSAL', ...)
    date_files : str, optional - date filter string (currently unused)

    Returns
    -------
    ListZst : list of str - paths to successfully opened ZST files
    ZSTfile : dict        - per-file metadata with keys:
        platform_number, wmo_inst_type,
        latitude_min, latitude_max,
        longitude_min, longitude_max,
        juld_min, juld_max, nbfiles
    """
    pattern = os.path.join(DirZst, f"ZST_CO_DMQCGL01_ISAS*{VarName}.nc")
    files = glob.glob(pattern)

    ListZst = []
    ZSTfile = {
        "platform_number": [],
        "wmo_inst_type": [],
        "latitude_min": [],
        "latitude_max": [],
        "longitude_min": [],
        "longitude_max": [],
        "juld_min": [],
        "juld_max": [],
    }

    for fname in files:
        try:
            print(fname)
            ds = xr.open_dataset(fname)

            ref_date = np.datetime64('1950-01-01T00:00:00')

            juld = ds['JULD'].data  # datetime64[ns] array

            # Remove NaT values
            juld_valid = juld[~np.isnat(juld)]

            # Convert to days since 1950-01-01
            juld_days = (juld_valid - ref_date) / np.timedelta64(1, 'D')

            ZSTfile["juld_min"].append(float(np.min(juld_days)))
            ZSTfile["juld_max"].append(float(np.max(juld_days)))
           
            ListZst.append(fname)

            ZSTfile["platform_number"].append(ds['PLATFORM_NUMBER'].data[0])
            ZSTfile["wmo_inst_type"].append(ds['WMO_INST_TYPE'].data[0])

            ZSTfile["latitude_min"].append(float(ds['LATITUDE'].min()))
            ZSTfile["latitude_max"].append(float(ds['LATITUDE'].max()))

            ZSTfile["longitude_min"].append(float(ds['LONGITUDE'].min()))
            ZSTfile["longitude_max"].append(float(ds['LONGITUDE'].max()))
        except Exception as e:
            print(f"Error processing {fname}: {e}")
            continue

    ZSTfile["nbfiles"] = len(ListZst)

    return ListZst, ZSTfile

def UT_cell_list(txt_list):
    """
    Split a comma-separated string into a list of trimmed substrings.

    Equivalent to parsing a MATLAB cell-array literal of strings.

    Parameters
    ----------
    txt_list : str - comma-separated values, e.g. "DM, CO, RT"

    Returns
    -------
    cell_list : list of str
    """
    # Ensure trailing comma
    if not txt_list.endswith(','):
        txt_list = txt_list + ','

    # Find indices of commas
    isep = [i for i, ch in enumerate(txt_list) if ch == ',']
    nb_dir = len(isep)

    cell_list = []
    i1 = 0

    for i in range(nb_dir):
        i2 = isep[i] - 1
        cell_list.append(txt_list[i1:i2+1].strip())
        i1 = i2 + 2

    return cell_list

def flddat_init():
    """
    Return the standard field-name lists for STD/dat NetCDF files.

    Returns
    -------
    fldHDR_pos : list - position header fields: juld, latitude, longitude, deph
    fldHDR_dat : list - data header fields: dc_reference, wmo_inst_type, platform_number
    fldVAR_dat : list - variable fields: proc, qc, data, clmn, clsd, erme, erur, resi
    """

    fldHDR_pos = ['juld', 'latitude','longitude', 'deph']

    fldHDR_dat = ['cycle_number', 'data_centre','data_mode', 'dc_reference', 'pi_name', 'wmo_inst_type', 'platform_number']


    fldHDR_dat = ['dc_reference', 'wmo_inst_type', 'platform_number' ]

    fldVAR_dat = ['proc', 'qc', 'data', 'clmn','clsd','erme','erur','resi']
    return fldHDR_pos, fldHDR_dat, fldVAR_dat


def parse_num_string(val):
    """
    Parse a MATLAB-style numeric string (e.g. "1 2 3", "[1 2 3]", "42.5") into a value.

    Returns a numpy array for multiple values, a float for a scalar, or the
    original string if parsing fails.
    """
    if not val or not val.strip():
        return np.array([])
    
    # Supprime les crochets MATLAB et caractères parasites
    cleaned = val.strip().strip('[]').replace(',', ' ').replace(';', ' ')
    
    try:
        arr = np.fromstring(cleaned, sep=' ', dtype=float)
        if arr.size == 1:
            return float(arr[0])
        return arr
    except ValueError:
        return val  # garde la chaîne si vraiment non parseable

def get_XML_databloc(root, tag_name):
    """
    Find a named XML block and return its children as a {field: value} dict.

    Equivalent to MATLAB's get_XML_databloc.
    """
    block = root.find('.//' + tag_name)
    if block is None:
        return {}
    result = {}
    for child in block:
        result[child.tag] = (child.text or '').strip()
    return result


def ANA_ini(conf_file_name):
    """
    Build and return the INIT configuration dict from an XML config file.

    Parameters
    ----------
    conf_file_name : str - full path to the XML configuration file

    Returns
    -------
    INIT : dict - configuration structure with directories, filenames,
                  spatial/temporal parameters, and derived paths; None if file not found
    """

    if not os.path.isfile(conf_file_name):
        print(f' MSG -ANA_ini :\n      File {conf_file_name} not found')
        return None

    tree = ET.parse(conf_file_name)
    root = tree.getroot()

    INIT = {}

    # ------------------------------------------------------------------
    #  Métadonnées (NCmeta)
    # ------------------------------------------------------------------
    block = get_XML_databloc(root, 'NCmeta')
    for key, val in block.items():
        INIT[key] = val

    # ------------------------------------------------------------------
    #  Répertoires (directories) — s'assure que chaque chemin finit par '/'
    # ------------------------------------------------------------------
    block = get_XML_databloc(root, 'directories')
    for key, val in block.items():
        if not val.endswith('/'):
            val += '/'
        INIT[key] = val

    # ------------------------------------------------------------------
    #  Noms de fichiers (filenames) — s'assure que chaque nom finit par '_'
    # ------------------------------------------------------------------
    block = get_XML_databloc(root, 'filenames')
    for key, val in block.items():
        if not val.endswith('_'):
            val += '_'
        INIT[key] = val

    # ------------------------------------------------------------------
    #  Noms de fichiers 2 (filenames2) — pas de suffixe forcé
    # ------------------------------------------------------------------
    block = get_XML_databloc(root, 'filenames2')
    for key, val in block.items():
        INIT[key] = val

    # ------------------------------------------------------------------
    #  Paramètres spatio-temporels (TimeSpaceAna)
    # ------------------------------------------------------------------
    block = get_XML_databloc(root, 'TimeSpaceAna')
    # area_limits est un vecteur numérique (str2num en MATLAB)
    # TimeSpaceAna
    INIT['area_limits'] = parse_num_string(block.get('area_limits', ''))
    # Mparameters
    block = get_XML_databloc(root, 'Parameters')
    for key, val in block.items():
        INIT[key] = parse_num_string(val)  # garde la chaîne si la conversion échoue

    # ------------------------------------------------------------------
    #  Paramètres vectoriels (Mparameters) — convertis en array (str2num)
    # ------------------------------------------------------------------
    block = get_XML_databloc(root, 'Parameters')
    for key, val in block.items():
        try:
            INIT[key] = np.fromstring(val, sep=' ', dtype=float)
        except Exception:
            INIT[key] = val

    # Mparameters
    block = get_XML_databloc(root, 'Mparameters')
    for key, val in block.items():
        INIT[key] = parse_num_string(val)  # garde la chaîne si la conversion échoue

    # ------------------------------------------------------------------
    #  Paramètres vectoriels (Mparameters) — convertis en array (str2num)
    # ------------------------------------------------------------------
    block = get_XML_databloc(root, 'Mparameters')
    for key, val in block.items():
        try:
            INIT[key] = np.fromstring(val, sep=' ', dtype=float)
        except Exception:
            INIT[key] = val

    # ------------------------------------------------------------------
    #  Paramètres chaînes multiples (MparamStr) — gardés comme strings
    # ------------------------------------------------------------------
    block = get_XML_databloc(root, 'MparamStr')
    for key, val in block.items():
        INIT[key] = val

    # ------------------------------------------------------------------
    #  Chemins dérivés
    # ------------------------------------------------------------------
    INIT['DIR_CLIMREF'] = INIT['DirIsasConfStd']
    INIT['DIR_BATHY']   = INIT['DirIsasConfStd']
    INIT['file_bathy']  = INIT['DIR_BATHY']  + INIT['FilenameBathy']
    INIT['file_climOA'] = INIT['DIR_CLIMREF'] + INIT['FilenameClimTS']
    INIT['file_presOA'] = INIT['file_climOA'] + 'ann_PRES.nc'
    INIT['file_StdOA']  = INIT['DIR_CLIMREF'] + INIT['FilenameStdTS']

    # Affichage de contrôle (équivalent au `INIT` nu en fin de fonction MATLAB)
    print(f'\nINIT loaded from: {conf_file_name}')
    for k, v in INIT.items():
        print(f'  {k}: {v}')

    return INIT

SOFTWARE = "PREOA - 7.0"


def globatt_init(CONFIG, SOFTWARE):
    """
    Build and return the global-attribute dictionary for NetCDF output files.

    Parameters
    ----------
    CONFIG   : dict - configuration dict (from ANA_ini)
    SOFTWARE : str  - software version string

    Returns
    -------
    GLOB_ATT : dict - global attributes including institution, project, spatial bounds, etc.
    """
    GLOB_ATT=dict()
    GLOB_ATT['reference_date_time'] = '1950-01-01T00:00:00Z'
    GLOB_ATT['TITLE']  = ''

    GLOB_ATT['SOFTWARE_VERSION'] = SOFTWARE;

    GLOB_ATT['INSTITUTION']  = CONFIG['INSTITUTION']
    GLOB_ATT['PROJECT_NAME'] = CONFIG['PROJECT_NAME']
    GLOB_ATT['DATA_MANAGER'] = CONFIG['DATA_MANAGER']

    if  'area_limit' in CONFIG.keys(): 
        GLOB_ATT['SOUTH_LAT'] = CONFIG.area_limits(1)
        GLOB_ATT['NORTH_LAT'] = CONFIG.area_limits(2)
        GLOB_ATT['WEST_LONG'] = CONFIG.area_limits(3)
        GLOB_ATT['EAST_LONG'] = CONFIG.area_limits(4)
    else:
        GLOB_ATT['SOUTH_LAT'] = -999
        GLOB_ATT['NORTH_LAT'] = -999
        GLOB_ATT['WEST_LONG'] = -999
        GLOB_ATT['EAST_LONG'] = -999

    GLOB_ATT['LAT_RESO'] = -999
    GLOB_ATT['LON_RESO'] = -999

    if  'REFERENCES' in CONFIG.keys():
        GLOB_ATT['REFERENCES'] = CONFIG['REFERENCES']
    else:
        GLOB_ATT['REFERENCES'] = 'none'
    

    GLOB_ATT['COMMENT'] = '...'
    GLOB_ATT['DATE_START'] = ''
    GLOB_ATT['DATE_STOP'] = ''
    GLOB_ATT['ANALYSIS_NAME'] = ''

    return GLOB_ATT

def convqc(qc):
    """Convert QC values to int8, replacing NaN or masked entries with 0."""
    qc_arr = np.asarray(qc, dtype=float)
    out = np.where(np.isfinite(qc_arr), qc_arr, 0).astype(np.int8)
    return out


def NCW_data_var(ncfile_name, PARAM, STDVAR):
    """
    Écrit les données du paramètre PARAM dans un fichier _dat_ existant
    (créé préalablement par NCW_data_hdr_PREOA).

    Parameters
    ----------
    ncfile_name : str  - chemin complet du fichier NetCDF à compléter
    PARAM       : str  - nom du paramètre (TEMP, PSAL, DOXY)
    STDVAR      : dict - champs : proc, qc, data, clmn, clsd, erme,
                                  erur (optionnel), resi (optionnel)

    Returns
    -------
    msg_error : str
    """

    msg_error = 'Ok: NCW_data_var'

    # ------------------------------------------------------------------
    #  Vérifie les champs obligatoires
    # ------------------------------------------------------------------
    _, _, fldVAR_dat = flddat_init()
    nbfield_min = 6
    fldVAR_PAR_min = fldVAR_dat[:nbfield_min]  # ['proc','qc','data','clmn','clsd','erme']

    present = set(STDVAR.keys())
    missing = [f for f in fldVAR_PAR_min if f not in present]
    if len(present) < nbfield_min or missing:
        return 'Error: NCW_data_var, fields are missing'

    nbfield = len(present)
    has_optional = ('erur' in present) and ('resi' in present)

    # ------------------------------------------------------------------
    #  Paramètres selon la variable
    # ------------------------------------------------------------------
    FILLVALUE = np.float32(99999.0)

    PARAM_DEFS = {
        'TEMP': dict(
            units      = 'degree_Celsius',
            long_name  = 'Temperature (T90) (interpolated on Z_levels)',
            std_name   = 'sea_water_temperature',
            par_valid  = [-3.0,   40.0],
            err_valid  = [ 0.001,  10.0],
            std_valid  = [ 0.001,  10.0],
        ),
        'PSAL': dict(
            units      = 'none',
            long_name  = 'Salinity (S78 - PSS) (interpolated on Z_levels)',
            std_name   = 'sea_water_salinity',
            par_valid  = [ 3.0,   60.0],
            err_valid  = [ 0.001,  10.0],
            std_valid  = [ 0.001,  10.0],
        ),
        'DOXY': dict(
            units      = 'micromole/kg',
            long_name  = 'Disolved oxygen (interpolated on Z_levels)',
            std_name   = 'Disolved oxygen',
            par_valid  = [   0.0, 1400.0],
            err_valid  = [   0.01,  100.0],
            std_valid  = [   0.01, 1400.0],
        ),
    }

    if PARAM not in PARAM_DEFS:
        return f'W: NCW_data_var: Parameter {PARAM} not defined'

    p = PARAM_DEFS[PARAM]
    UNITS     = p['units']
    LONG_NAME = p['long_name']
    STND_NAME = p['std_name']
    PAR_VALID = p['par_valid']
    ERR_VALID = p['err_valid']
    STD_VALID = p['std_valid']

    # Noms des variables NetCDF
    PAR_PROC = f'{PARAM}_PROC'
    PAR_QC   = f'{PARAM}_QC'
    PAR_CLMN = f'{PARAM}_CLMN'
    PAR_CLSD = f'{PARAM}_CLSD'
    PAR_ERME = f'{PARAM}_ERME'
    PAR_ERUR = f'{PARAM}_ERUR'
    PAR_RESI = f'{PARAM}_RESI'

    PAR_QC_flag_values   = np.array([1,2,3,4,5,6,7,8,9], dtype=np.int8)
    PAR_QC_flag_meanings = ('1--4:good_to_acceptable;'
                            '5:bad_quality_interpolation;'
                            '6:bad-manual_flag;'
                            '7--8:not_used;'
                            '9:not_interpolated;')

    COMPRESS = dict(zlib=True, complevel=9, shuffle=True)
    COMPRESS = {}

    # ------------------------------------------------------------------
    #  Ouvre le fichier en écriture (append)
    # ------------------------------------------------------------------
    with Dataset(ncfile_name, 'a') as nc:

            # Récupère les dimensions existantes
            n_prof = len(nc.dimensions['N_PROF'])
            n_lev  = len(nc.dimensions['N_LEVELS'])

            # ----------------------------------------------------------
            #  Définition des variables
            # ----------------------------------------------------------
            # OBLIGATOIRE — ordre exigé par OA_ncreaddata.f90

            v = nc.createVariable(PAR_PROC, 'i1', ('N_PROF',), **COMPRESS)
            v.long_name   = 'profile processing level'
            v.conventions = 'processing status: 1=raw, 2=adjusted'
            v.valid_min   = np.int8(0)
            v.valid_max   = np.int8(5)
            v = nc.createVariable(PAR_QC,   'i1', ('N_LEVELS', 'N_PROF'), fill_value=np.int8(0), **COMPRESS)
            v.long_name      = 'Quality flag on interpolated variable'
            v.flag_values    = PAR_QC_flag_values
            v.flag_meanings  = PAR_QC_flag_meanings

            v = nc.createVariable(PARAM,    'f4',  ('N_LEVELS', 'N_PROF'), fill_value=FILLVALUE,  **COMPRESS)
            v.long_name     = LONG_NAME
            v.units         = UNITS
            v.standard_name = STND_NAME
            v.valid_min     = np.float32(PAR_VALID[0])
            v.valid_max     = np.float32(PAR_VALID[1])

            v = nc.createVariable(PAR_CLMN, 'f4',  ('N_LEVELS', 'N_PROF'), fill_value=FILLVALUE,  **COMPRESS)
            v.long_name = 'Climatology mean for profile'
            v.units     = UNITS
            v.valid_min = np.float32(PAR_VALID[0])
            v.valid_max = np.float32(PAR_VALID[1])

            v = nc.createVariable(PAR_CLSD, 'f4',  ('N_LEVELS', 'N_PROF'), fill_value=FILLVALUE,  **COMPRESS)
            v.long_name = 'Climatology standard deviation for profile'
            v.units     = UNITS
            v.valid_min = np.float32(STD_VALID[0])
            v.valid_max = np.float32(STD_VALID[1])

            v = nc.createVariable(PAR_ERME, 'f4', ('N_LEVELS', 'N_PROF'), fill_value=FILLVALUE,  **COMPRESS)
            v.long_name = 'Measurement error'
            v.units     = UNITS
            v.valid_min = np.float32(ERR_VALID[0])
            v.valid_max = np.float32(ERR_VALID[1])

            if has_optional:
                v = nc.createVariable(PAR_ERUR, 'f4',  ('N_LEVELS', 'N_PROF'), fill_value=FILLVALUE,  **COMPRESS)
                v.long_name = 'Error from unresolved scales'
                v.units     = UNITS
                v.valid_min = np.float32(ERR_VALID[0])
                v.valid_max = np.float32(ERR_VALID[1])

                v = nc.createVariable(PAR_RESI, 'f4',  ('N_LEVELS', 'N_PROF'), fill_value=FILLVALUE,  **COMPRESS)
                v.long_name = 'Residual'
                v.units     = UNITS
                v.valid_min = np.float32(-2 * ERR_VALID[1])
                v.valid_max = np.float32(+2 * ERR_VALID[1])

            # ----------------------------------------------------------
            #  Prépare les tableaux
            # ----------------------------------------------------------
            tab_data = np.asarray(STDVAR['data'], dtype=float)
            tab_erme = np.asarray(STDVAR['erme'], dtype=float)
            tab_clmn = np.asarray(STDVAR['clmn'], dtype=float)
            tab_clsd = np.asarray(STDVAR['clsd'], dtype=float)
            tab_qcc  = convqc(STDVAR['qc'])

            isok_dat = np.isfinite(tab_data)
            isok_cli = np.isfinite(tab_clmn)
            isok_csd = np.isfinite(tab_clsd)
            inok = ~(isok_dat & isok_cli & isok_csd)

            if np.any(inok):
                tab_qcc[inok]  = 9
                tab_data[inok] = FILLVALUE
                tab_erme[inok] = FILLVALUE
                tab_clmn[inok] = FILLVALUE
                tab_clsd[inok] = FILLVALUE
            # ----------------------------------------------------------
            #  Écriture — transpose (N_LEVELS, N_PROF) → (N_PROF, N_LEVELS)
            # ----------------------------------------------------------
            def to_fort(arr, n_prof):
                """Garantit shape (N_LEVELS, N_PROF) pour écriture."""
                a = np.asarray(arr)
                if a.ndim == 2:
                    if a.shape[0] == n_prof and a.shape[1] != n_prof:
                        return a.T   
                return a

            proc = np.asarray(STDVAR['proc']).flatten()
            nc.variables[PAR_PROC][:] = proc.astype(np.int8)
            nc.variables[PAR_QC][:]   = to_fort(tab_qcc,  n_prof).astype(np.int8)
            nc.variables[PARAM][:]    = to_fort(tab_data, n_prof).astype(np.float32)
            nc.variables[PAR_ERME][:] = to_fort(tab_erme, n_prof).astype(np.float32)
            nc.variables[PAR_CLMN][:] = to_fort(tab_clmn, n_prof).astype(np.float32)
            nc.variables[PAR_CLSD][:] = to_fort(tab_clsd, n_prof).astype(np.float32)

            if has_optional:
                tab_erur = np.asarray(STDVAR['erur'], dtype=float)
                tab_resi = np.asarray(STDVAR['resi'], dtype=float)
                if np.any(inok):
                    tab_erur[inok] = FILLVALUE
                    tab_resi[inok] = FILLVALUE
                nc.variables[PAR_ERUR][:] = to_fort(tab_erur, n_prof).astype(np.float32)
                nc.variables[PAR_RESI][:] = to_fort(tab_resi, n_prof).astype(np.float32)

    return msg_error

def datenum(year, month, day):
    """Mimics MATLAB's datenum: days since 0000-01-00 (approx via reference)."""
    from datetime import date
    return (date(year, month, day) - date(1, 1, 1)).days + 367  # MATLAB epoch offset


def datenum_from_str(s, fmt):
    """Parse a date string and return a datenum."""
    from datetime import date, datetime as dt
    d = dt.strptime(s, fmt).date()
    return (d - date(1, 1, 1)).days + 367


def PREOA_main(config_fname, date_ana, VarName, PLOT_DISP, FilenameClimTS=None,
               input_dir=None, output_dir=None):
    """
    Prepares the dataset for each area.

    Parameters
    ----------
    config_fname    : str   - full path to the configuration XML file
    date_ana        : list  - [dd, mm, yyyy]
    VarName         : str   - variable name (e.g. 'TEMP', 'PSAL', 'DOXY')
    PLOT_DISP       : int   - 0 = no plots
    FilenameClimTS  : str, optional
    """

    # -----------------------------------------------------------------------
    #  Constants
    # -----------------------------------------------------------------------
    icheck = 0
    mx_platform = 10000

    # -----------------------------------------------------------------------
    #  Read configuration
    # -----------------------------------------------------------------------
    INIT = ANA_ini(config_fname)
    INIT['PLOT_DISP'] = PLOT_DISP

    if input_dir is not None:
        INIT['DirZstDm'] = input_dir if input_dir.endswith('/') else input_dir + '/'
    if output_dir is not None:
        INIT['DirPreoa'] = output_dir if output_dir.endswith('/') else output_dir + '/'

    if VarName == 'DOXY':
        INIT.nam_clim = INIT.NamClimOxy
        INIT.nam_std  = INIT.NamStdOxy
    else:
        INIT['nam_clim'] = INIT['FilenameClimTS']
        INIT['nam_std']  = INIT['FilenameStdTS']
    INIT['PreoaWmoInstrOut']=INIT['PreoaWmoInstrOut'][1:-1]
    if INIT['PreoaWmoInstrOut']=='':
        INIT['PreoaWmoInstrOut']=False
    # -----------------------------------------------------------------------
    #  Define directories
    # -----------------------------------------------------------------------
    dir_clim     = INIT['DirIsasConfStd']
    dirpreoa_log  = os.path.join(INIT['DirPreoa'], 'PREOA_LOG', '')
    dirpreoa_plot = os.path.join(INIT['DirPreoa'], 'PREOA_PLOT', '')
    dirpreoa_data = os.path.join(INIT['DirPreoa'], 'PREOA_DATA', '')

    Path(dirpreoa_log).mkdir(parents=True, exist_ok=True)
    Path(dirpreoa_plot).mkdir(parents=True, exist_ok=True)

    id_date = datetime.now().strftime('%Y%m%d%H%M%S')
    preoa_err_file = os.path.join(dirpreoa_log, f'PREOA_ERR_{id_date}.txt')
    preoa_err = open(preoa_err_file, 'w')
    preoa_err.write(f'***   std_err: {id_date}\n')

    # Error bounds
    err_max = INIT[f"{VarName[0].upper()}{VarName[1:].lower()}ErrMax"]
    err_min = INIT[f"{VarName[0].upper()}{VarName[1:].lower()}ErrMin"]

    # -----------------------------------------------------------------------
    #  Field/variable descriptors
    # -----------------------------------------------------------------------
    _, _, fldVAR_dat = flddat_init()
    nbfield_min = 6
    nb_fldvar   = nbfield_min

    # -----------------------------------------------------------------------
    #  Date handling
    # -----------------------------------------------------------------------
    dd, mm, yy = date_ana
    jest_abs = datenum(yy, mm, dd)
    from datetime import date as _date
    _d = _date(yy, mm, dd)
    DATE_EST = _d.strftime('%Y%m%d')   # yyyymmdd string

    dir_zst_co = INIT['DirZstRtCorio']
    dir_zst_dm = os.path.join(INIT['DirZstDm'], str(yy), f'{mm:02d}', '')

    dirpreoa_data_tmp = os.path.join(dirpreoa_data, VarName, '')
    Path(dirpreoa_data_tmp).mkdir(parents=True, exist_ok=True)

    OA_FNAM_TMP = f'OA_{DATE_EST}'

    # -----------------------------------------------------------------------
    #  Messages (English)
    # -----------------------------------------------------------------------
    lang = getattr(INIT, 'LANG', 'En')
    if lang == 'En':
        msg_run    = '\n  nb area: %4i, nb_zst: %3i, VarName: %s  \n'
        msg_areast = '\n ***** Starting area %4i ***** \n'
        msg_fileok = ' File: %s \n       %4i profiles \n'
        msg_remove = ' %8i profiles removed on %s criteria \n'
        msg_remdat = '\n %8i data removed on %s criteria \n'
        msg_areaok = ('\n   area %4i, nb_prof: %5i, nb_prof_nored:  %5i, '
                      'processing time: %6.2f sec\n'
                      '  ************************************************\n\n')
        msg_misserr = '  Missing error data \n'
        msg_end    = f'\n  {SOFTWARE} processing done  \n \n'
    else:
        msg_run    = '\n  nb area: %4i, nb_zst: %3i, VarName: %s   \n'
        msg_areast = '\n *****  Commence la zone %4i *****\n'
        msg_fileok = ' Fichier: %s \n       %4i profils \n'
        msg_remove = ' %8i profiles removed on %s criteria \n'
        msg_remdat = '\n %8i data removed on %s criteria \n'
        msg_areaok = ('\n   zone %4i, nb_prof: %5i, nb_prof_nored:  %5i, '
                      'temps de traitement: %6.2f sec\n'
                      '  ************************************************\n\n')
        msg_misserr = '  Missing error data \n'
        msg_end    = f'\n  {SOFTWARE} traitement termine  \n\n'

    # -----------------------------------------------------------------------
    #  Open log file
    # -----------------------------------------------------------------------
    id_date2 = datetime.now().strftime('%Y%m%d%H%M%S')
    preoa_mess_file = os.path.join(dirpreoa_log,
                                   f'preoa_mess_{VarName}_{DATE_EST}_{id_date2}.txt')
    try:
        preoa_mess = open(preoa_mess_file, 'w')
    except OSError:
        print(f'\n  pb opening log file {preoa_mess_file}\n')
        return

    def log(msg, *args):
        """Write to stdout and log file."""
        line = msg % args if args else msg
        print(line, end='')
        preoa_mess.write(line)

    log(f'\n\n >>>>>>>  Running {SOFTWARE}\n')
    log(f'\n Last update: {datetime.now()}\n')
    log(f'\n ZST_CORIO:   {dir_zst_co}\n')
    log(f'\n ZST_DM   :   {dir_zst_dm}\n')
    log(f'\n PREOA: {dirpreoa_data}\n')

    # -----------------------------------------------------------------------
    #  Read climatology coordinates
    # -----------------------------------------------------------------------
    nam_clim = os.path.join(dir_clim, f'ISAS15_DMFD_ann_{VarName}.nc')
    with Dataset(nam_clim, 'r') as nc:
        lon_clim = nc.variables['longitude'][:]
        lat_clim = nc.variables['latitude'][:]
        dep_clim = nc.variables['depth'][:]
        str_ref_date_time = getattr(nc, 'reference_date')

    if isinstance(str_ref_date_time, (list, np.ndarray)):
        str_ref_date_time = ''.join(str_ref_date_time).strip()

    if len(str_ref_date_time) > 16:
        jref = datenum_from_str(str_ref_date_time[:10], '%Y-%m-%d')
    else:
        jref = datenum_from_str(str_ref_date_time[:8], '%Y%m%d')

    jest = jest_abs - jref

    # -----------------------------------------------------------------------
    #  Depth levels to analyse
    # -----------------------------------------------------------------------
    dep_ana = INIT['PreoaDephAna']
    dep_ana=int(dep_ana[1:-1])
    if dep_ana < 0:
        dep_max   = -dep_ana
        nb_depana = int(np.sum(dep_clim <= dep_max))
        dep_ana   = dep_clim[:nb_depana]
        idep_ana  = np.arange(nb_depana)
    else:
        nb_depana = len(dep_ana)
        idep_ana  = np.zeros(nb_depana, dtype=int)
        for k in range(nb_depana):
            iz = int(np.argmin(np.abs(dep_clim - dep_ana[k])))
            dep_ana[k]  = dep_clim[iz]
            idep_ana[k] = iz
        dep_ana = np.array(dep_ana)

    # -----------------------------------------------------------------------
    #  Read bathymetry / area definition
    # -----------------------------------------------------------------------
    file_bathy = os.path.join(INIT['DirIsasConfStd'], f'{INIT['FilenameBathy']}.nc')
    with Dataset(file_bathy, 'r') as nc:
        def_area = nc.variables['basin_area'][:]
        lat_b = nc.variables['latitude'][:]
        lon_b = nc.variables['longitude'][:]
        print(f'def_area shape: {def_area.shape}')
        print(f'lat_b: {lat_b[[0,-1]]}  ({len(lat_b)} pts)')
        print(f'lon_b: {lon_b[[0,-1]]}  ({len(lon_b)} pts)')
        print(f'lat_clim: {lat_clim[[0,-1]]}  ({len(lat_clim)} pts)')
        print(f'lon_clim: {lon_clim[[0,-1]]}  ({len(lon_clim)} pts)')


    nb_ana_area = len(INIT['PreoaAreaList'])

    # -----------------------------------------------------------------------
    #  Build data-type list  (DM first)
    # -----------------------------------------------------------------------
    type_list = UT_cell_list(INIT['PreoaTypList'])
    type_list = sorted(type_list, reverse=True)  # DM first

    # -----------------------------------------------------------------------
    #  Date window for file selection
    # -----------------------------------------------------------------------
    jlim = np.array([jest - INIT['PreoaTimeInt'], jest + INIT['PreoaTimeInt']])

    date_files = f'{yy}{mm:02d}'
    t0 = time.time()
    list_zst, ZSTfile = list_ZSTFiles_coriolis(dir_zst_dm, VarName, date_files)
    print(f'list_ZSTFiles_coriolis elapsed: {time.time()-t0:.2f} s')
    NB_ZST_FILES = len(ZSTfile['juld_max'])

    log(msg_run, nb_ana_area, NB_ZST_FILES, VarName)

    log('\n File types: \n')
    for t in type_list:
        log(f'  {t:5s}\n')

    # -----------------------------------------------------------------------
    #  Loop over areas
    # -----------------------------------------------------------------------
    fname_prefx = f'preoa_{OA_FNAM_TMP}_{VarName}'
    INIT['PreoaAreaList'] = expand_ranges(INIT['PreoaAreaList'])
    nb_ana_area=len(INIT['PreoaAreaList'])
    for ii in range(nb_ana_area):
        t_area_start = time.time()
        iarea   = INIT['PreoaAreaList'][ii]
        str_area = f'{iarea:03d}'

        # Load area mask
        area_def_name = os.path.join(INIT['DirAreaDef'], f'msk_def{str_area}')
        area_data =  loadmat(f"{area_def_name}.mat")
        tab_msk = area_data['tab_msk'].T
        lat_msk = area_data['lat_msk']
        lon_msk = area_data['lon_msk']

        # def_area shape NetCDF Python : (lat, lon)
        # i_ok = indices lat, j_ok = indices lon


        # Après (correct) — utilise lat_msk/lon_msk comme référence
        lat_msk_flat = lat_msk.flatten()
        lon_msk_flat = lon_msk.flatten()

        # Trouve les indices dans lat_clim/lon_clim correspondant au masque
        i1_area = int(np.argmin(np.abs(lat_clim - lat_msk_flat[0])))
        i2_area = int(np.argmin(np.abs(lat_clim - lat_msk_flat[-1])))
        j1_area = int(np.argmin(np.abs(lon_clim - lon_msk_flat[0])))
        j2_area = int(np.argmin(np.abs(lon_clim - lon_msk_flat[-1])))

        lat_area = lat_clim[i1_area:i2_area+1]
        lon_area = lon_clim[j1_area:j2_area+1]
        print(f'Zone {iarea}: lon_area={lon_area[[0,-1]]}, lat_area={lat_area[[0,-1]]}')
        lat_lim = np.array([lat_msk[0], lat_msk[-1]])
        lon_lim = np.array([lon_msk[0], lon_msk[-1]])

        log(msg_areast, iarea)

        # Initialise area accumulators
        OAHDR = {'nbdat': 0}
        OAVAR = {}
        nb_platform = 0
        DM_platform = [None] * mx_platform

        # -------
        #  Loop over ZST files
        # -------
        for i_zst in range(NB_ZST_FILES):
            fnam_zst_all = list_zst[i_zst]

            # Parse TYPE from filename (between last two underscores)
            parts = fnam_zst_all.split('_')
            TYPE  = parts[-2]   # equivalent to MATLAB regex extraction
            typeok = TYPE in type_list

            dateok = not (jlim[0] > ZSTfile['juld_max'][i_zst] or
                          jlim[1] < ZSTfile['juld_min'][i_zst])
            latok  = not (lat_lim[0] > ZSTfile['latitude_max'][i_zst] or
                          lat_lim[1] < ZSTfile['latitude_min'][i_zst])

            ZSTLonLim = np.array([ZSTfile['longitude_min'][i_zst],
                                   ZSTfile['longitude_max'][i_zst]])
            if np.any(lon_lim > 180):
                ineg = ZSTLonLim < 0
                ZSTLonLim[ineg] += 360
            elif np.any(lon_lim < -180):
                ipos = ZSTLonLim > 0
                ZSTLonLim[ipos] -= 360

            lonok = not (lon_lim[0] > ZSTLonLim[1] or
                         lon_lim[1] < ZSTLonLim[0])

            if not (typeok and dateok and latok and lonok):
                continue

            #try:
            ZSTHDR, ZSTVAR = NCR_data(fnam_zst_all, VarName,
                                        0, '0', jlim, tab_msk,
                                        lat_msk, lon_msk)
            if 1:
                if ZSTHDR['nbdat'] <= 0:
                    continue

                zmax   = ZSTHDR['deph'][-1]
                iz_max = int(np.sum(dep_ana <= zmax))
                idep_sel = idep_ana[:iz_max]
                ZSTHDR['deph'] = dep_ana[idep_sel]

                if len(ZSTHDR['deph']) == 0:
                    continue

                for ifi in range(nb_fldvar):
                    nam_fld = fldVAR_dat[ifi]
                    if nam_fld != 'proc':
                        ZSTVAR[nam_fld] = ZSTVAR[nam_fld][idep_sel, :]

                # Adjust longitudes
                if np.max(lon_msk) > 180:
                    ineg = ZSTHDR['longitude'] < 0
                    ZSTHDR['longitude'][ineg] += 360
                elif np.min(lon_msk) < -180:
                    ipos = ZSTHDR['longitude'] > 0
                    ZSTHDR['longitude'][ipos] -= 360

                # Patch missing error
                is_miss = (~np.isfinite(ZSTVAR['erme'])) & np.isfinite(ZSTVAR['data'])
                if np.any(is_miss):
                    ZSTVAR['erme'][is_miss] = 0.01

                ZSTVAR['erme'][ZSTVAR['erme'] < err_min] = err_min

                ZSTVAR['erur'] = np.full((iz_max, ZSTHDR['nbdat']), np.nan)
                ZSTVAR['resi'] = np.full((iz_max, ZSTHDR['nbdat']), np.nan)

                sz = ZSTVAR['data'].shape
                if sz[1] == len(ZSTHDR['juld']) and sz[0] != 0:
                        OAHDR, OAVAR, DM_platform, nb_platform = PREOA_newdata(
                            TYPE, preoa_mess, msg_remove,
                            ZSTHDR, ZSTVAR, OAHDR, OAVAR,
                            DM_platform, nb_platform)

        # End ZST loop
        n_prof_nored = OAHDR['nbdat']

        # -------
        #  Remove excluded instrument types
        # -------
        if OAHDR['nbdat'] > 0 and INIT['PreoaWmoInstrOut']:

            list_ok = np.ones(OAHDR['nbdat'], dtype=bool)
            for i in range(OAHDR['nbdat']):
                num_i = float(OAHDR['wmo_inst_type'][i].strip())
                if num_i in INIT['PreoaWmoInstrOut']:
                    list_ok[i] = False
            nb_removed = OAHDR['nbdat'] - int(list_ok.sum())
            log(msg_remove, nb_removed, 'type')

            if nb_removed > 0:
                OAHDR_new, OAVAR_new = PREOA_select(OAHDR, OAVAR, list_ok)
            else:
                OAHDR_new, OAVAR_new = OAHDR, OAVAR
        elif OAHDR['nbdat'] > 0:
            nb_removed = 0
            OAHDR_new, OAVAR_new = OAHDR, OAVAR
        else:
            OAHDR_new = {'nbdat': 0}

        # -------
        #  Prepare final tables
        # -------
        if OAHDR_new['nbdat'] > 0:
            OAHDR = OAHDR_new
            nb_depana_iarea = len(OAHDR_new['deph'])

            if nb_depana_iarea < nb_depana:
                OAHDR['deph'] = dep_ana
                nb_add  = nb_depana - nb_depana_iarea
                tab_add = np.full((nb_add, OAHDR['nbdat']), np.nan)
                int_add = np.zeros((nb_add, OAHDR['nbdat']), dtype=np.int8)
                for ifi in range(nb_fldvar):
                    nam_fld = fldVAR_dat[ifi]
                    if nam_fld == 'proc':
                        OAVAR['proc'] = OAVAR_new['proc']
                    elif nam_fld == 'qc':
                        OAVAR[nam_fld] = np.vstack([OAVAR_new[nam_fld], int_add])
                    else:
                        OAVAR[nam_fld] = np.vstack([OAVAR_new[nam_fld], tab_add])
            else:
                OAVAR = OAVAR_new

            # Remove data where QC > QC_MAX
            if OAVAR.get('qc') is not None and OAVAR['qc'].size > 0:
                qc     = OAVAR['qc']
                nb_removed = 0
                for j in range(qc.shape[1]):
                    idx = qc[:, j] > INIT['PreoaQcMax']

                    nb_out = int(np.asarray(idx).sum())
                    if nb_out > 0:
                        OAVAR['data'][idx, j] = np.nan
                        nb_removed += nb_out
                log(msg_remdat, nb_removed, 'QC')
        else:
            OAHDR['nbdat'] = 0

        # -------
        #  Write output files
        # -------
        if OAHDR['nbdat'] > 0:
            if iarea == '209':
                99
            if iarea == 209:
                99
            print(f'Zone {iarea}: lon_area={lon_area[[0,-1]]}, lat_area={lat_area[[0,-1]]}')
            print(f'  profils lon={OAHDR["longitude"].min():.2f}..{OAHDR["longitude"].max():.2f}')
            print(f'  profils lat={OAHDR["latitude"].min():.2f}..{OAHDR["latitude"].max():.2f}')
            # --- dat file ---
            file_name   = f'{OA_FNAM_TMP}_{str_area}_dat_{VarName}.nc'
            if os.path.exists(f"{dirpreoa_data_tmp}/{date_ana[2]}/")==False : 
                os.system(f"mkdir {dirpreoa_data_tmp}/{date_ana[2]}/")
            if os.path.exists(f"{dirpreoa_data_tmp}/{date_ana[2]}/{'%02d'%date_ana[1]}")==False : 
                os.system(f"mkdir {dirpreoa_data_tmp}/{date_ana[2]}/{'%02d'%date_ana[1]}")  
            fname_OA_dat = os.path.join(f"{dirpreoa_data_tmp}/{date_ana[2]}/{'%02d'%date_ana[1]}", file_name)

            GLOB_ATT = globatt_init(INIT, SOFTWARE)
            GLOB_ATT['reference_date_time'] = str_ref_date_time
            GLOB_ATT['rANALYSIS_NAME'] = 'area'
            OAHDR['reference_date_time'] = str_ref_date_time
            OAHDR['data_type'] = 'ISAS-areaDataSet'
            '''
            msg_error1 = NCW_data_hdr_PREOA(fname_OA_dat, GLOB_ATT, OAHDR)
            log(f'\n{msg_error1}\n')

            if msg_error1[0].upper() in ('O', 'W'):
                OAVAR['erur'] = np.full(OAVAR['data'].shape, np.nan)
                OAVAR['resi'] = np.full(OAVAR['data'].shape, np.nan)
                NCW_data_var(fname_OA_dat, VarName, OAVAR)
            '''
            msg_error1 = NCW_data_hdr_and_var(fname_OA_dat, GLOB_ATT, OAHDR, VarName, OAVAR)

            # --- fld file ---

            fname_OA_fld = os.path.join(f"{dirpreoa_data_tmp}/{date_ana[2]}/{'%02d'%date_ana[1]}",
                f'{OA_FNAM_TMP}_{str_area}_fld_{VarName}.nc')
            PARANO = f'{VarName}_ANO'

            msg_error = NCW_OA_field(fname_OA_fld, PARANO, GLOB_ATT,
                                      lon_area, lat_area, dep_ana,
                                      DATE_EST, None, None, None)
            log(f'{msg_error}\n')

        elapsed = time.time() - t_area_start
        log(msg_areaok, iarea, OAHDR['nbdat'], n_prof_nored, elapsed)

    # End area loop
    log(msg_end)
    preoa_mess.close()
    preoa_err.close()


def PREOA_launcher(year, month, input_dir=None, output_dir=None,
                   config_fname=None, varname='PSAL'):
    """
    Entry point: run PREOA_main for the 15th of the given year/month.

    Parameters
    ----------
    year       : int - analysis year (e.g. 2010)
    month      : int - analysis month (1-12)
    input_dir  : str, optional - CORA input directory, overrides DirZstDm in XML
    output_dir : str, optional - zones output directory, overrides DirPreoa in XML
    config_fname : str, optional - XML config file (default: ./configfile/conf_isasana_MY_OA.xml)
    varname    : str - variable name (default: 'PSAL')
    """
    if config_fname is None:
        config_fname = './configfile/conf_isasana_MY_OA.xml'

    PREOA_main(config_fname, [15, month, year], varname, 0,
               input_dir=input_dir, output_dir=output_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='PREOA launcher — prepare observation areas for ISAS OA')
    parser.add_argument('--input',  required=True,
                        help='Path to the CORA input directory (overrides DirZstDm in XML)')
    parser.add_argument('--output', required=True,
                        help='Path to the zones output directory (overrides DirPreoa in XML)')
    parser.add_argument('--year',   required=True, type=int,
                        help='Analysis year (e.g. 2010)')
    parser.add_argument('--month',  required=True, type=int,
                        help='Analysis month (1-12, or zero-padded e.g. 01)')
    parser.add_argument('--config',
                        default='./configfile/conf_isasana_MY_OA.xml',
                        help='XML configuration file (default: ./configfile/conf_isasana_MY_OA.xml)')
    parser.add_argument('--var', default='PSAL',
                        help='Variable name to process (default: PSAL)')
    args = parser.parse_args()

    PREOA_launcher(
        year=args.year,
        month=args.month,
        input_dir=args.input,
        output_dir=args.output,
        config_fname=args.config,
        varname=args.var,
    )




