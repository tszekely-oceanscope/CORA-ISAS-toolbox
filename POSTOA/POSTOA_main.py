import os
import time
import logging
import configparser
import numpy as np
from pathlib import Path
from datetime import datetime, date
from netCDF4 import Dataset
import sys
from xml.etree import ElementTree as ET

VARDESC = {
    'TEMP': dict(
        long_name='Temperature', standard_name='sea_water_temperature',
        units='degree_Celsius', nc_type='NC_FLOAT',
        valid_min=-3.0, valid_max=40.0),
    'PSAL': dict(
        long_name='Salinity', standard_name='sea_water_salinity',
        units='PSS-78', nc_type='NC_FLOAT',
        valid_min=0, valid_max=40.0),
    'DOXY': dict(
        long_name='Dissolved oxygen', standard_name='dissolved_oxygen',
        units='micromole/kg', nc_type='NC_FLOAT',
        valid_min=0.0, valid_max=1400.0),
    'TEMP_ANO': dict(
        long_name='Temperature anomaly',
        standard_name='sea_water_temperature_anomaly',
        units='degree_Celsius', nc_type='NC_FLOAT',
        valid_min=-10.0, valid_max=10.0),
    'PSAL_ANO': dict(
        long_name='Salinity anomaly',
        standard_name='sea_water_salinity_anomaly',
        units='PSS-78', nc_type='NC_FLOAT',
        valid_min=-5.0, valid_max=5.0),
    'TEMP_ERR': dict(
        long_name='Temperature error',
        standard_name='sea_water_temperature_error',
        units='degree_Celsius', nc_type='NC_FLOAT',
        valid_min=0.0, valid_max=10.0),
    'PSAL_ERR': dict(
        long_name='Salinity error',
        standard_name='sea_water_salinity_error',
        units='PSS-78', nc_type='NC_FLOAT',
        valid_min=0.0, valid_max=5.0),
    'TEMP_ANO_ERR': dict(
        long_name='Temperature anomaly error',
        standard_name='',
        units='degree_Celsius', nc_type='NC_FLOAT',
        valid_min=0.0, valid_max=10.0),
    'PSAL_ANO_ERR': dict(
        long_name='Salinity anomaly error',
        standard_name='',
        units='PSS-78', nc_type='NC_FLOAT',
        valid_min=0.0, valid_max=5.0),
    'TEMP_PCTVAR': dict(
        long_name='Temperature percent variance',
        standard_name='', units='%', nc_type='NC_FLOAT',
        valid_min=0.0, valid_max=100.0),
    'PSAL_PCTVAR': dict(
        long_name='Salinity percent variance',
        standard_name='', units='%', nc_type='NC_FLOAT',
        valid_min=0.0, valid_max=100.0),
    'TEMP_ANO_PCTVAR': dict(
        long_name='Temperature anomaly percent variance',
        standard_name='', units='%', nc_type='NC_FLOAT',
        valid_min=0.0, valid_max=100.0),
    'PSAL_ANO_PCTVAR': dict(
        long_name='Salinity anomaly percent variance',
        standard_name='', units='%', nc_type='NC_FLOAT',
        valid_min=0.0, valid_max=100.0),
}

NC_TYPE_MAP = {
    'NC_FLOAT':  ('f4', np.float32,  np.float32(9999999.0)),
    'NC_DOUBLE': ('f8', np.float64,  np.float64(9999999.0)),
    'NC_INT':    ('i4', np.int32,    np.int32(2**31 - 1)),
    'NC_SHORT':  ('i2', np.int16,    np.int16(2**15 - 1)),
    'NC_BYTE':   ('i1', np.int8,     np.int8(2**7 - 1)),
}

COMPRESS = dict(zlib=True, complevel=9, shuffle=True)

def convqc(qc):
    """Convertit les QC en entiers int8, remplace NaN/masqué par 0."""
    qc_arr = np.asarray(qc, dtype=float)
    out = np.where(np.isfinite(qc_arr), qc_arr, 0).astype(np.int8)
    return out
def NCW_data_hdr_and_var(ncfile_name, GLOB_ATT, STDHDR, PARAM, STDVAR):
    """
    Crée le fichier et écrit header + variables en une seule ouverture.
    N_LEVELS et N_PROF sont les seules dimensions.
    Les variables 2D sont déclarées en (N_PROF, N_LEVELS) pour compenser
    l'inversion des dimids par ifx/netcdf-fortran.
    Les données sont transposées en conséquence.
    """
    msg_err = 'ERROR NCW_data_hdr_and_var'
    msg_ok  = 'OK: NCW_data_hdr_and_var'
    STDHDR['reference_date_time']='19500101T000000Z'
    ref_date_time_iso, msg_iso = iso_date_time(
        STDHDR.get('reference_date_time', ''))
    if msg_iso.startswith('er'):
        return f'{msg_err}, \n   {msg_iso}'

    for ch in ('-', ':'):
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
    #  Prépare les tableaux
    # ------------------------------------------------------------------
    tab_data = np.asarray(STDVAR['data'], dtype=float)
    tab_erme = np.asarray(STDVAR['erme'], dtype=float)
    tab_clmn = np.asarray(STDVAR['clmn'], dtype=float)
    tab_clsd = np.asarray(STDVAR['clsd'], dtype=float)
    tab_qcc  = convqc(STDVAR['qc'])

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
        tab_erur = np.asarray(STDVAR['erur'], dtype=float)
        tab_resi = np.asarray(STDVAR['resi'], dtype=float)
        if np.any(inok):
            tab_erur[inok] = FILLVALUE
            tab_resi[inok] = FILLVALUE

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

            if has_optional:
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

def NCW_4D_variable(ncfile_name, VARIABLE_XLS, VARIABLE_SAV=None):
    """
    Ajoute une variable 4D (longitude, latitude, depth, time) dans un
    fichier NetCDF4 existant créé par NCW_4D_create.

    Parameters
    ----------
    ncfile_name   : str - chemin du fichier NetCDF existant
    VARIABLE_XLS  : str - nom de variable dans VARDESC (ex: 'TEMP_ANO')
    VARIABLE_SAV  : str optionnel - nom à sauvegarder si différent

    Returns
    -------
    msg_error : str
    """

    fct_name = 'NCW_4D_variable'
    msg_err  = f'ERROR:{fct_name}'
    msg_ok   = f'OK:{fct_name}'

    vd = VARDESC.get(VARIABLE_XLS)
    if vd is None:
        return (f'{msg_err}\n  on file {ncfile_name}\n'
                f'  Variable {VARIABLE_XLS} not in VARDESC\n')

    # Nom et long_name
    if VARIABLE_SAV is not None:
        PARAM = VARIABLE_SAV
        special = {
            'ML_TEMP': 'Mixed layer temperature',
            'ML_PSAL': 'Mixed layer salinity',
        }
        LONG_NAME = special.get(VARIABLE_SAV, vd['long_name'])
    else:
        PARAM     = VARIABLE_XLS
        LONG_NAME = vd['long_name']

    nc_str, np_t, fillvalue = NC_TYPE_MAP.get(
        vd['nc_type'], ('f4', np.float32, np.float32(9999999.0)))

    try:
        with Dataset(ncfile_name, 'a') as nc:
            v = nc.createVariable(
                PARAM, nc_str,
                ('time', 'depth', 'latitude', 'longitude'),
                fill_value=fillvalue,
                **COMPRESS)

            v.long_name = LONG_NAME
            if vd.get('standard_name'):
                v.standard_name = vd['standard_name']
            if vd.get('units'):
                v.units = vd['units']
            v.valid_min = np_t(vd['valid_min'])
            v.valid_max = np_t(vd['valid_max'])

    except Exception as e:
        return f'{msg_err} on {ncfile_name} — {e}'

    return msg_ok
def NCW_4D_create(ncfile_name, GLOB_ATT, longitude, latitude, depth, time=None):
    """
    Create a NetCDF4 file with 4D grid structure (longitude, latitude, depth, time).

    Removes any existing file at `ncfile_name`, creates parent directories as
    needed, writes coordinate variables and global attributes.  Data variables
    must be added afterwards with NCW_4D_variable.

    Parameters
    ----------
    ncfile_name : str             - path of the NetCDF file to create
    GLOB_ATT    : namespace/dict  - global attributes (TITLE, INSTITUTION, …)
    longitude   : array-like      - 1-D longitude values (degree_east)
    latitude    : array-like      - 1-D latitude values  (degree_north)
    depth       : array-like      - 1-D depth values     (m, positive down)
    time        : array-like or None - 1-D time values in days since 1950-01-01;
                  if None, the time dimension is unlimited

    Returns
    -------
    'OK:NCW_4D_create'
    """

    from pathlib import Path
    from datetime import datetime, timezone

    if Path(ncfile_name).exists():
        print(f'Removing: {ncfile_name}')
        Path(ncfile_name).unlink()

    Path(ncfile_name).parent.mkdir(parents=True, exist_ok=True)

    lon = np.asarray(longitude).flatten().astype(np.float32)
    lat = np.asarray(latitude).flatten().astype(np.float32)
    dep = np.asarray(depth).flatten().astype(np.float32)

    COMPRESS = dict(zlib=True, complevel=9, shuffle=True)

    # Sans try/except pour voir l'erreur exacte
    with Dataset(ncfile_name, 'w', format='NETCDF4') as nc:

        # Dimensions en premier — AVANT toute createVariable
        nc.createDimension('longitude', len(lon))
        nc.createDimension('latitude',  len(lat))
        nc.createDimension('depth',     len(dep))
        nc.createDimension('time',
                           None if time is None
                           else len(np.asarray(time).flatten()))

        print(f'Dimensions created: {list(nc.dimensions.keys())}')

        v = nc.createVariable('longitude', 'f4', ('longitude',), **COMPRESS)
        v.standard_name = 'longitude'
        v.units         = 'degree_east'
        v.valid_min     = np.float32(-180.)
        v.valid_max     = np.float32(180.)
        v.axis          = 'X'
        v[:] = lon

        v = nc.createVariable('latitude', 'f4', ('latitude',), **COMPRESS)
        v.standard_name = 'latitude'
        v.units         = 'degree_north'
        v.valid_min     = np.float32(-90.)
        v.valid_max     = np.float32(90.)
        v.axis          = 'Y'
        v[:] = lat

        v = nc.createVariable('depth', 'f4', ('depth',), **COMPRESS)
        v.standard_name = 'depth'
        v.units         = 'm'
        v.positive      = 'down'
        v.valid_min     = np.float32(0.)
        v.valid_max     = np.float32(12000.)
        v.axis          = 'Z'
        v[:] = dep

        v = nc.createVariable('time', 'f4', ('time',), **COMPRESS)
        v.standard_name = 'time'
        v.units         = 'days since 1950-01-01T00:00:00Z'
        v.axis          = 'T'
        v.calendar      = 'standard'
        if time is not None:
            v[:] = np.asarray(time).flatten().astype(np.float32)

        def gatt(name, default=''):
            if hasattr(GLOB_ATT, name):
                return getattr(GLOB_ATT, name)
            if isinstance(GLOB_ATT, dict):
                return GLOB_ATT.get(name, default)
            return default

        str_now = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')
        nc.Conventions           = gatt('CONVENTIONS', 'CF-1.4')
        nc.title                 = gatt('TITLE')
        nc.history               = f'{str_now} : Creation'
        nc.institution           = gatt('INSTITUTION')
        nc.project_name          = gatt('PROJECT_NAME')
        nc.source                = gatt('SOURCE', '')
        nc.analysis_name         = gatt('ANALYSIS_NAME')
        nc.data_manager          = gatt('DATA_MANAGER')
        nc.software_version      = gatt('SOFTWARE_VERSION')
        nc.southernmost_latitude = gatt('SOUTH_LAT')
        nc.northernmost_latitude = gatt('NORTH_LAT')
        nc.latitude_resolution   = gatt('LAT_RESO')
        nc.westernmost_longitude = gatt('WEST_LONG')
        nc.easternmost_longitude = gatt('EAST_LONG')
        nc.longitude_resolution  = gatt('LON_RESO')
        nc.start_date            = gatt('DATE_START')
        nc.stop_date             = gatt('DATE_STOP')
        nc.creation_date         = str_now
        nc.reference_date        = gatt('reference_date_time')
        nc.references            = gatt('REFERENCES')
        nc.comment               = gatt('COMMENT', '')

    return 'OK:NCW_4D_create'



def data_hdr_format(val, str_dim, n_prof):
    """Retourne array (str_dim, n_prof) en bytes 'S1'."""
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

SOFTWARE = 'POSTOA_main - 7.0'
def POSTOA_filter(CONFIG, idx_lim, FIELD_TOT, PCT_VAR_TOT, filt_errmax):
    """
    Filtre le champ et l'erreur aux limites de zones où l'erreur
    dépasse le seuil filt_errmax.

    Parameters
    ----------
    CONFIG      : namespace/dict - configuration
    idx_lim     : list - [i1clim, i2clim, j1clim, j2clim]
    FIELD_TOT   : array (NX, NY, NZ)
    PCT_VAR_TOT : array (NX, NY, NZ)
    filt_errmax : float - seuil d'erreur

    Returns
    -------
    FIELD_TOT, PCT_VAR_TOT : arrays filtrés
    """

    smth_di = 2
    smth_dj = 2

    # ------------------------------------------------------------------
    #  Lecture bathymétrie
    # ------------------------------------------------------------------
    if hasattr(CONFIG, 'DIR_BATHY'):
        dir_bathy = CONFIG.DIR_BATHY
        fname_bathy = CONFIG.FilenameBathy
    else:
        dir_bathy = CONFIG['DIR_BATHY']
        fname_bathy = CONFIG['FilenameBathy']

    file_bathy = f'{dir_bathy}{fname_bathy}.nc'
    with Dataset(file_bathy, 'r') as nc:
        def_area = nc.variables['basin_area'][:]

    # NetCDF Python lit en (lat, lon), MATLAB lit en (lon, lat) — transposez
    def_area = np.array(def_area).T  # (lat, lon) → (lon, lat)

    i1, i2, j1, j2 = idx_lim
    def_area0 = def_area[j1:j2+1, i1:i2+1].astype(float)
    def_area0[def_area0 < 0] = 0

    NX_GLO, NY_GLO, NZ_GLO = FIELD_TOT.shape

    # ------------------------------------------------------------------
    #  Filtre discontinuités zonales (axe X)
    # ------------------------------------------------------------------
    dx = np.diff(def_area0, axis=0)
    dx = np.vstack([dx, np.zeros((1, NY_GLO))])

    tab_discont = np.zeros((NX_GLO, NY_GLO))
    tab_discont[np.abs(dx) > 0.5] = 1

    # Supprime les transitions terre-mer
    for j in range(NY_GLO):
        isout = np.where(tab_discont[:, j] == 1)[0]
        if len(isout) > 0:
            if isout[-1] + 1 >= NX_GLO:
                isout = isout[:-1]
            for i in isout:
                if def_area0[i, j] * def_area0[i+1, j] <= 0:
                    tab_discont[i, j] = 0

    for iz in range(NZ_GLO):
        fld_filt = FIELD_TOT[:, :, iz].copy()
        err_filt = PCT_VAR_TOT[:, :, iz].copy()
        tab_filt = tab_discont.copy()

        tab_filt[err_filt > filt_errmax] += 1
        tab_filt[err_filt < 100]         += 1

        for j in range(NY_GLO):
            ilist = np.where(tab_filt[:, j] == 3)[0]
            for i in ilist:
                for k in range(-smth_di, smth_di + 1):
                    ii = i + k
                    if 0 <= ii < NX_GLO:
                        i1s = max(0, ii - smth_di)
                        i2s = min(NX_GLO, ii + smth_di + 1)

                        tab_mean_err = PCT_VAR_TOT[i1s:i2s, j, iz]
                        tab_mean_fld = FIELD_TOT[i1s:i2s, j, iz]
                        isok = tab_mean_err < 100

                        if np.any(isok):
                            err_filt[ii, j] = np.mean(tab_mean_err[isok])
                            fld_filt[ii, j] = np.mean(tab_mean_fld[isok])

        FIELD_TOT[:, :, iz]   = fld_filt
        PCT_VAR_TOT[:, :, iz] = err_filt

    # ------------------------------------------------------------------
    #  Filtre discontinuités méridiennes (axe Y)
    # ------------------------------------------------------------------
    dy = np.diff(def_area0, axis=1)
    dy = np.hstack([dy, np.zeros((NX_GLO, 1))])

    tab_discont = np.zeros((NX_GLO, NY_GLO))
    tab_discont[np.abs(dy) > 0.5] = 1

    # Supprime les transitions terre-mer
    for i in range(NX_GLO):
        isout = np.where(tab_discont[i, :] == 1)[0]
        if len(isout) > 0:
            if isout[-1] + 1 >= NY_GLO:
                isout = isout[:-1]
            for j in isout:
                if def_area0[i, j] * def_area0[i, j+1] <= 0:
                    tab_discont[i, j] = 0

    for iz in range(NZ_GLO):
        fld_filt = FIELD_TOT[:, :, iz].copy()
        err_filt = PCT_VAR_TOT[:, :, iz].copy()
        tab_filt = tab_discont.copy()

        tab_filt[err_filt > filt_errmax] += 1
        tab_filt[err_filt < 100]         += 1

        for i in range(NX_GLO):
            jlist = np.where(tab_filt[i, :] == 3)[0]
            for j in jlist:
                for k in range(-smth_dj, smth_dj + 1):
                    jj = j + k
                    if 0 <= jj < NY_GLO:
                        j1s = max(0, jj - smth_dj)
                        j2s = min(NY_GLO, jj + smth_dj + 1)

                        tab_mean_err = PCT_VAR_TOT[i, j1s:j2s, iz]
                        tab_mean_fld = FIELD_TOT[i, j1s:j2s, iz]
                        isok = tab_mean_err < 100

                        if np.any(isok):
                            err_filt[i, jj] = np.mean(tab_mean_err[isok])
                            fld_filt[i, jj] = np.mean(tab_mean_fld[isok])

        FIELD_TOT[:, :, iz]   = fld_filt
        PCT_VAR_TOT[:, :, iz] = err_filt

    return FIELD_TOT, PCT_VAR_TOT
def iso_date_time(s):
    """
    Retourne la date au format ARGO/ISO : 'yyyymmddTHHMMSSZ' (16 chars)
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
def _iso_date_time(s):
    """Retourne une date ISO 8601 'yyyy-mm-ddTHH:MM:SSZ' depuis diverses entrées."""
    s = str(s).strip().replace('-', '').replace(':', '').replace('T', '') \
               .replace('Z', '').replace(' ', '')[:14].ljust(14, '0')
    try:
        datetime.strptime(s, '%Y%m%d%H%M%S')
        return f'{s[0:4]}-{s[4:6]}-{s[6:8]}T{s[8:10]}:{s[10:12]}:{s[12:14]}Z', 'ok'
    except ValueError as e:
        return '', f'error: {e}'
def NCW_data_hdr_PREOA(ncfile_name, GLOB_ATT, STDHDR):
    """
    Crée et écrit le header (métadonnées) pour les fichiers de données
    (types STD et dat) avec compression maximale.

    Parameters
    ----------
    ncfile_name : str          - chemin complet du fichier à créer
    GLOB_ATT    : namespace/dict - attributs globaux
    STDHDR      : dict           - header et métadonnées

    Returns
    -------
    msg_error : str - 'OK: ...' ou 'ERROR: ...'
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
        with Dataset(ncfile_name, 'w', format='NETCDF3_64BIT') as nc:

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
                STDHDR['platform_number'], 8, n_prof).T

            nc.variables['WMO_INST_TYPE'][:] = data_hdr_format(
                STDHDR['wmo_inst_type'], 4, n_prof).T

            nc.variables['DC_REFERENCE'][:] = data_hdr_format(
                STDHDR['dc_reference'], 32, n_prof).T

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


def NCR_OA_field(file_nc_OA, PARAM, list_level=-1, ilim=None, jlim=None):
    """
    Lit le champ paramètre et les coordonnées depuis un fichier NetCDF OA.

    Parameters
    ----------
    file_nc_OA  : str   - chemin du fichier NetCDF
    PARAM       : str   - nom du paramètre (PSAL, TEMP, ...)
    list_level  : int ou list
                  -1 : tous les niveaux (défaut)
                   0 : coordonnées seulement
                  -n : niveaux 0 à n-1
                  list : niveaux spécifiés
    ilim        : [i1, i2] - limites latitude (indices, inclusifs)
    jlim        : [j1, j2] - limites longitude (indices, inclusifs)

    Returns
    -------
    lon_ana, lat_ana, dep_ana, param_OA, pct_var
    """

    with Dataset(file_nc_OA, 'r') as nc:
        var_list = list(nc.variables.keys())

        # Détecte version V6
        nam_dep = 'depth' if 'depth' in var_list else 'DEPTH'

        # Coordonnée verticale
        if PARAM.startswith('Z'):
            ZCOO   = PARAM[1:]
            PARSAV = nam_dep
        else:
            ZCOO   = nam_dep
            PARSAV = PARAM

        dep_ana = nc.variables[ZCOO][:]
        lat_ana = nc.variables['latitude'][:]
        lon_ana = nc.variables['longitude'][:]

        nb_dep = len(dep_ana)

        # Construit list_level
        if np.isscalar(list_level):
            if list_level == 0:
                return lon_ana, lat_ana, dep_ana, None, None
            elif list_level == -1:
                list_level = list(range(nb_dep))
            elif list_level < -1:
                list_level = list(range(0, -int(list_level)))
            else:
                list_level = [int(list_level) - 1]  # 1-based → 0-based
        else:
            list_level = [int(l) - 1 for l in list_level]  # 1-based → 0-based

        # Limites lat/lon
        if ilim is None:
            ilim = [0, len(lat_ana) - 1]
        else:
            ilim = [int(ilim[0]), int(ilim[1])]

        if jlim is None:
            jlim = [0, len(lon_ana) - 1]
        else:
            jlim = [int(jlim[0]), int(jlim[1])]

        lat_ana = lat_ana[ilim[0]:ilim[1]+1]
        lon_ana = lon_ana[jlim[0]:jlim[1]+1]
        dep_ana = dep_ana[list_level]

        nb_lon   = jlim[1] - jlim[0] + 1
        nb_lat   = ilim[1] - ilim[0] + 1
        nb_level = len(list_level)

        # Scaling
        var_nc   = nc.variables[PARSAV]
        fill_val = getattr(var_nc, '_FillValue', None)

        if hasattr(var_nc, 'scale_factor'):
            scaling      = True
            scale_factor = float(var_nc.scale_factor)
            add_offset   = float(getattr(var_nc, 'add_offset', 0))
        else:
            scaling      = False
            scale_factor = 1.0
            add_offset   = 0.0

        # Lecture niveau par niveau
        param_OA = np.full((nb_lon, nb_lat, nb_level), np.nan)

        for ilev, iz in enumerate(list_level):
            # Fichier stocké en (time, depth, lat, lon) ou (lon, lat, depth, time)
            raw = var_nc[:]
            if raw.ndim == 4:
                # (time, depth, lat, lon) → extrait (lat, lon) pour ce niveau
                tab = np.array(raw[0, iz,
                                   ilim[0]:ilim[1]+1,
                                   jlim[0]:jlim[1]+1], dtype=float)
                # (lat, lon) → (lon, lat)
                tab = tab.T
            elif raw.ndim == 3:
                tab = np.array(raw[iz,
                                   ilim[0]:ilim[1]+1,
                                   jlim[0]:jlim[1]+1], dtype=float).T
            else:
                tab = np.array(raw[ilev], dtype=float)

            if fill_val is not None:
                tab[tab == fill_val] = np.nan

            tab = tab * scale_factor + add_offset
            param_OA[:, :, ilev] = tab

        # Lecture pct_var si disponible
        pct_var = None
        for vname in var_list:
            if vname in ('PCTVAR', 'pct_variance'):
                pv_nc    = nc.variables[vname]
                fv_pct   = getattr(pv_nc, '_FillValue', None)
                pct_var  = np.full((nb_lon, nb_lat, nb_level), np.nan)
                raw_pv   = pv_nc[:]

                for ilev, iz in enumerate(list_level):
                    if raw_pv.ndim == 4:
                        tab = np.array(raw_pv[0, iz,
                                              ilim[0]:ilim[1]+1,
                                              jlim[0]:jlim[1]+1],
                                       dtype=float).T
                    elif raw_pv.ndim == 3:
                        tab = np.array(raw_pv[iz,
                                              ilim[0]:ilim[1]+1,
                                              jlim[0]:jlim[1]+1],
                                       dtype=float).T
                    else:
                        tab = np.array(raw_pv[ilev], dtype=float)

                    if fv_pct is not None:
                        tab[tab == fv_pct] = np.nan
                    pct_var[:, :, ilev] = tab
                break

    return lon_ana, lat_ana, dep_ana, param_OA, pct_var

    
def globatt_init(CONFIG, SOFTWARE):
    """
    Initialise and return the global-attributes dictionary for NetCDF output.

    Parameters
    ----------
    CONFIG   : dict - configuration mapping (must contain INSTITUTION,
               PROJECT_NAME, DATA_MANAGER; optionally REFERENCES, area_limits)
    SOFTWARE : str  - software version string written to GLOB_ATT['SOFTWARE_VERSION']

    Returns
    -------
    GLOB_ATT : dict - global attributes ready to pass to NCW_4D_create or
               NCW_data_hdr_and_var
    """
    GLOB_ATT=dict()
    GLOB_ATT['reference_date_time'] = '1950-01-01T00:00:00Z'
    GLOB_ATT['TITLE']  = ''

    GLOB_ATT['SOFTWARE_VERSION'] = SOFTWARE

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

def normalize_tab(tab, n_prof):
    """Garantit shape (nlevels, nprofs)."""
    tab = np.ma.atleast_2d(tab)
    if tab.shape[0] == n_prof and tab.shape[1] != n_prof:
        tab = tab.T
    elif tab.shape[0] != n_prof and tab.shape[1] == n_prof:
        pass  # déjà correct
    return tab

def flddat_init_copernicus_PREOA():
    """
    Initialise les listes de champs pour les fichiers STD et _dat_.

    Returns
    -------
    fldHDR_pos : list - champs de position du header
    fldHDR_dat : list - champs de données du header
    fldVAR_dat : list - champs de variables
    """

    fldHDR_pos = ['juld', 'latitude', 'longitude', 'deph']

    fldHDR_dat = ['dc_reference', 'wmo_inst_type', 'platform_number']

    fldVAR_dat = ['proc', 'qc', 'data', 'clmn', 'clsd', 'erme', 'erur', 'resi']

    return fldHDR_pos, fldHDR_dat, fldVAR_dat

def NCR_data(fnam_data, PARAM, iopt_nan=0, pltnum='0',
             jlim=None, tab_msk=None, lat_msk=None, lon_msk=None):
    """
    Lit les données depuis un fichier STD ou dat et retourne les profils
    sélectionnés sur critères plateforme, temps, position.

    Parameters
    ----------
    fnam_data : str   - chemin complet du fichier NetCDF
    PARAM     : str   - nom du paramètre (TEMP, PSAL, DOXY, ...)
    iopt_nan  : int   - si 1, remplace fill_value par NaN
    pltnum    : str   - numéro de plateforme ('0' = pas de sélection)
    jlim      : array - [jmin, jmax] limites temporelles en jours juliens
    tab_msk   : array - masque géographique 2D
    lat_msk   : array - latitudes du masque
    lon_msk   : array - longitudes du masque

    Returns
    -------
    STDHDR : dict - métadonnées et attributs globaux
    STDVAR : dict - données
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
def parse_matlab_range_list(s):
    """
    Parse une chaîne MATLAB comme '[101:161,201:238,301:382]'
    en liste d'entiers Python.
    """
    s = s.strip().strip('[]')
    result = []
    for part in s.split(','):
        part = part.strip()
        if ':' in part:
            bounds = part.split(':')
            start, end = int(bounds[0]), int(bounds[1])
            result.extend(range(start, end + 1))
        else:
            result.append(int(part))
    return result

def flddat_init():
    """
    Return the field-name lists used to read/write STD and dat NetCDF files.

    Returns
    -------
    fldHDR_pos : list - positional header fields  ['juld', 'latitude', 'longitude', 'deph']
    fldHDR_dat : list - string header fields      ['dc_reference', 'wmo_inst_type', 'platform_number']
    fldVAR_dat : list - data variable fields      ['proc', 'qc', 'data', 'clmn', 'clsd', 'erme', 'erur', 'resi']
    """

    fldHDR_pos = ['juld', 'latitude','longitude', 'deph']

    fldHDR_dat = ['cycle_number', 'data_centre','data_mode', 'dc_reference', 'pi_name', 'wmo_inst_type', 'platform_number']


    fldHDR_dat = ['dc_reference', 'wmo_inst_type', 'platform_number' ]

    fldVAR_dat = ['proc', 'qc', 'data', 'clmn','clsd','erme','erur','resi']
    return fldHDR_pos, fldHDR_dat, fldVAR_dat


def datenum(year, month, day):
    """Return a MATLAB-compatible serial date number for the given calendar date."""
    return (date(year, month, day) - date(1, 1, 1)).days + 367


def datenum_from_str(s, fmt):
    """
    Parse a date string and return a MATLAB-compatible serial date number.

    Parameters
    ----------
    s   : str - date string (e.g. '2004-01-15')
    fmt : str - strptime format string (e.g. '%Y-%m-%d')

    Returns
    -------
    int - serial date number (days since 0001-01-01 + 366)
    """
    d = datetime.strptime(s, fmt).date()
    return (d - date(1, 1, 1)).days + 367

def get_XML_databloc(root, tag_name):
    """
    Reproduit get_XML_databloc de MATLAB :
    cherche le bloc <tag_name> dans le XML et retourne un dict {field: value}.
    """
    block = root.find('.//' + tag_name)
    if block is None:
        return {}
    result = {}
    for child in block:
        result[child.tag] = (child.text or '').strip()
    return result



def PREOA_select(OAHDR, OAVAR, liste_ok):
    """
    Réduit les structures OAHDR/OAVAR selon la liste de profils sélectionnés.

    Parameters
    ----------
    OAHDR    : dict - header et métadonnées
    OAVAR    : dict - données
    liste_ok : array booléen ou d'indices - profils à conserver

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

def parse_num_string(val):
    """
    Équivalent de str2num MATLAB : parse une chaîne comme "1 2 3", "[1 2 3]",
    "1,2,3", ou un scalaire "42.5".
    Retourne un np.ndarray, ou un float si scalaire, ou la chaîne originale si échec.
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

def ANA_ini(conf_file_name):
    """
    Crée et retourne le dictionnaire INIT depuis le fichier de configuration XML.

    Parameters
    ----------
    conf_file_name : str - chemin complet vers le fichier de configuration XML

    Returns
    -------
    INIT : dict - structure de configuration
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

def POSTOA_main(config_fname, DATE_EST, VarName, PLOT_DISP,
                nam_atlas_pres, nam_clim):
    """
    Post-traitement OA : agrège les fichiers par zone, fusionne les champs,
    ajoute la climatologie et sauvegarde les fichiers de résultats.

    Parameters
    ----------
    config_fname   : str  - chemin du fichier de configuration XML
    DATE_EST       : list - [dd, mm, yyyy]
    VarName        : str  - 'TEMP', 'PSAL', 'DOXY'
    PLOT_DISP      : int  - 0 = pas de plots
    nam_atlas_pres : str  - fichier atlas pression
    nam_clim       : str  - fichier climatologie
    """

    # ------------------------------------------------------------------
    #  Configuration
    # ------------------------------------------------------------------
    INIT = ANA_ini(config_fname)
    INIT['PLOT_DISP'] = PLOT_DISP
    lev_plt = -1 if PLOT_DISP == 0 else INIT.PostoaPlotLev

    if VarName == 'DOXY':
        INIT['nam_clim'] = INIT['FilenameClimOxy']
        INIT['nam_std']  = INIT['FilenameStdOxy']
        ClimRefOA     = str(INIT.ClimRefOAOxy).strip()
    else:
        INIT['nam_clim'] = INIT['FilenameClimTS']
        INIT['nam_std']  = INIT['FilenameStdTS']
        ClimRefOA     = str(INIT['ClimRefOATS']).strip()

    area_lim    = INIT['area_limits']
    filt_errmax = INIT['PostoaFiltErrMax']
    LANG        = getattr(INIT, 'LANG', 'En')

    # ------------------------------------------------------------------
    #  Date
    # ------------------------------------------------------------------
    dd, mm, yy   = DATE_EST
    nam_year     = str(yy)
    nam_month    = f'{mm:02d}'
    nam_day      = f'{dd:02d}'
    jest_abs     = datenum(yy, mm, dd)
    dstr         = f'{yy}{mm:02d}{dd:02d}'

    # ------------------------------------------------------------------
    #  Messages
    # ------------------------------------------------------------------
    def log(msg, *args):
        line = msg % args if args else msg
        print(line, end='')
        postoa_mess.write(line)

    # ------------------------------------------------------------------
    #  Champs
    # ------------------------------------------------------------------
    fldHDR_pos, fldHDR_dat, fldVAR_dat = flddat_init()
    nb_field  = len(fldHDR_dat)
    nb_fldvar = len(fldVAR_dat)

    # ------------------------------------------------------------------
    #  Répertoires
    # ------------------------------------------------------------------
    dir_POSTOA_log   = os.path.join(INIT['DirPostoa'], 'POSTOA_LOG')
    dir_POSTOA_field = os.path.join(INIT['DirPostoa'], 'field', str(yy))
    dir_POSTOA_data  = os.path.join(INIT['DirPostoa'], 'data',  str(yy))
    dir_OA_DATA      = os.path.join(INIT['DirPostoa'], 'POSTOA_DATA', VarName)
    dir_clim         = INIT['DirIsasConfStd']

    for d in [dir_POSTOA_log, dir_POSTOA_field,
              dir_POSTOA_data, dir_OA_DATA]:
        Path(d).mkdir(parents=True, exist_ok=True)

    if lev_plt >= 0:
        dir_POSTOA_plot = os.path.join(INIT['DirPostoa'], 'POSTOA_PLOT')
        Path(dir_POSTOA_plot).mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    #  Log
    # ------------------------------------------------------------------
    id_date = f'{dstr}_{datetime.now().strftime("%Y%m%d%H%M%S")}'
    postoa_mess_file = os.path.join(
        dir_POSTOA_log, f'postoa_mess_{VarName}_{id_date}.txt')
    try:
        postoa_mess = open(postoa_mess_file, 'w')
    except OSError:
        print(f'\n  pb opening log file {postoa_mess_file}\n')
        return

    log(f'\n\n >>>>>>>  Running {SOFTWARE}\n')
    log(f'\n {INIT['FilenameAna']} {VarName}  {dstr}\n')
    log(f' dir_OA_DATA       : {dir_OA_DATA}\n')
    log(f' dir_POSTOA_field  : {dir_POSTOA_field}\n')
    log(f' dir_POSTOA_data   : {dir_POSTOA_data}\n')
    log(f' dir_clim          : {dir_clim}\n')

    # ------------------------------------------------------------------
    #  Noms de fichiers
    # ------------------------------------------------------------------
    sep = '' if str(INIT['FilenameAna']).endswith('_') else '_'
    POSTOA_FILENAME = f'{INIT['FilenameAna']}{sep}{dstr}'
    OA_FNAM         = f'OA_{dstr}'

    # ------------------------------------------------------------------
    #  Coordonnées climatologie
    # ------------------------------------------------------------------
    with Dataset(nam_atlas_pres, 'r') as nc:
        str_ref_date_time = getattr(nc, 'reference_date', '19500101T000000Z')
        lat_clim = nc.variables['latitude'][:]
        lon_clim = nc.variables['longitude'][:]
        dep_clim = nc.variables['depth'][:]

    if isinstance(str_ref_date_time, (list, np.ndarray)):
        str_ref_date_time = ''.join(str_ref_date_time).strip()

    if len(str_ref_date_time) > 16:
        jref = datenum_from_str(str_ref_date_time[:10], '%Y-%m-%d')
    else:
        jref = datenum_from_str(str_ref_date_time[:8],  '%Y%m%d')

    lat_lim = area_lim[:2]
    lon_lim = area_lim[2:4]

    def find_idx(arr, val, side='last'):
        idx = np.where(arr <= val)[0] if side == 'last' else np.where(arr >= val)[0]
        return int(idx[-1]) if side == 'last' and len(idx) else \
               int(idx[0])  if side == 'first' and len(idx) else \
               (0 if side == 'last' else len(arr) - 1)

    i1clim = find_idx(lat_clim, lat_lim[0], 'last')
    i2clim = find_idx(lat_clim, lat_lim[1], 'first')
    j1clim = find_idx(lon_clim, lon_lim[0], 'last')
    j2clim = find_idx(lon_clim, lon_lim[1], 'first')
    idx_lim = [i1clim, i2clim, j1clim, j2clim]

    lat_ana = lat_clim[i1clim:i2clim+1]
    lon_ana = lon_clim[j1clim:j2clim+1]

    # ------------------------------------------------------------------
    #  Comptage des fichiers dat disponibles
    # ------------------------------------------------------------------
    t_start = time.time()
    
    nb_datamax  = 0
    nb_levmax   = 0
    fnam_ok     = None

    INIT['PreoaAreaList'] = parse_matlab_range_list(INIT['PreoaAreaList'])
    nb_ana_area = len(INIT['PreoaAreaList'])
    list_ok     = np.zeros(nb_ana_area, dtype=bool)
    for ia, iarea in enumerate(INIT['PreoaAreaList']):
        str_area     = f'{iarea:03d}'
        fname_OA_dat = os.path.join(f"{dir_OA_DATA}/{nam_year}/{nam_month}",
                                    f'{OA_FNAM}_{str_area}_dat_{VarName}.nc')
        if os.path.isfile(fname_OA_dat):
            log(f'Reading data area {iarea:4d}\n')
            fnam_ok = fname_OA_dat
            with Dataset(fname_OA_dat, 'r') as nc:
                nb_prof_i = len(nc.dimensions['N_PROF'])
                nb_lev_i  = len(nc.dimensions['N_LEVELS'])
            list_ok[ia]  = True
            nb_datamax  += nb_prof_i
            nb_levmax    = max(nb_levmax, nb_lev_i)
        else:
            log(f'Area {iarea:4d} missing\n')

    isok_area  = np.where(list_ok)[0]
    num_areaok = len(isok_area)
    log(f'{num_areaok} areas found\n')

    # ------------------------------------------------------------------
    #  Lecture et agrégation des fichiers dat
    # ------------------------------------------------------------------
    OAHDR = {'nbdat': 0}
    OAVAR = {}

    for fld in fldHDR_pos[:3]:
        OAHDR[fld] = np.full(nb_datamax, np.nan)
    OAHDR['deph'] = np.full(nb_levmax, np.nan)

    # Lit un fichier pour obtenir les dimensions des champs HDR
    STDHDR, _ = NCR_data(fnam_ok, VarName)

    for fld in fldHDR_dat:
        val = STDHDR.get(fld)
        if val is None:
            OAHDR[fld] = np.zeros(nb_datamax, dtype=object)
            continue
        arr = np.asarray(val)
        if arr.ndim == 2:
            n2 = arr.shape[1]
            OAHDR[fld] = np.zeros((nb_datamax, n2), dtype=arr.dtype)
        else:
            OAHDR[fld] = np.zeros(nb_datamax, dtype=arr.dtype)
    for fld in fldVAR_dat:
        if fld == 'proc':
            OAVAR[fld] = np.zeros(nb_datamax)
        else:
            OAVAR[fld] = np.zeros((nb_levmax, nb_datamax))

    i_cur = 0

    for ii in range(num_areaok):
        ia    = isok_area[ii]
        iarea = INIT['PreoaAreaList'][ia]
        str_area     = f'{iarea:03d}'
        fname_OA_dat = os.path.join(f"{dir_OA_DATA}/{nam_year}/{nam_month}",
                                    f'{OA_FNAM}_{str_area}_dat_{VarName}.nc')
        print(fname_OA_dat)
        STDHDR, STDVAR = NCR_data(fname_OA_dat, VarName)

        if STDHDR['nbdat'] > 0:
            i_new  = i_cur + STDHDR['nbdat']
            nz_max = len(np.asarray(STDHDR['deph']).flatten())

            for fld in fldHDR_pos[:3]:
                OAHDR[fld][i_cur:i_new] = np.asarray(STDHDR[fld]).flatten()

            OAHDR['deph'][:nz_max] = np.asarray(STDHDR['deph']).flatten()

            for fld in fldHDR_dat:
                if fld in STDHDR and STDHDR[fld] is not None:
                    OAHDR[fld][i_cur:i_new] = np.asarray(STDHDR[fld])

            for fld in fldVAR_dat:
                if fld not in STDVAR or STDVAR[fld] is None:
                    continue
                arr = np.asarray(STDVAR[fld])
                if fld == 'proc':
                    OAVAR[fld][i_cur:i_new] = arr.flatten()
                else:
                    if arr.ndim == 2:
                        OAVAR[fld][:nz_max, i_cur:i_new] = arr[:nz_max, :]
                    elif arr.ndim == 1:
                        OAVAR[fld][:nz_max, i_cur:i_new] = arr[:nz_max, np.newaxis]

            i_cur          = i_new
            OAHDR['nbdat'] = i_cur

    if num_areaok == 0:
        log(f'\n No area file found on {f"{dir_OA_DATA}/{nam_year}/{nam_month}"}, POSTOA stopped\n')
        postoa_mess.close()
        return

    # ------------------------------------------------------------------
    #  Suppression des doublons sur dc_reference
    # ------------------------------------------------------------------
    dc_ref = OAHDR.get('dc_reference')
    if dc_ref is not None:
        arr = np.asarray(dc_ref[:OAHDR['nbdat']])
        if arr.ndim == 2:
            rows    = [arr[i].tobytes() for i in range(arr.shape[0])]
            _, idx  = np.unique(rows, return_index=True)
        else:
            _, idx = np.unique(arr, return_index=True)

        if len(idx) < OAHDR['nbdat']:
            mask = np.zeros(OAHDR['nbdat'], dtype=bool)
            mask[idx] = True
            OAHDR_FIN, OAVAR_FIN = PREOA_select(OAHDR, OAVAR, mask)
        else:
            OAHDR_FIN = OAHDR
            OAVAR_FIN = OAVAR
    else:
        OAHDR_FIN = OAHDR
        OAVAR_FIN = OAVAR

    # Normalise les longitudes
    lon = np.asarray(OAHDR_FIN.get('longitude', []))
    lon[lon < -180] += 360
    lon[lon >  180] -= 360
    OAHDR_FIN['longitude'] = lon

    # ------------------------------------------------------------------
    #  Sauvegarde fichier dat
    # ------------------------------------------------------------------
    fname_POSTOA_dat = os.path.join(
        dir_POSTOA_data, f'{POSTOA_FILENAME}_dat_{VarName}.nc')

    GLOB_ATT = globatt_init(INIT, SOFTWARE)
    GLOB_ATT['COMMENT']       = 'V8.0 T and S fields'
    GLOB_ATT['TITLE']         = 'Monthly analysis'
    GLOB_ATT['reference_date_time'] = str_ref_date_time
    GLOB_ATT['ANALYSIS_NAME'] = INIT['FilenameAna']
    GLOB_ATT['SOUTH_LAT']     = str(lat_ana[0])
    GLOB_ATT['NORTH_LAT']     = str(lat_ana[-1])
    GLOB_ATT['WEST_LONG']     = str(lon_ana[0])
    GLOB_ATT['EAST_LONG']     = str(lon_ana[-1])
    GLOB_ATT['LAT_RESO']      = float(np.max(np.abs(np.diff(lat_ana))))
    GLOB_ATT['LON_RESO']      = float(np.max(np.abs(np.diff(lon_ana))))
    GLOB_ATT['DATE_START']    = f'{nam_year}-{nam_month}-{nam_day}'
    GLOB_ATT['DATE_STOP']     = GLOB_ATT['DATE_START']

    if OAHDR_FIN['nbdat'] > 0:
        log(f'\n Number of OA data: {OAHDR_FIN["nbdat"]}\n')
        OAHDR_FIN['reference_date_time'] = str_ref_date_time
        OAHDR_FIN['data_type']           = 'ISAS-DataSet    '
        msg_error1 = NCW_data_hdr_and_var(fname_POSTOA_dat, GLOB_ATT, OAHDR, VarName, OAVAR)
        #msg_error1 = NCW_data_hdr_PREOA(fname_POSTOA_dat, GLOB_ATT, OAHDR_FIN)
        log(f'\n{msg_error1}\n')

        #if msg_error1[0].upper() in ('O', 'W'):
        #    msg_error = NCW_data_var(fname_POSTOA_dat, VarName, OAVAR_FIN)
        #    log(f'{msg_error}\n')

    tps1 = time.time() - t_start
    log(f' >>>> Processing time for data: {tps1:5.2f}\n\n')

    # ------------------------------------------------------------------
    #  Lecture du premier fichier fld pour les coordonnées
    # ------------------------------------------------------------------
    ia_first  = isok_area[0]
    iarea_first = INIT['PreoaAreaList'][ia_first]
    str_area_first = f'{iarea_first:03d}'
    file_nc_fld_first = os.path.join(
        f"{dir_OA_DATA}/{nam_year}/{nam_month}", f'{OA_FNAM}_{str_area_first}_fld_{VarName}.nc')

    with Dataset(file_nc_fld_first, 'r') as nc:
        dep_ana = nc.variables['depth'][:]

    NX_GLO = len(lon_ana)
    NY_GLO = len(lat_ana)
    NZ_GLO = len(dep_ana)

    # Indices des niveaux dans la climatologie
    idep_clim = np.array([int(np.where(dep_clim == d)[0][0])
                           for d in dep_ana])

    PCT_VAR_TOT = 100 * np.ones((NX_GLO, NY_GLO, NZ_GLO))
    FIELD_TOT   = np.zeros((NX_GLO, NY_GLO, NZ_GLO))

    PARANO = f'{VarName}_ANO'

    # ------------------------------------------------------------------
    #  Agrégation des champs anomalie
    # ------------------------------------------------------------------
    t_ano = time.time()

    for ii in range(num_areaok):
        ia    = isok_area[ii]
        iarea = INIT['PreoaAreaList'][ia]
        str_area    = f'{iarea:03d}'
        file_nc_fld = os.path.join(
            f"{dir_OA_DATA}/{nam_year}/{nam_month}", f'{OA_FNAM}_{str_area}_fld_{VarName}.nc')
        print(file_nc_fld)

        with Dataset(file_nc_fld, 'r') as nc:
            try:
                lon_est       = nc.variables['longitude'][:]
                lat_est       = nc.variables['latitude'][:]
                
                # Skip si fichier vide
                if lon_est.size == 0 or lat_est.size == 0:
                    log(f'  Skipping area {iarea} — empty fld file\n')
                    continue
                anomaly_field = nc.variables[PARANO][:]
                pct_var_name  = f'{PARANO}_PCTVAR'
                if pct_var_name in nc.variables:
                    pct_variance = nc.variables[pct_var_name][:]
                else:
                    pct_variance = np.zeros_like(anomaly_field)
            except:
                log(f'  Skipping area {iarea} — error fld file\n')
                continue
        # Indices dans la grille globale
        i1 = np.where(lat_ana == lat_est[0])[0]
        i2 = np.where(lat_ana == lat_est[-1])[0]
        j1 = np.where(lon_ana == lon_est[0])[0]
        j2 = np.where(lon_ana == lon_est[-1])[0]

        if not (len(i1) and len(i2) and len(j1) and len(j2)):
            log(f'\n  AREA_LIMITS incompatible with selected area\n')
            postoa_mess.close()
            return

        i1, i2, j1, j2 = int(i1[0]), int(i2[0]), int(j1[0]), int(j2[0])

        # Squeeze dimension time si présente
        if anomaly_field.ndim == 4:
            anomaly_field = anomaly_field[0]   # (depth, lat, lon)
            pct_variance  = pct_variance[0]

        # Réorganise en (lon, lat, depth) si nécessaire
        if anomaly_field.shape == (NZ_GLO, i2-i1+1, j2-j1+1):
            anomaly_field = anomaly_field.transpose(2, 1, 0)
            pct_variance  = pct_variance.transpose(2, 1, 0)

        pct_variance  = np.where(np.isnan(pct_variance), 100, pct_variance)
        anomaly_field = np.where(pct_variance <= 0,   np.nan, anomaly_field)
        anomaly_field = np.where(pct_variance >= 100, np.nan, anomaly_field)

        isok = np.isfinite(anomaly_field)
        if np.any(isok):
            tab = FIELD_TOT[j1:j2+1, i1:i2+1, :]
            tab[isok] = anomaly_field[isok]
            FIELD_TOT[j1:j2+1, i1:i2+1, :] = tab

            tab = PCT_VAR_TOT[j1:j2+1, i1:i2+1, :]
            tab[isok] = pct_variance[isok]
            PCT_VAR_TOT[j1:j2+1, i1:i2+1, :] = tab

    tps2 = time.time() - t_ano
    log(f' ... Processing time for field anomaly: {tps2:5.2f}\n')

    # ------------------------------------------------------------------
    #  Filtrage aux limites de zones
    # ------------------------------------------------------------------
    t_smooth = time.time()
    if filt_errmax < 100:
        FIELD_TOT, PCT_VAR_TOT = POSTOA_filter(
            INIT, idx_lim, FIELD_TOT, PCT_VAR_TOT, filt_errmax)

    tps3 = time.time() - t_smooth
    log(f' ... Processing time for smoothing: {tps3:5.2f}\n')

# ------------------------------------------------------------------
    #  Ajout climatologie → champ absolu (version optimisée)
    # ------------------------------------------------------------------
    t_abs = time.time()

    # Lecture directe de la climatologie — une seule ouverture de fichier
    with Dataset(nam_clim, 'r') as nc:
        var_list = list(nc.variables.keys())
        nam_dep  = 'depth' if 'depth' in var_list else 'DEPTH'

        # Lecture des niveaux sélectionnés en une seule opération
        clim_var = nc.variables[VarName]

        # Détecte fill_value et scaling
        fill_val    = getattr(clim_var, '_FillValue', None)
        scale_fact  = float(getattr(clim_var, 'scale_factor', 1.0))
        add_offset  = float(getattr(clim_var, 'add_offset',   0.0))
        scaling     = not (scale_fact == 1.0 and add_offset == 0.0)

        # Lit le bloc (NZ_dep_clim, i1:i2, j1:j2) directement
        # Shape NetCDF : (time, depth, lat, lon) ou (depth, lat, lon)
        raw = clim_var[:]
        raw = np.squeeze(raw)   # supprime dim time si présente

        # Transpose en (lon, lat, depth) = (NX, NY, NZ)
        if raw.ndim == 3:
            raw = raw.transpose(2, 1, 0)  # (depth, lat, lon) → (lon, lat, depth)

        # Extrait la zone et les niveaux utiles
        clim_zone = raw[j1clim:j2clim+1, i1clim:i2clim+1, :]  # (NX, NY, NZ_clim)

        if fill_val is not None:
            clim_zone = np.where(clim_zone == fill_val, np.nan,
                                 clim_zone.astype(float))
        if scaling:
            clim_zone = clim_zone * scale_fact + add_offset

    # Sélectionne les niveaux d'analyse et ajoute en une seule opération vectorisée
    clim_sel = clim_zone[:, :, idep_clim]  # (NX, NY, NZ_GLO)
    FIELD_TOT = FIELD_TOT + clim_sel

    PCT_VAR_TOT = np.where(np.isnan(PCT_VAR_TOT), 100, PCT_VAR_TOT)

    # ------------------------------------------------------------------
    #  Calcul erreur absolue — vectorisé
    # ------------------------------------------------------------------
    ConvVarAna = np.sqrt(
        INIT['FactVar'] * (INIT['W1'] + INIT['W2'])
        / (INIT['W1'] + INIT['W2'] + INIT['W3']))

    nam_atlas_std = os.path.join(
        INIT['DirIsasConfStd'], f"{INIT['nam_std']}STD_{VarName}.nc")

    with Dataset(nam_atlas_std, 'r') as nc:
        std_var  = nc.variables[f'{VarName}_STD'][:]
        std_fill = getattr(nc.variables[f'{VarName}_STD'], '_FillValue', None)

    var_std = np.squeeze(std_var)
    if var_std.ndim == 3:
        var_std = var_std.transpose(2, 1, 0)  # (lon, lat, depth)
    elif var_std.ndim == 4:
        var_std = var_std.transpose(3, 2, 1, 0)[..., 0]

    if std_fill is not None:
        var_std = np.where(var_std == std_fill, np.nan, var_std.astype(float))

    i, j, k  = PCT_VAR_TOT.shape
    var_std   = ConvVarAna * var_std[:i, :j, :k]
    var_error = var_std * np.sqrt(PCT_VAR_TOT * 0.01)

    tps4 = time.time() - t_abs
    log(f' ... Processing time for absolute field: {tps4:5.2f}\n')


    # ------------------------------------------------------------------
    #  Sauvegarde fichier fld
    # ------------------------------------------------------------------
    fname_POSTOA_field = os.path.join(
        dir_POSTOA_field, f'{POSTOA_FILENAME}_fld_{VarName}.nc')

    jest = jest_abs - jref

    NCW_4D_create(fname_POSTOA_field, GLOB_ATT, lon_ana, lat_ana, dep_ana)
    NCW_4D_variable(fname_POSTOA_field, VarName)
    NCW_4D_variable(fname_POSTOA_field, f'{VarName}_ERR')
    msg_error = NCW_4D_variable(fname_POSTOA_field, f'{VarName}_PCTVAR')

    with Dataset(fname_POSTOA_field, 'a') as nc:
        nc.variables['time'][0] = float(jest)

        # FIELD_TOT shape: (lon, lat, depth) → (time, depth, lat, lon)
        nc.variables[VarName][0, :, :, :]          = FIELD_TOT.transpose(2, 1, 0).astype(np.float32)
        nc.variables[f'{VarName}_ERR'][0, :, :, :] = var_error.transpose(2, 1, 0).astype(np.float32)
        nc.variables[f'{VarName}_PCTVAR'][0, :, :, :] = PCT_VAR_TOT.transpose(2, 1, 0).astype(np.float32)
        log(f'\n{msg_error}\n')

    tps5 = time.time() - t_ano
    log(f' >>>> Processing time for field, total {tps5:5.2f}\n')

    log(f'\n  {SOFTWARE} processing done  \n \n')
    postoa_mess.close()

def POSTOA_launcher(inpt):
    """
    Launch the OA post-processing for a given month, reading all paths from
    config.ini located next to this script.

    Parameters
    ----------
    inpt : str - date in 'YYYYMM' format (e.g. '196001')
    """

    cfg = configparser.ConfigParser(interpolation=None)
    cfg.read(Path(__file__).parent / 'config.ini')

    yy = int(inpt[0:4])
    mm = int(inpt[4:6])
    month_str = f'{mm:02d}'

    VarName    = cfg['analysis']['var_name']
    config_fname  = cfg['paths']['config_fname']
    inpath        = cfg['paths']['dir_oa_data']
    POSTOAPATH    = cfg['paths']['dir_postoa']
    clim          = cfg['paths']['clim_template'].format(
                        year=inpt[0:4], month=inpt[4:6], var=VarName)
    nam_atlas_pres = cfg['paths']['atlas_pres_template'].format(month=month_str)

    Path(f'{POSTOAPATH}{VarName}/{yy}/{month_str}').mkdir(parents=True, exist_ok=True)
    os.system(f"cp {inpath}{VarName}/{yy}/{month_str}/*.nc "
              f"{POSTOAPATH}{VarName}/{yy}/{month_str}")

    PLOT_DISP = 0

    print(f'Processing {inpt}')

    date_ana = [15, mm, yy]

    POSTOA_main(config_fname, date_ana, VarName, PLOT_DISP,
                nam_atlas_pres, clim)

#POSTOA_launcher('200401')
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python POSTOA_launcher.py <YYYYMM>')
        print('Example: python POSTOA_launcher.py 196001')
        sys.exit(1)

    POSTOA_launcher(sys.argv[1])
