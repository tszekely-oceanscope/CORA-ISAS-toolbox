# CORA-ISAS-toolbox

Python and Fortran 90 toolbox for producing global monthly ocean temperature and salinity fields from the CORA dataset using the ISAS objective analysis method.

This code was developed to produce the dataset described in:

> Szekely, T. and Miadana, J.R.  (2026).The consistent in-situ global gridded temperature and salinity dataset, CORA OA 1960-2024 .submitted to Earth System Science Data. 

Dataset :Global Ocean in situ - Delayed Mode temperature and salinity CORA -objective analysis. https://doi.org/10.48670/mds-00383 


---

## Overview

The toolbox is split into three sequential steps:

```
CORA (NetCDF) --> [PREOA] --> Ocean zones --> [OA] --> Analyzed zones --> [POSTOA] --> Global monthly field (NetCDF)
```

| Module | Language | Role |
|--------|----------|------|
| `PREOA` | Python | Reads and formats CORA NetCDF input files, splits the ocean into zones to reduce computational load |
| `OA` | Fortran 90 | Performs the ISAS objective analysis over each ocean zone |
| `POSTOA` | Python | Assembles the analyzed ocean zones into a single global monthly NetCDF file |

---

## Requirements

### Python
- Python 3.8+
- `numpy`
- `netCDF4`
- `scipy`

Install dependencies with:
```bash
pip install numpy netCDF4 scipy
```

### Fortran
- gfortran (or equivalent Fortran 90 compiler)

Compile the OA module with:
```bash
cd OA/
make
```

---

## Input data

This toolbox expects **CORA 5.2** input files (global in-situ temperature and salinity measurements):

> Szekely, T., Gourrion, J., Pouliquen, S., and Reverdin, G.: The CORA 5.2 dataset: global in-situ temperature and salinity measurements, SEANOE [data set], https://doi.org/10.17882/46219, 2019.

---

## Usage

Run the three steps sequentially:

**Step 1 — Pre-processing**
```bash
python PREOA/PREOA_launcher.py --input /path/to/CORA/ --output /path/to/zones/ --year 2010 --month 01
```

**Step 2 — Objective analysis**
```bash

cd OA/
python launcher.py
Runn all the lines of commandlist.lst
```

**Step 3 — Post-processing**
```bash
python POSTOA/POSTOA_launcher.py --input /path/to/zones/ --output /path/to/output/ --year 2010 --month 01
```

---

## Output

The final output is a global monthly NetCDF file containing gridded temperature and salinity fields, following the ISAS format convention.

The dataset produced using this toolbox is available at:
> https://doi.org/10.48670/mds-00383

---

## Reference for the objective analysis method

> Gaillard, F., Reynaud, T., Thierry, V., Kolodziejczyk, N., and von Schuckmann, K.: In-situ based reanalysis of the global ocean temperature and salinity with ISAS: variability of the heat content and steric height, J. Climate, 29, 1305–1323, https://doi.org/10.1175/JCLI-D-15-0028.1, 2016.

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Contact

Tanguy Szekely — OceanScope - tanguy.szekely@ocean-scope.com

