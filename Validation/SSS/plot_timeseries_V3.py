import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

import matplotlib
matplotlib.use("Agg")  # Pas besoin de X11

plt.figure(figsize=(10, 2.5))
records=[]


col = dict()
col['CORA']="#D55E00"
col['IAP'] = "#009E73"
col['EN4'] = "#0072B2"
col['ISAS 20'] = "#000000"
col['ERSST'] = "#E69F00"
col['SSS L4'] = "#56B4E9"
for reg in ['Global','Pacific','Atlantic','Indian','South']:
 ind = -1
 for suf in ['CORA','EN4','IAP','SMOS','ISAS20']:
  ind +=1 
  df = pd.read_csv('Mean_SSS_' + suf + '.csv')  
        # Parse dates
  df['date'] = pd.to_datetime(df['date'])
  df = df[(df['date'] >= '1960-01-01') & (df['date'] <= '2024-12-31')]
  df['year'] = df['date'].dt.year
  annual_mean = df.groupby('year')[reg].mean()
  if suf == 'SMOS':
    suf = 'SSS L4'
  if suf == 'ISAS20':
    suf = 'ISAS 20'
  plt.plot(annual_mean.index, annual_mean.values, color=col[suf], label=f'{suf}')

 plt.title(reg + ' ocean')
 plt.ylabel("PSS-78")
 plt.grid()
 plt.legend()
 plt.savefig('./Mean_SSS_' + reg + '.png')
 plt.clf()

