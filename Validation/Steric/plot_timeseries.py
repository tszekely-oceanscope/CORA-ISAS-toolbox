import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

import matplotlib
matplotlib.use("Agg")  # Pas besoin de X11


def moving_average(vector, window_size=6):
    return np.convolve(vector, np.ones(window_size)/window_size, mode='valid')
# File pattern (adjust if filenames change)


col = dict()
col['CORA2026']="#D55E00"
col['IAP'] = "#009E73"
col['EN4'] = "#0072B2"
col['ISAS 20'] = "#000000"
col['ERSST'] = "#E69F00"
col['SSS L4'] = "#56B4E9"

plt.figure(figsize=(10, 2.5))
records=[]
ind = -1
for suf in ['CORA2026','IAP','EN4','ISAS20']:
 for lev in [1]:
  ind +=1 
  df = pd.read_csv('Mean_HSH_' + suf +'.csv')  
  df['date'] = pd.to_datetime(df['date'])

  df = df[(df['date'] < '1990-01-01') | (df['date'] >= '1992-12-31')]

  df['year'] = df['date'].dt.year
  annual_mean = df.groupby('year')['HSH'].mean()

  if suf == 'ISAS20':
    suf = 'ISAS 20'
  plt.plot(annual_mean.index, -1*(annual_mean.values-np.mean(annual_mean.values))*1000., color=col[suf], label=f'{suf} ')



plt.ylabel("Mean HSSL anomaly (mm)")
plt.grid()
plt.legend()
plt.savefig('./HSSL_Global2026.png')


