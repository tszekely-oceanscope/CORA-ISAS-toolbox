import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.ar_model import AutoReg
from statsmodels.nonparametric.smoothers_lowess import lowess

import matplotlib
matplotlib.use('Agg')
# Charger les données


col = dict()
col['CORA']="#D55E00"
col['IAP'] = "#009E73"
col['EN4'] = "#0072B2"
col['ISAS 20'] = "#000000"
col['ERSST'] = "#E69F00"
col['SSS'] = "#56B4E9"
plt.figure(figsize=(10, 2.5))
fid = open('./trand_list.csv','w')
fid.write('zone,dataset,trend,ci_90\n')
for zone in ["Global","Pacific","Atlantic","Indian","South"]:
 print(zone)
 records=[]
 ind = -1
 
 for suf in ['CORA','ERSST','IAP','EN4','ISAS20']:
    print(f"    {suf}")
    ind+=1
    df = pd.read_csv('Mean_SST-' + suf +'.csv')
    df['date'] = pd.to_datetime(df['date'])
    df['year'] = df['date'].dt.year
    annual_mean = df.groupby('year')[zone].mean()


    df = df.sort_values("date")
    df["date"] = pd.to_datetime(df["date"])
    years = df["date"].dt.year + (df["date"].dt.dayofyear - 1) / 365.25
    years = years.to_numpy()
    values = df[zone].values


  # LOWESS principal
    frac = 10 / (years[-1] - years[0] + 1)
    smoothed = lowess(values, years, frac=frac, return_sorted=False)

  # Résidus et fit AR(1)
    residuals = values - smoothed
    model_ar1 = AutoReg(residuals, lags=1, old_names=False).fit()
    phi = model_ar1.params[1]
    sigma = model_ar1.resid.std()

  # Monte Carlo
    n_sim = 500
    trend_diffs = []

    for _ in range(n_sim):
      ar_sim = [0]
      for _ in range(1, len(years)):
        ar_sim.append(phi * ar_sim[-1] + np.random.normal(0, sigma))
      ar_sim = np.array(ar_sim)
      surrogate_series = smoothed + ar_sim
      surrogate_lowess = lowess(surrogate_series, years, frac=frac, return_sorted=False)
      trend_diffs.append(surrogate_lowess[-1] - surrogate_lowess[0])

    trend_diffs = np.array(trend_diffs)
    mean_trend = np.mean(trend_diffs)
    std_trend = np.std(trend_diffs)
    ci90 = 1.65 * std_trend


    mean_trend_century = mean_trend / (years[-1] - years[0] + 1)* 100.
    ci90_century = ci90 * 100./ (years[-1] - years[0] + 1)

    print(f"Trend: {mean_trend_century:.3f} ± {ci90_century:.3f} {df.columns[1]}/century "
      f"over {years[-1]-years[0]} years at 90% CI")
    if "CORA" in suf:
      suf = "CORA"
    if suf == 'ISAS20':
      suf = 'ISAS 20' 
    fid.write(f"{zone},{suf},{mean_trend_century:.3f},{ci90_century:.3f}\n")
    plt.plot(annual_mean.index, annual_mean.values-np.mean(annual_mean.values), color=col[suf], label=f'{suf} ' + f"  {mean_trend_century:.3f} ± {ci90_century:.3f} ({90}% CI)")


 plt.title(zone + ' ocean')
 plt.xlabel("Year")
 plt.ylabel("Mean SST (°C)")
 plt.legend()
 plt.grid(True)
 plt.savefig(f"SST_trend_with_CI_{zone}.png", dpi=150)
 plt.clf()
fid.close()


