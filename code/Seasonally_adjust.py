## This code seasonally adjusts the timeseries of airfoil (A) and flag (F) motion
import scipy.io
import numpy as np
import pandas as pd
from statsmodels.tsa.seasonal import MSTL 
import matplotlib.pyplot as plt

files_dir = "./../data/"
filename = "03_19hz_ct_s2no_d008"

stidx=2000; endidx=3001 # Plot
period_foil = 19 ###### Get freq from fft
period_flag = 19 ###### Get freq from fft
period_flagfront = 19 ###### Get freq from fft

## Foil
data = scipy.io.loadmat(files_dir+filename+'_Foil.mat')

time = np.transpose(np.array(data['time']))
time = np.reshape(time,(time.shape[0],))
foilang=np.array(data['foilang'])
data1 = pd.DataFrame(data=foilang, index=time)

res1 = MSTL(data1, periods=(period_foil)).fit() 

fig, ax = plt.subplots(nrows=5, figsize=[7,7])
data1[0].iloc[stidx:endidx].plot(ax=ax[0])
ax[0].set_ylabel("Time series")
res1.trend.iloc[stidx:endidx].plot(ax=ax[1])
ax[1].set_ylabel("Trend")
res1.seasonal.iloc[stidx:endidx].plot(ax=ax[2])
ax[2].set_ylabel("Seasonal")
res1.resid.iloc[stidx:endidx].plot(ax=ax[3])
ax[3].set_ylabel("Residual")
res1.trendres = res1.trend + res1.resid
res1.trendres.iloc[stidx:endidx].plot(ax=ax[4])
ax[4].set_ylabel("Trend+Resid")
plt.tight_layout()

## Flag
data = scipy.io.loadmat(files_dir+filename+'_Flag2.mat')

time = np.transpose(np.array(data['time']))
time = np.reshape(time,(time.shape[0],))
flagtipY=np.array(data['flagtipY'])
data1 = pd.DataFrame(data=flagtipY, index=time)

res2 = MSTL(data1, periods=(period_flag)).fit() 

fig, ax = plt.subplots(nrows=5, figsize=[7,7])
data1[0].iloc[stidx:endidx].plot(ax=ax[0])
ax[0].set_ylabel("Time series")
res2.trend.iloc[stidx:endidx].plot(ax=ax[1])
ax[1].set_ylabel("Trend")
res2.seasonal.iloc[stidx:endidx].plot(ax=ax[2])
ax[2].set_ylabel("Seasonal")
res2.resid.iloc[stidx:endidx].plot(ax=ax[3])
ax[3].set_ylabel("Residual")
res2.trendres = res2.trend + res2.resid
res2.trendres.iloc[stidx:endidx].plot(ax=ax[4])
ax[4].set_ylabel("Trend+Resid")
plt.tight_layout()

print(np.array(res2.trend))


## Flag - front pitching

flagfrontY=np.array(data['flagfrontY'])
data1 = pd.DataFrame(data=flagfrontY, index=time)

res3 = MSTL(data1, periods=(period_flagfront)).fit() 

fig, ax = plt.subplots(nrows=5, figsize=[7,7])
data1[0].iloc[stidx:endidx].plot(ax=ax[0])
ax[0].set_ylabel("Time series")
res3.trend.iloc[stidx:endidx].plot(ax=ax[1])
ax[1].set_ylabel("Trend")
res3.seasonal.iloc[stidx:endidx].plot(ax=ax[2])
ax[2].set_ylabel("Seasonal")
res3.resid.iloc[stidx:endidx].plot(ax=ax[3])
ax[3].set_ylabel("Residual")
res3.trendres = res3.trend + res3.resid
res3.trendres.iloc[stidx:endidx].plot(ax=ax[4])
ax[4].set_ylabel("Trend+Resid")
plt.tight_layout()

print(np.array(res3.trend))

## Save
mdic = {"foilang_trend": np.array(res1.trend), 
        "foilang_seas": np.array(res1.seasonal),
        "foilang_resid": np.array(res1.resid)}
scipy.io.savemat(files_dir+filename+'_Foil_SA.mat',mdic)

mdic = {"flagtipY_trend": np.array(res2.trend), 
        "flagtipY_seas": np.array(res2.seasonal),
        "flagtipY_resid": np.array(res2.resid),
        "flagfrontY_trend": np.array(res3.trend), 
        "flagfrontY_seas": np.array(res3.seasonal),
        "flagfrontY_resid": np.array(res3.resid)}
scipy.io.savemat(files_dir+filename+'_Flag2_SA.mat',mdic)


plt.show()

