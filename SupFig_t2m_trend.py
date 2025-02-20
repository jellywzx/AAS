import numpy as np
import matplotlib.pyplot as plt
import glob
import matplotlib.dates as mdates
from pylab import *
from matplotlib import dates, ticker
from scipy import stats
import xarray as xr
import matplotlib.path as mpath
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import pandas as pd
import scipy.signal as signal
import cmaps
from mpl_toolkits.basemap import Basemap, shiftgrid
from scipy.stats import t
from global_land_mask import globe

def independent_ttest(mean1,mean2,std1,std2,n1,n2, alpha):
    # calculate means
    # calculate standard errors
    #se1, se2 = sem(data1), sem(data2)
    se1, se2 = std1/sqrt(n1), std2/sqrt(n2)
    # standard error on the difference between the samples
    sed = sqrt(se1**2.0 + se2**2.0)
    # calculate the t statistic
    t_stat = (mean1 - mean2) / sed
    # degrees of freedom
    df = n1 + n2 - 2
    # calculate the critical value
    cv = t.ppf(1.0 - alpha, df)
    # calculate the p-value
    p = (1.0 - t.cdf(abs(t_stat), df)) * 2.0
    # return everything
    return p
        
t2m = xr.open_dataset('D:/ERA5/ERA5-single-level-1x1-resolution-2meter-airtemp-1979-202208.nc')['t2m'].sel(time=slice('1992','2021'))
lons = t2m.longitude
lats = t2m.latitude

t2m_ON = t2m[(t2m.time.dt.month==10)|(t2m.time.dt.month==11)].values.reshape(30,2,181,360).mean(axis=1)
# t2m_ON = t2m[(t2m.time.dt.month==10)|(t2m.time.dt.month==11)].values.reshape(2,30,181,360).mean(axis=0)
years = np.linspace(1992,2021,30,dtype=int)
aa1 = np.zeros((181,360),'float')
pa1 = np.zeros((181,360),'float')

for i in range(181):
    for j in range(360):
        if (np.isnan(t2m_ON[:,i,j]).any()):
            aa1[i,j] = np.nan
            pa1[i,j] = np.nan
        else:
            #aa1[i,j]= np.corrcoef(va1[0:39],temp1[0:39,i,j])[0,1]
    #        aa1[i,j], intercept, r_value, pa1[i,j], std_err = stats.linregress(year[20:44],data[20:44,i,j])
            aa1[i,j], intercept, r_value, pa1[i,j], std_err = stats.linregress(years,t2m_ON[:,i,j])

levels1=np.linspace(-0.8,0.8,11)

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 18
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 20  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6
plt.close()
fig=plt.figure(figsize=(8,10))
ax = fig.add_subplot(111)
figname = '(c) ON t2m trend'
ax.set_title(figname)
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
x1,y1 = m(*np.meshgrid(lons,lats))
im1 = m.contourf(x1,y1,aa1*10,
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)
# m.contour(x1,y1,data_clm,levels=[15],colors='black',linewidth=4)
csp = m.contourf(x1,y1, 1 - pa1, levels=[0.9,1] ,colors='none',hatches=['..', None],alpha=0)
m.drawcoastlines(color='blue')
# m.fillcontinents(color='lightgray')
m.drawmeridians([90,120,150,180,210,240,270,300,330], linewidth=1.2,labels=[True,False,False,True],fmt='%g',latmax=80,)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,True,False,True],latmax=80,)
cb = plt.colorbar(im1, location='bottom', ticks =levels1[::2], extend='both', extendfrac='auto',aspect=30,shrink=0.6, pad=0.05)
cb.set_label('Trend in t2m per decade (K)',)
cb.ax.tick_params(width=0)
fig.tight_layout(rect=[0,0.08,1,0.95])
plt.savefig('C:/Users/fzjxw/python/code/Figures/t2m_trend.png',dpi=300,bbox_inches='tight')
# plt.show()