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
from scipy.signal import detrend


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

def draw_latlon_polygon(bmap, lons, lats, *args, **kwargs):
    '''Plot a polygon in lat/lon coordinates.

    bmap -- Basemap object.
    lons, lats -- Sequences of polygon vertices.
    *args, **kwargs -- Aditional arguments to pyplot.plot().

    You should use 'k-' in *args to draw the lines in color black.

    '''
    if len(lons) != len(lats):
        raise IndexError('lons and lats have different lenghts')
    if lons[-1] != lons[0] or lats[-1] != lats[0]:
        lons = np.concatenate((lons, lons[:1]))
        lats = np.concatenate((lats, lats[:1]))
    n = len(lons) - 1
    res = 10000
    for i in range(n):
        x = np.linspace(lons[i], lons[i + 1], res)
        y = np.linspace(lats[i], lats[i + 1], res)
        x, y = bmap(x, y)
        bmap.plot(x, y, *args, **kwargs)

def sh_project(ax,lon,lat):
    #sps = Basemap(lon_0=180, boundinglat=0, projection='spstere',round=True,resolution='l')
    sps = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=60,\
                llcrnrlon=0,urcrnrlon=360,resolution='l')
    sps.ax = ax
    mer = np.arange(0, 360, 60.)
    #par = np.arange(-90, -50, 10.)
    x, y = sps(*np.meshgrid(lon, lat))
    sps.drawmapboundary(fill_color='AntiqueWhite')
    sps.drawparallels([-90,-45,0,45,90], linewidth=1,labels=[True,True,False, True],alpha=0.6,latmax=80)
    sps.drawmeridians([0,60,120,180,240,300,360], linewidth=1,labels=[True,True,False, True],fmt='%g',alpha=0.6,latmax=80)
    sps.drawcoastlines(color='k',linewidth=0.8)
    sps.fillcontinents(color='lightgray', lake_color='aqua')
    return sps,x,y

#region  SST 观测
years = np.linspace(1982,2021,40)
year_min = [1994., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2016., 2018., 2021.]
num_y = 5

#对读取的SST数据进行处理
fname = 'sst.mnmean.nc'
sst = xr.open_dataset('D:/SST/'+fname)['sst'].loc['1982-01-01':'2021-12-31']
sst_detrended = xr.apply_ufunc(
    detrend, sst, kwargs={'axis': 0},  # axis=0 表示时间维度
    input_core_dims=[['time']],
    output_core_dims=[['time']]
)
sst_transposed = sst_detrended.transpose("time",'lon','lat')

#sic1 = xr.open_dataset('/Users/wangshaoyin/data-archive/ERA5-1x1grid/1979-2021winter/'+fname)['t2m'].sel(latitude=slice(90,-90)).loc['1979-01-01':'2021-12-31'][:,::-1,:]
dtime1 = pd.date_range(start='1982-01-01', end='2021-12-01',freq='MS')
#春季平均
# temp1 = sst[(dtime1.month==9)|(dtime1.month==10)].values.reshape(40,2,180,360).mean(axis=1)
lat1 = sst.lat.values
lon1 = sst.lon.values

t2mn1 = np.zeros((num_y,12,180,360),'float')
t2mn2 = np.zeros((num_y,12,180,360),'float')
for i in range(num_y):
    t2mn1[i,:,:,:] = sst_transposed[(sst_transposed['time.year']==year_max[i])].values.reshape(12,180,360)
for i in range(num_y):
    t2mn2[i,:,:,:] = sst_transposed[(sst_transposed['time.year']==year_min[i])].values.reshape(12,180,360)

#[8:11]表示10月和11月，[8:11]表示整个春季SON
pre_clm_son = t2mn1.mean(axis=0)[8:11].mean(axis=0)
pre_std_son = t2mn1.std(axis=0)[8:11].mean(axis=0)
post_clm_son  = t2mn2.mean(axis=0)[8:11].mean(axis=0)
post_std_son  = t2mn2.std(axis=0)[8:11].mean(axis=0)

comp3 = t2mn1.mean(axis=0)  - t2mn2.mean(axis=0)
corSON3 = comp3[8:11].mean(axis=0)

pval3=np.zeros((180,360),'float')

for i in range(180):
    for j in range(360):
        pval3[i,j] = independent_ttest(pre_clm_son[i,j], post_clm_son[i,j], pre_std_son[i,j],post_std_son[i,j], num_y,num_y,0.1)


file2 =xr.open_dataset('D:/ERA5/ERA5-pressure-level-1x1-resolution-geopotential-1979-2023.nc')
gph = file2['z'].sel(level=500).loc['1979-01-01':'2023-12-31']/9.8
dtime2 = pd.date_range(start='1979-01-01', end='2023-12-01',freq='MS')
# temp2 = gph[(dtime2.month==9)|(dtime2.month==10)|(dtime2.month==11)].values.reshape(45,3,181,360).mean(axis=1)
lat2 = file2.latitude.values
lon2 = file2.longitude.values

t2mn1 = np.zeros((num_y,12,181,360),'float')
t2mn2 = np.zeros((num_y,12,181,360),'float')
# gph_clm = gph.values.reshape(30,12,181,360).mean(axis=0)*1e-2
# gph_std = gph.values.reshape(30,12,181,360).std(axis=0)*1e-2
  
for i in range(num_y):
    t2mn1[i,:,:,:] = gph[(gph['time.year']==year_max[i])].values.reshape(12,181,360)
for i in range(num_y):
    t2mn2[i,:,:,:] = gph[(gph['time.year']==year_min[i])].values.reshape(12,181,360)
pre_clm_son = t2mn1.mean(axis=0)[8:11].mean(axis=0)
pre_std_son = t2mn1.std(axis=0)[8:11].mean(axis=0)
post_clm_son  = t2mn2.mean(axis=0)[8:11].mean(axis=0)
post_std_son  = t2mn2.std(axis=0)[8:11].mean(axis=0)

comp3 = t2mn1.mean(axis=0)  - t2mn2.mean(axis=0)
corSON3_gph = comp3[8:11].mean(axis=0)


#设置全局变量
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 9
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.4


#画图
plt.close()
fig=plt.figure(figsize=(6,4))
ax1 = fig.add_subplot(111)
#m = Basemap(lon_0=180, boundinglat=-50, projection='spstere',round=True,resolution='l')
m1_ax1, x1, y1 = sh_project(ax1, lon1, lat1)
m2_ax1, x2, y2 = sh_project(ax1, lon2, lat2)
levels1=np.linspace(-1,1,21)
im1 = m1_ax1.contourf(x1,y1,corSON3,levels=levels1,extend='both',shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)
csp = m1_ax1.contourf(x1,y1,1-pval3, levels=[0.95,1] ,colors='none',hatches=['...', None],alpha=0)
levels2=np.linspace(-100,100,11)
im2 = m2_ax1.contour(x2,y2,corSON3_gph,levels=levels2,extend='both',colors=['blue'])
ax1.clabel(im2, im2.levels, inline=True)
m1_ax1.fillcontinents(color='lightgray', lake_color='aqua')
ax1.set_title('(a) Observation', position=(0.5, 1.3),pad=15)
#cax_ax1 = fig.add_axes([0.1, 0.51, 0.5,0.12], aspect=0.02)
cb1 = fig.colorbar(im1,orientation='horizontal',ticks=levels1[::2],extend='both', extendfrac=0.03,drawedges=False)
cb1.outline.set_edgecolor('black')
cb1.set_label("SST (°C)")
fig.savefig('C:/users/fzjxw/python/code/Figures/test_0109.png')
#fig.suptitle('sea level pressure anomaly (2002)', fontsize=14)
# fig.tight_layout(rect=[0.05,0.1,0.95,0.95])
