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
sst_detrended = detrend(sst,type='linear',axis=0)
sst_rsp = sst_detrended.reshape(40,12,180,360)
t2mn1 = sst_rsp[np.isin(years,year_max)]
t2mn2 = sst_rsp[np.isin(years,year_min)]
#sic1 = xr.open_dataset('/Users/wangshaoyin/data-archive/ERA5-1x1grid/1979-2021winter/'+fname)['t2m'].sel(latitude=slice(90,-90)).loc['1979-01-01':'2021-12-31'][:,::-1,:]
#dtime1 = pd.date_range(start='1982-01-01', end='2021-12-01',freq='MS')
#春季平均
# temp1 = sst[(dtime1.month==9)|(dtime1.month==10)].values.reshape(40,2,180,360).mean(axis=1)
lat1 = sst.lat.values
lon1 = sst.lon.values

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


#GPH500hPa
file2 =xr.open_dataset('D:/ERA5/ERA5-pressure-level-1x1-resolution-geopotential-1979-2023.nc')
gph500 = file2['z'].sel(level=500).loc['1979-01-01':'2023-12-31']/9.8
dtime2 = pd.date_range(start='1979-01-01', end='2023-12-01',freq='MS')
# temp2 = gph[(dtime2.month==9)|(dtime2.month==10)|(dtime2.month==11)].values.reshape(45,3,181,360).mean(axis=1)
lat2 = file2.latitude.values
lon2 = file2.longitude.values
gph500_detrended = detrend(gph500,type='linear',axis=0)
gph500_rsp = gph500_detrended.reshape(45,12,181,360)
years = np.linspace(1979,2023,45)
t2mn1 = gph500_rsp[np.isin(years,year_max)]
t2mn2 = gph500_rsp[np.isin(years,year_min)]
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
fig=plt.figure(figsize=(160/25.4, 6))
ax1 = fig.add_subplot(2,3,(1,2))
#m = Basemap(lon_0=180, boundinglat=-50, projection='spstere',round=True,resolution='l')
m1_ax1, x1, y1 = sh_project(ax1, lon1, lat1)
m2_ax1, x2, y2 = sh_project(ax1, lon2, lat2)
levels1=np.linspace(-1,1,21)
im1 = m1_ax1.contourf(x1,y1,corSON3,levels=levels1,extend='both',shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)
csp = m1_ax1.contourf(x1,y1,1-pval3, levels=[0.95,1] ,colors='none',hatches=['...', None],alpha=0)
#levels2=np.linspace(-100,100,11)
im2 = m2_ax1.contour(x2,y2,corSON3_gph,
                     #levels=levels2,
                     extend='both',colors=['blue'])
ax1.clabel(im2, im2.levels, inline=True)
m1_ax1.fillcontinents(color='lightgray', lake_color='aqua')
ax1.set_title('(a) Observation', position=(0.5, 1.3),pad=15)
cax_ax1 = fig.add_axes([0.1, 0.51, 0.5,0.12], aspect=0.02)
cb1 = fig.colorbar(im1,orientation='horizontal',ticks=levels1[::2],extend='both', extendfrac=0.03,drawedges=False,cax = cax_ax1)
cb1.outline.set_edgecolor('black')
cb1.set_label("SST (°C)")
#fig.suptitle('sea level pressure anomaly (2002)', fontsize=14)
# fig.tight_layout(rect=[0.05,0.1,0.95,0.95])


#region 模式
def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} \%" if plt.rcParams["text.usetex"] else f"{s} %"
tos = xr.open_dataset('D:/CMIP6/CESM2-WACCM-FV2/piControl/Omon/tos/tos_Omon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all_remap.nc')['tos']
lon1 = tos.lon
lat1 = tos.lat
CMIP6_opwa = np.load('D:/npz/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']*1e-12
ross_sie_cmip6 = np.nanmean(CMIP6_opwa,axis=1)
Q1 = np.quantile(ross_sie_cmip6,0.05)
Q3 = np.quantile(ross_sie_cmip6,0.95)
sic_Feb = xr.open_dataset('D:/CMIP6/CESM2-WACCM-FV2/piControl/SImon/siconc/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_Feb.nc')['siconc']
year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year
#仅需要10月和11月的数据结果 
# year_min = tos.time[ross_sie_cmip6<Q1].dt.year

# year_max = tos.time[ross_sie_cmip6<Q1].dt.year
tos_max_lst = []
for i in range(len(year_max)):
    y = str(year_max[i].values).zfill(4)
#     print(y)
    tos_max_lst.append(tos.sel(time=y))
    
tos_min_lst = []
for i in range(len(year_min)):
    y = str(year_min[i].values).zfill(4)
#     print(y)
    tos_min_lst.append(tos.sel(time=y))
    
tos_max_mean = np.array(tos_max_lst).mean(axis=0)
tos_max_std = np.array(tos_max_lst).std(axis=0)
tos_min_mean = np.array(tos_min_lst).mean(axis=0)
tos_min_std = np.array(tos_min_lst).std(axis=0)

num_y = 25

pre_clm_son = tos_max_mean[9:11].mean(axis=0)
pre_std_son = tos_max_std[9:11].mean(axis=0)

# tos_clm_son = tos_clm
# tos_std_son = tos_std

post_clm_son = tos_min_mean[9:11].mean(axis=0)
post_std_son = tos_min_std[9:11].mean(axis=0)


# pval1=np.zeros((6,332,316),'float')
# pval2=np.zeros((6,332,316),'float')
pval3=np.zeros((180,360),'float')


for i in range(96):
    for j in range(144):
#             pval1[k,i,j] = independent_ttest(pre_clm_son[k,i,j], tos_clm_son[k,i,j], pre_std_son[k,i,j],tos_std_son[k,i,j], num_y, 30, 0.1)
#             pval2[k,i,j] = independent_ttest(post_clm_son[k,i,j], tos_clm_son[k,i,j], post_std_son[k,i,j],tos_std_son[k,i,j], num_y, 30, 0.1)
        pval3[i,j] = independent_ttest(pre_clm_son[i,j], post_clm_son[i,j], pre_std_son[i,j],post_std_son[i,j], num_y, num_y, 0.1)
tos_dif = (tos_max_mean-tos_min_mean)[9:11].mean(axis=0)
zg = xr.open_dataset('D:/CMIP6/CESM2-WACCM-FV2/piControl/Amon/zg/zg_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-050012_all.nc')['zg'].sel(plev = 50000)
lon2 = zg.lon
lat2 = zg.lat
 # year_min = zg.time[ross_sie_cmip6<Q1].dt.year
zg_min_lst = []
for i in range(len(year_min)):
    y = str(year_min[i].values).zfill(4)
#     print(y)
    zg_min_lst.append(zg.sel(time=y))
zg_min_mean = np.array(zg_min_lst).squeeze().mean(axis=0)
zg_min_std = np.array(zg_min_lst).squeeze().std(axis=0)

# year_max = zg.time[ross_sie_cmip6<Q1].dt.year
zg_max_lst = []
for i in range(len(year_max)):
    y = str(year_max[i].values).zfill(4)
#     print(y)
    zg_max_lst.append(zg.sel(time=y))
zg_max_mean = np.array(zg_max_lst).squeeze().mean(axis=0)
zg_max_std = np.array(zg_max_lst).squeeze().std(axis=0)


num_y = 25

pre_clm_son = zg_max_mean[9:11].mean(axis=0)
pre_std_son = zg_max_std[9:11].mean(axis=0)

# zg_clm_son = zg_clm
# zg_std_son = zg_std

post_clm_son = zg_min_mean[9:11].mean(axis=0)
post_std_son = zg_min_std[9:11].mean(axis=0)


zg_dif = (zg_max_mean-zg_min_mean)[9:11].mean(axis=0)
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

ax3 = fig.add_subplot(2,3,(4,5))
#m = Basemap(lon_0=180, boundinglat=-50, projection='spstere',round=True,resolution='l')
m1_ax3, x1, y1 = sh_project(ax3, lon1, lat1)
m2_ax3, x2, y2 = sh_project(ax3, lon2, lat2)

levels1=np.linspace(-1,1,21)
im3 = m1_ax3.contourf(x1,y1,tos_dif,
                 levels=levels1,extend='both',
                 shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed
                )
csp = m1_ax3.contourf(x1,y1,1-pval3, levels=[0.95,1] ,colors='none',hatches=['...', None],alpha=0)

levels2=np.linspace(-100,100,11)
im4 = m2_ax3.contour(x2,y2,zg_dif,levels=levels2,extend='both',colors=['blue'])
ax3.clabel(im4, im4.levels, inline=True)

m1_ax3.fillcontinents(color='lightgray', lake_color='aqua')
# cs = m.contourf(x,y, p1, levels=[0.95, 1] ,colors='none',hatches=['///', None],alpha=0)
#cs = m.contourf(x,y, p1, levels=[0.9, 0.95] ,colors='none',hatches=['++', None],alpha=0)
ax3.set_title('(c) CESM2-WACCM-FV2', position=(0.5, 1.12),pad=15)
#m.drawcoastlines(color='b')
#ax3.clabel(cf, fontsize=16, inline=1,fmt='%1.f')
cax_ax3 = fig.add_axes([0.1, 0.03, 0.5, 0.12], aspect=0.02)
cb3 = fig.colorbar(im3,orientation='horizontal',ticks=levels1[::2],extend='both', extendfrac=0.03,drawedges=False,cax = cax_ax3)
# cb2.ax.tick_params(labelsize=14)
cb3.outline.set_edgecolor('black')
#cb.set_label('GPH500 [meter]',fontsize=16)
#cb.set_label("Sea Level Pressure [hPa]",fontsize=12)
cb3.set_label("SST (°C)")
#fig.suptitle('sea level pressure anomaly (2002)', fontsize=14)
# fig.tight_layout(rect=[0.05,0.1,0.95,0.95])
# plt.subplots_adjust(hspace=0.25)
#plt.savefig('seasonal_std_sea_level_pressure_beafort_high.png',dpi=600,bbox_inches='tight')
# plt.savefig('/stu02/weizx24/figures/0924/Figure10_all.png' ,dpi=300,bbox_inches='tight')

#endregion

#region 极涡 观测
year_min = [1994., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2016., 2018., 2021.]
num_y = 5
file2 =xr.open_dataset('D:/ERA5/ERA5-pressure-level-1x1-resolution-geopotential-1979-2023.nc')
gph50 = file2['z'].sel(level=50).loc['1979-01-01':'2023-12-31']/9.8
dtime2 = pd.date_range(start='1979-01-01', end='2023-12-01',freq='MS')
# temp2 = gph[(dtime2.month==9)|(dtime2.month==10)|(dtime2.month==11)].values.reshape(45,3,181,360).mean(axis=1)
lat2 = file2.latitude.values
lon2 = file2.longitude.values
years = np.linspace(1979,2023,45)
gph50_detrended = detrend(gph50,type='linear',axis=0)
gph50_rsp = gph50_detrended.reshape(45,12,181,360)
t2mn1 = gph50_rsp[np.isin(years,year_max)]
t2mn2 = gph50_rsp[np.isin(years,year_min)]
pre_clm_son = t2mn1.mean(axis=0)[9:11].mean(axis=0)
pre_std_son = t2mn1.std(axis=0)[9:11].mean(axis=0)
post_clm_son  = t2mn2.mean(axis=0)[9:11].mean(axis=0)
post_std_son  = t2mn2.std(axis=0)[9:11].mean(axis=0)

gph_clm = gph50_rsp.mean(axis=0)
gph_std = gph50_rsp.std(axis=0)
gph_clm_son   = gph_clm[9:11].mean(axis=0)
gph_std_son   = gph_std[9:11].mean(axis=0)
comp3 = t2mn1.mean(axis=0)  - t2mn2.mean(axis=0)
corSON3_gph = comp3[9:11].mean(axis=0)

# pval1=np.zeros((181,360),'float')
# pval2=np.zeros((181,360),'float')
pval3=np.zeros((181,360),'float')

for i in range(181):
    for j in range(360):
# #         pval1[i,j] = independent_ttest(pre_clm_son[i,j], slp_clm_son[i,j], pre_std_son[i,j],slp_std_son[i,j],num_y, 30, 0.1)
# #         pval2[i,j] = independent_ttest(post_clm_son[i,j], slp_clm_son[i,j], post_std_son[i,j],slp_std_son[i,j], num_y, 30, 0.1)
        pval3[i,j] = independent_ttest(pre_clm_son[i,j], post_clm_son[i,j], pre_std_son[i,j],post_std_son[i,j], num_y,num_y,0.1)
gph_dif = pre_clm_son-post_clm_son

def sh_project(ax,lon,lat):
    #sps = Basemap(lon_0=180, boundinglat=0, projection='spstere',round=True,resolution='l')
    sps = Basemap(projection='spstere',
                  lon_0=180, boundinglat=-40,round=True,resolution='l')
    sps.ax = ax
    mer = np.arange(0, 360, 60.)
    #par = np.arange(-90, -50, 10.)
    x, y = sps(*np.meshgrid(lon, lat))
    sps.drawmapboundary(fill_color='AntiqueWhite')

    sps.drawparallels([-90,-70,-50,], labels=[True,True,True, True],alpha=0.6,latmax=80)
    sps.drawmeridians([0,60,120,180,240,300], labels=[True,True,True, True],fmt='%g',alpha=0.6,latmax=80)
    sps.drawcoastlines(color='k')
    #sps.fillcontinents(color='',lake_color='aqua')
#     sps.fillcontinents(color='lightgray', lake_color='aqua')
    return sps,x,y

#plt.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['hatch.color'] ='springgreen'
# mpl.rcParams['hatch.linewidth'] = 0.4
# plt.close()
# fig=plt.figure(figsize=(10, 8))
ax2 = fig.add_subplot(2,3,3)
#m = Basemap(lon_0=180, boundinglat=-50, projection='spstere',round=True,resolution='l')
m1_ax2, x1, y1 = sh_project(ax2, lon2, lat2)

levels1=np.linspace(-200,200,9)
im5 = m1_ax2.contourf(x1,y1,gph_dif,
                 levels=levels1,
                 extend='both',
                 shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)
csp = m1_ax2.contourf(x1,y1,1-pval3, levels=[0.95,1] ,colors='none',hatches=['.', None],alpha=0)

#levels2=np.linspace(19250,20250,3)
im6 = m1_ax2.contour(x1,y1,gph_clm_son,
                 #levels=levels2,
                 extend='both',colors=['black'])
ax2.clabel(im6, im6.levels, inline=True)

# m.fillcontinents(color='lightgray', lake_color='aqua')
# cs = m.contourf(x,y, p1, levels=[0.95, 1] ,colors='none',hatches=['///', None],alpha=0)
#cs = m.contourf(x,y, p1, levels=[0.9, 0.95] ,colors='none',hatches=['++', None],alpha=0)
ax2.set_title('(b) Observation', position=(0.5, 1.12),pad=20)
#m.drawcoastlines(color='b')
#ax2.clabel(cf, fontsize=16, inline=1,fmt='%1.f')
cax_ax2 = fig.add_axes([0.68, 0.52, 0.3, 0.1],aspect=0.03)
cb2 = fig.colorbar(im5,orientation='horizontal',
                  ticks=levels1[::2],
                  extend='both', extendfrac=0.03,drawedges=False,
                  cax = cax_ax2
                 )
# cb.ax.tick_params(labelsize=14)
cb2.outline.set_edgecolor('black')
#cb.set_label('GPH500 [meter]',fontsize=16)
#cb.set_label("Sea Level Pressure [hPa]",fontsize=12)
cb2.set_label("GPH (m)")
#fig.suptitle('sea level pressure anomaly (2002)', fontsize=14)
# fig.tight_layout(rect=[0.05,0.1,0.95,0.95])
#plt.savefig('seasonal_std_sea_level_pressure_beafort_high.png',dpi=600,bbox_inches='tight')
# plt.savefig('/stu02/weizx24/figures/0924/Figure10_obs_polarvortex.png' ,dpi=300,bbox_inches='tight')

#region 模式
CMIP6_opwa = np.load('D:/npz/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']*1e-12
ross_sie_cmip6 = np.nanmean(CMIP6_opwa,axis=1)
Q1 = np.quantile(ross_sie_cmip6,0.05)
Q3 = np.quantile(ross_sie_cmip6,0.95)
sic_Feb = xr.open_dataset('D:/CMIP6/CESM2-WACCM-FV2/piControl/SImon/siconc/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_Feb.nc')['siconc']
year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year
zg = xr.open_dataset('D:/CMIP6/CESM2-WACCM-FV2/piControl/Amon/zg/zg_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-050012_all.nc')['zg'].sel(plev = 5000)
lon2 = zg.lon
lat2 = zg.lat
zg_min_lst = []
for i in range(len(year_min)):
    y = str(year_min[i].values).zfill(4)
#     print(y)
    zg_min_lst.append(zg.sel(time=y))
zg_min_mean = np.array(zg_min_lst).squeeze().mean(axis=0)
zg_min_std = np.array(zg_min_lst).squeeze().std(axis=0)

zg_max_lst = []
for i in range(len(year_max)):
    y = str(year_max[i].values).zfill(4)
#     print(y)
    zg_max_lst.append(zg.sel(time=y))
zg_max_mean = np.array(zg_max_lst).squeeze().mean(axis=0)
zg_max_std = np.array(zg_max_lst).squeeze().std(axis=0)
num_y = 25

#选择Oct-Nov mean
pre_clm_son = zg_max_mean[9:11].mean(axis=0)
pre_std_son = zg_max_std[9:11].mean(axis=0)

zg_clm_son = zg.values.reshape(500,12,96,144).mean(axis=0)[9:11].mean(axis=0)
zg_std_son = zg.values.reshape(500,12,96,144).std(axis=0)[9:11].mean(axis=0)

post_clm_son = zg_min_mean[9:11].mean(axis=0)
post_std_son = zg_min_std[9:11].mean(axis=0)

pval3=np.zeros((96,144),'float')


for i in range(96):
    for j in range(144):
        pval3[i,j] = independent_ttest(pre_clm_son[i,j], post_clm_son[i,j], pre_std_son[i,j],post_std_son[i,j], num_y, num_y, 0.1)

zg_dif = pre_clm_son - post_clm_son
def sh_project(ax,lon,lat):
    #sps = Basemap(lon_0=180, boundinglat=0, projection='spstere',round=True,resolution='l')
    sps = Basemap(projection='spstere',
                  lon_0=180, boundinglat=-40,round=True,resolution='l')
    sps.ax = ax
    mer = np.arange(0, 360, 60.)
    #par = np.arange(-90, -50, 10.)
    x, y = sps(*np.meshgrid(lon, lat))
    sps.drawmapboundary(fill_color='AntiqueWhite')
    sps.drawparallels([-90,-70,-50,], labels=[True,True,True, True],alpha=0.6,latmax=80)
    sps.drawmeridians([0,60,120,180,240,300], labels=[True,True,True, True],fmt='%g',alpha=0.6,latmax=80)
    sps.drawcoastlines(color='k')
    #sps.fillcontinents(color='',lake_color='aqua')
#     sps.fillcontinents(color='lightgray', lake_color='aqua')
    return sps,x,y

#plt.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['hatch.color'] ='springgreen'
# mpl.rcParams['hatch.linewidth'] = 0.4
# plt.close()
# fig=plt.figure(figsize=(10, 8))
ax4 = fig.add_subplot(2,3,6)
#m = Basemap(lon_0=180, boundinglat=-50, projection='spstere',round=True,resolution='l')
m1_ax4, x1, y1 = sh_project(ax4, lon2, lat2)
# m2, x2, y2 = sh_project(ax4, lon2, lat2)
levels1=np.linspace(-200,200,9)
im7 = m1_ax4.contourf(x1,y1,zg_dif,
                 levels=levels1,
                 extend='both',
                 shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)
# csp = m.contourf(x1,y1,1-pval3, levels=[0.95,1] ,colors='none',hatches=['.', None],alpha=0)
levels2=np.linspace(19250,20250,3)
im8 = m1_ax4.contour(x1,y1,zg_clm_son,
                 levels=levels2,
                 extend='both',colors=['black'])
ax4.clabel(im8, im8.levels, inline=True)

# m.fillcontinents(color='lightgray', lake_color='aqua')
# cs = m.contourf(x,y, p1, levels=[0.95, 1] ,colors='none',hatches=['///', None],alpha=0)
#cs = m.contourf(x,y, p1, levels=[0.9, 0.95] ,colors='none',hatches=['++', None],alpha=0)
ax4.set_title('(d) CESM2-WACCM-FV2', position=(0.5, 1.5),pad=20)
#m.drawcoastlines(color='b')
#ax4.clabel(cf, fontsize=16, inline=1,fmt='%1.f')
cax_ax4 = fig.add_axes([0.68, 0.04, 0.3, 0.1],aspect=0.03)
cb4 = fig.colorbar(im7,orientation='horizontal',
                  ticks=levels1[::2],
                  extend='both', extendfrac=0.03,drawedges=False,
                  cax = cax_ax4
                 )
# cb.ax.tick_params(labelsize=14)
cb4.outline.set_edgecolor('black')
#cb.set_label('GPH500 [meter]',fontsize=16)
#cb.set_label("Sea Level Pressure [hPa]",fontsize=12)
cb4.set_label("GPH (m)")
#fig.suptitle('sea level pressure anomaly (2002)', fontsize=14)
# fig.tight_layout(rect=[0.05,0.1,0.95,0.95])
#plt.savefig('seasonal_std_sea_level_pressure_beafort_high.png',dpi=600,bbox_inches='tight')
fig.tight_layout()
# plt.subplots_adjust(hspace=0.2, wspace=0.5)
plt.savefig('C:/Users/fzjxw/python/code/Figures/Figure10_all.png' ,dpi=300,bbox_inches='tight')
plt.show()
