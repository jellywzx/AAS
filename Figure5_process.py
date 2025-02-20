'''sea ice concentration longterm evolution
first start with CESM2
'''
import numpy as np
import matplotlib.pyplot as plt
import glob
from pylab import *
from scipy import stats
import xarray as xr
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import pandas as pd
import scipy.signal as signal
import cmaps
import matplotlib.dates as mdates
import matplotlib.path as mpath
from matplotlib import dates, ticker
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.gridspec as gridspec
from scipy.stats import t
# from pyproj import Proj, transform
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import AxesGrid
'''
#region SST
tos = xr.open_dataset('/stu02/weizx24/data/tos_Omon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['tos']
# tos
lons = tos.lon
lats = tos.lat
CMIP6_opwa = np.load('/stu02/weizx24/data/npz/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']*1e-12
ross_sie_cmip6 = np.nanmean(CMIP6_opwa,axis=1)
Q1 = np.quantile(ross_sie_cmip6,0.05)
Q3 = np.quantile(ross_sie_cmip6,0.95)
sic_Feb = xr.open_dataset('/stu02/weizx24/data/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_Feb.nc')['siconc']
year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year
# year_min = tos.time[ross_sie_cmip6<Q1].dt.year
tos_min_lst = []
for i in range(len(year_min)):
    y = str(year_min[i].values).zfill(4)
#     print(y)
    tos_min_lst.append(tos.sel(time=y))
tos_min_mean = np.array(tos_min_lst).squeeze().mean(axis=0)
tos_min_std = np.array(tos_min_lst).squeeze().std(axis=0)
# year_max = tos.time[ross_sie_cmip6<Q1].dt.year
tos_max_lst = []
for i in range(len(year_max)):
    y = str(year_max[i].values).zfill(4)
#     print(y)
    tos_max_lst.append(tos.sel(time=y))
tos_max_mean = np.array(tos_max_lst).squeeze().mean(axis=0)
tos_max_std = np.array(tos_max_lst).squeeze().std(axis=0)
year_min_1 = sic_Feb.time[ross_sie_cmip6<Q1].dt.year+1
year_max_1 = sic_Feb.time[ross_sie_cmip6>Q3].dt.year+1
# year_min_1 = tos.time[ross_sie_cmip6<Q1].dt.year
tos_min_1_lst = []
for i in range(len(year_min_1)):
    y = str(year_min_1[i].values).zfill(4)
#     print(y)
    tos_min_1_lst.append(tos.sel(time=y))
tos_min_1_mean = np.array(tos_min_1_lst).squeeze().mean(axis=0)
tos_min_1_std = np.array(tos_min_1_lst).squeeze().std(axis=0)
# year_max_1 = tos.time[ross_sie_cmip6<Q1].dt.year
tos_max_1_lst = []
for i in range(len(year_max_1)):
    y = str(year_max_1[i].values).zfill(4)
#     print(y)
    tos_max_1_lst.append(tos.sel(time=y))
tos_max_1_mean = np.array(tos_max_1_lst).squeeze().mean(axis=0)
tos_max_1_std = np.array(tos_max_1_lst).squeeze().std(axis=0)
tos_min_mean_all = np.concatenate([tos_min_mean,tos_min_1_mean],axis=0)[9:15]
tos_max_mean_all = np.concatenate([tos_max_mean,tos_max_1_mean],axis=0)[9:15]
tos_min_std_all = np.concatenate([tos_min_std, tos_min_1_std],axis=0)[9:15]
tos_max_std_all = np.concatenate([tos_max_std, tos_max_1_std],axis=0)[9:15]

tos_dif = tos_max_mean_all - tos_min_mean_all
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

def sh_project(ax,lon,lat,i):
    sps = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-70)
    sps.ax = ax
    mer = np.arange(-180, 30, 30.)
    par = np.arange(-90, -50, 10.)
    x,y = sps(lon,lat)
#     x, y = sps(*np.meshgrid(lon, lat))
    if i==0:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[True,True,True,True])
    if (i==1)|(i==2):
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,True,True])
    if i==3:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,True,True])
    if i==4:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,True,True])
    if (i==5):
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,True,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,True,True])
#     if i==7:
#         sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[True,True,True,True])
#         sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,False, True])
#     sps.drawmapboundary(fill_color='AntiqueWhite')
#     draw_latlon_polygon(sps, [170,170, 298, 298], [-80, -60, -60, -80], 'k-')
    sps.fillcontinents(color='gray',lake_color='aqua')
#     
    return sps,x,y
num_y = 25

pre_clm_son = tos_max_mean_all
pre_std_son = tos_max_std_all

# tos_clm_son = tos_clm
# tos_std_son = tos_std

post_clm_son = tos_min_mean_all
post_std_son = tos_min_std_all


# pval1=np.zeros((6,332,316),'float')
# pval2=np.zeros((6,332,316),'float')
pval3=np.zeros((6,384,320),'float')

for k in range(6):
    for i in range(384):
        for j in range(320):
#             pval1[k,i,j] = independent_ttest(pre_clm_son[k,i,j], tos_clm_son[k,i,j], pre_std_son[k,i,j],tos_std_son[k,i,j], num_y, 30, 0.1)
#             pval2[k,i,j] = independent_ttest(post_clm_son[k,i,j], tos_clm_son[k,i,j], post_std_son[k,i,j],tos_std_son[k,i,j], num_y, 30, 0.1)
            pval3[k,i,j] = independent_ttest(pre_clm_son[k,i,j], post_clm_son[k,i,j], pre_std_son[k,i,j],post_std_son[k,i,j], num_y, num_y, 0.1)

#endregion
np.savez('/stu02/weizx24/data/npz/Figure5_SST.npz',tos_dif=tos_dif,pval3=pval3,lons=lons,lats=lats,)

#region ocean temperature

thetao = xr.open_dataset('/stu02/weizx24/data/thetao_Omon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['thetao'].sel(lev=slice(0,5000)).mean('lev')
# .sel(lev=slice(0,30000)).sel(nlat=slice(-80,-50))

# year_min = thetao.time[ross_sie_cmip6<Q1].dt.year
thetao_min_lst = []
for i in range(len(year_min)):
    y = str(year_min[i].values).zfill(4)
#     print(y)
    thetao_min_lst.append(thetao.sel(time=y))
thetao_min_mean = np.array(thetao_min_lst).squeeze().mean(axis=0)
thetao_min_std = np.array(thetao_min_lst).squeeze().std(axis=0)

# year_max = thetao.time[ross_sie_cmip6<Q1].dt.year
thetao_max_lst = []
for i in range(len(year_max)):
    y = str(year_max[i].values).zfill(4)
#     print(y)
    thetao_max_lst.append(thetao.sel(time=y))
thetao_max_mean = np.array(thetao_max_lst).squeeze().mean(axis=0)
thetao_max_std = np.array(thetao_max_lst).squeeze().std(axis=0)

year_min_1 = sic_Feb.time[ross_sie_cmip6<Q1].dt.year+1
year_max_1 = sic_Feb.time[ross_sie_cmip6>Q3].dt.year+1

# year_min_1 = thetao.time[ross_sie_cmip6<Q1].dt.year
thetao_min_1_lst = []
for i in range(len(year_min_1)):
    y = str(year_min_1[i].values).zfill(4)
#     print(y)
    thetao_min_1_lst.append(thetao.sel(time=y))
thetao_min_1_mean = np.array(thetao_min_1_lst).squeeze().mean(axis=0)
thetao_min_1_std = np.array(thetao_min_1_lst).squeeze().std(axis=0)

# year_max_1 = thetao.time[ross_sie_cmip6<Q1].dt.year
thetao_max_1_lst = []
for i in range(len(year_max_1)):
    y = str(year_max_1[i].values).zfill(4)
#     print(y)
    thetao_max_1_lst.append(thetao.sel(time=y))
thetao_max_1_mean = np.array(thetao_max_1_lst).squeeze().mean(axis=0)
thetao_max_1_std = np.array(thetao_max_1_lst).squeeze().std(axis=0)

thetao_min_mean_all = np.concatenate([thetao_min_mean,thetao_min_1_mean],axis=0)[9:15]
thetao_max_mean_all = np.concatenate([thetao_max_mean,thetao_max_1_mean],axis=0)[9:15]
thetao_min_std_all = np.concatenate([thetao_min_std, thetao_min_1_std],axis=0)[9:15]
thetao_max_std_all = np.concatenate([thetao_max_std, thetao_max_1_std],axis=0)[9:15]

thetao_dif = thetao_max_mean_all - thetao_min_mean_all

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

def sh_project(ax,lon,lat,i):
    sps = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-70)
    sps.ax = ax
    mer = np.arange(-180, 30, 30.)
    par = np.arange(-90, -50, 10.)
    x,y = sps(lon,lat)
#     x, y = sps(*np.meshgrid(lon, lat))
    if i==0:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[True,True,True,True])
    if (i==1)|(i==2):
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,True,True])
    if i==3:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,True,True])
    if i==4:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,True,True])
    if (i==5):
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,True,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,True,True])
#     if i==7:
#         sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[True,True,True,True])
#         sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,False, True])
#     sps.drawmapboundary(fill_color='AntiqueWhite')
#     draw_latlon_polygon(sps, [170,170, 298, 298], [-80, -60, -60, -80], 'k-')
    sps.fillcontinents(color='gray',lake_color='aqua')
#     
    return sps,x,y
num_y = 25

pre_clm_son = thetao_max_mean_all
pre_std_son = thetao_max_std_all

# thetao_clm_son = thetao_clm
# thetao_std_son = thetao_std

post_clm_son = thetao_min_mean_all
post_std_son = thetao_min_std_all


# pval1=np.zeros((6,332,316),'float')
# pval2=np.zeros((6,332,316),'float')
pval3=np.zeros((6,384,320),'float')

for k in range(6):
    for i in range(384):
        for j in range(320):
#             pval1[k,i,j] = independent_ttest(pre_clm_son[k,i,j], thetao_clm_son[k,i,j], pre_std_son[k,i,j],thetao_std_son[k,i,j], num_y, 30, 0.1)
#             pval2[k,i,j] = independent_ttest(post_clm_son[k,i,j], thetao_clm_son[k,i,j], post_std_son[k,i,j],thetao_std_son[k,i,j], num_y, 30, 0.1)
            pval3[k,i,j] = independent_ttest(pre_clm_son[k,i,j], post_clm_son[k,i,j], pre_std_son[k,i,j],post_std_son[k,i,j], num_y, num_y, 0.1)

lons = thetao.lon
lats = thetao.lat

np.savez('/stu02/weizx24/data/npz/Figure5_UOT.npz',thetao_dif=thetao_dif,lons=lons,lats=lats,pval3=pval3)
#endregion
'''


#region SIC
CMIP6_opwa = np.load('/stu02/weizx24/data/npz/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']*1e-12
ross_sie_cmip6 = np.nanmean(CMIP6_opwa,axis=1)
Q1 = np.quantile(ross_sie_cmip6,0.05)
Q3 = np.quantile(ross_sie_cmip6,0.95)
sic_Feb = xr.open_dataset('/stu02/weizx24/data/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_Feb.nc')['siconc']
year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year

sic = xr.open_dataset('/stu02/weizx24/data/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912.nc')['siconc']
# sic
lons = sic.lon
lats = sic.lat
sic_clm = sic.values.reshape(499,12,384,320).mean(axis=0)
sic_clm_sel = np.concatenate([sic_clm[9:12],sic_clm[0:3]])
sic_Feb = xr.open_dataset('/stu02/weizx24/data/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_Feb.nc')['siconc']
# year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
# year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year
# year_min = sic.time[ross_sie_cmip6<Q1].dt.year
sic_min_lst = []
for i in range(len(year_min)):
    y = str(year_min[i].values).zfill(4)
#     print(y)
    sic_min_lst.append(sic.sel(time=y))
sic_min_mean = np.array(sic_min_lst).squeeze().mean(axis=0)
sic_min_std = np.array(sic_min_lst).squeeze().std(axis=0)
# year_max = sic.time[ross_sie_cmip6<Q1].dt.year
sic_max_lst = []
for i in range(len(year_max)):
    y = str(year_max[i].values).zfill(4)
#     print(y)
    sic_max_lst.append(sic.sel(time=y))
sic_max_mean = np.array(sic_max_lst).squeeze().mean(axis=0)
sic_max_std = np.array(sic_max_lst).squeeze().std(axis=0)
year_min_1 = sic_Feb.time[ross_sie_cmip6<Q1].dt.year+1
year_max_1 = sic_Feb.time[ross_sie_cmip6>Q3].dt.year+1
# year_min_1 = sic.time[ross_sie_cmip6<Q1].dt.year
sic_min_1_lst = []
for i in range(len(year_min_1)):
    y = str(year_min_1[i].values).zfill(4)
#     print(y)
    sic_min_1_lst.append(sic.sel(time=y))
sic_min_1_mean = np.array(sic_min_1_lst).squeeze().mean(axis=0)
sic_min_1_std = np.array(sic_min_1_lst).squeeze().std(axis=0)
# year_max_1 = sic.time[ross_sie_cmip6<Q1].dt.year
sic_max_1_lst = []
for i in range(len(year_max_1)):
    y = str(year_max_1[i].values).zfill(4)
#     print(y)
    sic_max_1_lst.append(sic.sel(time=y))
sic_max_1_mean = np.array(sic_max_1_lst).squeeze().mean(axis=0)
sic_max_1_std = np.array(sic_max_1_lst).squeeze().std(axis=0)
sic_min_mean_all = np.concatenate([sic_min_mean,sic_min_1_mean],axis=0)[9:15]
sic_max_mean_all = np.concatenate([sic_max_mean,sic_max_1_mean],axis=0)[9:15]
sic_min_std_all = np.concatenate([sic_min_std, sic_min_1_std],axis=0)[9:15]
sic_max_std_all = np.concatenate([sic_max_std, sic_max_1_std],axis=0)[9:15]

sic_dif = sic_max_mean_all - sic_min_mean_all
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
num_y = 25

pre_clm_son = sic_max_mean_all
pre_std_son = sic_max_std_all

# sic_clm_son = sic_clm
# sic_std_son = sic_std

post_clm_son = sic_min_mean_all
post_std_son = sic_min_std_all


# pval1=np.zeros((6,332,316),'float')
# pval2=np.zeros((6,332,316),'float')
pval3=np.zeros((6,384,320),'float')

for k in range(6):
    for i in range(384):
        for j in range(320):
#             pval1[k,i,j] = independent_ttest(pre_clm_son[k,i,j], sic_clm_son[k,i,j], pre_std_son[k,i,j],sic_std_son[k,i,j], num_y, 30, 0.1)
#             pval2[k,i,j] = independent_ttest(post_clm_son[k,i,j], sic_clm_son[k,i,j], post_std_son[k,i,j],sic_std_son[k,i,j], num_y, 30, 0.1)
            pval3[k,i,j] = independent_ttest(pre_clm_son[k,i,j], post_clm_son[k,i,j], pre_std_son[k,i,j],post_std_son[k,i,j], num_y, num_y, 0.1)

def sh_project(ax,lon,lat,i):
    sps = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-70)
    sps.ax = ax
    mer = np.arange(-180, 30, 30.)
    par = np.arange(-90, -50, 10.)
    x,y = sps(lon,lat)
#     x, y = sps(*np.meshgrid(lon, lat))
    if i==0:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[True,True,True,True])
    if (i==1)|(i==2):
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,True,True])
    if i==3:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,True,True])
    if i==4:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,True,True])
    if (i==5):
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,True,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,True,True])
#     if i==7:
#         sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[True,True,True,True])
#         sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,False, True])
#     sps.drawmapboundary(fill_color='AntiqueWhite')
#     draw_latlon_polygon(sps, [170,170, 298, 298], [-80, -60, -60, -80], 'k-')
    sps.fillcontinents(color='gray',lake_color='aqua')
    return sps,x,y

#endregion
np.savez('/stu02/weizx24/data/Figure5_SIC.npz',sic_dif=sic_dif,lons=lons,lats=lats,pval3=pval3,sic_clm_sel=sic_clm_sel)