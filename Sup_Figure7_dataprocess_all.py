#这个计算的是5个月份...............................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................


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
import cmocean as cm
from scipy.stats import t

# 筛选最大最小年份
CMIP6_opwa = np.load('/stu02/weizx24/data/npz/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']*1e-12
ross_sie_cmip6 = np.nanmean(CMIP6_opwa,axis=1)
Q1 = np.quantile(ross_sie_cmip6,0.05)
Q3 = np.quantile(ross_sie_cmip6,0.95)
sic_Feb = xr.open_dataset('/stu02/weizx24/data/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_Feb.nc')['siconc']
year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year

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

lw = xr.open_dataset('/stu02/weizx24/data/ERA5/ERA5_79-22_lwsw.nc')['msnlwrf'].sel(time=slice('1992','2022'))
sw = xr.open_dataset('/stu02/weizx24/data/ERA5/ERA5_79-22_lwsw.nc')['msnswrf'].sel(time=slice('1992','2022'))

lhf = xr.open_dataset('/stu02/weizx24/data/ERA5/ERA5_79-22_lhfshf.nc')['mslhf'].sel(time=slice('1992','2022'))
shf = xr.open_dataset('/stu02/weizx24/data/ERA5/ERA5_79-22_lhfshf.nc')['msshf'].sel(time=slice('1992','2022'))

Qnet = lhf+shf+lw+sw

votemp = xr.open_dataset('/stu02/weizx24/data/ORAS5/votemp_79-22_50m.nc')['votemper'].mean(dim = 'LEV').sel(time=slice('1992','2022'))

# Qnet_all = np.array(Qnet)
lats = Qnet.latitude
lons = Qnet.longitude
time_len = len(Qnet)
ross_Qnet_all = np.zeros(372,'float')
for i in range(time_len):
    mask_ross = (lats<-70.)&(lons>180.)&(lons<230.)
    Qnet_anom_ross = np.where(mask_ross,Qnet[i],np.nan)
    ross_Qnet_all[i] = np.nanmean(Qnet_anom_ross)

lats = votemp.latitude
lons = votemp.longitude
time_len = len(votemp)
ross_votemp_all = np.zeros(372,'float')
for i in range(time_len):
    mask_ross = (lats<-70.)&(lons>180.)&(lons<230.)
    votemp_anom_ross = np.where(mask_ross,votemp[i],np.nan)
    ross_votemp_all[i] = np.nanmean(votemp_anom_ross)

np.savez('/stu02/weizx24/data/npz/Figure7_obs_allmon.npz',ross_Qnet_all=ross_Qnet_all,ross_votemp_all=ross_votemp_all)



CMIP6_opwa = np.load('/stu02/weizx24/data/npz/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']*1e-12
ross_sie_cmip6 = np.nanmean(CMIP6_opwa,axis=1)
Q1 = np.quantile(ross_sie_cmip6,0.05)
Q3 = np.quantile(ross_sie_cmip6,0.95)
sic_Feb = xr.open_dataset('/stu02/weizx24/data/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_Feb.nc')['siconc']
# year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
# year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year

# surface_upward_latent_heat_flux
hfls = xr.open_dataset('/stu02/weizx24/data/hfls_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['hfls']
# surface_upward_sensible_heat_flux
hfss = xr.open_dataset('/stu02/weizx24/data/hfss_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['hfss']
# Surface Downwelling Longwave Radiation
rlds = xr.open_dataset('/stu02/weizx24/data/rlds_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['rlds']
# Surface Upwelling Longwave Radiation
rlus = xr.open_dataset('/stu02/weizx24/data/rlus_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['rlus']
# Surface Downwelling Shortwave Radiation
rsds = xr.open_dataset('/stu02/weizx24/data/rsds_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['rsds']
# Surface Upwelling Shortwave Radiation
rsus = xr.open_dataset('/stu02/weizx24/data/rsus_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['rsus']
Qnet1 = (rlds-rlus)+(rsds-rsus)+hfls+hfss
thetao = xr.open_dataset('/stu02/weizx24/data/thetao_Omon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['thetao'].sel(lev=slice(0,5000)).mean('lev')

lats = rsus.lat
lons = rsus.lon
time_len = len(Qnet1)
ross_Qnet1_all = np.zeros(time_len,'float')
for j in range(time_len):
        mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
        Qnet1_anom_ross = np.where(mask_ross,Qnet1[j],np.nan)
        ross_Qnet1_all[j] = np.nanmean(Qnet1_anom_ross)

lats = thetao.lat
lons = thetao.lon
time_len = len(thetao)
ross_thetao_all = np.zeros(time_len,'float')
for i in range(time_len):
        mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
        thetao_anom_ross = np.where(mask_ross,thetao[i],np.nan)
        ross_thetao_all[i] = np.nanmean(thetao_anom_ross)

np.savez('/stu02/weizx24/data/npz/Figure7_mdl_allmon.npz',ross_Qnet1_all=ross_Qnet1_all,ross_thetao_all=ross_thetao_all)

'''
#region 观测

#气候平均态取1979-2018年就好
lw = xr.open_dataset('/stu02/weizx24/data/ERA5/ERA5_79-22_lwsw.nc')['msnlwrf'].sel(time=slice('1992','2022'))
sw = xr.open_dataset('/stu02/weizx24/data/ERA5/ERA5_79-22_lwsw.nc')['msnswrf'].sel(time=slice('1992','2022'))

lhf = xr.open_dataset('/stu02/weizx24/data/ERA5/ERA5_79-22_lhfshf.nc')['mslhf'].sel(time=slice('1992','2022'))
shf = xr.open_dataset('/stu02/weizx24/data/ERA5/ERA5_79-22_lhfshf.nc')['msshf'].sel(time=slice('1992','2022'))

Qnet = lhf+shf+lw+sw

votemp = xr.open_dataset('/stu02/weizx24/data/ORAS5/votemp_79-22_50m.nc')['votemper'].mean(dim = 'LEV').sel(time=slice('1992','2022'))


date = ['Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr']

#计算所有年份
# lw_all = list()
# sw_all = list()
# shf_all = list()
# lhf_all = list()
votemp_all = list()
Qnet_all = list()
year = np.linspace(1992,2021,30)

for i in range(len(year)):
    starttime = str(int(year[i]))+'-12-01'
    endtime = str(int(year[i])+1)+'-04-30'
    # lw_all.append(lw.sel(time=slice(starttime,endtime)))
    # sw_all.append(sw.sel(time=slice(starttime,endtime)))
    # shf_all.append(shf.sel(time=slice(starttime,endtime)))
    # lhf_all.append(lhf.sel(time=slice(starttime,endtime)))
    votemp_all.append(votemp.sel(time=slice(starttime,endtime)))
    Qnet_all.append(Qnet.sel(time=slice(starttime,endtime)))

# Qnet
Qnet_all = np.array(Qnet_all)
lats = Qnet.latitude
lons = Qnet.longitude
# time_len = len(Qnet_all)
ross_Qnet_all = np.zeros((30,5),'float')
for j in range(len(year)):
    for i in range(5):
        mask_ross = (lats<-70.)&(lons>180.)&(lons<230.)
        Qnet_anom_ross = np.where(mask_ross,Qnet_all[j,i],np.nan)
        ross_Qnet_all[j,i] = np.nanmean(Qnet_anom_ross)

# votemp
votemp_all = np.array(votemp_all)
lats = votemp.latitude
lons = votemp.longitude
# time_len = len(votemp_all)
ross_votemp_all = np.zeros((30,5),'float')
for j in range(len(year)):
    for i in range(5):
        mask_ross = (lats<-70.)&(lons>180.)&(lons<230.)
        votemp_anom_ross = np.where(mask_ross,votemp_all[j,i],np.nan)
        ross_votemp_all[j,i] = np.nanmean(votemp_anom_ross)
#
# np.savez('/stu02/weizx24/data/npz/Figure7_obs_0924.npz',ross_lhf_dif=ross_lhf_dif,ross_lw_dif=ross_lw_dif,ross_shf_dif=ross_shf_dif,ross_sw_dif=ross_sw_dif,ross_Qnet_dif=ross_Qnet_dif,ross_votemp_dif=ross_votemp_dif)
# print('观测数据处理完毕')

#保存全部的数据
np.savez('/stu02/weizx24/data/npz/Figure7_obs_all.npz',ross_Qnet_all=ross_Qnet_all,ross_votemp_all=ross_votemp_all)
#endregion


#region 模式
# 筛选最大最小年份
CMIP6_opwa = np.load('/stu02/weizx24/data/npz/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']*1e-12
ross_sie_cmip6 = np.nanmean(CMIP6_opwa,axis=1)
Q1 = np.quantile(ross_sie_cmip6,0.05)
Q3 = np.quantile(ross_sie_cmip6,0.95)
sic_Feb = xr.open_dataset('/stu02/weizx24/data/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_Feb.nc')['siconc']
# year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
# year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year

# surface_upward_latent_heat_flux
hfls = xr.open_dataset('/stu02/weizx24/data/hfls_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['hfls']
# surface_upward_sensible_heat_flux
hfss = xr.open_dataset('/stu02/weizx24/data/hfss_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['hfss']
# Surface Downwelling Longwave Radiation
rlds = xr.open_dataset('/stu02/weizx24/data/rlds_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['rlds']
# Surface Upwelling Longwave Radiation
rlus = xr.open_dataset('/stu02/weizx24/data/rlus_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['rlus']
# Surface Downwelling Shortwave Radiation
rsds = xr.open_dataset('/stu02/weizx24/data/rsds_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['rsds']
# Surface Upwelling Shortwave Radiation
rsus = xr.open_dataset('/stu02/weizx24/data/rsus_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['rsus']

#downward-upward
Qnet1 = (rlds-rlus)+(rsds-rsus)+hfls+hfss
years = np.linspace(1,499,499)
#目标是计算所有年份的1月Qnet和2月海洋温度，看看两者之间的滞后系数
Qnet1_lst = []
for i in range(498):
    starttime = str(int(years[i])).zfill(4)+'-12-01'
    endtime = str(int(years[i])+1).zfill(4)+'-04-30'
    # y = str(int(years[i])).zfill(4)
    Qnet1_lst.append(Qnet1.sel(time=slice(starttime,endtime)))
    # print(Qnet1.sel(time=slice(starttime,endtime)))
Qnet1_arr = np.array(Qnet1_lst)

lats = rsus.lat
lons = rsus.lon
# time_len = len(Qnet1_arr)
ross_Qnet_all = np.zeros((498,5),'float')
for j in range(498):
    for i in range(5):
        mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
        Qnet_anom_ross = np.where(mask_ross,Qnet1_arr[j,i],np.nan)
        ross_Qnet_all[j,i] = np.nanmean(Qnet_anom_ross)

thetao = xr.open_dataset('/stu02/weizx24/data/thetao_Omon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['thetao'].sel(lev=slice(0,5000)).mean('lev')
thetao_lst = []
for i in range(498):
    starttime = str(int(years[i])).zfill(4)+'-12-01'
    endtime = str(int(years[i])+1).zfill(4)+'-04-30'
    # y = str(int(years[i])).zfill(4)
    thetao_lst.append(thetao.sel(time=slice(starttime,endtime)))
    # print(Qnet1.sel(time=slice(starttime,endtime)))
thetao_arr = np.array(thetao_lst)

lats = thetao.lat
lons = thetao.lon
# time_len = len(thetao_arr)
ross_thetao_all = np.zeros((498,5),'float')
for j in range(498):
    for i in range(5):
        mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
        thetao_anom_ross = np.where(mask_ross,thetao_arr[j,i],np.nan)
        ross_thetao_all[j,i] = np.nanmean(thetao_anom_ross)

# np.savez('/stu02/weizx24/data/npz/Figure7_model_0924.npz',ross_sw_dif=ross_sw_dif,ross_lw_dif=ross_lw_dif,ross_hfls_dif=ross_hfls_dif,ross_hfss_dif=ross_hfss_dif,ross_Qnet_dif=ross_Qnet_dif,ross_thetao_dif=ross_thetao_dif)
# print('观测数据处理完毕')
np.savez('/stu02/weizx24/data/npz/Figure7_mdl_all.npz',ross_Qnet_all=ross_Qnet_all,ross_thetao_all=ross_thetao_all)

#endregion
'''