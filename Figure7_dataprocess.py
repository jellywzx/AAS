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


#region 观测

#气候平均态取1979-2018年就好
lw = xr.open_dataset('/stu02/weizx24/data/ERA5/ERA5_79-22_lwsw.nc')['msnlwrf'].sel(time=slice('1992','2022'))
sw = xr.open_dataset('/stu02/weizx24/data/ERA5/ERA5_79-22_lwsw.nc')['msnswrf'].sel(time=slice('1992','2022'))

lw_clm = lw.sel(time=slice('1992','2021')).values.reshape(30,12,161,1440).mean(axis=0)
lw_std = lw.sel(time=slice('1992','2021')).values.reshape(30,12,161,1440).std(axis=0)
sw_clm = sw.sel(time=slice('1992','2021')).values.reshape(30,12,161,1440).mean(axis=0)
sw_std = sw.sel(time=slice('1992','2021')).values.reshape(30,12,161,1440).std(axis=0)

lw_clm_sel = np.concatenate([lw_clm[8:12],lw_clm[0:4]])
lw_std_sel = np.concatenate([lw_std[8:12],lw_std[0:4]])
sw_clm_sel = np.concatenate([sw_clm[8:12],sw_clm[0:4]])
sw_std_sel = np.concatenate([sw_std[8:12],sw_std[0:4]])


lhf = xr.open_dataset('/stu02/weizx24/data/ERA5/ERA5_79-22_lhfshf.nc')['mslhf'].sel(time=slice('1992','2022'))
shf = xr.open_dataset('/stu02/weizx24/data/ERA5/ERA5_79-22_lhfshf.nc')['msshf'].sel(time=slice('1992','2022'))

shf_clm = shf.sel(time=slice('1992','2021')).values.reshape(30,12,161,1440).mean(axis=0)
shf_std = shf.sel(time=slice('1992','2021')).values.reshape(30,12,161,1440).std(axis=0)
lhf_clm = lhf.sel(time=slice('1992','2021')).values.reshape(30,12,161,1440).mean(axis=0)
lhf_std = lhf.sel(time=slice('1992','2021')).values.reshape(30,12,161,1440).std(axis=0)

shf_clm_sel = np.concatenate([shf_clm[8:12],shf_clm[0:4]])
shf_std_sel = np.concatenate([shf_std[8:12],shf_std[0:4]])
lhf_clm_sel = np.concatenate([lhf_clm[8:12],lhf_clm[0:4]])
lhf_std_sel = np.concatenate([lhf_std[8:12],lhf_std[0:4]])

Qnet = lhf+shf+lw+sw
Qnet_clm = Qnet.sel(time=slice('1992','2021')).values.reshape(30,12,161,1440).mean(axis=0)
Qnet_std = Qnet.sel(time=slice('1992','2021')).values.reshape(30,12,161,1440).std(axis=0)
Qnet_clm_sel = np.concatenate([Qnet_clm[8:12],Qnet_clm[0:4]])
Qnet_std_sel = np.concatenate([Qnet_std[8:12],Qnet_std[0:4]])

votemp = xr.open_dataset('/stu02/weizx24/data/ORAS5/votemp_79-22_50m.nc')['votemper'].mean(dim = 'LEV').sel(time=slice('1992','2022'))
votemp_clm = votemp.sel(time=slice('1992','2021')).values.reshape((30,12,40,360)).mean(axis=0)
votemp_std = votemp.sel(time=slice('1992','2021')).values.reshape((30,12,40,360)).std(axis=0)
votemp_clm_sel = np.concatenate([votemp_clm[8:12],votemp_clm[0:4]])
votemp_std_sel = np.concatenate([votemp_std[8:12],votemp_std[0:4]])


# votemp_anom_f = [votemp_anom[1],votemp_anom[3],votemp_anom[5],votemp_anom[7]]
# votemp_std_f = [votemp_std_sel[1],votemp_std_sel[3],votemp_std_sel[5],votemp_std_sel[7]]

date = ['Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr']

#--------试一下5年的------------
year_min = [1994., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2016., 2018., 2021.]
num_y = 5

# year_min = [2002., 2003., 2007., 2014.]
# year_max = [2004.,2005.,2016.,2021.]

#最大年份
lw_sel_max = list()
sw_sel_max = list()
shf_sel_max = list()
lhf_sel_max = list()
votemp_sel_max = list()
Qnet_sel_max = list()

for i in range(num_y):
    starttime = str(int(year_max[i]))+'-09-01'
    endtime = str(int(year_max[i])+1)+'-04-30'

    lw_sel_max.append(lw.sel(time=slice(starttime,endtime)))
    sw_sel_max.append(sw.sel(time=slice(starttime,endtime)))
    shf_sel_max.append(shf.sel(time=slice(starttime,endtime)))
    lhf_sel_max.append(lhf.sel(time=slice(starttime,endtime)))
    votemp_sel_max.append(votemp.sel(time=slice(starttime,endtime)))
    Qnet_sel_max.append(Qnet.sel(time=slice(starttime,endtime)))

#最小年份
lw_sel_min = list()
sw_sel_min = list()
shf_sel_min = list()
lhf_sel_min = list()
votemp_sel_min = list()
Qnet_sel_min = list()

for i in range(num_y):
    starttime = str(int(year_min[i]))+'-09-01'
    endtime = str(int(year_min[i])+1)+'-04-30'

    lw_sel_min.append(lw.sel(time=slice(starttime,endtime)))
    sw_sel_min.append(sw.sel(time=slice(starttime,endtime)))
    shf_sel_min.append(shf.sel(time=slice(starttime,endtime)))
    lhf_sel_min.append(lhf.sel(time=slice(starttime,endtime)))
    votemp_sel_min.append(votemp.sel(time=slice(starttime,endtime)))
    Qnet_sel_min.append(Qnet.sel(time=slice(starttime,endtime)))

lw_dif = np.nanmean(np.array(lw_sel_max)-np.array(lw_sel_min),axis=0)
sw_dif = np.nanmean(np.array(sw_sel_max)-np.array(sw_sel_min),axis=0)
shf_dif = np.nanmean(np.array(shf_sel_max)-np.array(shf_sel_min),axis=0)
lhf_dif = np.nanmean(np.array(lhf_sel_max)-np.array(lhf_sel_min),axis=0)
votemp_dif = np.nanmean(np.array(votemp_sel_max)-np.array(votemp_sel_min),axis=0)
Qnet_dif = np.nanmean(np.array(Qnet_sel_max)-np.array(Qnet_sel_min),axis=0)


# lw_dif = np.array(lw_sel).mean(axis=0)-lw_clm_sel
# sw_dif = np.array(sw_sel).mean(axis=0)-sw_clm_sel
# shf_dif =np.array(shf_sel).mean(axis=0)-shf_clm_sel
# lhf_dif = np.array(lhf_sel).mean(axis=0)-lhf_clm_sel
# Qnet_dif = np.array(Qnet_sel).mean(axis=0)-Qnet_clm_sel
# votemp_dif = np.array(votemp_sel).mean(axis=0)-votemp_clm_sel

# ross_mask =  (lats<-70.)&(lons>180.)&(lons<230.)
# all_mask = (lats<-70.)&(lons>160.)&(lons<300.)

# lw
lats = lw.latitude
lons = lw.longitude
time_len = len(lw_dif)
ross_lw_dif = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
    lw_dif_ross = np.where(mask_ross,lw_dif[i],np.nan)
    ross_lw_dif[i] = np.nanmean(lw_dif_ross)
    
# sw
lats = sw.latitude
lons = sw.longitude
time_len = len(sw_dif)
ross_sw_dif = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
    sw_anom_ross = np.where(mask_ross,sw_dif[i],np.nan)
    ross_sw_dif[i] = np.nanmean(sw_anom_ross)

# shf
lats = shf.latitude
lons = shf.longitude
time_len = len(shf_dif)
ross_shf_dif = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross = (lats<-70.)&(lons>180.)&(lons<230.)
    shf_anom_ross = np.where(mask_ross,shf_dif[i],np.nan)
    ross_shf_dif[i] = np.nanmean(shf_anom_ross)

# lhf
lats = lhf.latitude
lons = lhf.longitude
time_len = len(lhf_dif)
ross_lhf_dif = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
    lhf_anom_ross = np.where(mask_ross,lhf_dif[i],np.nan)
    ross_lhf_dif[i] = np.nanmean(lhf_anom_ross)

# Qnet
lats = Qnet.latitude
lons = Qnet.longitude
time_len = len(Qnet_dif)
ross_Qnet_dif = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross = (lats<-70.)&(lons>180.)&(lons<230.)
    Qnet_anom_ross = np.where(mask_ross,Qnet_dif[i],np.nan)
    ross_Qnet_dif[i] = np.nanmean(Qnet_anom_ross)

# votemp
lats = votemp.latitude
lons = votemp.longitude
time_len = len(votemp_dif)
ross_votemp_dif = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross = (lats<-70.)&(lons>180.)&(lons<230.)
    votemp_anom_ross = np.where(mask_ross,votemp_dif[i],np.nan)
    ross_votemp_dif[i] = np.nanmean(votemp_anom_ross)

np.savez('/stu02/weizx24/data/npz/Figure7_obs_0924.npz',ross_lhf_dif=ross_lhf_dif,ross_lw_dif=ross_lw_dif,ross_shf_dif=ross_shf_dif,ross_sw_dif=ross_sw_dif,ross_Qnet_dif=ross_Qnet_dif,ross_votemp_dif=ross_votemp_dif)
print('观测数据处理完毕')
#endregion

#region 模式
# 筛选最大最小年份
CMIP6_opwa = np.load('/stu02/weizx24/data/npz/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']*1e-12
ross_sie_cmip6 = np.nanmean(CMIP6_opwa,axis=1)
Q1 = np.quantile(ross_sie_cmip6,0.05)
Q3 = np.quantile(ross_sie_cmip6,0.95)
sic_Feb = xr.open_dataset('/stu02/weizx24/data/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_Feb.nc')['siconc']
year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year

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

year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year

# year_min = Qnet1.time[ross_sie_cmip6<Q1].dt.year
Qnet1_min_lst = []
for i in range(len(year_min)):
    y = str(year_min[i].values).zfill(4)
#     print(y)
    Qnet1_min_lst.append(Qnet1.sel(time=y))
Qnet1_min_mean = np.array(Qnet1_min_lst).squeeze().mean(axis=0)
Qnet1_min_std = np.array(Qnet1_min_lst).squeeze().std(axis=0)

# year_max = Qnet1.time[ross_sie_cmip6<Q1].dt.year
Qnet1_max_lst = []
for i in range(len(year_max)):
    y = str(year_max[i].values).zfill(4)
#     print(y)
    Qnet1_max_lst.append(Qnet1.sel(time=y))
Qnet1_max_mean = np.array(Qnet1_max_lst).squeeze().mean(axis=0)
Qnet1_max_std = np.array(Qnet1_max_lst).squeeze().std(axis=0)

year_min_1 = sic_Feb.time[ross_sie_cmip6<Q1].dt.year+1
year_max_1 = sic_Feb.time[ross_sie_cmip6>Q3].dt.year+1

# year_min_1 = Qnet1.time[ross_sie_cmip6<Q1].dt.year
Qnet1_min_1_lst = []
for i in range(len(year_min_1)):
    y = str(year_min_1[i].values).zfill(4)
#     print(y)
    Qnet1_min_1_lst.append(Qnet1.sel(time=y))
Qnet1_min_1_mean = np.array(Qnet1_min_1_lst).squeeze().mean(axis=0)
Qnet1_min_1_std = np.array(Qnet1_min_1_lst).squeeze().std(axis=0)

# year_max_1 = Qnet1.time[ross_sie_cmip6<Q1].dt.year
Qnet1_max_1_lst = []
for i in range(len(year_max_1)):
    y = str(year_max_1[i].values).zfill(4)
#     print(y)
    Qnet1_max_1_lst.append(Qnet1.sel(time=y))
Qnet1_max_1_mean = np.array(Qnet1_max_1_lst).squeeze().mean(axis=0)
Qnet1_max_1_std = np.array(Qnet1_max_1_lst).squeeze().std(axis=0)

Qnet1_min_mean_all = np.concatenate([Qnet1_min_mean,Qnet1_min_1_mean],axis=0)[8:16]
Qnet1_max_mean_all = np.concatenate([Qnet1_max_mean,Qnet1_max_1_mean],axis=0)[8:16]
Qnet1_min_std_all = np.concatenate([Qnet1_min_std, Qnet1_min_1_std],axis=0)[8:16]
Qnet1_max_std_all = np.concatenate([Qnet1_max_std, Qnet1_max_1_std],axis=0)[8:16]

Qnet1_dif = Qnet1_max_mean_all - Qnet1_min_mean_all

year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year

# year_min = rsus.time[ross_sie_cmip6<Q1].dt.year
rsus_min_lst = []
for i in range(len(year_min)):
    y = str(year_min[i].values).zfill(4)
#     print(y)
    rsus_min_lst.append(rsus.sel(time=y))
rsus_min_mean = np.array(rsus_min_lst).squeeze().mean(axis=0)
rsus_min_std = np.array(rsus_min_lst).squeeze().std(axis=0)

# year_max = rsus.time[ross_sie_cmip6<Q1].dt.year
rsus_max_lst = []
for i in range(len(year_max)):
    y = str(year_max[i].values).zfill(4)
#     print(y)
    rsus_max_lst.append(rsus.sel(time=y))
rsus_max_mean = np.array(rsus_max_lst).squeeze().mean(axis=0)
rsus_max_std = np.array(rsus_max_lst).squeeze().std(axis=0)

year_min_1 = sic_Feb.time[ross_sie_cmip6<Q1].dt.year+1
year_max_1 = sic_Feb.time[ross_sie_cmip6>Q3].dt.year+1

# year_min_1 = rsus.time[ross_sie_cmip6<Q1].dt.year
rsus_min_1_lst = []
for i in range(len(year_min_1)):
    y = str(year_min_1[i].values).zfill(4)
#     print(y)
    rsus_min_1_lst.append(rsus.sel(time=y))
rsus_min_1_mean = np.array(rsus_min_1_lst).squeeze().mean(axis=0)
rsus_min_1_std = np.array(rsus_min_1_lst).squeeze().std(axis=0)

# year_max_1 = rsus.time[ross_sie_cmip6<Q1].dt.year
rsus_max_1_lst = []
for i in range(len(year_max_1)):
    y = str(year_max_1[i].values).zfill(4)
#     print(y)
    rsus_max_1_lst.append(rsus.sel(time=y))
rsus_max_1_mean = np.array(rsus_max_1_lst).squeeze().mean(axis=0)
rsus_max_1_std = np.array(rsus_max_1_lst).squeeze().std(axis=0)

rsus_min_mean_all = np.concatenate([rsus_min_mean,rsus_min_1_mean],axis=0)[8:16]
rsus_max_mean_all = np.concatenate([rsus_max_mean,rsus_max_1_mean],axis=0)[8:16]
rsus_min_std_all = np.concatenate([rsus_min_std, rsus_min_1_std],axis=0)[8:16]
rsus_max_std_all = np.concatenate([rsus_max_std, rsus_max_1_std],axis=0)[8:16]

rsus_dif = rsus_max_mean_all - rsus_min_mean_all

year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year

# year_min = rsds.time[ross_sie_cmip6<Q1].dt.year
rsds_min_lst = []
for i in range(len(year_min)):
    y = str(year_min[i].values).zfill(4)
#     print(y)
    rsds_min_lst.append(rsds.sel(time=y))
rsds_min_mean = np.array(rsds_min_lst).squeeze().mean(axis=0)
rsds_min_std = np.array(rsds_min_lst).squeeze().std(axis=0)

# year_max = rsds.time[ross_sie_cmip6<Q1].dt.year
rsds_max_lst = []
for i in range(len(year_max)):
    y = str(year_max[i].values).zfill(4)
#     print(y)
    rsds_max_lst.append(rsds.sel(time=y))
rsds_max_mean = np.array(rsds_max_lst).squeeze().mean(axis=0)
rsds_max_std = np.array(rsds_max_lst).squeeze().std(axis=0)

year_min_1 = sic_Feb.time[ross_sie_cmip6<Q1].dt.year+1
year_max_1 = sic_Feb.time[ross_sie_cmip6>Q3].dt.year+1

# year_min_1 = rsds.time[ross_sie_cmip6<Q1].dt.year
rsds_min_1_lst = []
for i in range(len(year_min_1)):
    y = str(year_min_1[i].values).zfill(4)
#     print(y)
    rsds_min_1_lst.append(rsds.sel(time=y))
rsds_min_1_mean = np.array(rsds_min_1_lst).squeeze().mean(axis=0)
rsds_min_1_std = np.array(rsds_min_1_lst).squeeze().std(axis=0)

# year_max_1 = rsds.time[ross_sie_cmip6<Q1].dt.year
rsds_max_1_lst = []
for i in range(len(year_max_1)):
    y = str(year_max_1[i].values).zfill(4)
#     print(y)
    rsds_max_1_lst.append(rsds.sel(time=y))
rsds_max_1_mean = np.array(rsds_max_1_lst).squeeze().mean(axis=0)
rsds_max_1_std = np.array(rsds_max_1_lst).squeeze().std(axis=0)

rsds_min_mean_all = np.concatenate([rsds_min_mean,rsds_min_1_mean],axis=0)[8:16]
rsds_max_mean_all = np.concatenate([rsds_max_mean,rsds_max_1_mean],axis=0)[8:16]
rsds_min_std_all = np.concatenate([rsds_min_std, rsds_min_1_std],axis=0)[8:16]
rsds_max_std_all = np.concatenate([rsds_max_std, rsds_max_1_std],axis=0)[8:16]

rsds_dif = rsds_max_mean_all - rsds_min_mean_all

sw_dif = rsds_dif-rsus_dif
#如果认为向上为正
# sw_dif = rsus_dif-rsds_dif

year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year

# year_min = rlus.time[ross_sie_cmip6<Q1].dt.year
rlus_min_lst = []
for i in range(len(year_min)):
    y = str(year_min[i].values).zfill(4)
#     print(y)
    rlus_min_lst.append(rlus.sel(time=y))
rlus_min_mean = np.array(rlus_min_lst).squeeze().mean(axis=0)
rlus_min_std = np.array(rlus_min_lst).squeeze().std(axis=0)

# year_max = rlus.time[ross_sie_cmip6<Q1].dt.year
rlus_max_lst = []
for i in range(len(year_max)):
    y = str(year_max[i].values).zfill(4)
#     print(y)
    rlus_max_lst.append(rlus.sel(time=y))
rlus_max_mean = np.array(rlus_max_lst).squeeze().mean(axis=0)
rlus_max_std = np.array(rlus_max_lst).squeeze().std(axis=0)

year_min_1 = sic_Feb.time[ross_sie_cmip6<Q1].dt.year+1
year_max_1 = sic_Feb.time[ross_sie_cmip6>Q3].dt.year+1

# year_min_1 = rlus.time[ross_sie_cmip6<Q1].dt.year
rlus_min_1_lst = []
for i in range(len(year_min_1)):
    y = str(year_min_1[i].values).zfill(4)
#     print(y)
    rlus_min_1_lst.append(rlus.sel(time=y))
rlus_min_1_mean = np.array(rlus_min_1_lst).squeeze().mean(axis=0)
rlus_min_1_std = np.array(rlus_min_1_lst).squeeze().std(axis=0)

# year_max_1 = rlus.time[ross_sie_cmip6<Q1].dt.year
rlus_max_1_lst = []
for i in range(len(year_max_1)):
    y = str(year_max_1[i].values).zfill(4)
#     print(y)
    rlus_max_1_lst.append(rlus.sel(time=y))
rlus_max_1_mean = np.array(rlus_max_1_lst).squeeze().mean(axis=0)
rlus_max_1_std = np.array(rlus_max_1_lst).squeeze().std(axis=0)

rlus_min_mean_all = np.concatenate([rlus_min_mean,rlus_min_1_mean],axis=0)[8:16]
rlus_max_mean_all = np.concatenate([rlus_max_mean,rlus_max_1_mean],axis=0)[8:16]
rlus_min_std_all = np.concatenate([rlus_min_std, rlus_min_1_std],axis=0)[8:16]
rlus_max_std_all = np.concatenate([rlus_max_std, rlus_max_1_std],axis=0)[8:16]

rlus_dif = rlus_max_mean_all - rlus_min_mean_all

year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year

# year_min = rlds.time[ross_sie_cmip6<Q1].dt.year
rlds_min_lst = []
for i in range(len(year_min)):
    y = str(year_min[i].values).zfill(4)
#     print(y)
    rlds_min_lst.append(rlds.sel(time=y))
rlds_min_mean = np.array(rlds_min_lst).squeeze().mean(axis=0)
rlds_min_std = np.array(rlds_min_lst).squeeze().std(axis=0)

# year_max = rlds.time[ross_sie_cmip6<Q1].dt.year
rlds_max_lst = []
for i in range(len(year_max)):
    y = str(year_max[i].values).zfill(4)
#     print(y)
    rlds_max_lst.append(rlds.sel(time=y))
rlds_max_mean = np.array(rlds_max_lst).squeeze().mean(axis=0)
rlds_max_std = np.array(rlds_max_lst).squeeze().std(axis=0)

year_min_1 = sic_Feb.time[ross_sie_cmip6<Q1].dt.year+1
year_max_1 = sic_Feb.time[ross_sie_cmip6>Q3].dt.year+1

# year_min_1 = rlds.time[ross_sie_cmip6<Q1].dt.year
rlds_min_1_lst = []
for i in range(len(year_min_1)):
    y = str(year_min_1[i].values).zfill(4)
#     print(y)
    rlds_min_1_lst.append(rlds.sel(time=y))
rlds_min_1_mean = np.array(rlds_min_1_lst).squeeze().mean(axis=0)
rlds_min_1_std = np.array(rlds_min_1_lst).squeeze().std(axis=0)

# year_max_1 = rlds.time[ross_sie_cmip6<Q1].dt.year
rlds_max_1_lst = []
for i in range(len(year_max_1)):
    y = str(year_max_1[i].values).zfill(4)
#     print(y)
    rlds_max_1_lst.append(rlds.sel(time=y))
rlds_max_1_mean = np.array(rlds_max_1_lst).squeeze().mean(axis=0)
rlds_max_1_std = np.array(rlds_max_1_lst).squeeze().std(axis=0)

rlds_min_mean_all = np.concatenate([rlds_min_mean,rlds_min_1_mean],axis=0)[8:16]
rlds_max_mean_all = np.concatenate([rlds_max_mean,rlds_max_1_mean],axis=0)[8:16]
rlds_min_std_all = np.concatenate([rlds_min_std, rlds_min_1_std],axis=0)[8:16]
rlds_max_std_all = np.concatenate([rlds_max_std, rlds_max_1_std],axis=0)[8:16]

rlds_dif = rlds_max_mean_all - rlds_min_mean_all

lw_dif = rlds_dif-rlus_dif
#如果认为向上为正的话，用向上的减去想下的
# lw_dif = rlus_dif-rlds_dif

year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year

# year_min = hfls.time[ross_sie_cmip6<Q1].dt.year
hfls_min_lst = []
for i in range(len(year_min)):
    y = str(year_min[i].values).zfill(4)
#     print(y)
    hfls_min_lst.append(hfls.sel(time=y))
hfls_min_mean = np.array(hfls_min_lst).squeeze().mean(axis=0)
hfls_min_std = np.array(hfls_min_lst).squeeze().std(axis=0)

# year_max = hfls.time[ross_sie_cmip6<Q1].dt.year
hfls_max_lst = []
for i in range(len(year_max)):
    y = str(year_max[i].values).zfill(4)
#     print(y)
    hfls_max_lst.append(hfls.sel(time=y))
hfls_max_mean = np.array(hfls_max_lst).squeeze().mean(axis=0)
hfls_max_std = np.array(hfls_max_lst).squeeze().std(axis=0)

year_min_1 = sic_Feb.time[ross_sie_cmip6<Q1].dt.year+1
year_max_1 = sic_Feb.time[ross_sie_cmip6>Q3].dt.year+1

# year_min_1 = hfls.time[ross_sie_cmip6<Q1].dt.year
hfls_min_1_lst = []
for i in range(len(year_min_1)):
    y = str(year_min_1[i].values).zfill(4)
#     print(y)
    hfls_min_1_lst.append(hfls.sel(time=y))
hfls_min_1_mean = np.array(hfls_min_1_lst).squeeze().mean(axis=0)
hfls_min_1_std = np.array(hfls_min_1_lst).squeeze().std(axis=0)

# year_max_1 = hfls.time[ross_sie_cmip6<Q1].dt.year
hfls_max_1_lst = []
for i in range(len(year_max_1)):
    y = str(year_max_1[i].values).zfill(4)
#     print(y)
    hfls_max_1_lst.append(hfls.sel(time=y))
hfls_max_1_mean = np.array(hfls_max_1_lst).squeeze().mean(axis=0)
hfls_max_1_std = np.array(hfls_max_1_lst).squeeze().std(axis=0)

hfls_min_mean_all = np.concatenate([hfls_min_mean,hfls_min_1_mean],axis=0)[8:16]
hfls_max_mean_all = np.concatenate([hfls_max_mean,hfls_max_1_mean],axis=0)[8:16]
hfls_min_std_all = np.concatenate([hfls_min_std, hfls_min_1_std],axis=0)[8:16]
hfls_max_std_all = np.concatenate([hfls_max_std, hfls_max_1_std],axis=0)[8:16]

hfls_dif = hfls_max_mean_all - hfls_min_mean_all

year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year

# year_min = hfss.time[ross_sie_cmip6<Q1].dt.year
hfss_min_lst = []
for i in range(len(year_min)):
    y = str(year_min[i].values).zfill(4)
#     print(y)
    hfss_min_lst.append(hfss.sel(time=y))
hfss_min_mean = np.array(hfss_min_lst).squeeze().mean(axis=0)
hfss_min_std = np.array(hfss_min_lst).squeeze().std(axis=0)

# year_max = hfss.time[ross_sie_cmip6<Q1].dt.year
hfss_max_lst = []
for i in range(len(year_max)):
    y = str(year_max[i].values).zfill(4)
#     print(y)
    hfss_max_lst.append(hfss.sel(time=y))
hfss_max_mean = np.array(hfss_max_lst).squeeze().mean(axis=0)
hfss_max_std = np.array(hfss_max_lst).squeeze().std(axis=0)

year_min_1 = sic_Feb.time[ross_sie_cmip6<Q1].dt.year+1
year_max_1 = sic_Feb.time[ross_sie_cmip6>Q3].dt.year+1

# year_min_1 = hfss.time[ross_sie_cmip6<Q1].dt.year
hfss_min_1_lst = []
for i in range(len(year_min_1)):
    y = str(year_min_1[i].values).zfill(4)
#     print(y)
    hfss_min_1_lst.append(hfss.sel(time=y))
hfss_min_1_mean = np.array(hfss_min_1_lst).squeeze().mean(axis=0)
hfss_min_1_std = np.array(hfss_min_1_lst).squeeze().std(axis=0)

# year_max_1 = hfss.time[ross_sie_cmip6<Q1].dt.year
hfss_max_1_lst = []
for i in range(len(year_max_1)):
    y = str(year_max_1[i].values).zfill(4)
#     print(y)
    hfss_max_1_lst.append(hfss.sel(time=y))
hfss_max_1_mean = np.array(hfss_max_1_lst).squeeze().mean(axis=0)
hfss_max_1_std = np.array(hfss_max_1_lst).squeeze().std(axis=0)

hfss_min_mean_all = np.concatenate([hfss_min_mean,hfss_min_1_mean],axis=0)[8:16]
hfss_max_mean_all = np.concatenate([hfss_max_mean,hfss_max_1_mean],axis=0)[8:16]
hfss_min_std_all = np.concatenate([hfss_min_std, hfss_min_1_std],axis=0)[8:16]
hfss_max_std_all = np.concatenate([hfss_max_std, hfss_max_1_std],axis=0)[8:16]

hfss_dif = hfss_max_mean_all - hfss_min_mean_all

# lw
lats = rlds.lat
lons = rlds.lon
time_len = len(lw_dif)
ross_lw_dif = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
    lw_anom_ross = np.where(mask_ross,lw_dif[i],np.nan)
    ross_lw_dif[i] = np.nanmean(lw_anom_ross)
# sw
lats = rsds.lat
lons = rsds.lon
time_len = len(sw_dif)
ross_sw_dif = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
    sw_anom_ross = np.where(mask_ross,sw_dif[i],np.nan)
    ross_sw_dif[i] = np.nanmean(sw_anom_ross)
# lhf
lats = hfls.lat
lons = hfls.lon
time_len = len(hfls_dif)
ross_hfls_dif = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
    hfls_anom_ross = np.where(mask_ross,hfls_dif[i],np.nan)
    ross_hfls_dif[i] = np.nanmean(hfls_anom_ross)
# shf
lats = hfss.lat
lons = hfss.lon
time_len = len(hfss_dif)
ross_hfss_dif = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
    hfss_anom_ross = np.where(mask_ross,hfss_dif[i],np.nan)
    ross_hfss_dif[i] = np.nanmean(hfss_anom_ross)
Qnet_dif = hfls_dif+hfss_dif+sw_dif+lw_dif
# Qnet
lats = hfss.lat
lons = hfss.lon
time_len = len(Qnet_dif)
ross_Qnet_dif = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
    Qnet_anom_ross = np.where(mask_ross,Qnet_dif[i],np.nan)
    ross_Qnet_dif[i] = np.nanmean(Qnet_anom_ross)

thetao = xr.open_dataset('/stu02/weizx24/data/thetao_Omon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['thetao'].sel(lev=slice(0,5000)).mean('lev')
year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year

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

thetao_min_mean_all = np.concatenate([thetao_min_mean,thetao_min_1_mean],axis=0)[8:16]
thetao_max_mean_all = np.concatenate([thetao_max_mean,thetao_max_1_mean],axis=0)[8:16]
thetao_min_std_all = np.concatenate([thetao_min_std, thetao_min_1_std],axis=0)[8:16]
thetao_max_std_all = np.concatenate([thetao_max_std, thetao_max_1_std],axis=0)[8:16]

thetao_dif = thetao_max_mean_all - thetao_min_mean_all
# lhf
lats = thetao.lat
lons = thetao.lon
time_len = len(thetao_dif)
ross_thetao_dif = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
    thetao_anom_ross = np.where(mask_ross,thetao_dif[i],np.nan)
    ross_thetao_dif[i] = np.nanmean(thetao_anom_ross)

np.savez('/stu02/weizx24/data/npz/Figure7_model_0924.npz',ross_sw_dif=ross_sw_dif,ross_lw_dif=ross_lw_dif,ross_hfls_dif=ross_hfls_dif,ross_hfss_dif=ross_hfss_dif,ross_Qnet_dif=ross_Qnet_dif,ross_thetao_dif=ross_thetao_dif)
print('观测数据处理完毕')
#endregion