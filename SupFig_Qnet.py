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

#气候平均态取1979-2018年就好
lw = xr.open_dataset('D:/ERA5/ERA5_79-22_lwsw.nc')['msnlwrf'].sel(time=slice('1992','2022'))
sw = xr.open_dataset('D:/ERA5/ERA5_79-22_lwsw.nc')['msnswrf'].sel(time=slice('1992','2022'))

lw_clm = lw.sel(time=slice('1992','2021')).values.reshape(30,12,161,1440).mean(axis=0)
lw_std = lw.sel(time=slice('1992','2021')).values.reshape(30,12,161,1440).std(axis=0)
sw_clm = sw.sel(time=slice('1992','2021')).values.reshape(30,12,161,1440).mean(axis=0)
sw_std = sw.sel(time=slice('1992','2021')).values.reshape(30,12,161,1440).std(axis=0)

lw_clm_sel = np.concatenate([lw_clm[8:12],lw_clm[0:4]])
lw_std_sel = np.concatenate([lw_std[8:12],lw_std[0:4]])
sw_clm_sel = np.concatenate([sw_clm[8:12],sw_clm[0:4]])
sw_std_sel = np.concatenate([sw_std[8:12],sw_std[0:4]])


lhf = xr.open_dataset('D:/ERA5/ERA5_79-22_lhfshf.nc')['mslhf'].sel(time=slice('1992','2022'))
shf = xr.open_dataset('D:/ERA5/ERA5_79-22_lhfshf.nc')['msshf'].sel(time=slice('1992','2022'))

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

votemp = xr.open_dataset('D:/ORAS5/votemp_79-22_50m.nc')['votemper'].mean(dim = 'LEV').sel(time=slice('1992','2022'))
votemp_clm = votemp.sel(time=slice('1992','2021')).values.reshape((30,12,40,360)).mean(axis=0)
votemp_std = votemp.sel(time=slice('1992','2021')).values.reshape((30,12,40,360)).std(axis=0)
votemp_clm_sel = np.concatenate([votemp_clm[8:12],votemp_clm[0:4]])
votemp_std_sel = np.concatenate([votemp_std[8:12],votemp_std[0:4]])


# votemp_anom_f = [votemp_anom[1],votemp_anom[3],votemp_anom[5],votemp_anom[7]]
# votemp_std_f = [votemp_std_sel[1],votemp_std_sel[3],votemp_std_sel[5],votemp_std_sel[7]]
date = ['Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr']

years = np.linspace(1992,2021,30,dtype=int)

lw_sel = list()
sw_sel = list()
shf_sel = list()
lhf_sel = list()
votemp_sel = list()
Qnet_sel = list()

for i in range(len(years)):
    starttime = str(int(years[i]))+'-09-01'
    endtime = str(int(years[i])+1)+'-04-30'

    lw_sel.append(lw.sel(time=slice(starttime,endtime)))
    sw_sel.append(sw.sel(time=slice(starttime,endtime)))
    shf_sel.append(shf.sel(time=slice(starttime,endtime)))
    lhf_sel.append(lhf.sel(time=slice(starttime,endtime)))
    votemp_sel.append(votemp.sel(time=slice(starttime,endtime)))
    Qnet_sel.append(Qnet.sel(time=slice(starttime,endtime)))

lons = lw.longitude
lats = lw.latitude

# lw
def lm_trend(data):
    years = np.linspace(1992,2021,30,dtype=int)
    aa, intercept, r_value, p1, std_err = stats.linregress(years, data)
    return (aa, intercept, r_value, p1, std_err)

ross_lw_trend_list = []
for i in range(8):
    data_sel = np.array(lw_sel)[:,i]
#     print(data_sel.shape)
    time_len = len(data_sel)
    ross_lw = np.zeros((time_len),'float')
    for j in range(time_len):
        mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
        lw_anom_ross = np.where(mask_ross,data_sel[j],np.nan)
        ross_lw[j] = np.nanmean(lw_anom_ross)
#     print(ross_lw)
    ross_lw_trend = xr.apply_ufunc(
        lm_trend,
        np.array(ross_lw),
        input_core_dims=[["time"]],
        output_core_dims=[[], [], [], [], []],
        vectorize=True)
    ross_lw_trend_list.append(ross_lw_trend[0])

    # sw
def lm_trend(data):
    years = np.linspace(1992,2021,30,dtype=int)
    aa, intercept, r_value, p1, std_err = stats.linregress(years, data)
    return (aa, intercept, r_value, p1, std_err)

ross_sw_trend_list = []
for i in range(8):
    data_sel = np.array(sw_sel)[:,i]
#     print(data_sel.shape)
    time_len = len(data_sel)
    ross_sw = np.zeros((time_len),'float')
    for j in range(time_len):
        mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
        sw_anom_ross = np.where(mask_ross,data_sel[j],np.nan)
        ross_sw[j] = np.nanmean(sw_anom_ross)
#     print(ross_sw)
    ross_sw_trend = xr.apply_ufunc(
        lm_trend,
        np.array(ross_sw),
        input_core_dims=[["time"]],
        output_core_dims=[[], [], [], [], []],
        vectorize=True)
    ross_sw_trend_list.append(ross_sw_trend[0])

# shf
def lm_trend(data):
    years = np.linspace(1992,2021,30,dtype=int)
    aa, intercept, r_value, p1, std_err = stats.linregress(years, data)
    return (aa, intercept, r_value, p1, std_err)

ross_shf_trend_list = []
for i in range(8):
    data_sel = np.array(shf_sel)[:,i]
#     print(data_sel.shape)
    time_len = len(data_sel)
    ross_shf = np.zeros((time_len),'float')
    for j in range(time_len):
        mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
        shf_anom_ross = np.where(mask_ross,data_sel[j],np.nan)
        ross_shf[j] = np.nanmean(shf_anom_ross)
#     print(ross_shf)
    ross_shf_trend = xr.apply_ufunc(
        lm_trend,
        np.array(ross_shf),
        input_core_dims=[["time"]],
        output_core_dims=[[], [], [], [], []],
        vectorize=True)
    ross_shf_trend_list.append(ross_shf_trend[0])

# lhf
def lm_trend(data):
    years = np.linspace(1992,2021,30,dtype=int)
    aa, intercept, r_value, p1, std_err = stats.linregress(years, data)
    return (aa, intercept, r_value, p1, std_err)

ross_lhf_trend_list = []
for i in range(8):
    data_sel = np.array(lhf_sel)[:,i]
#     print(data_sel.shape)
    time_len = len(data_sel)
    ross_lhf = np.zeros((time_len),'float')
    for j in range(time_len):
        mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
        lhf_anom_ross = np.where(mask_ross,data_sel[j],np.nan)
        ross_lhf[j] = np.nanmean(lhf_anom_ross)
#     print(ross_lhf)
    ross_lhf_trend = xr.apply_ufunc(
        lm_trend,
        np.array(ross_lhf),
        input_core_dims=[["time"]],
        output_core_dims=[[], [], [], [], []],
        vectorize=True)
    ross_lhf_trend_list.append(ross_lhf_trend[0])

# votemp
def lm_trend(data):
    years = np.linspace(1992,2021,30,dtype=int)
    aa, intercept, r_value, p1, std_err = stats.linregress(years, data)
    return (aa, intercept, r_value, p1, std_err)

lons = votemp.longitude
lats = votemp.latitude

ross_votemp_trend_list = []
for i in range(8):
    data_sel = np.array(votemp_sel)[:,i]
#     print(data_sel.shape)
    time_len = len(data_sel)
    ross_votemp = np.zeros((time_len),'float')
    for j in range(time_len):
        mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
        votemp_anom_ross = np.where(mask_ross,data_sel[j],np.nan)
        ross_votemp[j] = np.nanmean(votemp_anom_ross)
#     print(ross_votemp)
    ross_votemp_trend = xr.apply_ufunc(
        lm_trend,
        np.array(ross_votemp),
        input_core_dims=[["time"]],
        output_core_dims=[[], [], [], [], []],
        vectorize=True)
    ross_votemp_trend_list.append(ross_votemp_trend[0])

# Qnet
def lm_trend(data):
    years = np.linspace(1992,2021,30,dtype=int)
    aa, intercept, r_value, p1, std_err = stats.linregress(years, data)
    return (aa, intercept, r_value, p1, std_err)

lons = Qnet.longitude
lats = Qnet.latitude

ross_Qnet_trend_list = []
for i in range(8):
    data_sel = np.array(Qnet_sel)[:,i]
#     print(data_sel.shape)
    time_len = len(data_sel)
    ross_Qnet = np.zeros((time_len),'float')
    for j in range(time_len):
        mask_ross =  (lats<-70.)&(lons>180.)&(lons<230.)
        Qnet_anom_ross = np.where(mask_ross,data_sel[j],np.nan)
        ross_Qnet[j] = np.nanmean(Qnet_anom_ross)
#     print(ross_Qnet)
    ross_Qnet_trend = xr.apply_ufunc(
        lm_trend,
        np.array(ross_Qnet),
        input_core_dims=[["time"]],
        output_core_dims=[[], [], [], [], []],
        vectorize=True)
    ross_Qnet_trend_list.append(ross_Qnet_trend[0])

plt.close()

fig = plt.figure(1, figsize=(12,9))
ax = fig.add_subplot(111)
ax.plot(date,ross_lw_trend_list,color='C2',markersize=15,linewidth=4,label='LW',alpha=0.5)
ax.plot(date,ross_sw_trend_list,color='y',markersize=15,linewidth=4,label='SW',alpha=0.5)
ax.plot(date,ross_lhf_trend_list,color='r',markersize=13,linewidth=4,label='LHF',alpha=0.5)
ax.plot(date,ross_shf_trend_list,color='k',markersize=13,linewidth=4,label='SHF',alpha=0.5)
ax.plot(date,ross_Qnet_trend_list,color='b',marker='^',markersize=13,linewidth=4,label='Qnet')
ax1 = ax.twinx()
# ax1.plot(date,ross_sst_anom,color='m',marker='8',markersize=13,linewidth=4,label='SST_anom')
ax1.plot(date,ross_votemp_trend_list,color='C3',marker='8',markersize=13,linewidth=4,label='UOT')
lines = ax.get_lines()+ax1.get_lines()
plt.legend(lines, [l.get_label() for l in lines],edgecolor='k',fontsize=20,loc='upper left')
# plt.legend(lines, [l.get_label() for l in lines],edgecolor='k',fontsize=20,loc='lower left')
ax1.set_ylabel('Trend in UOT (℃)',fontsize=22,color='red')
ax1.set_yticks(np.linspace(0,0.02,9))
ax1.yaxis.set_tick_params(labelsize=20)
ax1.tick_params('both', length=8, width=1.2, which='major')
ax1.tick_params('both', length=4, width=0.5, which='minor')
ax1.tick_params(colors='red')
# ax1.set_ylim([-0.01,0.14])
ax.set_xlabel('Month',fontsize=20)
ax.set_ylabel('Trend in Qnet ' + r'(W m$^{\text{-2}}$)',fontsize=20)
ax.tick_params('both', length=8, width=1.2, which='major')
ax.tick_params('both', length=4, width=0.5, which='minor')
# ax.set_ylim([0,1.1])
#ax.set_xlim([0,151])
# ax.set_title('1979-2022 Feb Eastern Ross Sea SIE time series',fontsize=24)
#x.set_title('Equatorial-zonal-wind-50-hPa')
#ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.yaxis.set_tick_params(labelsize=20)
ax.xaxis.set_tick_params(labelsize=20)
ax.set_title('Trend, 1992-2021',fontsize=22)
#ax.axhline(y=0.4,color='C0',linestylae='--')
#ax.set_xlim([time_s2s[0], time_s2s[-1]+ pd.Timedelta(days=32)])
#ax.axhline(y=0,color='C0',linestyle='--')
#ax.axhline(y=ant_ref,color='C0',linestyle='--',linewidth=5,label='piControl')
ax.grid(color='lightgray',linestyle='--',linewidth=3)
#给两个y轴的线条添加Legend
# lines = ax.get_lines()+ax1.get_lines()
plt.savefig('C:/Users/fzjxw/python/code/Figures/Qnet_timeseries_trend.png',dpi=300,bbox_inches='tight')
plt.show()