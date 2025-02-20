import numpy as np
import matplotlib.pyplot as plt
import glob
import matplotlib.dates as mdates
from pylab import *
from matplotlib import dates, ticker
from scipy import stats
import xarray as xr
import matplotlib.path as mpath
from datetime import datetime, timedelta, date
from dateutil.relativedelta import relativedelta
import pandas as pd
import scipy.signal as signal
import cmaps
from mpl_toolkits.basemap import Basemap
from scipy.stats import t
import os 
from scipy.stats import pearsonr

ftools='D:/tools/'
with open(ftools+'pss25lats_v3.dat','rb') as flat:
    lats = np.fromfile(flat,dtype='<i4').reshape(332, 316)/100000.

with open(ftools+'pss25lons_v3.dat','rb') as flon:
    lons = np.fromfile(flon,dtype='<i4').reshape(332, 316)/100000.

with open(ftools+'pss25area_v3.dat','rb') as flon:
    area = np.fromfile(flon,dtype='<i4').reshape(332, 316)*1e-9

'''
#合成分析年份
year_min = [1994., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2016., 2018., 2021.]
num_y = 5

#region SIC
sic = xr.open_dataset('D:/SIC/monthly_sic/seaice_conc_monthly_sh_197811_202212_sub_latlon.nc')['nsidc_nt_seaice_conc_monthly'].sel(time=slice('1992','2022'))
year = np.linspace(1992,2021,30)
sic_list = []
for i in range(30):
    start_time = str(int(year[i]))+'-09-01'
    end_time = str(int(year[i]+1))+'-04-30'
    sic_list.append(sic.sel(time=slice(start_time,end_time)))
#     print(sic_sel)
sic_arr = np.array(sic_list)

#计算罗斯海区域平均SIC
#计算所有月份的罗斯海平均SIC，顺序从9月到次年4月
ross_sic = np.zeros((30,8),'float')
for j in range(8):
    for i in range(30):
        #仅包含东罗斯海
        mask_ross_nt = (lons>-180.)& (lons<-130.) & (lats<-70.) & (sic_arr[i,j]>=0)
        sic_ross_nt = np.where(mask_ross_nt,sic_arr[i,j],np.nan)
        ross_sic[i,j] = np.nanmean(sic_ross_nt)


#region 海洋热含量
date = ['Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr']
#计算每一年的罗斯海上层50m海洋温度
votemp = xr.open_dataset('D:/ORAS5/votemp_79-22_50m.nc')['votemper'].mean(dim = 'LEV').sel(time=slice('1992','2022'))
year = np.linspace(1992,2021,30)
votemp_list = []
for i in range(30):
    start_time = str(int(year[i]))+'-09-01'
    end_time = str(int(year[i]+1))+'-04-30'
    votemp_list.append(votemp.sel(time=slice(start_time,end_time)))
votemp_arr = np.array(votemp_list)

lats = votemp.latitude
lons = votemp.longitude
ross_votemp = np.zeros((30,8),'float')
for j in range(8):
    for i in range(30):
        #包含整个罗斯海
        mask_ross = (lats<-70.)&(lons>160.)&(lons<230.)
        votemp_anom_ross = np.where(mask_ross,votemp_arr[i,j],np.nan)
        ross_votemp[i,j] = np.nanmean(votemp_anom_ross)
#endregion

#region 作图
#设置全局变量并开始作图
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 14  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
# plt.rc('font', family='Arial', size=9)

#————————作图————————————
#次年3月SIC和次年2月UOT的相关性图
fig = plt.figure(figsize=(5.5,4.5))
ax1 = fig.add_subplot(111)
#3月SIC和2月votemp
ax1.scatter(ross_sic[:,6],ross_votemp[:,5],color='black',label='Observation',s=50)
corre,p = pearsonr(ross_sic[:,6],ross_votemp[:,5])
print(f'相关系数是{corre}')
print(f'p值是{p}')
aa,intercept,r_value,p1,std_err = stats.linregress(ross_sic[:,6],ross_votemp[:,5])
y = aa*np.array(ross_sic[:,6])+intercept
#ax1.plot(ross_sic[:,6],y,color='lightgray',linestyle='--')
plt.xlabel('Mar SIC')
plt.ylabel('Feb UOT (°C) ')
plt.grid(True,linestyle='--',color='lightgray')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax1.text(0.7, 0.9,'pearson r: '+str(round(corre,2))+'\n'+'p value:<0.05',
        transform=ax1.transAxes, 
        verticalalignment='top', bbox=props,)
ax1.text(0, 1.2, '(d)', fontsize=12, transform=ax1.transAxes, va='top', ha='right')
#plt.legend(markerscale=1)
plt.subplots_adjust(wspace=0.3)
plt.savefig('C:/Users/fzjxw/python/code/Figures/3SIC&2UOT_0104.png',dpi=300)

#次年4月SIC和次年2月UOT的相关性图
fig = plt.figure(figsize=(5.5,4.5))
ax1 = fig.add_subplot(111)
#3月SIC和2月votemp
ax1.scatter(ross_sic[:,7],ross_votemp[:,5],color='black',label='Observation',s=50)
corre,p = pearsonr(ross_sic[:,7],ross_votemp[:,5])
print(f'相关系数是{corre}')
print(f'p值是{p}')
aa,intercept,r_value,p1,std_err = stats.linregress(ross_sic[:,6],ross_votemp[:,5])
y = aa*np.array(ross_sic[:,6])+intercept
#ax1.plot(ross_sic[:,6],y,color='lightgray',linestyle='--')
plt.xlabel('Apr SIC')
plt.ylabel('Feb UOT (°C) ')
plt.grid(True,linestyle='--',color='lightgray')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax1.text(0.7, 0.9,'pearson r: '+str(round(corre,2))+'\n'+'p value:<0.05',
        transform=ax1.transAxes, 
        verticalalignment='top', bbox=props,)
ax1.text(0, 1.2, '(d)', fontsize=12, transform=ax1.transAxes, va='top', ha='right')
#plt.legend(markerscale=1)
plt.subplots_adjust(wspace=0.3)
plt.savefig('C:/Users/fzjxw/python/code/Figures/4SIC&2UOT_0104.png',dpi=300)
#endregion
'''

#region 其他：计算5月份，要重新算SIC
#次年5月SIC和次年2月UOT的相关性图
#合成分析年份
year_min = [1994., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2016., 2018., 2021.]
num_y = 5
sic = xr.open_dataset('D:/SIC/monthly_sic/seaice_conc_monthly_sh_197811_202212_sub_latlon.nc')['nsidc_nt_seaice_conc_monthly'].sel(time=slice('1992','2022'))
year = np.linspace(1992,2021,30)
sic_list = []
for i in range(30):
    start_time = str(int(year[i]))+'-10-01'
    end_time = str(int(year[i]+1))+'-05-31'
    sic_list.append(sic.sel(time=slice(start_time,end_time)))
#     print(sic_sel)
sic_arr = np.array(sic_list)
#计算罗斯海区域平均SIC
#计算所有月份的罗斯海平均SIC，顺序从9月到次年4月
ross_sic = np.zeros((30,8),'float')
for j in range(8):
    for i in range(30):
        #仅包含东罗斯海
        mask_ross_nt = (lons>-180.)& (lons<-130.) & (lats<-70.) & (sic_arr[i,j]>=0)
        sic_ross_nt = np.where(mask_ross_nt,sic_arr[i,j],np.nan)
        ross_sic[i,j] = np.nanmean(sic_ross_nt)

date = ['Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr']
#计算每一年的罗斯海上层50m海洋温度
votemp = xr.open_dataset('D:/ORAS5/votemp_79-22_50m.nc')['votemper'].mean(dim = 'LEV').sel(time=slice('1992','2022'))
year = np.linspace(1992,2021,30)
votemp_list = []
for i in range(30):
    start_time = str(int(year[i]))+'-09-01'
    end_time = str(int(year[i]+1))+'-04-30'
    votemp_list.append(votemp.sel(time=slice(start_time,end_time)))
votemp_arr = np.array(votemp_list)

lats = votemp.latitude
lons = votemp.longitude
ross_votemp = np.zeros((30,8),'float')
for j in range(8):
    for i in range(30):
        #包含整个罗斯海
        mask_ross = (lats<-70.)&(lons>160.)&(lons<230.)
        votemp_anom_ross = np.where(mask_ross,votemp_arr[i,j],np.nan)
        ross_votemp[i,j] = np.nanmean(votemp_anom_ross)


#region 作图
#设置全局变量并开始作图
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 14  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
# plt.rc('font', family='Arial', size=9)

fig = plt.figure(figsize=(5.5,4.5))
ax1 = fig.add_subplot(111)
#3月SIC和2月votemp
ax1.scatter(ross_sic[:,7],ross_votemp[:,5],color='black',label='Observation',s=50)
corre,p = pearsonr(ross_sic[:,7],ross_votemp[:,5])
print(f'相关系数是{corre}')
print(f'p值是{p}')
aa,intercept,r_value,p1,std_err = stats.linregress(ross_sic[:,6],ross_votemp[:,5])
y = aa*np.array(ross_sic[:,6])+intercept
#ax1.plot(ross_sic[:,6],y,color='lightgray',linestyle='--')
plt.xlabel('May SIC')
plt.ylabel('Feb UOT (°C) ')
plt.grid(True,linestyle='--',color='lightgray')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax1.text(0.7, 0.9,'pearson r: '+str(round(corre,2))+'\n'+'p value:<0.05',
        transform=ax1.transAxes, 
        verticalalignment='top', bbox=props,)
ax1.text(0, 1.2, '(d)', fontsize=12, transform=ax1.transAxes, va='top', ha='right')
#plt.legend(markerscale=1)
plt.subplots_adjust(wspace=0.3)
plt.savefig('C:/Users/fzjxw/python/code/Figures/5SIC&2UOT_0104.png',dpi=300)
#endregion
