'''sea ice concentration longterm evolution
first start with CESM2
'''
import numpy as np
import matplotlib.pyplot as plt
import glob
from pylab import *
import xarray as xr
import pandas as pd
from scipy import stats
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import scipy.signal as signal
import cmaps
import matplotlib.dates as mdates
import matplotlib.path as mpath
from matplotlib import dates, ticker
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.gridspec as gridspec
from scipy.stats import t
from pyproj import Proj, transform
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import AxesGrid
from scipy.stats import pearsonr
from scipy.stats import pearsonr

#————————————————观测————————————————
import os
folder_path = 'D:/opw_rec/nsidc/new/'
files = os.listdir(folder_path)
files.sort()
files
a = list()
for i in files:
#     prini(i)
    a.append(pd.read_csv(folder_path+i)['eastross']+pd.read_csv(folder_path+i)['westross'])
a = np.array(a)
from scipy.stats import pearsonr

opwa_12_1 = a[:,27:34].mean(axis=1)

#春季SON平均的ASL Relative Center Pressure
ASL_SON = pd.read_excel('D:/ASL/ASL_SON.xlsx')['ActCenPres']
ASL_SON_y = list()
n = 2
for i in range(43):
    ASL_SON_y.append((ASL_SON[3*i+1]+ASL_SON[3*i+2])/n)
ASL_SON_Act_arr = np.array(ASL_SON_y)[13:]

#春季SON平均的ASL Actual Center Pressure
ASL_SON = pd.read_excel('D:/ASL/ASL_SON.xlsx')['SectorPres']
ASL_SON_y = list()
n = 2
for i in range(43):
    ASL_SON_y.append((ASL_SON[3*i+1]+ASL_SON[3*i+2])/n)
ASL_SON_Sec_arr = np.array(ASL_SON_y)[13:]

yy = np.linspace(1992,2021,30)

flux_sea= np.load('D:/npz/Ross_adv_flux_0414.npz')['flux_sea'][13:]/1000

# #如果不做去趋势处理
opwa = opwa_12_1
asl_act = ASL_SON_Act_arr
asl_sec = ASL_SON_Sec_arr
flux_sea = flux_sea
#去趋势
# opwa = signal.detrend(opwa_12_1)
# asl_act = signal.detrend(ASL_SON_Act_arr)
# asl_sec = signal.detrend(ASL_SON_Sec_arr)
# flux_sea = signal.detrend(flux_sea)

#设置全局变量并开始作图
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 9
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
# plt.rc('font', family='Arial', size=9)

#——————————图1————————————————
fig = plt.figure(figsize=(5.9,5.2))
ax1 = fig.add_subplot(221)
#画出散点图分布，计算变量的相关性
ax1.scatter(opwa,asl_sec,marker='o',color='k')
p = pearsonr(opwa,asl_sec)
aa,intercept,r_value,p1,std_err = stats.linregress(opwa,asl_sec)
print(p)
y = aa*np.array(opwa)+intercept
ax1.plot(opwa,y,color='k')
ax1.set_xticks([0,0.05,0.1,0.15,0.2,0.25,0.3])
ax1.set_yticks(np.linspace(970,1000,5))
# ax1.set_xticks([-0.1,-0.05,0,0.05,0.1,])
# ax1.set_yticks(np.linspace(-5,15,5))
ax1.set_ylabel('Sector pressure (hPa)')
ax1.set_xlabel('Areas'+ r' (10$^{6}$ km$^{2}$)')
ax1.grid(linestyle='--',alpha=0.5)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax1.text(0.4, 0.9,'Pearson\'s r: '+str(round(p[0],2))+'\n'+'p value <0.05',
        transform=ax1.transAxes, 
        verticalalignment='top', bbox=props,)
ax1.text(0, 1.13, '(a)', fontsize=12, transform=ax1.transAxes, va='top', ha='right')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

#——————————图2————————————————
ax2= fig.add_subplot(222)
ax2.scatter(opwa,flux_sea,marker='o',color='k')
p = pearsonr(opwa,flux_sea)
#画出散点的拟合直线
aa,intercept,r_value,p1,std_err = stats.linregress(opwa,flux_sea)
print(p)
y = aa*np.array(opwa)+intercept
ax2.plot(opwa,y,color='k')
ax2.set_xticks([0,0.05,0.1,0.15,0.2,0.25,0.3])
# ax.set_yticks([0,0.5,1.0,1.5,2.0])
ax2.set_ylabel('Area flux'+r' (10$^{\text{6}}$ km$^{\text{2}}$ day$^{\text{-1}}$)')
ax2.set_xlabel('Areas'+ r' (10$^{6}$ km$^{2}$)')
# ax.yaxis.set_tick_params(labelsize=20)
# ax.xaxis.set_tick_params(labelsize=20)
ax2.grid(linestyle='--',alpha=0.5)
# figname = 'ASL & sea ice area flux'
# ax.set_title(figname,fontsize=22)
# ax.grid(color='lightgray',linestyle='--',linewidth=3)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax2.text(0.5, 0.2,'Pearson\'s r: '+str(round(p[0],2))+'\n'+'p value < 0.05',
        transform=ax2.transAxes,
        verticalalignment='top', bbox=props,)
ax2.text(0, 1.13, '(b)', fontsize=12, transform=ax2.transAxes, va='top', ha='right')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# plt.tight_layout(rect=[0,0,0.9,1])
# plt.savefig('/Users/guizijin/DATA/FIGURES/Figure7_obs_20240603.pdf' ,dpi=600,bbox_inches = 'tight')



#————————————————模式————————————————
# 筛选最大最小年份
CMIP6_opwa = np.load('D:/npz/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']*1e-12
ross_sie_cmip6 = np.nanmean(CMIP6_opwa,axis=1)
Q1 = np.quantile(ross_sie_cmip6,0.05)
Q3 = np.quantile(ross_sie_cmip6,0.95)
sic_Feb = xr.open_dataset('D:/CMIP6/CESM2-WACCM-FV2/piControl/SImon/siconc/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_Feb.nc')['siconc']
year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year

psl_regionmean = xr.open_dataset('D:/CMIP6/CESM2-WACCM-FV2/piControl/Amon/psl/psl_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all_sellonlatbox_fldmean.nc')['psl']*1e-2
psl_ON = psl_regionmean[(psl_regionmean.time.dt.month==10)|(psl_regionmean.time.dt.month==11)].values.reshape(499,2).mean(axis=1)

file =np.load('D:/npz/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']*1e-12
ross_sie_cmip6 = np.nanmean(file,axis=1)

ua = xr.open_dataset('D:/CMIP6/CESM2-WACCM-FV2/piControl/Amon/ua/ua_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['ua'].sel(plev=92500)
va = xr.open_dataset('D:/CMIP6/CESM2-WACCM-FV2/piControl/Amon/va/va_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['va'].sel(plev=92500)
lons = va.lon
lats = va.lat
va_ON = va[(va.time.dt.month==10)|(va.time.dt.month==11)].values.reshape(499,2,96,144).mean(axis=1)
# va_ON = va[(va.time.dt.month==11)].values.reshape(499,1,96,144).mean(axis=1)
time_len = len(va_ON)
ross_v = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-60.)&(lats>-80.)&(lons>160.)&(lons<230.)
    va_ross = np.where(mask_ross,va_ON[i],np.nan)
    ross_v[i] = np.nanmean(va_ross)

flux_sea = np.load('D:/npz/CMIP6_seaiceareaflux_0615.npz')['flux_sea']

#去趋势处理
# ross_sie_cmip6 = signal.detrend(ross_sie_cmip6)
# psl_ON = signal.detrend(psl_ON)
# flux_sea = signal.detrend(flux_sea)


#——————————图3————————————————
# fig = plt.figure(figsize=(20,6))
ax3 = fig.add_subplot(223)
#画出散点图分布，计算变量的相关性
ax3.scatter(ross_sie_cmip6,psl_ON,marker='o',color='k')
p = pearsonr(ross_sie_cmip6,psl_ON)
print(p)
#画出散点的拟合直线
aa,intercept,r_value,p1,std_err = stats.linregress(ross_sie_cmip6,psl_ON,)
y = aa*np.array(ross_sie_cmip6)+intercept
ax3.plot(ross_sie_cmip6,y,color='k')
ax3.set_ylabel('Sector pressure (hPa)')
ax3.set_xlabel('Areas '+ r' (10$^{6}$ km$^{2}$)')
ax3.set_xticks([0,0.05,0.1,0.15,0.2,0.25,0.3])
ax3.set_yticks(np.linspace(970,1000,5))
# ax3.set_yticks(np.linspace(-5,15,5))
# ax3.set_xticks([-0.05,0,0.05,0.1,0.15,0.2])
# ax.yaxis.set_tick_params(labelsize=20)
# ax.xaxis.set_tick_params(labelsize=20)
ax3.grid(linestyle='--',alpha=0.5)
# figname = 'Ross Sea'
# ax.set_title(figname,fontsize=22)
# ax.grid(color='lightgray',linestyle='--',linewidth=3)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax3.text(0.4, 0.9,'Pearson\'s r: '+str(round(p[0],2))+'\n'+'p value:<0.05',
        transform=ax3.transAxes, 
        verticalalignment='top', bbox=props,)
ax3.text(0, 1.13, '(c)', fontsize=12, transform=ax3.transAxes, va='top', ha='right')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

#——————————图4————————————————
from scipy.stats import pearsonr
ax4 = fig.add_subplot(224)
#画出散点图分布，计算变量的相关性
ax4.scatter(ross_sie_cmip6,flux_sea,marker='o',color='k')
p = pearsonr(ross_sie_cmip6,flux_sea)
print(p)
#画出散点的拟合直线
aa,intercept,r_value,p1,std_err = stats.linregress(ross_sie_cmip6,flux_sea)
y = aa*np.array(ross_sie_cmip6)+intercept
ax4.plot(ross_sie_cmip6,y,color='k')
ax4.set_ylabel('Area flux '+r' (10$^{\text{6}}$ km$^{\text{2}}$ day$^{\text{-1}}$)')
ax4.set_xlabel('Areas '+ r' (10$^{6}$ km$^{2}$)')
ax4.set_xticks([0,0.05,0.1,0.15,0.2,0.25,0.3])
# ax3.set_yticks(np.linspace(-5,15,5))
# ax4.set_xticks([-0.05,0,0.05,0.1,0.15,0.2])
# ax.yaxis.set_tick_params(labelsize=20)
# ax.xaxis.set_tick_params(labelsize=20)
ax4.grid(linestyle='--',alpha=0.5)
# figname = ''
# ax.set_title(figname,fontsize=22)
# ax.grid(color='lightgray',linestyle='--',linewidth=3)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax4.text(0.5, 0.2,'Pearson\'s r: '+str(round(p[0],2))+'\n'+'p value:<0.05',
        transform=ax4.transAxes, 
        verticalalignment='top', bbox=props,)
ax4.text(0, 1.13, '(d)', fontsize=12, transform=ax4.transAxes, va='top', ha='right')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
plt.subplots_adjust(hspace=0.35, wspace=0.35)
plt.savefig('C:/Users/fzjxw/python/code/Figures/Figure9.pdf' ,dpi=300,bbox_inches = 'tight')
plt.show()
