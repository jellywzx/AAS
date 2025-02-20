
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
import os

#设置全局变量
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 9
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6
'''
#观测数据
ASL_lon = pd.read_excel('/stu02/weizx24/data/ASL/ASL_SON.xlsx')['lon']
ASL_lon_lst = list()
n = 2
for i in range(43):
    ASL_lon_lst.append((ASL_lon[3*i+1]+ASL_lon[3*i+2])/n)
ASL_lon_arr = np.array(ASL_lon_lst)[13:]

folder_path = '/stu02/weizx24/data/opw_rec/nsidc/new/'
files = os.listdir(folder_path)
files.sort()
files
a = list()
for i in files:
#     prini(i)
    a.append(pd.read_csv(folder_path+i)['eastross']+pd.read_csv(folder_path+i)['westross'])
a = np.array(a)

opwa_12_1 = a[:,27:34].mean(axis=1)
# opwa = signal.detrend(opwa_12_1)
# asl_lon = signal.detrend(ASL_lon_arr)

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
ax.scatter(opwa_12_1,ASL_lon_arr,marker='o',color='k',s=20)
p = pearsonr(opwa_12_1,ASL_lon_arr)
aa,intercept,r_value,p1,std_err = stats.linregress(opwa_12_1,ASL_lon_arr)
print(p)
y = aa*np.array(opwa_12_1)+intercept
ax.plot(opwa_12_1,y,color='k')
ax.set_ylabel('Longitude of ASL ',fontsize=20)
ax.set_xlabel('Areas'+ r' (10$^{6}$ km$^{2}$)',fontsize=20)
ax.yaxis.set_tick_params(labelsize=20)
ax.xaxis.set_tick_params(labelsize=20)
ax.grid(linestyle='--',linewidth=2,alpha=0.5)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.6, 0.2,'pearson r: '+str(round(p[0],2))+'\n'+'p value < 0.05',
        transform=ax.transAxes, fontsize=20,
        verticalalignment='top', bbox=props,)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
plt.tight_layout(rect=[0,0,0.9,1])
plt.savefig('/stu02/weizx24/figures/0924/majorrevison_ASL.png' ,dpi=300,bbox_inches = 'tight')
# plt.show()
'''
#————————————————————————————————————————————————————————
#计算latitude的相关性
ASL_lat = pd.read_excel('/stu02/weizx24/data/ASL/ASL_SON.xlsx')['lat']
ASL_lat_lst = list()
n = 2
for i in range(43):
    ASL_lat_lst.append((ASL_lat[3*i+1]+ASL_lat[3*i+2])/n)
ASL_lat_arr = np.array(ASL_lat_lst)[13:]

folder_path = '/stu02/weizx24/data/opw_rec/nsidc/new/'
files = os.listdir(folder_path)
files.sort()
files
a = list()
for i in files:
#     prini(i)
    a.append(pd.read_csv(folder_path+i)['eastross']+pd.read_csv(folder_path+i)['westross'])
a = np.array(a)

opwa_12_1 = a[:,27:34].mean(axis=1)
# opwa = signal.detrend(opwa_12_1)
# asl_lat = signal.detrend(ASL_lat_arr)

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
ax.scatter(opwa_12_1,ASL_lat_arr,marker='o',color='k',s=20)
p = pearsonr(opwa_12_1,ASL_lat_arr)
aa,intercept,r_value,p1,std_err = stats.linregress(opwa_12_1,ASL_lat_arr)
print(p)
y = aa*np.array(opwa_12_1)+intercept
ax.plot(opwa_12_1,y,color='k')
ax.set_ylabel('Latitude of ASL ',fontsize=20)
ax.set_xlabel('Areas'+ r' (10$^{6}$ km$^{2}$)',fontsize=20)
ax.yaxis.set_tick_params(labelsize=20)
ax.xaxis.set_tick_params(labelsize=20)
ax.grid(linestyle='--',linewidth=2,alpha=0.5)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.6, 0.2,'pearson r: '+str(round(p[0],2))+'\n'+'p value < 0.05',
        transform=ax.transAxes, fontsize=20,
        verticalalignment='top', bbox=props,)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
plt.tight_layout(rect=[0,0,0.9,1])
plt.savefig('/stu02/weizx24/figures/0924/majorrevison_ASL_lat.png' ,dpi=300,bbox_inches = 'tight')
plt.show()