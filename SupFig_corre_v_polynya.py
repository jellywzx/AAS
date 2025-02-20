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
from scipy.stats import pearsonr

# #region 实验一
# #实验一：图（a）表示东罗斯海经向风 和 东罗斯海冰间湖面积的相关性；（b）表示东罗斯海冰间湖面积 与 次年二月东罗斯海SIE的相关性。
# wind_file = xr.open_dataset('D:/ERA5/ERA5-pressure-level-1x1-resolution-vwnd-850-197901-202211.nc')
# v = wind_file['v'].squeeze()

# lons = v.longitude
# lats = v.latitude

# v_SON = v[(v.time.dt.month==10)|(v.time.dt.month==11)].values.reshape(44,2,181,360).mean(axis=1)

# time_len = len(v_SON)
# ross_v = np.zeros((time_len),'float')
# for i in range(time_len):
#     mask_ross =  (lats<-60.)&(lats>-80.)&(lons>180.)&(lons<230.)
#     v_ross = np.where(mask_ross,v_SON[i],np.nan)
#     ross_v[i] = np.nanmean(v_ross)
    

# import os
# folder_path = 'D:/opw_rec/nsidc/new/'
# files = os.listdir(folder_path)
# files.sort()
# b = list()
# for i in files:
#     b.append(pd.read_csv(folder_path+i)['eastross']+pd.read_csv(folder_path+i)['westross'])
# b = np.array(b)
# opwa_12_1 = b[:,27:34].mean(axis=1)

# file = np.load('D:/0324/92-21_Feb_eastross_SIE_0407.npz')

# ross_sie_nt = file['ross_sie_nt']
# #去趋势操作
# # ross_sie_nt = signal.detrend(file['ross_sie_nt'])
# # opwa_12_1 = signal.detrend(opwa_12_1)
# # ross_v = signal.detrend(ross_v)

# from scipy.stats import pearsonr

# #-----------------开始作图——------------
# #-----------第一张图 ----Actual Center Pressure---------------

# from scipy.stats import pearsonr
# plt.rcParams['font.family'] = 'Arial'
# # plt.rcParams['font.size'] = 9
# # plt.rcParams['lines.linewidth'] = 1
# # plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
# # plt.rcParams['lines.markersize'] = 1 
# # plt.rc('font', family='Arial', size=9)
# fig = plt.figure(figsize=(10,6))
# ax = fig.add_subplot(111)
# #画出散点图分布，计算变量的相关性
# ax.scatter(opwa_12_1,ross_v[13:-1],marker='o',color='k')
# # ax.scatter(ross_owa_lanina,ross_v_lanina,marker='o',color='r',label='La Nina years')
# p = pearsonr(opwa_12_1,ross_v[13:-1])
# print(p)
# # plt.legend(fontsize=16,loc='lower right')
# #画出散点的拟合直线
# aa,intercept,r_value,p1,std_err = stats.linregress(opwa_12_1,ross_v[13:-1])
# y = aa*np.array(opwa_12_1)+intercept
# ax.plot(opwa_12_1,y,color='k')
# # ax.set_ylim(-0.5,1.5)
# # ax.set_xlim(-5,5)
# # ax.axhline(0,color='r')
# ax.set_ylabel('ON meridional wind (m/s)',fontsize=22)
# ax.set_xlabel('Areas '+ r' (10$^{6}$ km$^{2}$)' ,fontsize=22)
# ax.yaxis.set_tick_params(labelsize=20)
# ax.xaxis.set_tick_params(labelsize=20)
# figname = 'eastern Ross'
# # ax.set_title(figname,fontsize=18)
# ax.grid(linestyle='--',linewidth=2,alpha=0.5)
# props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# ax.text(0.6, 0.23,'pearson r: '+str(round(p[0],2))+'\n'+'p value:<0.05',
#         transform=ax.transAxes, fontsize=20,
#         verticalalignment='top', bbox=props,)
# props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# plt.savefig('C:/Users/fzjxw/python/code/Figures/eastross_wind_polynya_corre.png' ,dpi=300,bbox_inches = 'tight')
# plt.show()

#-----------------开始作图——------------
#-----------第一张图 ----Actual Center Pressure---------------
# fig = plt.figure(figsize=(16,6))
# ax = fig.add_subplot(122)
# #画出散点图分布，计算变量的相关性
# ax.scatter(opwa_12_1[:-1],ross_sie_nt[1:],marker='o',color='k')
# p = pearsonr(opwa_12_1[:-1],ross_sie_nt[1:])
# print(p)
# #画出散点的拟合直线
# aa,intercept,r_value,p1,std_err = stats.linregress(opwa_12_1[:-1],ross_sie_nt[1:])
# y = aa*np.array(opwa_12_1[:-1])+intercept
# ax.plot(opwa_12_1[:-1],y,color='k')
# ax.set_ylabel('Feb eastross SIE '+r' (10$^{6}$ km$^{2}$)',fontsize=22)
# ax.set_xlabel('Areas '+ r' (10$^{6}$ km$^{2}$)',fontsize=22)
# ax.yaxis.set_tick_params(labelsize=20)
# ax.xaxis.set_tick_params(labelsize=20)
# ax.grid(linestyle='--',linewidth=2,alpha=0.5)
# figname = ' eastern Ross'
# # ax.set_title(figname,fontsize=22)
# # ax.grid(color='lightgray',linestyle='--',linewidth=3)
# props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# ax.text(0.6, 0.90,'pearson r: '+str(round(p[0],2))+'\n'+'p value:<0.05',
#         transform=ax.transAxes, fontsize=20,
#         verticalalignment='top', bbox=props,)
# props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# # ax.text(0, 1.05, '(b)',fontsize=24, transform=ax.transAxes, va='top', ha='right')
# plt.savefig('C:/Users/fzjxw/python/code/Figures/eastross_wind_polynya_corre.png' ,dpi=300,bbox_inches = 'tight')
# plt.show()
#endregion

'''
#region 实验二：罗斯海冰间湖面积 与 （a）东罗斯海经向风 和（b）次年2月SIE 的相关性。
wind_file = xr.open_dataset('D:/ERA5/ERA5-pressure-level-1x1-resolution-vwnd-850-197901-202211.nc')
v = wind_file['v'].squeeze()

lons = v.longitude
lats = v.latitude

v_SON = v[(v.time.dt.month==10)|(v.time.dt.month==11)].values.reshape(44,2,181,360).mean(axis=1)

time_len = len(v_SON)
ross_v = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-60.)&(lats>-80.)&(lons>180.)&(lons<230.)
    v_ross = np.where(mask_ross,v_SON[i],np.nan)
    ross_v[i] = np.nanmean(v_ross)

#纬向风
wind_file = xr.open_dataset('D:/ERA5/ERA5-1x1-resolution-uv10-197901-202212.nc')
u = wind_file['u10'].squeeze()

lons = u.longitude
lats = u.latitude

u_SON = u[(u.time.dt.month==10)|(u.time.dt.month==11)].values.reshape(44,2,181,360).mean(axis=1)

time_len = len(u_SON)
ross_u = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-60.)&(lats>-80.)&(lons>180.)&(lons<230.)
    u_ross = np.where(mask_ross,u_SON[i],np.nan)
    ross_u[i] = np.nanmean(u_ross)

#读取冰间湖面积数据
import os
folder_path = 'D:/opw_rec/nsidc/new/'
files = os.listdir(folder_path)
files.sort()
b = list()
for i in files:
    b.append(pd.read_csv(folder_path+i)['eastross']+pd.read_csv(folder_path+i)['westross'])
b = np.array(b)
opwa_12_1 = b[:,27:34].mean(axis=1)


file = np.load('D:/0324/92-21_Feb_eastross_SIE_0407.npz')
ross_sie_nt = file['ross_sie_nt']

fig = plt.figure(figsize=(20,6))
ax = fig.add_subplot(121)
#画出散点图分布，计算变量的相关性
ax.scatter(opwa_12_1,ross_v[13:-1],marker='o',color='k')
# ax.scatter(ross_owa_lanina,ross_v_lanina,marker='o',color='r',label='La Nina years')
p = pearsonr(opwa_12_1,ross_v[13:-1])
print(p)
# plt.legend(fontsize=16,loc='lower right')
#画出散点的拟合直线
aa,intercept,r_value,p1,std_err = stats.linregress(opwa_12_1,ross_v[13:-1])
y = aa*np.array(opwa_12_1)+intercept
ax.plot(opwa_12_1,y,color='k')
# ax.set_ylim(-0.5,1.5)
# ax.set_xlim(-5,5)
# ax.axhline(0,color='r')
ax.set_ylabel('ON meridional wind (m/s)',fontsize=22)
ax.set_xlabel('Areas '+ r' (10$^{6}$ km$^{2}$)' ,fontsize=22)
ax.yaxis.set_tick_params(labelsize=20)
ax.xaxis.set_tick_params(labelsize=20)
figname = 'Ross'
# ax.set_title(figname,fontsize=18)
ax.grid(linestyle='--',linewidth=2,alpha=0.5)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.6, 0.23,'pearson r: '+str(round(p[0],2))+'\n'+'p value:<0.05',
        transform=ax.transAxes, fontsize=20,
        verticalalignment='top', bbox=props,)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0, 1.05, '(a)', fontsize=24,transform=ax.transAxes, va='top', ha='right')
# plt.savefig('C:/Users/fzjxw/python/code/Figures/Ross_wind_polynya_corre.png' ,dpi=300,bbox_inches = 'tight')
# plt.show()

ax = fig.add_subplot(122)
#画出散点图分布，计算变量的相关性
ax.scatter(opwa_12_1[:-1],ross_sie_nt[1:],marker='o',color='k')
p = pearsonr(opwa_12_1[:-1],ross_sie_nt[1:])
print(p)
#画出散点的拟合直线
aa,intercept,r_value,p1,std_err = stats.linregress(opwa_12_1[:-1],ross_sie_nt[1:])
y = aa*np.array(opwa_12_1[:-1])+intercept
ax.plot(opwa_12_1[:-1],y,color='k')
ax.set_ylabel('Feb eastross SIE '+r' (10$^{6}$ km$^{2}$)',fontsize=22)
ax.set_xlabel('Areas '+ r' (10$^{6}$ km$^{2}$)',fontsize=22)
ax.yaxis.set_tick_params(labelsize=20)
ax.xaxis.set_tick_params(labelsize=20)
ax.grid(linestyle='--',linewidth=2,alpha=0.5)
figname = ' eastern Ross'
# ax.set_title(figname,fontsize=22)
# ax.grid(color='lightgray',linestyle='--',linewidth=3)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.6, 0.90,'pearson r: '+str(round(p[0],2))+'\n'+'p value:<0.05',
        transform=ax.transAxes, fontsize=20,
        verticalalignment='top', bbox=props,)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0, 1.05, '(b)',fontsize=24, transform=ax.transAxes, va='top', ha='right')
plt.savefig('C:/Users/fzjxw/python/code/Figures/SupFig3.png' ,dpi=300,bbox_inches = 'tight')
plt.show()

#endregion
'''

'''
#region 其他补充内容
#补充计算纬向风和冰间湖面积的相关性
wind_file = xr.open_dataset('D:/ERA5/ERA5-1x1-resolution-uv10-197901-202212.nc')
u = wind_file['u10'].squeeze()

lons = u.longitude
lats = u.latitude

u_SON = u[(u.time.dt.month==10)|(u.time.dt.month==11)].values.reshape(44,2,181,360).mean(axis=1)

time_len = len(u_SON)
ross_u = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-60.)&(lats>-80.)&(lons>180.)&(lons<230.)
    u_ross = np.where(mask_ross,u_SON[i],np.nan)
    ross_u[i] = np.nanmean(u_ross)

import os
folder_path = 'D:/opw_rec/nsidc/new/'
files = os.listdir(folder_path)
files.sort()
b = list()
for i in files:
    b.append(pd.read_csv(folder_path+i)['eastross']+pd.read_csv(folder_path+i)['westross'])
b = np.array(b)
opwa_12_1 = b[:,27:34].mean(axis=1)

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
#画出散点图分布，计算变量的相关性
ax.scatter(opwa_12_1,ross_u[13:-1],marker='o',color='k')
# ax.scatter(ross_owa_lanina,ross_v_lanina,marker='o',color='r',label='La Nina years')
p = pearsonr(opwa_12_1,ross_u[13:-1])
print(p)
# plt.legend(fontsize=16,loc='lower right')
#画出散点的拟合直线
aa,intercept,r_value,p1,std_err = stats.linregress(opwa_12_1,ross_u[13:-1])
y = aa*np.array(opwa_12_1)+intercept
ax.plot(opwa_12_1,y,color='k')
# ax.set_ylim(-0.5,1.5)
# ax.set_xlim(-5,5)
# ax.axhline(0,color='r')
ax.set_ylabel('ON zonal wind (m/s)',fontsize=22)
ax.set_xlabel('Areas '+ r' (10$^{6}$ km$^{2}$)' ,fontsize=22)
ax.yaxis.set_tick_params(labelsize=20)
ax.xaxis.set_tick_params(labelsize=20)
figname = 'Ross'
# ax.set_title(figname,fontsize=18)
ax.grid(linestyle='--',linewidth=2,alpha=0.5)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.6, 0.23,'pearson r: '+str(round(p[0],2))+'\n'+'p value:<0.05',
        transform=ax.transAxes, fontsize=20,
        verticalalignment='top', bbox=props,)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# ax.text(0, 1.05, '(a)', fontsize=24,transform=ax.transAxes, va='top', ha='right')
plt.savefig('C:/Users/fzjxw/python/code/Figures/uwind_polynya_corre.png' ,dpi=300,bbox_inches = 'tight')
plt.show()

#经向风的相关性单独出图
wind_file = xr.open_dataset('D:/ERA5/ERA5-pressure-level-1x1-resolution-vwnd-850-197901-202211.nc')
v = wind_file['v'].squeeze()

lons = v.longitude
lats = v.latitude

v_SON = v[(v.time.dt.month==10)|(v.time.dt.month==11)].values.reshape(44,2,181,360).mean(axis=1)

time_len = len(v_SON)
ross_v = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-60.)&(lats>-80.)&(lons>180.)&(lons<230.)
    v_ross = np.where(mask_ross,v_SON[i],np.nan)
    ross_v[i] = np.nanmean(v_ross)

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
#画出散点图分布，计算变量的相关性
ax.scatter(opwa_12_1,ross_v[13:-1],marker='o',color='k')
# ax.scatter(ross_owa_lanina,ross_v_lanina,marker='o',color='r',label='La Nina years')
p = pearsonr(opwa_12_1,ross_v[13:-1])
print(p)
# plt.legend(fontsize=16,loc='lower right')
#画出散点的拟合直线
aa,intercept,r_value,p1,std_err = stats.linregress(opwa_12_1,ross_v[13:-1])
y = aa*np.array(opwa_12_1)+intercept
ax.plot(opwa_12_1,y,color='k')
# ax.set_ylim(-0.5,1.5)
# ax.set_xlim(-5,5)
# ax.axhline(0,color='r')
ax.set_ylabel('ON meridional wind (m/s)',fontsize=22)
ax.set_xlabel('Areas '+ r' (10$^{6}$ km$^{2}$)' ,fontsize=22)
ax.yaxis.set_tick_params(labelsize=20)
ax.xaxis.set_tick_params(labelsize=20)
figname = 'Ross'
# ax.set_title(figname,fontsize=18)
ax.grid(linestyle='--',linewidth=2,alpha=0.5)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.6, 0.23,'pearson r: '+str(round(p[0],2))+'\n'+'p value:<0.05',
        transform=ax.transAxes, fontsize=20,
        verticalalignment='top', bbox=props,)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# ax.text(0, 1.05, '(a)', fontsize=24,transform=ax.transAxes, va='top', ha='right')
plt.savefig('C:/Users/fzjxw/python/code/Figures/vwind_polynya_corre.png' ,dpi=300,bbox_inches = 'tight')
plt.show()

#endregion 
'''
wind_file = xr.open_dataset('D:/ERA5/ERA5-pressure-level-1x1-resolution-vwnd-850-197901-202211.nc')
v = wind_file['v'].squeeze()

lons = v.longitude
lats = v.latitude

v_SON = v[(v.time.dt.month==10)|(v.time.dt.month==11)].values.reshape(44,2,181,360).mean(axis=1)

time_len = len(v_SON)
ross_v = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-60.)&(lats>-80.)&(lons>180.)&(lons<230.)
    v_ross = np.where(mask_ross,v_SON[i],np.nan)
    ross_v[i] = np.nanmean(v_ross)

#纬向风
wind_file = xr.open_dataset('D:/ERA5/ERA5-1x1-resolution-uv10-197901-202212.nc')
u = wind_file['u10'].squeeze()

lons = u.longitude
lats = u.latitude

u_SON = u[(u.time.dt.month==10)|(u.time.dt.month==11)].values.reshape(44,2,181,360).mean(axis=1)

time_len = len(u_SON)
ross_u = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-60.)&(lats>-80.)&(lons>180.)&(lons<230.)
    u_ross = np.where(mask_ross,u_SON[i],np.nan)
    ross_u[i] = np.nanmean(u_ross)

#读取冰间湖面积数据
import os
folder_path = 'D:/opw_rec/nsidc/new/'
files = os.listdir(folder_path)
files.sort()
b = list()
for i in files:
    b.append(pd.read_csv(folder_path+i)['eastross']+pd.read_csv(folder_path+i)['westross'])
b = np.array(b)
opwa_12_1 = b[:,27:34].mean(axis=1)


file = np.load('D:/0324/92-21_Feb_eastross_SIE_0407.npz')
ross_sie_nt = file['ross_sie_nt']

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
#画出散点图分布，计算变量的相关性
ax.scatter(opwa_12_1,ross_v[13:-1],marker='o',color='k')
# ax.scatter(ross_owa_lanina,ross_v_lanina,marker='o',color='r',label='La Nina years')
p = pearsonr(opwa_12_1,ross_v[13:-1])
print(p)
# plt.legend(fontsize=16,loc='lower right')
#画出散点的拟合直线
aa,intercept,r_value,p1,std_err = stats.linregress(opwa_12_1,ross_v[13:-1])
y = aa*np.array(opwa_12_1)+intercept
ax.plot(opwa_12_1,y,color='k')
# ax.set_ylim(-0.5,1.5)
# ax.set_xlim(-5,5)
# ax.axhline(0,color='r')
ax.set_ylabel('ON v wind in eastern Ross (m/s)',fontsize=22)
ax.set_xlabel('Areas '+ r' (10$^{6}$ km$^{2}$)' ,fontsize=22)
ax.yaxis.set_tick_params(labelsize=20)
ax.xaxis.set_tick_params(labelsize=20)
figname = 'Ross'
# ax.set_title(figname,fontsize=18)
ax.grid(linestyle='--',linewidth=2,alpha=0.5)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.6, 0.23,'pearson r: '+str(round(p[0],2))+'\n'+'p value:<0.05',
        transform=ax.transAxes, fontsize=20,
        verticalalignment='top', bbox=props,)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0, 1.05, '(a)', fontsize=24,transform=ax.transAxes, va='top', ha='right')
plt.savefig('C:/Users/fzjxw/python/code/Figures/minor_eastwind_polynya_corre.png' ,dpi=300,bbox_inches = 'tight')
# plt.show()

#region 补充西罗斯海经向风与冰间湖面积的相关性
wind_file = xr.open_dataset('D:/ERA5/ERA5-pressure-level-1x1-resolution-vwnd-850-197901-202211.nc')
v = wind_file['v'].squeeze()

lons = v.longitude
lats = v.latitude

v_SON = v[(v.time.dt.month==10)|(v.time.dt.month==11)].values.reshape(44,2,181,360).mean(axis=1)

#西罗斯海的范围是160-180
time_len = len(v_SON)
ross_v = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-60.)&(lats>-80.)&(lons>160.)&(lons<180.)
    v_ross = np.where(mask_ross,v_SON[i],np.nan)
    ross_v[i] = np.nanmean(v_ross)

#整个罗斯海冰间湖面积
import os
folder_path = 'D:/opw_rec/nsidc/new/'
files = os.listdir(folder_path)
files.sort()
b = list()
for i in files:
    b.append(pd.read_csv(folder_path+i)['eastross']+pd.read_csv(folder_path+i)['westross'])
b = np.array(b)
opwa_12_1 = b[:,27:34].mean(axis=1)


file = np.load('D:/0324/92-21_Feb_eastross_SIE_0407.npz')
ross_sie_nt = file['ross_sie_nt']

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
#画出散点图分布，计算变量的相关性
ax.scatter(opwa_12_1,ross_v[13:-1],marker='o',color='k')
# ax.scatter(ross_owa_lanina,ross_v_lanina,marker='o',color='r',label='La Nina years')
p = pearsonr(opwa_12_1,ross_v[13:-1])
print(p)
# plt.legend(fontsize=16,loc='lower right')
#画出散点的拟合直线
aa,intercept,r_value,p1,std_err = stats.linregress(opwa_12_1,ross_v[13:-1])
y = aa*np.array(opwa_12_1)+intercept
ax.plot(opwa_12_1,y,color='k')
# ax.set_ylim(-0.5,1.5)
# ax.set_xlim(-5,5)
# ax.axhline(0,color='r')
ax.set_ylabel('ON v in western Ross (m/s)',fontsize=22)
ax.set_xlabel('Areas '+ r' (10$^{6}$ km$^{2}$)' ,fontsize=22)
ax.yaxis.set_tick_params(labelsize=20)
ax.xaxis.set_tick_params(labelsize=20)
#figname = 'Ross'
# ax.set_title(figname,fontsize=18)
ax.grid(linestyle='--',linewidth=2,alpha=0.5)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.6, 0.23,'pearson r: '+str(round(p[0],2)),
        transform=ax.transAxes, fontsize=20,
        verticalalignment='top', bbox=props,)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0, 1.05, '(b)', fontsize=24,transform=ax.transAxes, va='top', ha='right')
plt.savefig('C:/Users/fzjxw/python/code/Figures/minor_westwind_polynya_corre.png' ,dpi=300,bbox_inches = 'tight')
plt.show()

#——————————————————————————————————————————
#整个罗斯海经向风和冰间湖面积的现相关性
wind_file = xr.open_dataset('D:/ERA5/ERA5-pressure-level-1x1-resolution-vwnd-850-197901-202211.nc')
v = wind_file['v'].squeeze()

lons = v.longitude
lats = v.latitude

v_SON = v[(v.time.dt.month==10)|(v.time.dt.month==11)].values.reshape(44,2,181,360).mean(axis=1)

#西罗斯海的范围是160-180
time_len = len(v_SON)
ross_v = np.zeros((time_len),'float')
for i in range(time_len):
    mask_ross =  (lats<-60.)&(lats>-80.)&(lons>160.)&(lons<230.)
    v_ross = np.where(mask_ross,v_SON[i],np.nan)
    ross_v[i] = np.nanmean(v_ross)

#整个罗斯海冰间湖面积
import os
folder_path = 'D:/opw_rec/nsidc/new/'
files = os.listdir(folder_path)
files.sort()
b = list()
for i in files:
    b.append(pd.read_csv(folder_path+i)['eastross']+pd.read_csv(folder_path+i)['westross'])
b = np.array(b)
opwa_12_1 = b[:,27:34].mean(axis=1)


file = np.load('D:/0324/92-21_Feb_eastross_SIE_0407.npz')
ross_sie_nt = file['ross_sie_nt']

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
#画出散点图分布，计算变量的相关性
ax.scatter(opwa_12_1,ross_v[13:-1],marker='o',color='k')
# ax.scatter(ross_owa_lanina,ross_v_lanina,marker='o',color='r',label='La Nina years')
p = pearsonr(opwa_12_1,ross_v[13:-1])
print(p)
# plt.legend(fontsize=16,loc='lower right')
#画出散点的拟合直线
aa,intercept,r_value,p1,std_err = stats.linregress(opwa_12_1,ross_v[13:-1])
y = aa*np.array(opwa_12_1)+intercept
ax.plot(opwa_12_1,y,color='k')
# ax.set_ylim(-0.5,1.5)
# ax.set_xlim(-5,5)
# ax.axhline(0,color='r')
ax.set_ylabel('ON v in Ross Sea (m/s)',fontsize=22)
ax.set_xlabel('Areas '+ r' (10$^{6}$ km$^{2}$)' ,fontsize=22)
ax.yaxis.set_tick_params(labelsize=20)
ax.xaxis.set_tick_params(labelsize=20)
#figname = 'Ross'
# ax.set_title(figname,fontsize=18)
ax.grid(linestyle='--',linewidth=2,alpha=0.5)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.6, 0.23,'pearson r: '+str(round(p[0],2))+'\n'+'p value:<0.05',
        transform=ax.transAxes, fontsize=20,
        verticalalignment='top', bbox=props,)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0, 1.05, '(c)', fontsize=24,transform=ax.transAxes, va='top', ha='right')
plt.savefig('C:/Users/fzjxw/python/code/Figures/minor_allwind_polynya_corre.png' ,dpi=300,bbox_inches = 'tight')
plt.show()
