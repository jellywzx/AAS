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
from pyproj import Proj, transform
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import AxesGrid

path_area_weight ='D:/CMIP6/CESM2-WACCM-FV2/areacello_Ofx_CESM2-WACCM_piControl_r1i1p1f1_gn.nc'
areacella= xr.open_dataset(path_area_weight)
aw_xr = areacella['areacello']

# CMIP6_opwa = np.load('/stu02/weizx24/data/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']
# CMIP6_opwa_7day_mean = np.nanmean(np.array(CMIP6_opwa),axis=1)

CMIP6_opwa = np.load('D:/npz/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']*1e-12
#时间维度上平均得到7天平均的冰间湖面积
ross_sie_cmip6 = np.nanmean(CMIP6_opwa,axis=1)
# print(ross_sie_cmip6.shape)
Q1 = np.quantile(ross_sie_cmip6,0.05)
Q3 = np.quantile(ross_sie_cmip6,0.95)
sic_Feb = xr.open_dataset('D:/CMIP6/CESM2-WACCM-FV2/piControl/SImon/siconc/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_Feb.nc')['siconc']
year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year
# print(f"最小年份: {year_min} ")
# print(f"最大年份: {year_max} ")

sic =  xr.open_dataset('D:/CMIP6/CESM2-WACCM-FV2/piControl/SIday/1128/siconc_SIday_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_1128-1204_all.nc')['siconc']
sic_rsp = sic.values.reshape(499,7,384,320).mean(axis=1)
sic_clm = np.nanmean(sic_rsp,axis=0)
lats = sic.lat
lons = sic.lon

sic_max = np.nanmean(sic_rsp[ross_sie_cmip6>Q3],axis=0)
sic_min = np.nanmean(sic_rsp[ross_sie_cmip6<Q1],axis=0)

#region 画图
levels1 = np.linspace(0,100,21)
# 使用basemap画图
# #----------图1----------------
fig=plt.figure(figsize=(20,10))
ax = fig.add_subplot(131)
levels1 = np.linspace(0,100,21)

figname = '(a) Climatology'
ax.set_title(figname,fontsize=24)
m = Basemap(projection='stere',resolution='c',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
x1,y1 = m(lons,lats)
im1 = m.contourf(x1,y1,sic_clm,
                     levels=levels1,
                     extend='both',
                    shading='faceted', antialiased=True,cmap='Blues')
# im2 = m.contourf(x1,y1,sic_sel,levels=[0.15,1],hatches=['///', None],alpha=0)
m.contour(x1,y1,sic_clm,levels=[15],colors='green',linewidths=3)
xpt1,ypt1 = m(218,-67)
ax.text(xpt1,ypt1,'Eastern Ross ',fontsize=20,color='white')
xpt2,ypt2 = m(182,-68)
ax.text(xpt2,ypt2,'Western \n Ross ',fontsize=19,color='white')
# xpt3,ypt3 = m(250,-63)
# ax.text(xpt3,ypt3,'Amundsen',fontsize=20)
# xpt4,ypt4 = m(273,-61.5)
# ax.text(xpt4,ypt4,'Bellingshausen ',fontsize=20)
m.drawcoastlines(color='black')
m.fillcontinents(color='gray')
m.drawmeridians([90,120,150,180,230,270,300,330], linewidth=1.2,labels=[True,False,False,True],fmt='%g',latmax=80,fontsize=18)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,True],latmax=80,fontsize=18)

# #----------图2------------------
ax = fig.add_subplot(132)
levels1 = np.linspace(0,100,21)
figname = '(b) Large'
ax.set_title(figname,fontsize=24)
# m = Basemap(projection='stere',resolution='h',
#             lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
#             llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
x1,y1 = m(lons,lats)
im1 = m.contourf(x1,y1,sic_max,
                     levels=levels1,
                     extend='both',
                    shading='faceted', antialiased=True,cmap='Blues')
# im2 = m.contourf(x1,y1,sic_sel,levels=[0.15,1],hatches=['///', None],alpha=0)
m.contour(x1,y1,sic_max,levels=[15],colors='green',linewidths=3)
xpt1,ypt1 = m(218,-67)
ax.text(xpt1,ypt1,'Eastern Ross ',fontsize=20,color='white')
xpt2,ypt2 = m(182,-68)
ax.text(xpt2,ypt2,'Western \n Ross ',fontsize=19,color='white')
# xpt3,ypt3 = m(250,-63)
# ax.text(xpt3,ypt3,'Amundsen',fontsize=20)
# xpt4,ypt4 = m(273,-61.5)
# ax.text(xpt4,ypt4,'Bellingshausen ',fontsize=20)
m.drawcoastlines(color='black')
m.fillcontinents(color='gray')

m.drawmeridians([90,120,150,180,230,270,300,330], linewidth=1.2,labels=[False,False,False,True],fmt='%g',latmax=80,fontsize=18)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,True],latmax=80,fontsize=18)

# #----------图3------------------
ax = fig.add_subplot(133)
levels1 = np.linspace(0,100,21)
figname = '(c) Small'
ax.set_title(figname,fontsize=24)
# m = Basemap(projection='stere',resolution='h',
#             lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
#             llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
x1,y1 = m(lons,lats)
im1 = m.contourf(x1,y1,sic_min,
                     levels=levels1,
                     extend='both',
                    shading='faceted', antialiased=True,cmap='Blues')
# im2 = m.contourf(x1,y1,sic_sel,levels=[0.15,1],hatches=['///', None],alpha=0)
m.contour(x1,y1,sic_min,levels=[15],colors='green',linewidths=3)
xpt1,ypt1 = m(218,-67)
ax.text(xpt1,ypt1,'Eastern Ross ',fontsize=20,color='white',)
xpt2,ypt2 = m(182,-68)
ax.text(xpt2,ypt2,'Western \n Ross ',fontsize=19,color='white')
# xpt3,ypt3 = m(250,-63)
# ax.text(xpt3,ypt3,'Amundsen',fontsize=20)
# xpt4,ypt4 = m(273,-61.5)
# ax.text(xpt4,ypt4,'Bellingshausen ',fontsize=20)
m.drawcoastlines(color='black')
m.fillcontinents(color='gray')
m.drawmeridians([90,120,150,180,230,270,300,330], linewidth=1.2,labels=[False,False,False,True],fmt='%g',latmax=80,fontsize=18)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,True,False,True],latmax=80,fontsize=18)



# cb = m.colorbar(im1, location='bottom', pad="5%",ticks =levels1,
# #                 boundaries=levels1[::2],
#                 extend='both', extendfrac='auto',drawedges=True, )
# cb.set_label('Trend in Sea ice concentration per decade(%)',fontsize=16)
# cb.ax.tick_params(labelsize=16, width=0)


cax = fig.add_axes([0.1, 0.02, 0.8, 0.1],aspect=0.02)
cb = fig.colorbar(im1,orientation='horizontal',
                  ticks =levels1[::2], 
                  boundaries=levels1[::2], 
                  extend='False', extendfrac=0.05, drawedges=False,
                cax = cax
                 )
cb.set_label('Sea ice concentration (%)',fontsize=20)
cb.ax.tick_params(labelsize=20, width=0)

fig.tight_layout()
plt.savefig('C:/Users/fzjxw/python/code/Figures/SupFIg2_CMIP6polynya.png',dpi=300,bbox_inches='tight')
# plt.show()
#endregion

