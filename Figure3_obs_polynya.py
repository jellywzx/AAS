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


#——————————空间分布图————————————————————————

sic = xr.open_dataset('C:/Users/fzjxw/python/data/nt_seaice_78-22.nc')['nsidc_nt_seaice_conc'].sel(time=slice('1992','2021'))
ftools='C:/Users/fzjxw/python/data/'
with open(ftools+'pss25lats_v3.dat','rb') as flat:
    lats = np.fromfile(flat,dtype='<i4').reshape(332, 316)/100000.

with open(ftools+'pss25lons_v3.dat','rb') as flon:
    lons = np.fromfile(flon,dtype='<i4').reshape(332, 316)/100000.

with open(ftools+'pss25area_v3.dat','rb') as flon:
    area = np.fromfile(flon,dtype='<i4').reshape(332, 316)*1e-9

east_ross_mask = ((lons>-180.)&(lons<-130.)&(lats<-70.))
west_ross_mask = ((lons>160.)&(lats<-68.))
amun_mask = ((lons>-130.)&(lons<-90.)&(lats<-70.))
bell_mask = ((lons>-90.)&(lons<-60.)&(lats<-70.))

from global_land_mask import globe
globe_land_mask = globe.is_land(lats, lons)

year = np.linspace(1992,2021,30)
sic_list = []
for i in range(30):
    start_time = str(int(year[i]))+'-11-28'
    end_time = str(int(year[i]))+'-12-04'
    sic_list.append(sic.sel(time=slice(start_time,end_time)).mean('time'))
sic_7days_avg = np.nanmean(np.array(sic_list),axis=0)

year_min = [1994., 1997., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2010., 2016., 2018., 2021.]
num_y = 6

sic_list = []
for i in range(6):
    start_time = str(int(year_max[i]))+'-11-28'
    end_time = str(int(year_max[i]))+'-12-04'
    sic_list.append(sic.sel(time=slice(start_time,end_time)).mean('time'))
sic_max = np.nanmean(np.array(sic_list),axis=0)

num_y = 6
sic_list = []
for i in range(6):
    start_time = str(int(year_min[i]))+'-11-28'
    end_time = str(int(year_min[i]))+'-12-04'
    sic_list.append(sic.sel(time=slice(start_time,end_time)).mean('time'))
sic_min = np.nanmean(np.array(sic_list),axis=0)


#设置全局变量
plt.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.8
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 9
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
levels1 = np.linspace(0,100,21)
#----------图1----------------
fig=plt.figure(figsize=(160/25.4,6.5))
ax1 = fig.add_subplot(331)
levels1 = np.linspace(0,100,21)
plt.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.8
figname = '(a) Climatology'
ax1.set_title(figname)
m = Basemap(projection='stere',resolution='c',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
x1,y1 = m(lons,lats)
im1 = m.contourf(x1,y1,sic_7days_avg*100,
                     levels=levels1,
                     extend='both',
                    shading='faceted', antialiased=True,cmap='Blues')
# im2 = m.contourf(x1,y1,sic_sel,levels=[0.15,1],hatches=['///', None],alpha=0)
m.contour(x1,y1,sic_7days_avg,levels=[0.15],colors='green')
# xpt1,ypt1 = m(218,-66)
# ax1.text(xpt1,ypt1,'Eastern\nRoss',color='white',fontsize=7)
# xpt2,ypt2 = m(191,-67)
# ax1.text(xpt2,ypt2,'Western\nRoss ',color='white',fontsize=7)
# xpt3,ypt3 = m(250,-63)
# ax1.text(xpt3,ypt3,'Amundsen',fontsize=20)
# xpt4,ypt4 = m(273,-61.5)
# ax1.text(xpt4,ypt4,'Bellingshausen ',fontsize=20)
m.drawcoastlines(color='black')
m.fillcontinents(color='gray')
m.drawmeridians([90,120,150,180,230,270,300,330],labels=[True,False,False,True],fmt='%g',latmax=80,)
m.drawparallels([-60,-70,-80],labels=[False,False,False,True],latmax=80,)

#----------图2------------------
ax2 = fig.add_subplot(332)
levels1 = np.linspace(0,100,21)
plt.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.8
figname = '(b) Large'
ax2.set_title(figname,)
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
x1,y1 = m(lons,lats)
im1 = m.contourf(x1,y1,sic_max*100,
                     levels=levels1,
                     extend='both',
                    shading='faceted', antialiased=True,cmap='Blues')
# im2 = m.contourf(x1,y1,sic_sel,levels=[0.15,1],hatches=['///', None],alpha=0)
m.contour(x1,y1,sic_max,levels=[0.15],colors='green',)
# xpt1,ypt1 = m(218,-67)
# ax2.text(xpt1,ypt1,'Eastern',color='white')
# xpt2,ypt2 = m(182,-68)
# ax2.text(xpt2,ypt2,'Western',color='white')
# xpt3,ypt3 = m(250,-63)
# ax2.text(xpt3,ypt3,'Amundsen',fontsize=20)
#  xpt4,ypt4 = m(273,-61.5)
# ax2.text(xpt4,ypt4,'Bellingshausen ',fontsize=20)
m.drawcoastlines(color='black')
m.fillcontinents(color='gray')

m.drawmeridians([90,120,150,180,230,270,300,330], linewidth=1.2,labels=[False,False,False,True],fmt='%g',latmax=80)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,True],latmax=80)

#----------图3------------------
ax3 = fig.add_subplot(333)
levels1 = np.linspace(0,100,21)
plt.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.8
figname = '(c) Small'
ax3.set_title(figname)
m = Basemap(projection='stere',resolution='c',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
x1,y1 = m(lons,lats)
im1 = m.contourf(x1,y1,sic_min*100,
                     levels=levels1,
                     extend='both',
                    shading='faceted', antialiased=True,cmap='Blues')
# im2 = m.contourf(x1,y1,sic_sel,levels=[0.15,1],hatches=['///', None],alpha=0)
m.contour(x1,y1,sic_min,levels=[0.15],colors='green')
# xpt1,ypt1 = m(218,-67)
# ax3.text(xpt1,ypt1,'Eastern',color='white',)
# xpt2,ypt2 = m(182,-68)
# ax3.text(xpt2,ypt2,'Western',color='white')
# xpt3,ypt3 = m(250,-63)
# ax3.text(xpt3,ypt3,'Amundsen',fontsize=20)
# xpt4,ypt4 = m(273,-61.5)
# ax3.text(xpt4,ypt4,'Bellingshausen ',fontsize=20)
m.drawcoastlines(color='black')
m.fillcontinents(color='gray')
m.drawmeridians([90,120,150,180,230,270,300,330], linewidth=1.2,labels=[False,False,False,True],fmt='%g',latmax=80)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,True,False,True],latmax=80)
# cb = m.colorbar(im1, location='bottom', pad="5%",ticks =levels1,
# #                 boundaries=levels1[::2],
#                 extend='both', extendfrac='auto',drawedges=True, )
# cb.set_label('Trend in Sea ice concentration per decade(%)',fontsize=16)
# cb.ax.tick_params(labelsize=16, width=0)
cax1 = fig.add_axes([0.1, 0.48, 0.8, 0.1],aspect=0.02)
cb1 = fig.colorbar(im1,orientation='horizontal',
                  ticks =levels1[::2], 
                  boundaries=levels1[::2], 
                  extend='False', extendfrac=0.05, drawedges=False,
                cax = cax1
                 )
cb1.set_label('Sea ice concentration (%)',)
# cb.ax.tick_params(width=0)
# fig.tight_layout()
# plt.savefig('/stu02/weizx24/figures/0924/Figure3_SIC.png',dpi=600,bbox_inches='tight')
# plt.show()


#——————————————线形图————————————————————
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 9
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
# 仅罗斯海冰间湖面积
import os
folder_path = 'C:/Users/fzjxw/python/data/opw_rec/nsidc/new/'
files = os.listdir(folder_path)
files.sort()
# files = np.delete(files,[-1])

a = list()
for i in files:
#     prini(i)
    a.append(pd.read_csv(folder_path+i)['eastross']+pd.read_csv(folder_path+i)['westross'])
a = np.array(a)

years = np.linspace(1992,2021,30)
year_min = [1994., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2016., 2018., 2021.]
bool_max = np.isclose(years[:, None], year_max).any(axis=1)
bool_min = np.isclose(years[:, None], year_min).any(axis=1)

a_max = a[bool_max].mean(axis=0)
a_min = a[bool_min].mean(axis=0)
a_clm = a.mean(axis=0)

# print('开始画图')
# plt.close()
date_name = pd.date_range('1979-11-01','1979-12-31')
date_list = list()
for i in range(len(date_name)):
    date_list.append(str(date_name[i])[5:10])

# fig = plt.figure(1, figsize=(160/25.4,3))
ax4 = fig.add_subplot(337)
ax4.plot(date_list[:-5],a_max[:-5],color='r',label='Large ')
ax4.plot(date_list[:-5],a_clm[:-5],color='k',label='Climatology ')
ax4.plot(date_list[:-5],a_min[:-5],color='b',label='Small ')
ax4.set_xlabel('Date')
ax4.set_ylabel('Areas'+ r' (10$^{6}$ km$^{2}$)',)
ax4.set_xticks(date_list[:-5][::5])
ax4.set_yticks([0,0.2,0.4,0.6])
ax4.grid(color='lightgray',linestyle='--',alpha=0.4)
ax4.legend(edgecolor='k',loc='upper left')
ax4.text(0, 1.05, '(d)', fontsize=10, transform=ax4.transAxes, va='top', ha='right')
plt.subplots_adjust(hspace=0.35)
plt.savefig('C:/Users/fzjxw/python/code/Figures/Figure3_all.png',dpi=300,bbox_inches='tight')
