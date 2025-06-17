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
from matplotlib.patches import Patch

#region 观测
sic = xr.open_dataset('/stu02/weizx24/data/Daily_SIC/process/nt_seaice_78-22.nc')['nsidc_nt_seaice_conc']
ftools='/stu02/weizx24/data//tools/'
with open(ftools+'pss25lats_v3.dat','rb') as flat:
    lats = np.fromfile(flat,dtype='<i4').reshape(332, 316)/100000.

with open(ftools+'pss25lons_v3.dat','rb') as flon:
    lons = np.fromfile(flon,dtype='<i4').reshape(332, 316)/100000.

with open(ftools+'pss25area_v3.dat','rb') as flon:
    area = np.fromfile(flon,dtype='<i4').reshape(332, 316)*1e-9
east_ross_mask = (((lons>-180.)&(lons<-155.)&(lats<-70.))|((lons>-155.)&(lons<-130.)&(lats<-73.5)))
west_ross_mask = ((lons>160.)&(lats<-68.))
amun_mask = ((lons>-130.)&(lons<-90.)&(lats<-72.))
bell_mask = ((lons>-90.)&(lons<-60.)&(lats<-72.))
from global_land_mask import globe
globe_land_mask = globe.is_land(lats, lons)

sic_sel = sic.sel(time='2021-12-01')
mask_ross_nt = (west_ross_mask|east_ross_mask|amun_mask|bell_mask)&(sic_sel<0.15)
sic_sel_ross_nt = np.where(mask_ross_nt,area,np.nan)


plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 9
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.8
plt.rcParams['lines.markersize']=10

fig=plt.figure(figsize=(150/25.4,6))
ax = fig.add_subplot(121)
levels1 = np.linspace(0,100,21)
figname = '(a) Observation, 2021-12-01'
ax.set_title(figname)
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72,)
x1,y1 = m(lons,lats)
im1 = m.contourf(x1,y1,sic_sel_ross_nt*100,
                     levels=levels1,
                     extend='both',
                    shading='faceted', antialiased=True,cmap='RdBu_r')
im2 = m.contourf(x1,y1,sic_sel,levels=[0.15,1],hatches=['///', None],alpha=0)
m.contour(x1,y1,sic_sel,levels=[0.15],colors='g')

xpt1,ypt1 = m(223,-69.5)
ax.text(xpt1,ypt1,'Eastern Ross ')
xpt2,ypt2 = m(182,-66)
ax.text(xpt2,ypt2,'Western \n Ross ',fontsize=8)
xpt3,ypt3 = m(250,-63)
ax.text(xpt3,ypt3,'Amundsen')
xpt4,ypt4 = m(273,-61.5)
ax.text(xpt4,ypt4,'Bellingshausen ')

m.drawcoastlines(color='black')
m.fillcontinents(color='gray')
m.drawmeridians([90,120,150,180,230,270,300,330], labels=[True,False,False,True],fmt='%g',latmax=80,)
m.drawparallels([-60,-70,-80], labels=[False,True,False,True],latmax=80,)

lonn, latt = 163.7,-74.9
xpt,ypt = m(lonn,latt)
m.plot(xpt,ypt,'r*',label='Qinling station')

hatch_dict = { 0:'', 1:'///', 2:'xx' ,3:'*'}
legend_elements = [ Line2D([0], [0], color='green', lw=1, label='SIE'),
                    Patch(edgecolor='springgreen',facecolor='white',hatch=hatch_dict[1],
                         label='Sea ice'),
                   Patch(facecolor='#0a3b70',hatch=hatch_dict[0],
                         label='Coastal polynya'),
                   Line2D([0], [0], marker='*', color='w', label='Qinling Station',
                          markerfacecolor='r',markersize=15),
                    ]
ax.legend(handles=legend_elements, loc='lower left',fontsize=8,bbox_to_anchor=(0,0))

#     cax = fig.add_axes([0.15, 0.02, 0.7, 0.1],aspect=0.04)
#     cb = fig.colorbar(im1,orientation='vertical',
#                       ticks =levels1[::4], 
#                       boundaries=levels1[::2], 
#                       extend='False', extendfrac=0.05, drawedges=False,
#                     # cax = cax
#                      )
#     cb.set_label('sea ice concentration(%)',fontsize=14)
#     cb.ax.tick_params(labelsize=14, width=0)
fig.tight_layout(rect=[0,0.08,1,0.95])
# plt.savefig('/stu02/weizx24/data//FIGURES/Figure1_obs_0618.pdf',dpi=600,bbox_inches='tight')
plt.show()
#endregion

#region 模式
path_area_weight ='/stu02/weizx24/data/areacello_Ofx_CESM2-WACCM_piControl_r1i1p1f1_gn.nc'
areacella= xr.open_dataset(path_area_weight)
aw_xr = areacella['areacello']
sic = xr.open_dataset('/stu02/weizx24/data/siconc_SIday_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_00010102-05000101_1201_all.nc')['siconc']
lons = sic.lon.values
lats = sic.lat.values

#     y = str(year_max[i].values).zfill(4)
#     print(y)
sic_sel = sic.sel(time='0085').squeeze()
mask_ross1 = (((lons>160.)&(lons<200.)&(lats<-69.))|((lons>200.)&(lons<230.)&(lats<-71.)))&(sic_sel<=15)
tempross1 = np.where(mask_ross1,aw_xr,np.nan)
#     ross_sie[i] = np.nansum(tempross1)

#-----------作图部分--------------------

mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.8
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 9
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 10 


ax = fig.add_subplot(122)
#     levels1 = np.linspace(0,30,21)

#     figname = str(int(year[i]))+'-12-01'
# figname = str(sic_sel.time.values)[0:5]+'12-01'
figname = '(b) CMIP6, 0085-12-01'
ax.set_title(figname,)
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
x1,y1 = m(lons,lats)
im1 = m.contourf(x1,y1,tempross1,
                     levels=levels1,
                     extend='both',
                    shading='faceted', antialiased=True,cmap='RdBu')
im2 = m.contourf(x1,y1,sic_sel,levels=[15,100],hatches=['///', None],alpha=0)
m.contour(x1,y1,sic_sel,levels=[15],colors='g',)

xpt1,ypt1 = m(224,-69)
ax.text(xpt1,ypt1,'Eastern Ross ')
xpt2,ypt2 = m(182,-66)
ax.text(xpt2,ypt2,'Western \n Ross ',fontsize=8)
xpt3,ypt3 = m(250,-63)
ax.text(xpt3,ypt3,'Amundsen')
xpt4,ypt4 = m(273,-61.5)
ax.text(xpt4,ypt4,'Bellingshausen ')

m.drawcoastlines(color='black')
m.fillcontinents(color='gray')
m.drawmeridians([90,120,150,180,230,270,300,330], labels=[True,False,False,True],fmt='%g',latmax=80,)
m.drawparallels([-60,-70,-80], labels=[False,True,False,True],latmax=80,)

lonn, latt = 163.7,-74.9
xpt,ypt = m(lonn,latt)
m.plot(xpt,ypt,'r*',label='Qinling station')


hatch_dict = { 0:'', 1:'///', 2:'xx' }
legend_elements = [ Line2D([0], [0], color='green', lw=1,label='SIE'),
                    Patch(edgecolor='springgreen',facecolor='white',hatch=hatch_dict[1],
                         label='Sea ice'),
                   Patch(facecolor='#0a3b70',hatch=hatch_dict[0],
                         label='Coastal polynya'),
                   Line2D([0], [0], marker='*', color='w', label='Qinling Station',
                          markerfacecolor='r',markersize=15,)
                    ]
ax.legend(handles=legend_elements, loc='lower left',fontsize=8,bbox_to_anchor=(0,0))

#     cax = fig.add_axes([0.15, 0.02, 0.7, 0.1],aspect=0.04)
#     cb = fig.colorbar(im1,orientation='vertical',
#                       ticks =levels1[::4], 
#                       boundaries=levels1[::2], 
#                       extend='False', extendfrac=0.05, drawedges=False,
#                     # cax = cax
#                      )
#     cb.set_label('sea ice concentration(%)',fontsize=14)
#     cb.ax.tick_params(labelsize=14, width=0)
fig.tight_layout(rect=[0,0.08,1,0.95])
plt.savefig('/stu02/weizx24/figures/0924/Figure1_all.pdf',dpi=300,bbox_inches='tight')
plt.show()
#endregion
