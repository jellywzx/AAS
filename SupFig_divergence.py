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
from scipy.stats import t
import cmocean as cm

def sh_project(ax,lon,lat):
    sps = Basemap(lon_0=180, boundinglat=-55, projection='spstere',round=True,resolution='l')

    #sps = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=-30, llcrnrlon=0,urcrnrlon=360,resolution='c')
    sps.ax = ax
    mer = np.arange(0, 360, 60.)
    #par = np.arange(-90, -50, 10.)
    x, y = sps(*np.meshgrid(lon, lat))
    #if i==0:
    #sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
    #sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[True,True,True,True])
    #sps.drawmapboundary(fill_color='AntiqueWhite')
    #sps.ax.annotate(r'$70\degree N$',fontsize=12,xy=sps(165,70),xycoords='data',xytext = sps(165,70),textcoords='data')
    #sps.ax.annotate(r'$60\degree N$',fontsize=12,xy=sps(170,60),xycoords='data',xytext = sps(170,60),textcoords='data')
    #sps.ax.annotate(r'$80\degree N$',fontsize=12,xy=sps(0,80),xycoords='data',xytext = sps(0,80),textcoords='data')
    #sps.ax.annotate(r'$70\degree N$',fontsize=12,xy=sps(0,70),xycoords='data',xytext = sps(0,70),textcoords='data')
    #sps.ax.annotate(r'$60\degree N$',fontsize=10,xy=sps(170,60),xycoords='data',xytext = sps(170,60),textcoords='data')
    #sps.ax.annotate(r'$80\degree N$',fontsize=10,xy=sps(160,80),xycoords='data',xytext = sps(160,80),textcoords='data')
    #draw_latlon_polygon(m, [0,360], [60, 60], '--',color='gray',linewidth=1.5,alpha=1)
    #draw_latlon_polygon(m, [0,360], [70, 70], '--',color='gray',linewidth=1.5,alpha=1)
    #draw_latlon_polygon(m, [0,360], [80, 80], '--',color='gray',linewidth=1.5,alpha=1)
    #draw_latlon_polygon(m, [0, 0], [90, 50], '--',color='gray',linewidth=1.5)
    #draw_latlon_polygon(m, [60, 60], [90, 50], '--',color='gray',linewidth=1.5)
    #draw_latlon_polygon(m, [300, 300], [90, 50], '--',color='gray',linewidth=1.5)
    #draw_latlon_polygon(m, [120, 120], [90, 50], '--',color='gray',linewidth=1.5)
    #draw_latlon_polygon(m, [180, 180], [90, 50], '--',color='gray',linewidth=1.5)
    #sps.drawparallels([-90,-60,-30], linewidth=1,labels=[True,True,True, True],alpha=0.6,latmax=80)
    #sps.drawmeridians([0,60,120,180,240,300], linewidth=1,labels=[True,True,True, True],fmt='%g',alpha=0.6,latmax=80)
    #sps.drawparallels([60,70,80], linewidth=1,labels=[True,True,True, True],alpha=0.6,latmax=80)
    #sps.drawcoastlines(color='gray',linewidth=0.6)
    #sps.drawcoastlines(color='blue',linewidth=0.8)
    #sps.fillcontinents(color='',lake_color='aqua')
    sps.fillcontinents(color='lightgray', lake_color='aqua')
    return sps,x,y

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


fpath='D:/SID/monthly_seaicemotion/'
icev = xr.open_dataset(fpath+'ice_drift_vice_r720x361_sh_197811-202112.nc')['v'].loc['1992-01-01':'2021-12-31']
iceu = xr.open_dataset(fpath+'ice_drift_uice_r720x361_sh_197811-202112.nc')['u'].loc['1992-01-01':'2021-12-31']
fpath2='D:/SID/monthly_seaicemotion/'
sic = xr.open_dataset(fpath2+'NSIDC_seaice_conc_monthly_0p5_0p5_197901_202112_sh.nc')['nsidc_nt_seaice_conc_monthly'].loc['1992-01-01':'2021-12-31']
lat1 = sic.lat.values
lon1 = sic.lon.values

phi = lat1*np.pi/180.0
cosphi = np.cos(phi)[np.newaxis,:,np.newaxis]
A = 6.37*1e3
dx = A * cosphi * 0.5/180*np.pi
dy = A * 0.5/180*np.pi

adv = (np.gradient(icev,axis=1,edge_order=2)/dy + np.gradient(iceu,axis=2,edge_order=2)/dx)*86.4

#var = np.nanmean(adv.reshape(42,12,81, 720)[:,8:11,:,:],axis=1)*100
var_nov = adv.reshape(30,12,81,720)[:,10,:,:]*100
var_nov_clm = var_nov.mean(axis=0)

# aa1 = var[-5:].mean(axis=0) - var[0:40].mean(axis=0)
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 15
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 16  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6

levels1=np.linspace(-1,1,21)
#aa1[np.where(np.abs(aa1)<0.2)] = np.nan
plt.close()
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.8
#x1, y1 = m(*np.meshgrid(lon,lat))
fig=plt.figure(figsize=(7,6.5))
ax = fig.add_subplot(111)
#ax.set_title('(a) ASIE (ERA5)',fontsize=12)
m, x1, y1 = sh_project(ax, lon1, lat1)

#m2, x2, y2 = sh_project(ax, lon1, lat1)
im1 = m.contourf(x1,y1,var_nov_clm,
                 levels=levels1,
                 extend='both',
                shading='faceted', antialiased=True,cmap=cm.cm.balance)
# ax.set_title('ON', position=(0.5, 0.48),bbox={ 'pad':5},fontsize=18)
m.drawmeridians([0,30,60,90,120,150,180,210,240,270,300,330], linewidth=1.2, dashes=[1, 5],labels=[True,True,True, True])
m.drawparallels([-80,-70,-60],linewidth=1.2, dashes=[1, 5],labels=[True,True,True, True])
m.drawcoastlines(color='gray')
m.drawcoastlines(color='gray')
ax.annotate('70° S',xy=m(10,-70),xycoords='data',xytext = m(10,-70),textcoords='data')
ax.annotate('60° S',xy=m(10,-60),xycoords='data',xytext = m(10,-60),textcoords='data')
ax.annotate('80° S',xy=m(10,-80),xycoords='data',xytext = m(10,-80),textcoords='data')
cax = fig.add_axes([0.16, 0.02, 0.7, 0.15], aspect=0.03)
draw_latlon_polygon(m, [180,180], [-80,-60], '-.',color='purple',linewidth=2,alpha=1)
draw_latlon_polygon(m, [230,230], [-80,-60], '-.',color='purple',linewidth=2,alpha=1)
draw_latlon_polygon(m, [180,230], [-80,-80], '-.',color='purple',linewidth=2,alpha=1)
draw_latlon_polygon(m, [180,230], [-60,-60], '-.',color='purple',linewidth=2,alpha=1)
fig.tight_layout(rect=[0,0.08,1,0.95])
cbar = plt.colorbar(im1,orientation='horizontal',cax=cax)
cbar.set_label('Sea ice divergence')
#plt.savefig('fig/Reg_asl_t2m.png' ,dpi=300,bbox_inches='tight')
plt.savefig('C:/Users/fzjxw/python/code/Figures/seaicedivergence.png' ,dpi=300,bbox_inches='tight')
plt.show()
