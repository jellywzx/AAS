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

def sh_project(ax,lon,lat):
#     sps = Basemap(lon_0=180, boundinglat=-55, projection='spstere',round=True,resolution='l')
    sps = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-70)
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
    #sps.drawcoastlines(color='blue',linewidth=0.6)
    #sps.drawcoastlines(color='blue',linewidth=0.8)
    #sps.fillcontinents(color='',lake_color='aqua')
    #sps.fillcontinents(color='lightgray', lake_color='aqua')
    return sps,x,y

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


#fpath2='/Users/wangshaoyin/data-archive/NSIDC-snow-ice-data-center/sea_ice_drift/icedrift/'
fpath='D:/SID/monthly_seaicemotion/'

#icev = xr.open_dataset(fpath2+'ice_drift_vice_r720x361_sh_197811-202112.nc')['v'].loc['1979-01-01':'2021-12-31']*86.4
icev = xr.open_dataset(fpath+'ice_drift_vice_r720x361_sh_197811-202112.nc')['v'].loc['1979-01-01':'2021-12-31'].values[:,:,:]*86.4

fpath2='D:/SID/monthly_seaicemotion/'
sic = xr.open_dataset(fpath2+'NSIDC_seaice_conc_monthly_0p5_0p5_197901_202112_sh.nc')['nsidc_nt_seaice_conc_monthly'].loc['1979-01-01':'2021-12-31']

adv = np.multiply(sic,icev)
#adv =  icev.values
lat1 = sic.lat.values
lon1 = sic.lon.values

adv_ON = np.nanmean(adv[(adv.time.dt.month==10)|(adv.time.dt.month==11)].values.reshape(43,2,81,720),axis=1)
# adv_ON[13:].shape

years = np.linspace(1992,2021,30)
adv_list = []
for i in range(30):
    start_time = str(int(years[i]))+'-10-01'
    end_time = str(int(years[i]))+'-11-30'
    adv_list.append(adv.sel(time=slice(start_time,end_time)))
    
adv_clm = np.nanmean(np.nanmean(np.array(adv_list),axis=1),axis=0)
adv_std = np.nanstd(np.nanmean(np.array(adv_list),axis=1),axis=0)

years = np.linspace(1992,2021,30)

aa1 = np.zeros((81,720),'float')
pa1 = np.zeros((81,720),'float')
for i in range(81):
    for j in range(720):
        aa1[i,j], intercept, r_value, pa1[i,j], std_err = stats.linregress(years, adv_ON[13:,i,j])

        
levels1 = np.linspace(-80,80,21)
tickmarks = np.linspace(-80,80,5)
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 18
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 20  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6
plt.close()
plt.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6

plt.close()
fig=plt.figure(figsize=(8,10))
ax = fig.add_subplot(111)
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
x1,y1 = m(*np.meshgrid(lon1,lat1))
im1 = m.contourf(x1,y1,aa1*100*10,levels=levels1,extend='both',
                shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)
cs = m.contourf(x1,y1, 1-pa1, levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
# Q1 = m.quiver(xx,yy,uproj*-1,vproj*-1,scale=11,color='purple',scale_units='inches')
# qk1 = plt.quiverkey(Q1, 0.2, -0.025, 2, '2 km/day', labelpos='W',fontproperties={'size':'large'})
m.drawmeridians([90,120,150,180,210,240,270,300,330], linewidth=1.2,labels=[True,False,False,True],fmt='%g',latmax=80,)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,False],latmax=80)
ax.set_title('(d) ON ice advection trend',)
# m.drawmeridians([160,230,300,20,90], linewidth=0.5, dashes=[1, 5],labels=[True,True,True, True])
m.drawcoastlines(color='blue')
m.fillcontinents(color='gray')
#cf = m.contour(x,y, fclm4,levels=np.array([980,986]),colors=['black'],alpha=0.5)
#m.drawparallels([-60],linewidth=1.)
#cs = m.contourf(x,y, p4, levels=[0.9, 0.95] ,colors='none',hatches=['++', None],alpha=0)
# m.drawcoastlines(color='blue')
# cax = fig.add_axes([0.2, 0.02, 0.6, 0.1], aspect=0.05)
# draw_latlon_polygon(m, [180,230], [-70,-70], '-.',color='black',linewidth=2,alpha=1)
cb = plt.colorbar(im1, location='bottom', ticks = tickmarks, extend='both', extendfrac='auto',aspect=30,shrink=0.6, pad=0.05)
cb.ax.tick_params()
cb.set_label('Trend in ice advection per decade' + ' (% km/day)')
fig.tight_layout(rect=[0,0.08,1,0.95])
plt.savefig('C:/Users/fzjxw/python/code/Figures/Trend_ice_advection.png' ,dpi=300,bbox_inches='tight')
plt.show()
