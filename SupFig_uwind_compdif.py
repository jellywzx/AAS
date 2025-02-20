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



#合成年份为5年
year_min = [1994., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2016., 2018., 2021.]
num_y = 5

u10_sel = xr.open_dataset('D:/ERA5/ERA5-1x1-resolution-uv10-197901-202212.nc')['u10'][:,::-1,:].sel(time=slice('1992','2021'))
v10_sel = xr.open_dataset('D:/ERA5/ERA5-1x1-resolution-uv10-197901-202212.nc')['v10'][:,::-1,:].sel(time=slice('1992','2021'))

lons = u10_sel.longitude.values
lats = u10_sel.latitude.values

u10_SON_clm = u10_sel[(u10_sel.time.dt.month==10)|(u10_sel.time.dt.month==11)].mean('time')
v10_SON_clm = v10_sel[(v10_sel.time.dt.month==10)|(v10_sel.time.dt.month==11)].mean('time')

u = list()
v = list()
for i in range(num_y):
    stime = str(int(year_max[i]))+'-10-01'
    etime = str(int(year_max[i]))+'-11-30'
    u.append(u10_sel.sel(time=slice(stime,etime)))
    v.append(v10_sel.sel(time=slice(stime,etime)))
# print(u)
pre_u = np.array(u).mean(axis=1)
pre_v = np.array(v).mean(axis=1)
pre_u_anom = pre_u - u10_SON_clm.values
pre_v_anom = pre_v - v10_SON_clm.values

u = list()
v = list()
for i in range(num_y): 
    stime = str(int(year_min[i]))+'-10-01'
    etime = str(int(year_min[i]))+'-11-30'
    u.append(u10_sel.sel(time=slice(stime,etime)))
    v.append(v10_sel.sel(time=slice(stime,etime)))
# print(u)
post_u = np.array(u).mean(axis=1)
post_v = np.array(v).mean(axis=1)
post_u_anom = post_u - u10_SON_clm.values
post_v_anom = post_v - v10_SON_clm.values

u_anom_diff = pre_u_anom - post_u_anom
v_anom_diff = pre_v_anom - post_v_anom


levels1=np.linspace(-3,3,50)
tickmarks = np.linspace(-3,3,5)
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 18
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 20  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6

m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=180,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-130,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
fig=plt.figure(figsize=(7,10))
# fig.subplots_adjust(top=0.99)
ax1 = fig.add_subplot(111)
x, y = m(*np.meshgrid(lons,lats))
im = m.contourf(x,y,u_anom_diff.mean(axis=0),
                levels=levels1 ,
                extend='both',
                shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)

# ugrid1,newlons = shiftgrid(180.,u_anom_diff.mean(axis=0),lons,start=False)
# vgrid1,newlons = shiftgrid(180.,v_anom_diff.mean(axis=0),lons,start=False) 
# uproj,vproj,xx,yy = m.transform_vector(ugrid1,vgrid1,newlons,lats,25,25,returnxy=True,masked=True)
# Q1 = m.quiver(xx,yy,uproj*-1,vproj*-1,scale=10,color='black',scale_units='inches')
# qk1 = plt.quiverkey(Q1, 0.16, -0.03, 2, '2 m/s', labelpos='W',color='k',)

m.drawmeridians([90,120,150,180,210,240,270,300,330],labels=[True,False,False,True],fmt='%g',latmax=80)
m.drawparallels([-60,-70,-80], labels=[False,False,False,False],latmax=80)
ax1.annotate(r'$70\degree S$',xy=m(190,-70),xycoords='data',xytext = m(190,-70),textcoords='data')
ax1.annotate(r'$60\degree S$',xy=m(190,-60),xycoords='data',xytext = m(190,-60),textcoords='data')
ax1.annotate(r'$80\degree S$',xy=m(190,-80),xycoords='data',xytext = m(190,-80),textcoords='data')
# cs = m.contourf(x,y, 1-pval3, levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
# ax1.set_title('Observation')
m.drawcoastlines(color='blue')
#------------------------------公共部分-------------------------
cax = fig.add_axes([0.15, 0.1, 0.7, 0.08],aspect=0.03)
cb = fig.colorbar(im,orientation='horizontal',
                  ticks=tickmarks,
                  extend='both', 
                  cax = cax
                  )
# cb.ax.tick_params(labelsize=16)
cb.set_label("U wind (m/s)")
# fig.suptitle('Composite SLP & UV10(ON)', fontsize=24)
# fig.tight_layout(rect=[0.02,0.1,0.98,0.95])
# ax1.text(0, 1.05, '(a)', transform=ax1.transAxes, va='top', ha='right')
plt.savefig('C:/Users/fzjxw/python/code/Figures/Uwind_compdif.png' ,dpi=300,bbox_inches='tight')
plt.show()
        