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
from global_land_mask import globe

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



#region SIC plot
def sh_project(ax,lon,lat,i):
    sps = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-70)
    sps.ax = ax
    mer = np.arange(-180, 30, 30.)
    par = np.arange(-90, -50, 10.)
    x,y = sps(lon,lat)
#     x, y = sps(*np.meshgrid(lon, lat))
    if i==0:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False], )
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[True,False,True,False], )
    if (i==1)|(i==2):
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False], )
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,False,True,False], )
    if i==3:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False], )
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,False,True,False], )
    if i==4:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False], )
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,False,True,False], )
    if (i==5):
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,True,False,False], )
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,True,False], )
#     if i==7:
#         sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[True,True,True,True])
#         sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,False, True])
#     sps.drawmapboundary(fill_color='AntiqueWhite')
#     draw_latlon_polygon(sps, [170,170, 298, 298], [-80, -60, -60, -80], 'k-')
    sps.fillcontinents(color='gray',lake_color='aqua')
    return sps,x,y

#设置全局变量
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 30
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.8


#----------图1--------------
sic_file = np.load('/stu02/weizx24/data/npz/Figure4_SIC.npz')
sic_dif = sic_file['sic_dif']
sic_clm = sic_file['sic_clm']
lons = sic_file['lons']
lats = sic_file['lats']
pval3 = sic_file['pval3']

from mpl_toolkits.axes_grid1 import AxesGrid

m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-45,urcrnrlat=-70)

clevs = np.linspace(-60, 60, 13)
levels = np.linspace(-60, 60, 13)

plt.close()
plt.clf()
#设置全局变量
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 30
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.8

fig = plt.figure(figsize=(20,9))
axgr = AxesGrid(fig, [0.05,0.01,0.9,0.9],
                    nrows_ncols=(1, 6),
                    axes_pad=0.3,
                    cbar_location='right',
                    cbar_mode='single',
                #调整colorbar距离图幅的位置
                    cbar_pad='38%',
                    cbar_size='5%',
                    label_mode='',
                    share_all=True)  # note the empty label_mode
case = ['Oct','Nov','Dec','Jan','Feb','Mar',]
# ab = list('abcdefghijklmn')
for i, ax in enumerate(axgr):
    #SIC有乘以100，数字整数比较好看
    axz = sic_dif[i]*100
    sps, x, y = sh_project(ax, lons, lats, i)
    cs1 = sps.contourf(x,y, axz, clevs, cmap='RdBu_r', extend='both')
    cs2 = sps.contourf(x,y, 1-pval3[i], levels=[0.95, 1] ,colors='none',hatches=['..', None],alpha=0)
    #气候平均态的海冰范围，还是0.15
    sps.contour(x,y,sic_clm[i],levels=[0.15],colors='k',linewidths=1)
    m.drawmapboundary(fill_color='AntiqueWhite')
    mer = np.arange(-180, 30, 30.)
    par = np.arange(-90, -50, 10.)
    # ax.set_title('%s'%((case[i])), position=(0.1, 0.9),fontsize=40)
    ax.text(-0.12, 1.1, case[i],
        transform=ax.transAxes,
        fontsize=40,
        # fontweight='bold',
        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.1'))

cb = axgr.cbar_axes[0].colorbar(cs1)
cb.set_ticks(np.linspace(-60, 60, 13)[::3]) 
cb.set_label('SIC (%)',fontsize=40)
cb.ax.tick_params(labelsize=40)
plt.savefig('/stu02/weizx24/figures/0924/Figure4/Figure4a.png',dpi=300,
            bbox_inches='tight'
            )
plt.show()
print('SIC出图完毕')


#region SST
file_SST = np.load('/stu02/weizx24/data/npz/Figure4_SST.npz')
sst_dif = file_SST['sst_dif']
# sst_clm = file_SST['sst_clm']
pval3 = file_SST['pval3']
lons = file_SST['lons']
lats = file_SST['lats']

def sh_project(ax,lon,lat,i):
    sps = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-70)
    sps.ax = ax
    mer = np.arange(-180, 30, 30.)
    # mer = [-120,-90,-60]
    par = np.arange(-90, -50, 10.)
#     x,y = sps(lon,lat)
    x, y = sps(*np.meshgrid(lon, lat))
    if i==0:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[True,False,False,False])
    if (i==1)|(i==2):
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
    if i==3:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
    if i==4:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
    if (i==5):
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,True,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,False,False])
    return sps,x,y

from mpl_toolkits.axes_grid1 import AxesGrid
plt.close()
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-45,urcrnrlat=-70)

clevs = np.linspace(-1, 1, 13)
levels = np.linspace(-1, 1, 13)
fig = plt.figure(figsize=(20, 9))
# plt.rcParams['font.family'] = 'Arial'
# plt.rcParams['font.size'] = 40
# plt.rcParams['lines.linewidth'] = 1
# plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
# plt.rcParams['lines.markersize'] = 1 
# # plt.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['hatch.color'] ='springgreen'
# mpl.rcParams['hatch.linewidth'] = 0.8

axgr = AxesGrid(fig, [0.05,0.01,0.9,0.9],
                    nrows_ncols=(1, 6),
                    axes_pad=0.3,
                    cbar_location='right',
                    cbar_mode='single',
                #调整colorbar距离图幅的位置
                    cbar_pad='38%',
                    cbar_size='5%',
                    label_mode='',
                    share_all=True)  # note the empty label_mode
case = ['Oct','Nov','Dec','Jan','Feb','Mar',]
# ab = list('abcdefghijklmn')
for i, ax in enumerate(axgr):
    #sst有乘以100，数字整数比较好看
    axz = sst_dif[i]
    sps, x, y = sh_project(ax, lons, lats, i)
    cs1 = sps.contourf(x,y, axz, clevs, cmap='RdBu_r', extend='both')
    cs2 = sps.contourf(x,y, 1-pval3[i], levels=[0.95, 1] ,colors='none',hatches=['..', None],alpha=0)
    sps.drawmapboundary(fill_color='AntiqueWhite')
    sps.fillcontinents(color='gray',lake_color='aqua')
cb = axgr.cbar_axes[0].colorbar(cs1)
cb.set_ticks(np.linspace(-1, 1, 13)[::3]) 
cb.set_label('SST (°C)',fontsize=40)
cb.ax.tick_params(labelsize=40)
plt.savefig('/stu02/weizx24/figures/0924/Figure4/Figure4b.png',dpi=300,bbox_inches='tight')
# plt.show()
print('SST出图完成')
#endregion



#region votemper
file_votemper = np.load('/stu02/weizx24/data/npz/Figure4_votemper.npz')
votemper_dif = file_votemper['votemper_dif']
pval3 = file_votemper['pval3']
lons = file_votemper['lons']
lats = file_votemper['lats']

def sh_project(ax,lon,lat,i):
    sps = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-70)
    sps.ax = ax
    mer = np.arange(-180, 30, 30.)
    par = np.arange(-90, -50, 10.)
#     x,y = sps(lon,lat)
    x, y = sps(*np.meshgrid(lon, lat))
    if i==0:
        # sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[True,True,False,True])
    if (i==1)|(i==2):
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,True])
    if i==3:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,True])
    if i==4:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,True])
    if (i==5):
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,True,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,True])
#     if i==7:
#         sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[True,True,True,True])
#         sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,False, True])
#     sps.drawmapboundary(fill_color='AntiqueWhite')
#     draw_latlon_polygon(sps, [170,170, 298, 298], [-80, -60, -60, -80], 'k-')
    sps.fillcontinents(color='gray',lake_color='aqua')
    return sps,x,y

# from mpl_toolkits.axes_grid1 import AxesGrid
plt.close()
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-45,urcrnrlat=-70)

clevs = np.linspace(-1, 1, 13)
levels = np.linspace(-1, 1, 13)

# plt.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['hatch.color'] ='springgreen'
# mpl.rcParams['hatch.linewidth'] = 0.8

fig = plt.figure(figsize=(20, 9))
# plt.rcParams['font.family'] = 'Arial'
# plt.rcParams['font.size'] = 20
# plt.rcParams['lines.linewidth'] = 1
# plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
# plt.rcParams['lines.markersize'] = 1 
# # plt.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['hatch.color'] ='springgreen'
# mpl.rcParams['hatch.linewidth'] = 0.8

axgr = AxesGrid(fig, [0.05,0.01,0.9,0.9],
                    nrows_ncols=(1, 6),
                    axes_pad=0.3,
                    cbar_location='right',
                    cbar_mode='single',
                #调整colorbar距离图幅的位置
                    cbar_pad='38%',
                    cbar_size='5%',
                    label_mode='',
                    share_all=True)  # note the empty label_mode
case = ['Oct','Nov','Dec','Jan','Feb','Mar',]
# ab = list('abcdefghijklmn')
for i, ax in enumerate(axgr):
    #votemper有乘以100，数字整数比较好看
    axz = votemper_dif[i]
    sps, x, y = sh_project(ax, lons, lats, i)
    cs1 = sps.contourf(x,y, axz, clevs, cmap='RdBu_r', extend='both')
    cs2 = sps.contourf(x,y, 1-pval3[i], levels=[0.95, 1] ,colors='none',hatches=['..', None],alpha=0)
    sps.drawmapboundary(fill_color='AntiqueWhite')
cb = axgr.cbar_axes[0].colorbar(cs1)
cb.set_ticks(np.linspace(-1, 1, 13)[::3]) 
cb.set_label('UOT (°C)',fontsize=40)
cb.ax.tick_params(labelsize=40)
plt.savefig('/stu02/weizx24/figures/0924/Figure4/Figure4c.png',dpi=300,bbox_inches='tight')
# plt.show()
print('UOT出图完成')
#endregion

