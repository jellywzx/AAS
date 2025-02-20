'''sea ice concentration longterm evolution
first start with CESM2
'''
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
# from pyproj import Proj, transform
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import AxesGrid



#region SST 
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
#     if i==7:
#         sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[True,True,True,True])
#         sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,False, True])
#     sps.drawmapboundary(fill_color='AntiqueWhite')
#     draw_latlon_polygon(sps, [170,170, 298, 298], [-80, -60, -60, -80], 'k-')
    sps.fillcontinents(color='gray',lake_color='aqua')
#     
    return sps,x,y

from mpl_toolkits.axes_grid1 import AxesGrid
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-45,urcrnrlat=-70)
clevs = np.linspace(-3,3,13)
# levels = np.linspace(-60, 60, 13)
plt.close()
plt.clf()
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 20
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.8

SST_file = np.load('/stu02/weizx24/data/npz/Figure5_SST.npz')
tos_dif = SST_file['tos_dif']
lons = SST_file['lons']
lats = SST_file['lats']
pval3 = SST_file['pval3']

fig = plt.figure(figsize=(20, 8))
axgr = AxesGrid(fig, [0.05,0.01,0.9,0.9],
                    nrows_ncols=(1, 6),
                    axes_pad=0.3,
                    cbar_location='right',
                    cbar_mode='single',
                #调整colorbar距离图幅的位置
                    cbar_pad='25%',
                    cbar_size='5%',
                    label_mode='',
                    share_all=True)  # note the empty label_mode
case = ['Oct','Nov','Dec','Jan','Feb','Mar',]
# ab = list('abcdefghijklmn')
for i, ax in enumerate(axgr):
    #tos有乘以100，数字整数比较好看
    axz = tos_dif[i]
#     pval1 = np.zeros((332,316),'float')
    #anomaly大于一倍标准差代表显著
#     pval1[np.where(abs(tos_2011_fig[i])>tos_std_sel[i]*1)]=0.99
    sps, x, y = sh_project(ax, lons, lats, i)
    cs1 = sps.contourf(x,y, axz,clevs,cmap='RdBu_r', extend='both')
    cs2 = sps.contourf(x,y, 1-pval3[i], levels=[0.95, 1] ,colors='none',hatches=['..', None],alpha=0)
    m.drawmapboundary(fill_color='AntiqueWhite')
    # ax.set_title('%s'%(case[i]), position=(0.1, 0.99),fontsize=20)
cb = axgr.cbar_axes[0].colorbar(cs1)
cb.set_label('SST (°C)',fontsize=20)
cb.ax.tick_params(labelsize=20)
plt.savefig('/stu02/weizx24/figures/0924/Figure5/Figure5_SST.png',dpi=300,bbox_inches='tight')
# plt.show()
# endregion
print('SST出图完毕')
'''
'''
#region ocean temperature
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
#     
    return sps,x,y

from mpl_toolkits.axes_grid1 import AxesGrid
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-45,urcrnrlat=-70)

clevs = np.linspace(-3,3,13)
# levels = np.linspace(-60, 60, 13)

UOT_file = np.load('/stu02/weizx24/data/npz/Figure5_UOT.npz')
thetao_dif = UOT_file['thetao_dif']
pval3 = UOT_file['pval3']
lons = UOT_file['lons']
lats = UOT_file['lats']



plt.close()
plt.clf()
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 20
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.8
fig = plt.figure(figsize=(20, 8))
axgr = AxesGrid(fig, [0.05,0.01,0.9,0.9],
                    nrows_ncols=(1, 6),
                    axes_pad=0.3,
                    cbar_location='right',
                    cbar_mode='single',
                #调整colorbar距离图幅的位置
                    cbar_pad='25%',
                    cbar_size='5%',
                    # label_mode='',
                    share_all=True)  # note the empty label_mode
case = ['Oct','Nov','Dec','Jan','Feb','Mar',]
# ab = list('abcdefghijklmn')
for i, ax in enumerate(axgr):
    #thetao有乘以100，数字整数比较好看
    axz = thetao_dif[i]
#     pval1 = np.zeros((332,316),'float')
    #anomaly大于一倍标准差代表显著
#     pval1[np.where(abs(thetao_2011_fig[i])>thetao_std_sel[i]*1)]=0.99
    sps, x, y = sh_project(ax, lons, lats, i)
    cs1 = sps.contourf(x,y, axz,clevs,cmap='RdBu_r', extend='both')
    cs2 = sps.contourf(x,y, 1-pval3[i], levels=[0.95, 1] ,colors='none',hatches=['..', None],alpha=0)
    m.drawmapboundary(fill_color='AntiqueWhite')
    # ax.set_title('%s'%(case[i]), position=(0.1, 0.99),fontsize=20)
cb = axgr.cbar_axes[0].colorbar(cs1)
cb.set_label('UOT (°C)',fontsize=20)
cb.ax.tick_params(labelsize=20)
plt.savefig('/stu02/weizx24/figures/0924/Figure5/Figure5_UOT.png',dpi=300,bbox_inches='tight')
plt.show()
print('UOT出图完毕')


#endregion


#region SIC


from mpl_toolkits.axes_grid1 import AxesGrid

file_SIC = np.load('/stu02/weizx24/data/npz/Figure5_SIC.npz')
sic_dif = file_SIC['sic_dif']
lons = file_SIC['lons']
lats = file_SIC['lats']
pval3 = file_SIC['pval3']
sic_clm_sel = file_SIC['sic_clm_sel']

m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=0,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-135,urcrnrlon=20,llcrnrlat=-45,urcrnrlat=-70)
clevs = np.linspace(-60, 60, 13)
levels = np.linspace(-60, 60, 13)

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 20
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.8

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
plt.close()
plt.clf()
fig = plt.figure(figsize=(20, 8))
axgr = AxesGrid(fig, [0.05,0.01,0.9,0.9],
                    nrows_ncols=(1, 6),
                    axes_pad=0.3,
                    cbar_location='right',
                    cbar_mode='single',
                #调整colorbar距离图幅的位置
                    cbar_pad='25%',
                    cbar_size='5%',
                    label_mode='',
                    share_all=True)  # note the empty label_mode
case = ['Oct','Nov','Dec','Jan','Feb','Mar',]
# ab = list('abcdefghijklmn')
for i, ax in enumerate(axgr):
    #SIC有乘以100，数字整数比较好看
    axz = sic_dif[i]
#     pval1 = np.zeros((332,316),'float')
    #anomaly大于一倍标准差代表显著
#     pval1[np.where(abs(sic_2011_fig[i])>sic_std_sel[i]*1)]=0.99
    sps, x, y = sh_project(ax, lons, lats, i)
    cs1 = sps.contourf(x,y, axz, clevs, cmap='RdBu_r', extend='both')
    cs2 = sps.contourf(x,y, 1-pval3[i], levels=[0.95, 1] ,colors='none',hatches=['..', None],alpha=0)
    #气候平均态的海冰范围，还是0.15
    sps.contour(x,y,sic_clm_sel[i],levels=[15],colors='k',linewidths=1)
#     sps.contour(x,y,sic_clm_sel[i],levels=[0.15],colors='r')
    m.drawmapboundary(fill_color='AntiqueWhite')
#     m.fillcontinents(color='lightgray')
#     draw_latlon_polygon(m, [180, 180, 230, 230], [-80, -60, -60, -80], 'k-')
    ax.set_title('%s'%(case[i]), position=(0.1, 0.99),fontsize=20)
cb = axgr.cbar_axes[0].colorbar(cs1)
cb.set_label('SIC (%)',fontsize=20)
cb.ax.tick_params(labelsize=20)
plt.savefig('/stu02/weizx24/figures/0924/Figure5/Figure5_SIC.png',dpi=300,bbox_inches='tight')
plt.show()

#endregion
