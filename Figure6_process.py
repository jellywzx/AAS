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
import cmocean as cmaps
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
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[True,True,True,True])
    if (i==1)|(i==2):
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,True,True])
    if i==3:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,True,True])
    if i==4:
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,False,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,True,True])
    if (i==5):
        sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[False,True,False,False])
        sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,True,True])
#     if i==7:
#         sps.drawparallels(par, linewidth=0.5, dashes=[1, 5],labels=[True,True,True,True])
#         sps.drawmeridians(mer, linewidth=0.5, dashes=[1, 5],labels=[False,True,False, True])
#     sps.drawmapboundary(fill_color='AntiqueWhite')
#     draw_latlon_polygon(sps, [170,170, 298, 298], [-80, -60, -60, -80], 'k-')
    sps.fillcontinents(color='#577777',lake_color='aqua')
    return sps,x,y
# labelfont=20
# tickfont=16
# colorbarfont=18
#设置全局变量
votemp = xr.open_dataset('/stu02/weizx24/data/ORAS5/79-22_votemp_1000m_30S.nc')['votemper'].sel(time=slice('1992','2022')).sel(lon=slice(160,230)).mean('lon').sel(lat=slice(-78,-50)).sel(LEV=slice(0,300))
year = np.linspace(1992,2021,30)

year_min = [1994., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2016., 2018., 2021.]
num_y = 5
votemp_list = list()
for i in range(30):
    start_time = str(int(year[i]))+'-10-01'
    end_time = str(int(year[i])+1)+'-03-31'
#     d_sel = str(int(year[i])+1)+'-02-15'
#     print(votemp.sel(time=d_sel))
    lat = votemp.sel(time=slice(start_time,end_time)).lat
    lev = votemp.sel(time=slice(start_time,end_time)).LEV
    votemp_list.append(votemp.sel(time=slice(start_time,end_time)))
votemp_mean = np.nanmean(np.array(votemp_list),axis=0)
votemp_std = np.nanstd(np.array(votemp_list),axis=0)
votemp_max = list()
for i in range(num_y):
    start_time = str(int(year_max[i]))+'-10-01'
    end_time = str(int(year_max[i])+1)+'-03-31'
#     d_sel = str(int(year_max[i])+1)+'-02-15'
#     print(votemp.sel(time=d_sel))
#     lat = votemp.sel(time=slice(start_time,end_time)).lat
#     lev = votemp.sel(time=slice(start_time,end_time)).LEV
    votemp_max.append(votemp.sel(time=slice(start_time,end_time)).squeeze())
votemp_max_mean = np.nanmean(np.array(votemp_max),axis=0)
votemp_max_std = np.nanstd(np.array(votemp_max),axis=0)
# votemp = xr.open_dataset('/stu02/weizx24/data/ORAS5/79-22_votemp_1000m_30S.nc')['votemper'].sel(time=slice('1992','2022')).sel(lon=slice(160,230)).mean('lon').sel(lat=slice(-78,-50)).sel(LEV=slice(0,550))
votemp_min = list()
for i in range(num_y):
    start_time = str(int(year_min[i]))+'-10-01'
    end_time = str(int(year_min[i])+1)+'-03-31'
#     d_sel = str(int(year_min[i])+1)+'-02-15'
#     print(votemp.sel(time=d_sel))
#     lat = votemp.sel(time=slice(start_time,end_time)).lat
#     lev = votemp.sel(time=slice(start_time,end_time)).LEV
    votemp_min.append(votemp.sel(time=slice(start_time,end_time)).squeeze())
votemp_min_mean = np.nanmean(np.array(votemp_min),axis=0)
votemp_min_std = np.nanstd(np.array(votemp_min),axis=0)
pre_clm_son = votemp_max_mean
pre_std_son = votemp_max_std

sic_clm_son = votemp_mean
sic_std_son = votemp_std

post_clm_son = votemp_min_mean
post_std_son = votemp_min_std

pval1=np.zeros((6,34,28),'float')
pval2=np.zeros((6,34,28),'float')
pval3=np.zeros((6,34,28),'float')

for k in range(6):
    for i in range(34):
        for j in range(28):
            pval1[k,i,j] = independent_ttest(pre_clm_son[k,i,j], sic_clm_son[k,i,j], pre_std_son[k,i,j],sic_std_son[k,i,j], num_y, 30, 0.1)
            pval2[k,i,j] = independent_ttest(post_clm_son[k,i,j], sic_clm_son[k,i,j], post_std_son[k,i,j],sic_std_son[k,i,j], num_y, 30, 0.1)
            pval3[k,i,j] = independent_ttest(pre_clm_son[k,i,j], post_clm_son[k,i,j], pre_std_son[k,i,j],post_std_son[k,i,j], num_y, num_y, 0.1)
votemp_dif = votemp_max_mean - votemp_min_mean
somx_data = xr.open_dataset('/stu02/weizx24/data/somx_79-23.nc')['somxl010'].sel(time=slice('1992','2022')).sel(lon=slice(160,230)).mean('lon')[:,12:40]
somx_max = list()
for i in range(num_y):
    start_time = str(int(year_max[i]))+'-10-01'
    end_time = str(int(year_max[i])+1)+'-03-31'
#     d_sel = str(int(year_max[i])+1)+'-02-15'
#     print(somx.sel(time=d_sel))
#     lat = somx.sel(time=slice(start_time,end_time)).lat
#     lev = somx.sel(time=slice(start_time,end_time)).LEV
    somx_max.append(somx_data.sel(time=slice(start_time,end_time)).squeeze())
somx_max_mean = np.nanmean(np.array(somx_max),axis=0)
somx_max_std = np.nanstd(np.array(somx_max),axis=0)
# somx = xr.open_dataset('/stu02/weizx24/data/ORAS5/79-22_somx_1000m_30S.nc')['somxer'].sel(time=slice('1992','2022')).sel(lon=slice(160,230)).mean('lon').sel(lat=slice(-78,-50)).sel(LEV=slice(0,550))
somx_min = list()
for i in range(num_y):
    start_time = str(int(year_min[i]))+'-10-01'
    end_time = str(int(year_min[i])+1)+'-03-31'
#     d_sel = str(int(year_min[i])+1)+'-02-15'
#     print(somx.sel(time=d_sel))
#     lat = somx.sel(time=slice(start_time,end_time)).lat
#     lev = somx.sel(time=slice(start_time,end_time)).LEV
    somx_min.append(somx_data.sel(time=slice(start_time,end_time)).squeeze())
somx_min_mean = np.nanmean(np.array(somx_min),axis=0)
somx_min_std = np.nanstd(np.array(somx_min),axis=0)

np.savez('/stu02/weizx24/data/npz/Figure6_obs.npz',votemp_dif=votemp_dif,somx_max_mean=somx_max_mean,somx_min_mean=somx_min_mean,lev=lev,lat=lat,pval3=pval3)