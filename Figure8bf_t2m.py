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




#region 气温
#--------试一下5年的------------
year_min = [1994., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2016., 2018., 2021.]
num_y = 5

sic = xr.open_dataset('/stu02/weizx24/data/ERA5/ERA5-single-level-1x1-resolution-2meter-airtemp-1979-202208.nc')['t2m'].sel(time=slice('1992','2021'))
t2mn1 = np.zeros((num_y,12,181,360),'float')
t2mn2 = np.zeros((num_y,12,181,360),'float')
sic_clm = np.nanmean(sic.values.reshape(30,12,181,360),axis=0)
sic_std = np.nanstd(sic.values.reshape(30,12,181,360),axis=0)
lons = sic.longitude
lats = sic.latitude

for i in range(num_y):
    t2mn1[i,:,:,:] = sic[(sic['time.year']==year_max[i])].values.reshape(12,181,360)
    t2mn2[i,:,:,:] = sic[(sic['time.year']==year_min[i])].values.reshape(12,181,360)
    
# t2mn1SON = np.nanmean(t2mn1,axis=0)[9:11]
# t2mn2SON = np.nanmean(t2mn2,axis=0)[9:11]
#2000年前春季海冰平均  
pre_clm_son = t2mn1.mean(axis=0)[9:11].mean(axis=0)
pre_std_son = t2mn1.std(axis=0)[9:11].mean(axis=0)
#2000年后春季海冰平均  
post_clm_son = t2mn2.mean(axis=0)[9:11].mean(axis=0)
post_std_son = t2mn2.std(axis=0)[9:11].mean(axis=0)
#春季海冰平均气候平均态  
sic_clm_son = sic_clm[9:11].mean(axis=0)
sic_std_son = sic_std[9:11].mean(axis=0)
#求合成差 
# comp1 = t2mn1.mean(axis=0)  - sic_clm
# comp2 = t2mn2.mean(axis=0) - sic_clm
comp3 = t2mn1.mean(axis=0) - t2mn2.mean(axis=0)

# corSON1 = comp1[9:11].mean(axis=0)
# corSON2 = comp2[9:11].mean(axis=0)
corSON3 =comp3[9:11].mean(axis=0)

# pval1=np.zeros((181,360),'float')
# pval2=np.zeros((181,360),'float')
pval3=np.zeros((181,360),'float')

for i in range(181):
    for j in range(360):
#         pval1[i,j] = independent_ttest(pre_clm_son[i,j], sic_clm_son[i,j], pre_std_son[i,j],sic_std_son[i,j], num_y, 21, 0.1)
#         pval2[i,j] = independent_ttest(post_clm_son[i,j], sic_clm_son[i,j], post_std_son[i,j],sic_std_son[i,j], num_y, 21, 0.1)
        pval3[i,j] = independent_ttest(pre_clm_son[i,j], post_clm_son[i,j], pre_std_son[i,j],post_std_son[i,j], num_y,num_y, 0.1)

levels1=np.linspace(-5,5,50)
tickmarks=np.linspace(-5,5,5)

#设置全局变量
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 25
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 20  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=180,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-130,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
#-----------------------------第一张图------------------------------------------
# plt.close()
fig=plt.figure(figsize=(7,10))
# fig.subplots_adjust(top=0.99)
ax3 = fig.add_subplot(211)
x, y = m(*np.meshgrid(lons,lats))
im = m.contourf(x,y,corSON3,levels=levels1 ,extend='both',
                shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)
m.drawmeridians([90,120,150,180,210,240,270,300,330], linewidth=1.2,labels=[False,False,False,False],fmt='%g',latmax=80)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,False],latmax=80)
# ax3.annotate(r'$70\degree\,\mathrm{S}$',xy=m(190,-70),xycoords='data')
# ax3.annotate(r'$80\degree\,\mathrm{S}$',xy=m(190,-80),xycoords='data')
cs = m.contourf(x,y, 1-pval3, levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
# ax3.set_title('Observation')
m.drawcoastlines(color='blue')

#------------------------------公共部分-------------------------
# cax = fig.add_axes([0.3, 0.49, 0.4, 0.08],aspect=0.03)
# cb = fig.colorbar(im,orientation='horizontal',ticks=tickmarks,extend='both',cax = cax)
# # cb.ax.tick_params(labelsize=16)
# cb.set_label("Air temperature (K)")
# fig.tight_layout(rect=[0.02,0.1,0.98,0.95])
ax3.text(0, 1.05, '(b)',  transform=ax3.transAxes, va='top', ha='right')
# plt.savefig('/stu02/weizx24/figures/0924/Figure8b_obs_t2m.png' ,dpi=300,bbox_inches='tight')
# plt.show()
#endregion

#region 模式
#————————模式——————————
# 筛选最大最小年份
CMIP6_opwa = np.load('/stu02/weizx24/data/npz/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']*1e-12
ross_sie_cmip6 = np.nanmean(CMIP6_opwa,axis=1)
Q1 = np.quantile(ross_sie_cmip6,0.05)
Q3 = np.quantile(ross_sie_cmip6,0.95)
sic_Feb = xr.open_dataset('/stu02/weizx24/data/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_Feb.nc')['siconc']
year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year

tas = xr.open_dataset('/stu02/weizx24/data/tas_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['tas']
lons = tas.lon
lats = tas.lat
tas_clm = np.nanmean(tas.values.reshape(499,12,96,144),axis=0)
tas_std = np.nanstd(tas.values.reshape(499,12,96,144),axis=0)
tas_min_lst = []
for i in range(len(year_min)):
    y = str(year_min[i].values).zfill(4)
#     print(tas.sel(time=y))
    tas_min_lst.append(tas.sel(time=y))
tas_min_mean = np.array(tas_min_lst).squeeze().mean(axis=0)
tas_min_std = np.array(tas_min_lst).squeeze().std(axis=0)

# year_max = tas.time[ross_sie_cmip6<Q1].dt.year
tas_max_lst = []
for i in range(len(year_max)):
    y = str(year_max[i].values).zfill(4)
#     print(y)
    tas_max_lst.append(tas.sel(time=y))
tas_max_mean = np.array(tas_max_lst).squeeze().mean(axis=0)
tas_max_std = np.array(tas_max_lst).squeeze().std(axis=0)
tas_dif = np.nanmean(tas_max_mean[9:11],axis=0)-np.nanmean(tas_min_mean[9:11],axis=0)
num_y = 25


pre_clm_son = np.nanmean(tas_max_mean[9:11],axis=0)
pre_std_son = np.nanmean(tas_max_std[9:11],axis=0)

slp_clm_son = np.nanmean(tas_clm[9:11],axis=0)
slp_std_son = np.nanmean(tas_std[9:11],axis=0)

post_clm_son = np.nanmean(tas_min_mean[9:11],axis=0)
post_std_son = np.nanmean(tas_min_std[9:11],axis=0)

pval1=np.zeros((96,144),'float')
pval2=np.zeros((96,144),'float')
pval3=np.zeros((96,144),'float')

for i in range(96):
    for j in range(144):
#         pval1[i,j] = independent_ttest(pre_clm_son[i,j], slp_clm_son[i,j], pre_std_son[i,j],slp_std_son[i,j],num_y, 30, 0.1)
#         pval2[i,j] = independent_ttest(post_clm_son[i,j], slp_clm_son[i,j], post_std_son[i,j],slp_std_son[i,j], num_y, 30, 0.1)
        pval3[i,j] = independent_ttest(pre_clm_son[i,j], post_clm_son[i,j], pre_std_son[i,j],post_std_son[i,j], num_y,num_y,0.1)

    
levels1=np.linspace(-5,5,50)
tick_marks = np.linspace(-5,5,5)
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=180,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-130,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
#-----------------------------第一张图------------------------------------------
# plt.close()
# fig=plt.figure(figsize=(11,7.5))
# fig.subplots_adjust(top=0.99)
ax4 = fig.add_subplot(212)
x, y = m(*np.meshgrid(lons,lats))
im = m.contourf(x,y,tas_dif,
                levels=levels1 ,
                extend='both',
                shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)
m.drawmeridians([90,120,150,180,210,240,270,300,330], labels=[False,False,False,True],fmt='%g',latmax=80)
m.drawparallels([-60,-70,-80],labels=[False,False,False,False],latmax=80,)
# ax4.annotate(r'$70\degree\,\mathrm{S}$',xy=m(190,-70),xycoords='data')
# ax4.annotate(r'$80\degree\,\mathrm{S}$',xy=m(190,-80),xycoords='data')
cs = m.contourf(x,y, 1-pval3, levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
# ax4.set_title('CESM2-WACCM-FV2')
m.drawcoastlines(color='blue')

#------------------------------公共部分-------------------------
cax = fig.add_axes([0.3, 0.02, 0.4, 0.08], aspect=0.03)
cb = fig.colorbar(im,orientation='horizontal',
                  ticks=tick_marks,
                  extend='both',cax = cax)
# cb.ax.tick_params(labelsize=22)
cb.set_label("T2M (K)")
# fig.suptitle('Composite SLP & UV10(ON)', fontsize=24)
# fig.tight_layout(rect=[0.02,0.1,0.98,0.95])
ax4.text(0, 1.05, '(f)', transform=ax4.transAxes, va='top', ha='right')
plt.subplots_adjust(hspace=0.05)
plt.savefig('/stu02/weizx24/figures/0924/Figure8/Figure8bf_t2m.png' ,dpi=300,bbox_inches='tight')
# plt.show()

#endregion
