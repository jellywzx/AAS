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



#region 观测
#fpath2='/Users/wangshaoyin/data-archive/NSIDC-snow-ice-data-center/sea_ice_drift/icedrift/'
fpath='/stu02/weizx24/data/monthly_seaicemotion/'

#icev = xr.open_dataset(fpath2+'ice_drift_vice_r720x361_sh_197811-202112.nc')['v'].loc['1979-01-01':'2021-12-31']*86.4
icev = xr.open_dataset(fpath+'ice_drift_vice_r720x361_sh_197811-202112.nc')['v'].loc['1979-01-01':'2021-12-31'].values[:,:,:]*86.4

fpath2='/stu02/weizx24/data/'
sic = xr.open_dataset(fpath2+'NSIDC_seaice_conc_monthly_0p5_0p5_197901_202112_sh.nc')['nsidc_nt_seaice_conc_monthly'].loc['1979-01-01':'2021-12-31']

adv = np.multiply(sic,icev)
#adv =  icev.values
lat1 = sic.lat.values
lon1 = sic.lon.values


year_min = [1994., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2016., 2018., 2021.]
num_y = 5

adv_ON = np.nanmean(adv[(adv.time.dt.month==10)|(adv.time.dt.month==11)].values.reshape(43,2,81,720),axis=1)

years = np.linspace(1992,2021,30)
adv_list = []
for i in range(30):
    start_time = str(int(years[i]))+'-10-01'
    end_time = str(int(years[i]))+'-11-30'
    adv_list.append(adv.sel(time=slice(start_time,end_time)))
    
adv_clm = np.nanmean(np.nanmean(np.array(adv_list),axis=1),axis=0)
adv_std = np.nanstd(np.nanmean(np.array(adv_list),axis=1),axis=0)

adv_max_list = []
for i in range(5):
    start_time = str(int(year_max[i]))+'-10-01'
    end_time = str(int(year_max[i]))+'-11-30'
    adv_max_list.append(adv.sel(time=slice(start_time,end_time)))
adv_max_arr = np.nanmean(np.array(adv_max_list),axis=1)

adv_min_list = []
for i in range(5):
    start_time = str(int(year_min[i]))+'-10-01'
    end_time = str(int(year_min[i]))+'-11-30'
    adv_min_list.append(adv.sel(time=slice(start_time,end_time)))
adv_min_arr = np.nanmean(np.array(adv_min_list),axis=1)

adv_dif = np.nanmean(adv_max_arr-adv_clm,axis=0) - np.nanmean(adv_min_arr-adv_clm,axis=0)

pre_clm_son = np.nanmean(adv_max_arr,axis=0)
pre_std_son = np.nanstd(adv_max_arr,axis=0)

# slp_clm_son = adv_clm
# slp_std_son = adv_std

post_clm_son = np.nanmean(adv_min_arr,axis=0)
post_std_son = np.nanstd(adv_min_arr,axis=0)

# pval1=np.zeros((81,720),'float')
# pval2=np.zeros((81,720),'float')
pval3=np.zeros((81,720),'float')

for i in range(81):
    for j in range(720):
#         pval1[i,j] = independent_ttest(pre_clm_son[i,j], slp_clm_son[i,j], pre_std_son[i,j],slp_std_son[i,j],num_y, 30, 0.1)
#         pval2[i,j] = independent_ttest(post_clm_son[i,j], slp_clm_son[i,j], post_std_son[i,j],slp_std_son[i,j], num_y, 30, 0.1)
        pval3[i,j] = independent_ttest(pre_clm_son[i,j], post_clm_son[i,j], pre_std_son[i,j],post_std_son[i,j], num_y,num_y,0.1)
        
        
levels1 = np.linspace(-500,500,50)
tickmarks = np.linspace(-500,500,5)

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 18
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 20  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6

plt.close()
#x1, y1 = m(*np.meshgrid(lon,lat))
fig=plt.figure(figsize=(7,10))
# fig.subplots_adjust(top=0.99)
ax1 = fig.add_subplot(211)
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=180,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-130,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
x1,y1 = m(*np.meshgrid(lon1,lat1))
im1 = m.contourf(x1,y1,adv_dif*100,levels=levels1,extend='both',
                shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)
cs = m.contourf(x1,y1, 1-pval3, levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
# Q1 = m.quiver(xx,yy,uproj*-1,vproj*-1,scale=11,color='purple',scale_units='inches')
# qk1 = plt.quiverkey(Q1, 0.2, -0.025, 2, '2 km/day', labelpos='W',fontproperties={'size':'large'},)
m.drawmeridians([90,120,150,180,210,240,270,300,330], linewidth=1.2,labels=[True,False,False,True],fmt='%g',latmax=80)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,False],latmax=80,)
ax1.annotate(r'$70\degree S$',xy=m(190,-70),xycoords='data',xytext = m(190,-70),textcoords='data')
ax1.annotate(r'$60\degree S$',xy=m(190,-60),xycoords='data',xytext = m(190,-60),textcoords='data')
ax1.annotate(r'$80\degree S$',xy=m(190,-80),xycoords='data',xytext = m(190,-80),textcoords='data')

#ax.annotate(r'$70\degree S$',fontsize=16,xy=m(10,-70),xycoords='data',xytext = m(10,-70),textcoords='data')
#ax.annotate(r'$60\degree S$',fontsize=16,xy=m(10,-60),xycoords='data',xytext = m(10,-60),textcoords='data')
#ax.annotate(r'$80\degree S$',fontsize=16,xy=m(10,-80),xycoords='data',xytext = m(10,-80),textcoords='data')
# ax1.set_title('Observation')

# m.drawmeridians([160,230,300,20,90], linewidth=0.5, dashes=[1, 5],labels=[True,True,True, True])
m.drawcoastlines(color='blue')
m.fillcontinents(color='gray')
ax1.text(0, 1.05, '(d)', fontsize=20, transform=ax1.transAxes, va='top', ha='right')
#cf = m.contour(x,y, fclm4,levels=np.array([980,986]),colors=['black'],alpha=0.5)
#m.drawparallels([-60],linewidth=1.)
#cs = m.contourf(x,y, p4, levels=[0.9, 0.95] ,colors='none',hatches=['++', None],alpha=0)
# m.drawcoastlines(color='blue')
# cax = fig.add_axes([0.2, 0.02, 0.6, 0.1], aspect=0.05)
# draw_latlon_polygon(m, [180,230], [-70,-70], '-.',color='black',linewidth=2,alpha=1)
# cax1 = fig.add_axes([0.3, 0.49, 0.4, 0.08],aspect=0.03)
# cb1 = fig.colorbar(im1,orientation='horizontal',ticks=tickmarks,extend='both', cax = cax1)
# cb1.ax.tick_params(labelsize=16)
#cb.dividers.set_color('k')
#cb.dividers.set_linewidth(1)
#cb.dividers.set_linestyle('dotted')
# cb1.set_label('Sea ice advection ' + r'(% km day$^{\text{-1}}$)')
#fig.suptitle('temperature [K]', fontsize=14)
#fig.tight_layout(rect=[0,0.15,1,0.95])
# fig.tight_layout(rect=[0,0.08,1,0.95])
# ax1.text(0, 1.05, '(c)', fontsize=12, transform=ax1.transAxes, va='top', ha='right')
#plt.savefig('fig/Reg_asl_t2m.png' ,dpi=300,bbox_inches='tight')
# plt.savefig('/stu02/weizx24/figures/0924/Figure8d_obs_ice_advection.png' ,dpi=300,bbox_inches='tight')
#endregion 

#region 模式
#————————模式————————————
# 筛选最大最小年份
CMIP6_opwa = np.load('/stu02/weizx24/data/npz/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']*1e-12
ross_sie_cmip6 = np.nanmean(CMIP6_opwa,axis=1)
Q1 = np.quantile(ross_sie_cmip6,0.05)
Q3 = np.quantile(ross_sie_cmip6,0.95)
sic_Feb = xr.open_dataset('/stu02/weizx24/data/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_Feb.nc')['siconc']
year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year
siconc_mon = xr.open_dataset('/stu02/weizx24/data/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912.nc')['siconc']
siv = xr.open_dataset('/stu02/weizx24/data/siv_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['siv']
adv = (siv * siconc_mon)
# .values.reshape(499,12,384,320)
adv_clm = np.nanmean(adv.values.reshape(499,12,384,320),axis=0)
adv_std = np.nanstd(adv.values.reshape(499,12,384,320),axis=0)
# year_min = adv.time[ross_sie_cmip6<Q1].dt.year
adv_min_lst = []
for i in range(len(year_min)):
    y = str(year_min[i].values).zfill(4)
#     print(adv.sel(time=y))
    adv_min_lst.append(adv.sel(time=y))
adv_min_mean = np.array(adv_min_lst).squeeze().mean(axis=0)
adv_min_std = np.array(adv_min_lst).squeeze().std(axis=0)

# year_max = adv.time[ross_sie_cmip6<Q1].dt.year
adv_max_lst = []
for i in range(len(year_max)):
    y = str(year_max[i].values).zfill(4)
#     print(y)
    adv_max_lst.append(adv.sel(time=y))
adv_max_mean = np.array(adv_max_lst).squeeze().mean(axis=0)
adv_max_std = np.array(adv_max_lst).squeeze().std(axis=0)
adv_dif = np.nanmean(adv_max_mean[9:11],axis=0)-np.nanmean(adv_min_mean[9:11],axis=0)
num_y = 25

pre_clm_son = np.nanmean(adv_max_mean[9:11],axis=0)
pre_std_son = np.nanmean(adv_max_std[9:11],axis=0)

slp_clm_son = np.nanmean(adv_clm[9:11],axis=0)
slp_std_son = np.nanmean(adv_std[9:11],axis=0)

post_clm_son = np.nanmean(adv_min_mean[9:11],axis=0)
post_std_son = np.nanmean(adv_min_std[9:11],axis=0)

pval1=np.zeros((384,320),'float')
pval2=np.zeros((384,320),'float')
pval3=np.zeros((384,320),'float')

for i in range(384):
    for j in range(320):
#         pval1[i,j] = independent_ttest(pre_clm_son[i,j], slp_clm_son[i,j], pre_std_son[i,j],slp_std_son[i,j],num_y, 30, 0.1)
#         pval2[i,j] = independent_ttest(post_clm_son[i,j], slp_clm_son[i,j], post_std_son[i,j],slp_std_son[i,j], num_y, 30, 0.1)
        pval3[i,j] = independent_ttest(pre_clm_son[i,j], post_clm_son[i,j], pre_std_son[i,j],post_std_son[i,j], num_y,num_y,0.1)
    
lons = siconc_mon['lon']
lats = siconc_mon['lat']


levels1=np.linspace(-500,500,50)
tickmarks = np.linspace(-500,500,5)
# plt.close()
# plt.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['hatch.color'] ='springgreen'
# mpl.rcParams['hatch.linewidth'] = 0.6
#x1, y1 = m(*np.meshgrid(lon,lat))
# fig=plt.figure(figsize=(11,7.5))
# fig.subplots_adjust(top=0.99)
ax = fig.add_subplot(212)
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=180,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-130,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
x1,y1 = m(lons,lats)
im1 = m.contourf(x1,y1,adv_dif*100,
                 levels=levels1,
                 extend='both',
                shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)
cs = m.contourf(x1,y1, 1-pval3, levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
# Q1 = m.quiver(xx,yy,uproj*-1,vproj*-1,scale=11,color='purple',scale_units='inches')
# qk1 = plt.quiverkey(Q1, 0.2, -0.025, 2, '2 km/day', labelpos='W',fontproperties={'size':'large'})
m.drawmeridians([90,120,150,180,210,240,270,300,330], linewidth=1.2,labels=[True,False,False,True],fmt='%g',latmax=80)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,False],latmax=80,)
ax.annotate(r'$70\degree S$',xy=m(190,-70),xycoords='data',xytext = m(190,-70),textcoords='data')
ax.annotate(r'$60\degree S$',xy=m(190,-60),xycoords='data',xytext = m(190,-60),textcoords='data')
ax.annotate(r'$80\degree S$',xy=m(190,-80),xycoords='data',xytext = m(190,-80),textcoords='data')

#ax.annotate(r'$70\degree S$',fontsize=16,xy=m(10,-70),xycoords='data',xytext = m(10,-70),textcoords='data')
#ax.annotate(r'$60\degree S$',fontsize=16,xy=m(10,-60),xycoords='data',xytext = m(10,-60),textcoords='data')
#ax.annotate(r'$80\degree S$',fontsize=16,xy=m(10,-80),xycoords='data',xytext = m(10,-80),textcoords='data')
# ax.set_title('CESM2-WACCM-FV2')

# m.drawmeridians([160,230,300,20,90], linewidth=0.5, dashes=[1, 5],labels=[True,True,True, True])
m.drawcoastlines(color='blue')
m.fillcontinents(color='gray')
#cf = m.contour(x,y, fclm4,levels=np.array([980,986]),colors=['black'],alpha=0.5)
#m.drawparallels([-60],linewidth=1.)
#cs = m.contourf(x,y, p4, levels=[0.9, 0.95] ,colors='none',hatches=['++', None],alpha=0)
# m.drawcoastlines(color='blue')
# cax = fig.add_axes([0.2, 0.02, 0.6, 0.1], aspect=0.05)
# draw_latlon_polygon(m, [180,230], [-70,-70], '-.',color='black',linewidth=2,alpha=1)
cax = fig.add_axes([0.3, 0.02, 0.4, 0.08], aspect=0.03)
cb = fig.colorbar(im1,orientation='horizontal',ticks=tickmarks,extend='both', cax = cax)
# cb.ax.tick_params(labelsize=16)
cb.set_label('Sea ice advection ' + r'(% km day$^{\text{-1}}$)')
# fig.suptitle('Composite SLP & UV10(ON)', fontsize=24)
# fig.tight_layout(rect=[0.02,0.1,0.98,0.95])
ax.text(0, 1.05, '(h)', fontsize=20, transform=ax.transAxes, va='top', ha='right')
# plt.savefig('/stu02/weizx24/figures/0924/Figure8a_obs_slp_uv.png' ,dpi=300,bbox_inches='tight')
# plt.show()
# cb.ax.tick_params(labelsize=16)
#cb.dividers.set_color('k')
#cb.dividers.set_linewidth(1)
#cb.dividers.set_linestyle('dotted')
#fig.suptitle('temperature [K]', fontsize=14)
#fig.tight_layout(rect=[0,0.15,1,0.95])
# fig.tight_layout(rect=[0,0.08,1,0.95])
#plt.savefig('fig/Reg_asl_t2m.png' ,dpi=300,bbox_inches='tight') 
# ax.text(0, 1.05, '(g)', fontsize=12, transform=ax1.transAxes, va='top', ha='right')
plt.subplots_adjust(hspace=0.15)
plt.savefig('/stu02/weizx24/figures/0924/Figure8/Figure8dh_seaiceadvection.png' ,dpi=300,bbox_inches='tight')
#endregion
 