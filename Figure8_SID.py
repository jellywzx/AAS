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


#region 海冰漂移
#——————————————————————海冰漂移——————————————————————

fpath='/stu02/weizx24/data/monthly_seaicemotion/'
lat2 = xr.open_dataset(fpath+'ice_drift_vice_r720x361_sh_197811-202112.nc')['lat'].values[:]
lon2 = xr.open_dataset(fpath+'ice_drift_vice_r720x361_sh_197811-202112.nc')['lon'].values[:]
vice = xr.open_dataset(fpath+'ice_drift_vice_r720x361_sh_197811-202112.nc')['v'].loc['1979-01-01':'2021-12-31'].values[:,:,:]*86.4
uice = xr.open_dataset(fpath+'ice_drift_uice_r720x361_sh_197811-202112.nc')['u'].loc['1979-01-01':'2021-12-31'].values[:,:,:]*86.4

uu = np.nanmean(uice.reshape(43,12,81,720)[:,9:11],axis=1)
vv = np.nanmean(vice.reshape(43,12,81,720)[:,9:11],axis=1)

uu_clm = np.nanmean(uu,axis=0)
vv_clm = np.nanmean(vv,axis=0)

year_min = [1994., 1997., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2011., 2016., 2018., 2021.]
num_y = 6
year = np.linspace(1979,2021,43)

uu_max_list = []
for i in range(6):
    index = np.where(year_max[i] == year)
#     print(index[0][0])
    uu_max_list.append(uu[index])
vv_max_list = []
for i in range(6):
    index = np.where(year_max[i] == year)
    vv_max_list.append(vv[index])
uu_min_list = []
for i in range(6):
    index = np.where(year_min[i] == year)
    uu_min_list.append(uu[index])
vv_min_list = []
for i in range(6):
    index = np.where(year_min[i] == year)
    vv_min_list.append(vv[index])
uu_min = np.nanmean(np.array(uu_min_list),axis=0)
vv_min = np.nanmean(np.array(vv_min_list),axis=0)
uu_max = np.nanmean(np.array(uu_max_list),axis=0)
vv_max = np.nanmean(np.array(vv_max_list),axis=0)

uu_dif =(uu_max-uu_min).squeeze()
vv_dif = (vv_max-vv_min).squeeze()

fpath='/stu02/weizx24/data/monthly_seaicemotion/'

#icev = xr.open_dataset(fpath2+'ice_drift_vice_r720x361_sh_197811-202112.nc')['v'].loc['1979-01-01':'2021-12-31']*86.4
icev = xr.open_dataset(fpath+'ice_drift_vice_r720x361_sh_197811-202112.nc')['v'].loc['1979-01-01':'2021-12-31'].values[:,:,:]*86.4

#------读取SIC数据----------
fpath2='/stu02/weizx24/data/'
sic = xr.open_dataset(fpath2+'NSIDC_seaice_conc_monthly_0p5_0p5_197901_202112_sh.nc')['nsidc_nt_seaice_conc_monthly'].loc['1979-01-01':'2021-12-31']

adv = np.multiply(sic,icev)
#adv =  icev.values
lat1 = sic.lat.values
lon1 = sic.lon.values

sic_ON = np.nanmean(sic[(sic.time.dt.month==10)|(sic.time.dt.month==11)].values.reshape(43,2,81,720),axis=1)
sic_ON_clm = np.nanmean(sic_ON,axis=0)

# sic_std = np.nanstd(sic[(sic.time.dt.month==10)|(sic.time.dt.month==11)].values.reshape(43,2,81,720),axis=0)
# sic_ON_std = np.nanmean(sic_std,axis=0)


sic_list = []
for i in range(6):
    start_time = str(int(year_max[i]))+'-10-01'
    end_time = str(int(year_max[i]))+'-11-30'
    sic_sel = sic.sel(time=slice(start_time,end_time)).mean('time')
    sic_list.append(sic_sel)
#     print(sic.sel(time=slice(start_time,end_time)))
sic_max = np.nanmean(np.array(sic_list),axis=0)
sic_std_max = np.nanstd(np.array(sic_list),axis=0)

sic_list = []
for i in range(6):
    start_time = str(int(year_min[i]))+'-10-01'
    end_time = str(int(year_min[i]))+'-11-30'
    sic_sel = sic.sel(time=slice(start_time,end_time)).mean('time')
    sic_list.append(sic_sel)
#     print(sic.sel(time=slice(start_time,end_time)))
sic_min = np.nanmean(np.array(sic_list),axis=0)
sic_std_min = np.nanstd(np.array(sic_list),axis=0)

sic_dif = sic_max-sic_min


#------海冰密集度的显著性计算------
num_y = 6
pre_clm_son = sic_max
pre_std_son = sic_std_max
# slp_clm_son = adv_clm
# slp_std_son = adv_std
post_clm_son = sic_min
post_std_son = sic_std_min

# pval1=np.zeros((81,720),'float')
# pval2=np.zeros((81,720),'float')
pval3=np.zeros((81,720),'float')

for i in range(81):
    for j in range(720):
#         pval1[i,j] = independent_ttest(pre_clm_son[i,j], slp_clm_son[i,j], pre_std_son[i,j],slp_std_son[i,j],num_y, 30, 0.1)
#         pval2[i,j] = independent_ttest(post_clm_son[i,j], slp_clm_son[i,j], post_std_son[i,j],slp_std_son[i,j], num_y, 30, 0.1)
        pval3[i,j] = independent_ttest(pre_clm_son[i,j], post_clm_son[i,j], pre_std_son[i,j],post_std_son[i,j], num_y,num_y,0.1)
        
#-------作图----------------

levels1=np.linspace(-40,40,50)
tickmarks = np.linspace(-40,40,5)
#aa1[np.where(np.abs(aa1)<0.2)] = np.nan
ugrid1,newlons = shiftgrid(180.,uu_dif,lon2,start=False)
vgrid1,newlons = shiftgrid(180.,vv_dif,lon2,start=False)

plt.close()
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 18
plt.rcParams['lines.linewidth'] = 2 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6
plt.rcParams['axes.titlesize'] = 20  # 轴标题字体大小


#x1, y1 = m(*np.meshgrid(lon,lat))
fig=plt.figure(figsize=(7,10))
# fig.subplots_adjust(top=0.99)
ax1 = fig.add_subplot(211)
m = Basemap(projection='stere',resolution='c',
            lat_0=-90, lon_0=180,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-130,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
x1,y1 = m(*np.meshgrid(lon1,lat1))
uproj,vproj,xx,yy = m.transform_vector(ugrid1,vgrid1,newlons,lat2,31,31,returnxy=True,masked=True)
im1 = m.contourf(x1,y1,sic_dif*100,levels=levels1,extend='both',
                shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)
Q1 = m.quiver(xx,yy,uproj*-1,vproj*-1,scale=15,color='purple',scale_units='inches')
qk1 = plt.quiverkey(Q1, 0.16, -0.03, 5, '5 km/day', labelpos='W')
cs = m.contourf(x1,y1, 1-pval3, levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
m.drawmeridians([90,120,150,180,210,240,270,300,330], labels=[True,False,False,True],fmt='%g',latmax=80,)
m.drawparallels([-60,-70,-80], labels=[False,False,False,False],latmax=80)
ax1.annotate(r'$70\degree S$',xy=m(190,-70),xycoords='data',xytext = m(190,-70),textcoords='data')
ax1.annotate(r'$60\degree S$',xy=m(190,-60),xycoords='data',xytext = m(190,-60),textcoords='data')
ax1.annotate(r'$80\degree S$',xy=m(190,-80),xycoords='data',xytext = m(190,-80),textcoords='data')
# ax1.set_title('Observation')
m.drawcoastlines(color='blue')
m.fillcontinents(color='gray')
# draw_latlon_polygon(m, [180,230], [-70,-70], '-.',color='black',linewidth=2,alpha=1)
ax1.text(0,1.05,'(c)', transform=ax1.transAxes, va='top', ha='right')
# cax = fig.add_axes([0.3, 0.49, 0.4, 0.08],aspect=0.03)
# cb = fig.colorbar(im1,orientation='horizontal',ticks=tickmarks,extend='both', cax = cax)
# cb.ax.tick_params(labelsize=16)
#cb.dividers.set_color('k')
#cb.dividers.set_linewidth(1)
#cb.dividers.set_linestyle('dotted')
# cb.set_label('Sea ice concentration ' + '(%)')
#fig.tight_layout(rect=[0,0.15,1,0.95])
# fig.tight_layout(rect=[0,0.08,1,0.95])
#plt.savefig('fig/Reg_asl_t2m.png' ,dpi=300,bbox_inches='tight')
# plt.savefig('/stu02/weizx24/figures/0924/Figure8_obs_seaicemotion.png' ,dpi=600,bbox_inches='tight')

# plt.show()
#endregion

#region 模式
#————————————————————模式———————————————————————————
# 筛选最大最小年份
CMIP6_opwa = np.load('/stu02/weizx24/data/npz/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']*1e-12
ross_sie_cmip6 = np.nanmean(CMIP6_opwa,axis=1)
Q1 = np.quantile(ross_sie_cmip6,0.05)
Q3 = np.quantile(ross_sie_cmip6,0.95)
sic_Feb = xr.open_dataset('/stu02/weizx24/data//siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_Feb.nc')['siconc']
year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year

#海冰漂移数据
siu = xr.open_dataset('/stu02/weizx24/data/siu_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['siu']*86.4
siv = xr.open_dataset('/stu02/weizx24/data/siv_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['siv']*86.4

lons = siu.lon
lats = siu.lat

siu_ON = np.nanmean(siu[(siu.time.dt.month==10)|(siu.time.dt.month==11)].values.reshape(499,2,384,320),axis=1)
siv_ON = np.nanmean(siv[(siv.time.dt.month==10)|(siv.time.dt.month==11)].values.reshape(499,2,384,320),axis=1)

siu_max = np.nanmean(siu_ON[ross_sie_cmip6>Q3],axis=0)
siu_min = np.nanmean(siu_ON[ross_sie_cmip6<Q1],axis=0)

siv_max = np.nanmean(siv_ON[ross_sie_cmip6>Q3],axis=0)
siv_min = np.nanmean(siv_ON[ross_sie_cmip6<Q1],axis=0)

siu_dif = siu_max-siu_min
siv_dif = siv_max-siv_min

#海冰密集度
sic = xr.open_dataset('/stu02/weizx24/data/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912.nc')['siconc']
# sic
lon2 = sic.lon
lat2 = sic.lat

sic_ON = np.nanmean(sic[(sic.time.dt.month==10)|(sic.time.dt.month==11)].values.reshape(499,2,384,320),axis=1)
sic_max = np.nanmean(sic_ON[ross_sie_cmip6>Q3],axis=0) 
sic_min = np.nanmean(sic_ON[ross_sie_cmip6<Q1],axis=0)

sic_dif = sic_max-sic_min

sic_max_mean = np.nanmean(sic_ON[ross_sie_cmip6>Q3],axis=0) 
sic_max_std = np.nanstd(sic_ON[ross_sie_cmip6>Q3],axis=0) 
sic_min_mean = np.nanmean(sic_ON[ross_sie_cmip6<Q1],axis=0)
sic_min_std = np.nanstd(sic_ON[ross_sie_cmip6<Q1],axis=0)

pre_clm_son = sic_max_mean
post_clm_son = sic_min_mean
pre_std_son = sic_max_std
post_std_son = sic_min_std
num_y = 25

pval3=np.zeros((384,320),'float')

for i in range(384):
    for j in range(320):
#         pval1[i,j] = independent_ttest(pre_clm_son[i,j], slp_clm_son[i,j], pre_std_son[i,j],slp_std_son[i,j],num_y, 30, 0.1)
#         pval2[i,j] = independent_ttest(post_clm_son[i,j], slp_clm_son[i,j], post_std_son[i,j],slp_std_son[i,j], num_y, 30, 0.1)
        pval3[i,j] = independent_ttest(pre_clm_son[i,j], post_clm_son[i,j], pre_std_son[i,j],post_std_son[i,j], num_y,num_y,0.1)


levels1=np.linspace(-40,40,50)
tick_marks = np.linspace(-40,40,5)
# plt.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['hatch.color'] ='springgreen'
# mpl.rcParams['hatch.linewidth'] = 0.6
m = Basemap(projection='stere',resolution='c',
            lat_0=-90, lon_0=180,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-130,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
#-----------------------------第一张图------------------------------------------
# plt.close()
# fig=plt.figure(figsize=(11,7.5))
# fig.subplots_adjust(top=0.99)
ax2 = fig.add_subplot(212)
x, y = m(lons,lats)
us_wdir,vs_wdir,xs,ys = m.rotate_vector(siu_dif,siv_dif,lons,lats,returnxy=True)
# ws_dif = np.sqrt(siu_dif**2+siv_dif**2)
im = m.contourf(x,y,sic_dif,
                levels=levels1 ,
                extend='both',
                shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)
# ugrid1,newlons = shiftgrid(180.,siu_dif,lons,start=False)
# vgrid1,newlons = shiftgrid(180.,siv_dif,lons,start=False) 
# uproj,vproj,xx,yy = m.transform_vector(ugrid1,vgrid1,newlons,lats,25,25,returnxy=True,masked=True)
# Q1 = m.quiver(xx,yy,uproj*-1,vproj*-1,scale=6,color='black',scale_units='inches')
q1 = m.quiver(xs[::5,::5],ys[::5,::5],us_wdir[::5,::5]*-1,vs_wdir[::5,::5]*-1,scale=15, color='purple',scale_units='inches')
# q1 = m.quiver(xs,ys,us_wdir*-1,vs_wdir*-1,scale=15, color='purple',scale_units='inches')
qk1 = plt.quiverkey(q1, 0.16, -0.03, 5, '5 km/day', labelpos='W',color='purple',)
m.drawmeridians([90,120,150,180,210,240,270,300,330], labels=[True,False,False,True],fmt='%g',latmax=80)
m.drawparallels([-60,-70,-80], labels=[False,False,False,False],latmax=80)
ax2.annotate(r'$70\degree S$',xy=m(190,-70),xycoords='data',xytext = m(190,-70),textcoords='data')
ax2.annotate(r'$60\degree S$',xy=m(190,-60),xycoords='data',xytext = m(190,-60),textcoords='data')
ax2.annotate(r'$80\degree S$',xy=m(190,-80),xycoords='data',xytext = m(190,-80),textcoords='data')
cs = m.contourf(x,y, 1-pval3, levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
# ax2.set_title('CESM2-WACCM-FV2')
m.drawcoastlines(color='blue')
m.fillcontinents(color='gray')
#------------------------------公共部分-------------------------
cax = fig.add_axes([0.3, 0.02, 0.4, 0.08], aspect=0.03)
cb = fig.colorbar(im,orientation='horizontal',ticks=tick_marks,extend='both',cax = cax)
# cb.ax.tick_params(labelsize=16)
cb.set_label("Sea ice concentration (%)")
ax2.text(0,1.05,'(g)',transform=ax2.transAxes, va='top', ha='right')
# fig.suptitle('Composite SLP & UV10(ON)', fontsize=24)
# fig.tight_layout(rect=[0.02,0.1,0.98,0.95])
plt.subplots_adjust(hspace=0.15)
plt.savefig('/stu02/weizx24/figures/0924/Figure8/Figure8cg_SID.png' ,dpi=300,bbox_inches='tight')
plt.show()
#endregion 
