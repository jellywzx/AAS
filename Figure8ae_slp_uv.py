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


#region 气压 no.1

#————————————————————观测数据——————————————————————————————
#------风场-------
u10_sel = xr.open_dataset('/stu02/weizx24/data/ERA5/ERA5-1x1-resolution-uv10-197901-202212.nc')['u10'][:,::-1,:].sel(time=slice('1992','2021'))
v10_sel = xr.open_dataset('/stu02/weizx24/data/ERA5/ERA5-1x1-resolution-uv10-197901-202212.nc')['v10'][:,::-1,:].sel(time=slice('1992','2021'))

lons = u10_sel.longitude.values
lats = u10_sel.latitude.values


#合成年份为5年
year_min = [1994., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2016., 2018., 2021.]
num_y = 5

u10_sel = xr.open_dataset('/stu02/weizx24/data/ERA5/ERA5-1x1-resolution-uv10-197901-202212.nc')['u10'][:,::-1,:].sel(time=slice('1992','2021'))
v10_sel = xr.open_dataset('/stu02/weizx24/data/ERA5/ERA5-1x1-resolution-uv10-197901-202212.nc')['v10'][:,::-1,:].sel(time=slice('1992','2021'))

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

#-------气压--------
slp = xr.open_dataset('/stu02/weizx24/data/ERA5/79-21_msl.nc')['msl'].sel(latitude=slice(-30,-90)).sel(time=slice('1992','2021'))
lat = slp.latitude.values
lon = slp.longitude.values
t2mn1 = np.zeros((num_y,12,161,1440),'float')
t2mn2 = np.zeros((num_y,12,161,1440),'float')
slp_clm = slp.values.reshape(30,12,161,1440).mean(axis=0)*1e-2
slp_std = slp.values.reshape(30,12,161,1440).std(axis=0)*1e-2
  
for i in range(num_y):
    t2mn1[i,:,:,:] = slp[(slp['time.year']==year_max[i])].values.reshape(12,161,1440)*1e-2
for i in range(num_y):
    t2mn2[i,:,:,:] = slp[(slp['time.year']==year_min[i])].values.reshape(12,161,1440)*1e-2
pre_clm_son = t2mn1.mean(axis=0)[9:11].mean(axis=0)
pre_std_son = t2mn1.std(axis=0)[9:11].mean(axis=0)
post_clm_son  = t2mn2.mean(axis=0)[9:11].mean(axis=0)
post_std_son  = t2mn2.std(axis=0)[9:11].mean(axis=0)
slp_clm_son   = slp_clm[9:11].mean(axis=0)
slp_std_son   = slp_std[9:11].mean(axis=0)

comp1 = t2mn1.mean(axis=0)  - slp_clm
comp2 = t2mn2.mean(axis=0)  - slp_clm
comp3 = t2mn1.mean(axis=0)  - t2mn2.mean(axis=0)
corSON1 = comp1[9:11].mean(axis=0)
corSON2 = comp2[9:11].mean(axis=0)
corSON3 = comp3[9:11].mean(axis=0)

pval1=np.zeros((161,1440),'float')
pval2=np.zeros((161,1440),'float')
pval3=np.zeros((161,1440),'float')

for i in range(161):
    for j in range(1440):
        pval1[i,j] = independent_ttest(pre_clm_son[i,j], slp_clm_son[i,j], pre_std_son[i,j],slp_std_son[i,j],num_y, 30, 0.1)
        pval2[i,j] = independent_ttest(post_clm_son[i,j], slp_clm_son[i,j], post_std_son[i,j],slp_std_son[i,j], num_y, 30, 0.1)
        pval3[i,j] = independent_ttest(pre_clm_son[i,j], post_clm_son[i,j], pre_std_son[i,j],post_std_son[i,j], num_y,num_y,0.1)


levels1=np.linspace(-15,15,50)
tickmarks = np.linspace(-15,15,5)
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
ax1 = fig.add_subplot(211)
x, y = m(*np.meshgrid(lon,lat))
im = m.contourf(x,y,corSON3,levels=levels1 ,extend='both',
                shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)

ugrid1,newlons = shiftgrid(180.,u_anom_diff.mean(axis=0),lons,start=False)
vgrid1,newlons = shiftgrid(180.,v_anom_diff.mean(axis=0),lons,start=False) 
uproj,vproj,xx,yy = m.transform_vector(ugrid1,vgrid1,newlons,lats,25,25,returnxy=True,masked=True)
Q1 = m.quiver(xx,yy,uproj*-1,vproj*-1,scale=10,color='black',scale_units='inches',width=0.005)
# qk1 = plt.quiverkey(Q1, 0.16, -0.05, 2, '2 m/s', labelpos='W',color='k',)

m.drawmeridians([90,120,150,180,210,240,270,300,330],labels=[True,False,False,False],fmt='%g',latmax=80)
m.drawparallels([-60,-70,-80], labels=[False,False,False,False],latmax=80)
ax1.annotate(r'$70\degree\,\mathrm{S}$',xy=m(190,-70),xycoords='data')
ax1.annotate(r'$80\degree\,\mathrm{S}$',xy=m(190,-80),xycoords='data')
cs = m.contourf(x,y, 1-pval3, levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
# ax1.set_title('Observation')
m.drawcoastlines(color='blue')
ax1.text(0, 1.05, '(a)', transform=ax1.transAxes, va='top', ha='right')
#endregion
        


#region 模式
#————————————————模式数据——————————————————————————
# 筛选最大最小年份
CMIP6_opwa = np.load('/stu02/weizx24/data/npz/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']*1e-12
ross_sie_cmip6 = np.nanmean(CMIP6_opwa,axis=1)
Q1 = np.quantile(ross_sie_cmip6,0.05)
Q3 = np.quantile(ross_sie_cmip6,0.95)
sic_Feb = xr.open_dataset('/stu02/weizx24/data/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_Feb.nc')['siconc']
year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year



file_path = '/stu02/weizx24/data/psl_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc'
psl = xr.open_dataset(file_path)['psl']*1e-2
lon = psl.lon
lat = psl.lat

#海平面气压
psl_ON = psl[(psl.time.dt.month==10)|(psl.time.dt.month==11)].values.reshape(499,2,96,144).mean(axis=1)
psl_clm = psl_ON.mean(axis=0)

psl_max = psl_ON[ross_sie_cmip6>Q3].mean(axis=0)
psl_min = psl_ON[ross_sie_cmip6<Q1].mean(axis=0)

psl_max_mean = psl_ON[ross_sie_cmip6>Q3].mean(axis=0)
psl_min_mean = psl_ON[ross_sie_cmip6<Q1].mean(axis=0)
psl_max_std = psl_ON[ross_sie_cmip6>Q3].std(axis=0)
psl_min_std = psl_ON[ross_sie_cmip6<Q1].std(axis=0)

#风
ua = xr.open_dataset('/stu02/weizx24/data/ua_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['ua'].sel(plev=92500)
va = xr.open_dataset('/stu02/weizx24/data/va_Amon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all.nc')['va'].sel(plev=92500)

lon2 = va.lon.values
lat2 = va.lat.values

ua_ON = ua[(ua.time.dt.month==10)|(ua.time.dt.month==11)].values.reshape(499,2,96,144).mean(axis=1)
ua_clm = ua_ON.mean(axis=0)

ua_max = ua_ON[ross_sie_cmip6>Q3].mean(axis=0)
ua_min = ua_ON[ross_sie_cmip6<Q1].mean(axis=0)

ua_max_mean = ua_ON[ross_sie_cmip6>Q3].mean(axis=0)
ua_min_mean = ua_ON[ross_sie_cmip6<Q1].mean(axis=0)
ua_max_std = ua_ON[ross_sie_cmip6>Q3].std(axis=0)
ua_min_std = ua_ON[ross_sie_cmip6<Q1].std(axis=0)

va_ON = va[(va.time.dt.month==10)|(va.time.dt.month==11)].values.reshape(499,2,96,144).mean(axis=1)
va_clm = va_ON.mean(axis=0)

va_max = va_ON[ross_sie_cmip6>Q3].mean(axis=0)
va_min = va_ON[ross_sie_cmip6<Q1].mean(axis=0)

va_max_mean = va_ON[ross_sie_cmip6>Q3].mean(axis=0)
va_min_mean = va_ON[ross_sie_cmip6<Q1].mean(axis=0)
va_max_std = va_ON[ross_sie_cmip6>Q3].std(axis=0)
va_min_std = va_ON[ross_sie_cmip6<Q1].std(axis=0)


#计算psl的显著性
pre_clm_son = psl_max_mean
post_clm_son = psl_min_mean
pre_std_son = psl_max_std
post_std_son = psl_min_std
num_y = 25

pval3=np.zeros((96,144),'float')

for i in range(96):
    for j in range(144):
        pval3[i,j] = independent_ttest(pre_clm_son[i,j], post_clm_son[i,j], pre_std_son[i,j],post_std_son[i,j], num_y,num_y,0.1)


#-------------开始作图----------------------
levels1=np.linspace(-15,15,50)
tick_marks = np.linspace(-15,15,5)
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=180,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-130,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
#-----------------------------第一张图------------------------------------------
ax2 = fig.add_subplot(212)

x, y = m(*np.meshgrid(lon,lat))
im = m.contourf(x,y,psl_max-psl_min,
                levels=levels1,
                extend='both',
                shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)
ugrid1,newlons = shiftgrid(180.,ua_max-ua_min,lon2,start=False)
vgrid1,newlons = shiftgrid(180.,va_max-va_min,lon2,start=False) 
uproj,vproj,xx,yy = m.transform_vector(ugrid1,vgrid1,newlons,lat2,25,25,returnxy=True,masked=True)
# Q1 = m.quiver(xx[::2,::2],yy[::2,::2],uproj[::2,::2]*-1,vproj[::2,::2]*-1,scale=7,color='black',scale_units='inches')
Q1 = m.quiver(xx,yy,uproj*-1,vproj*-1,scale=10,color='black',scale_units='inches', width=0.005)
qk1 = plt.quiverkey(Q1, 0.16, -0.05, 2, '2 m/s', labelpos='W',color='k',)
m.drawmeridians([90,120,150,180,210,240,270,300,330], linewidth=1.2,labels=[True,False,False,True],fmt='%g',latmax=80,)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,False],latmax=80,)
# ax2.annotate(r'$70\degree\,\mathrm{S}$',xy=m(190,-70),xycoords='data')
# ax2.annotate(r'$80\degree\,\mathrm{S}$',xy=m(190,-80),xycoords='data')
cs = m.contourf(x,y, 1-pval3, levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
# ax2.set_title('CESM2-WACCM-FV2')
m.drawcoastlines(color='blue')

cax = fig.add_axes([0.3, 0.02, 0.4, 0.08], aspect=0.03)
cb = fig.colorbar(im,orientation='horizontal',
                  ticks=tick_marks,
                  extend='both',
                  cax = cax,
                  )
cb.ax.tick_params(labelsize=22)
cb.set_label("SLP (hPa)")
ax2.text(0, 1.05, '(e)', transform=ax2.transAxes, va='top', ha='right')
plt.subplots_adjust(hspace=0.05)
plt.savefig('/stu02/weizx24/figures/0924/Figure8/Figure8ae_slp_uv.png' ,dpi=300,bbox_inches='tight')
# plt.show()
#endregion


