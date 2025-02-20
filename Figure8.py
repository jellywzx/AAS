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
plt.rcParams['font.size'] = '9'
plt.rcParams['lines.linewidth'] = 2 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6

#x1, y1 = m(*np.meshgrid(lon,lat))
fig=plt.figure(figsize=(6.3,4))
fig.subplots_adjust(top=0.99)
ax1 = fig.add_subplot(241)
m = Basemap(projection='stere',resolution='c',
            lat_0=-90, lon_0=180,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-130,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
x1,y1 = m(*np.meshgrid(lon1,lat1))
uproj,vproj,xx,yy = m.transform_vector(ugrid1,vgrid1,newlons,lat2,31,31,returnxy=True,masked=True)
im1 = m.contourf(x1,y1,sic_dif*100,levels=levels1,extend='both',
                shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)
Q1 = m.quiver(xx,yy,uproj*-1,vproj*-1,scale=15,color='purple',scale_units='inches')
qk1 = plt.quiverkey(Q1, 0.16, -0.03, 5, '5 km/day', labelpos='W',fontproperties={'size':'large'})
cs = m.contourf(x1,y1, 1-pval3, levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
m.drawmeridians([90,120,150,180,210,240,270,300,330], linewidth=1.2,labels=[True,False,False,True],fmt='%g',latmax=80,fontsize=15)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,False],latmax=80,fontsize=15)
ax1.annotate(r'$70\degree S$',fontsize=15,xy=m(190,-70),xycoords='data',xytext = m(190,-70),textcoords='data')
ax1.annotate(r'$60\degree S$',fontsize=15,xy=m(190,-60),xycoords='data',xytext = m(190,-60),textcoords='data')
ax1.annotate(r'$80\degree S$',fontsize=15,xy=m(190,-80),xycoords='data',xytext = m(190,-80),textcoords='data')
ax1.set_title('SID (ON)', fontsize=24)
m.drawcoastlines(color='blue')
m.fillcontinents(color='gray')
draw_latlon_polygon(m, [180,230], [-70,-70], '-.',color='black',linewidth=2,alpha=1)
cax = fig.add_axes([0.24, 0, 0.55, 0.1], aspect=0.05)
cb = fig.colorbar(im1,orientation='horizontal',ticks=tickmarks,extend='both', cax = cax)
cb.ax.tick_params(labelsize=16)
#cb.dividers.set_color('k')
#cb.dividers.set_linewidth(1)
#cb.dividers.set_linestyle('dotted')
cb.set_label('SIC ' + '(%)',fontsize=18)
#fig.suptitle('temperature [K]', fontsize=14)
#fig.tight_layout(rect=[0,0.15,1,0.95])
# fig.tight_layout(rect=[0,0.08,1,0.95])
#plt.savefig('fig/Reg_asl_t2m.png' ,dpi=300,bbox_inches='tight')
plt.savefig('/stu02/weizx24/figures/Figure8_obs_seaicemotion.png' ,dpi=600,bbox_inches='tight')

# plt.show()



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
plt.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6
m = Basemap(projection='stere',resolution='c',
            lat_0=-90, lon_0=180,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-130,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
#-----------------------------第一张图------------------------------------------
plt.close()
fig=plt.figure(figsize=(11,7.5))
fig.subplots_adjust(top=0.99)
ax1 = fig.add_subplot(111)
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
qk1 = plt.quiverkey(q1, 0.16, -0.03, 5, '5 km/day', labelpos='W',color='purple',fontproperties={'size':'large'})
m.drawmeridians([90,120,150,180,210,240,270,300,330], linewidth=1.2,labels=[True,False,False,True],fmt='%g',latmax=80,fontsize=15)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,False],latmax=80,fontsize=15)
ax1.annotate(r'$70\degree S$',fontsize=15,xy=m(190,-70),xycoords='data',xytext = m(190,-70),textcoords='data')
ax1.annotate(r'$60\degree S$',fontsize=15,xy=m(190,-60),xycoords='data',xytext = m(190,-60),textcoords='data')
ax1.annotate(r'$80\degree S$',fontsize=15,xy=m(190,-80),xycoords='data',xytext = m(190,-80),textcoords='data')
cs = m.contourf(x,y, 1-pval3, levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
ax1.set_title('SID (ON)', fontsize=24)
m.drawcoastlines(color='blue')
m.fillcontinents(color='gray')
#------------------------------公共部分-------------------------
cax = fig.add_axes([0.24, 0, 0.55, 0.1], aspect=0.05)
cb = fig.colorbar(im,orientation='horizontal',ticks=tick_marks,extend='both', extendfrac='auto',cax = cax)
cb.ax.tick_params(labelsize=16)
cb.set_label("SIC (%)",fontsize=18)
# fig.suptitle('Composite SLP & UV10(ON)', fontsize=24)
# fig.tight_layout(rect=[0.02,0.1,0.98,0.95])
plt.savefig('/stu02/weizx24/figures/Figure8_CMIP6_SID.png' ,dpi=600,bbox_inches='tight')
# plt.show()
#endregion 
'''
'''
#region 气压

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
plt.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=180,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-130,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
#-----------------------------第一张图------------------------------------------
plt.close()
fig=plt.figure(figsize=(11,7.5))
fig.subplots_adjust(top=0.99)
ax1 = fig.add_subplot(111)
x, y = m(*np.meshgrid(lon,lat))
im = m.contourf(x,y,corSON3,levels=levels1 ,extend='both',
                shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)

ugrid1,newlons = shiftgrid(180.,u_anom_diff.mean(axis=0),lons,start=False)
vgrid1,newlons = shiftgrid(180.,v_anom_diff.mean(axis=0),lons,start=False) 
uproj,vproj,xx,yy = m.transform_vector(ugrid1,vgrid1,newlons,lats,25,25,returnxy=True,masked=True)
Q1 = m.quiver(xx,yy,uproj*-1,vproj*-1,scale=7,color='black',scale_units='inches')
qk1 = plt.quiverkey(Q1, 0.16, -0.03, 2, '2 m/s', labelpos='W',color='k',fontproperties={'size':'large'})

m.drawmeridians([90,120,150,180,210,240,270,300,330], linewidth=1.2,labels=[True,False,False,True],fmt='%g',latmax=80,fontsize=15)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,False],latmax=80,fontsize=15)
ax1.annotate(r'$70\degree S$',fontsize=15,xy=m(190,-70),xycoords='data',xytext = m(190,-70),textcoords='data')
ax1.annotate(r'$60\degree S$',fontsize=15,xy=m(190,-60),xycoords='data',xytext = m(190,-60),textcoords='data')
ax1.annotate(r'$80\degree S$',fontsize=15,xy=m(190,-80),xycoords='data',xytext = m(190,-80),textcoords='data')
cs = m.contourf(x,y, 1-pval3, levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
ax1.set_title('SLP & UV10 (ON)', fontsize=24)
m.drawcoastlines(color='blue')

#------------------------------公共部分-------------------------
cax = fig.add_axes([0.24, 0, 0.55, 0.1], aspect=0.05)
cb = fig.colorbar(im,orientation='horizontal',ticks=tickmarks,extend='both', cax = cax)
cb.ax.tick_params(labelsize=16)
cb.set_label("Sea level pressure (hPa)",fontsize=18)
# fig.suptitle('Composite SLP & UV10(ON)', fontsize=24)
# fig.tight_layout(rect=[0.02,0.1,0.98,0.95])
plt.savefig('/stu02/weizx24/figures/Figure8a_obs_slp_uv.png' ,dpi=300,bbox_inches='tight')
# plt.show()
        



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
#         pval1[i,j] = independent_ttest(pre_clm_son[i,j], slp_clm_son[i,j], pre_std_son[i,j],slp_std_son[i,j],num_y, 30, 0.1)
#         pval2[i,j] = independent_ttest(post_clm_son[i,j], slp_clm_son[i,j], post_std_son[i,j],slp_std_son[i,j], num_y, 30, 0.1)
        pval3[i,j] = independent_ttest(pre_clm_son[i,j], post_clm_son[i,j], pre_std_son[i,j],post_std_son[i,j], num_y,num_y,0.1)


#-------------开始作图----------------------
levels1=np.linspace(-15,15,50)
tick_marks = np.linspace(-15,15,5)
plt.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=180,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-130,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
#-----------------------------第一张图------------------------------------------
plt.close()
fig=plt.figure(figsize=(11,7.5))
fig.subplots_adjust(top=0.99)
ax1 = fig.add_subplot(111)

x, y = m(*np.meshgrid(lon,lat))
im = m.contourf(x,y,psl_max-psl_min,
                levels=levels1,
                extend='both',
                shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)

ugrid1,newlons = shiftgrid(180.,ua_max-ua_min,lon2,start=False)
vgrid1,newlons = shiftgrid(180.,va_max-va_min,lon2,start=False) 
uproj,vproj,xx,yy = m.transform_vector(ugrid1,vgrid1,newlons,lat2,25,25,returnxy=True,masked=True)
# Q1 = m.quiver(xx[::2,::2],yy[::2,::2],uproj[::2,::2]*-1,vproj[::2,::2]*-1,scale=7,color='black',scale_units='inches')
Q1 = m.quiver(xx,yy,uproj*-1,vproj*-1,scale=7,color='black',scale_units='inches')
qk1 = plt.quiverkey(Q1, 0.16, -0.03, 2, '2 m/s', labelpos='W',color='k',fontproperties={'size':'large'})

m.drawmeridians([90,120,150,180,210,240,270,300,330], linewidth=1.2,labels=[True,False,False,True],fmt='%g',latmax=80,fontsize=15)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,False],latmax=80,fontsize=15)
ax1.annotate(r'$70\degree S$',fontsize=15,xy=m(190,-70),xycoords='data',xytext = m(190,-70),textcoords='data')
ax1.annotate(r'$60\degree S$',fontsize=15,xy=m(190,-60),xycoords='data',xytext = m(190,-60),textcoords='data')
ax1.annotate(r'$80\degree S$',fontsize=15,xy=m(190,-80),xycoords='data',xytext = m(190,-80),textcoords='data')
cs = m.contourf(x,y, 1-pval3, levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
ax1.set_title('SLP & UV (ON)', fontsize=24)
m.drawcoastlines(color='blue')


# cax = fig.add_axes([0.2, 0.02, 0.6, 0.04], aspect=0.03)
# # draw_latlon_polygon(m, [180,230], [-70,-70], '-.',color='purple',linewidth=2,alpha=1)

# cb = fig.colorbar(im,orientation='horizontal',ticks =levels1[::2], boundaries=levels1[::2], extend='False', extendfrac=0.05, drawedges=False,cax = cax)
# cb.ax.tick_params(labelsize=15, width=1)


# cb = m.colorbar(im, location='bottom', pad="8%",ticks =levels1[::2],
# #                 boundaries=levels1[::2],
#                 extend='both', extendfrac='auto',drawedges=True, )
#-----------------------------第二张图------------------------------------------
# ax1 = fig.add_subplot(122)
# im = m.contourf(x,y,corSON2,levels=levels1 ,extend='both',
#                 shading='faceted', antialiased=True,cmap='BlueWhiteOrangeRed')
# ugrid1,newlons = shiftgrid(180.,post_u_anom.mean(axis=0),lons,start=False)
# vgrid1,newlons = shiftgrid(180.,post_v_anom.mean(axis=0),lons,start=False) 
# uproj,vproj,xx,yy = m.transform_vector(ugrid1,vgrid1,newlons,lats,25,25,returnxy=True,masked=True)
# Q1 = m.quiver(xx,yy,uproj*-1,vproj*-1,scale=4,color='black',scale_units='inches')
# qk1 = plt.quiverkey(Q1, 0.16, -0.05, 1, '1 m/s', labelpos='W',color='k')
# m.drawmeridians([90,120,150,180,210,240,270,300,330], linewidth=1.2,labels=[True,False,False,True],fmt='%g',latmax=80,fontsize=15)
# m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,False],latmax=80,fontsize=15)
# # m.drawmeridians([160,230,300,20,90], linewidth=0.5, dashes=[1, 5],labels=[True,True,True, True])
# ax1.annotate(r'$70\degree S$',fontsize=10,xy=m(190,-70),xycoords='data',xytext = m(190,-70),textcoords='data')
# ax1.annotate(r'$60\degree S$',fontsize=10,xy=m(190,-60),xycoords='data',xytext = m(190,-60),textcoords='data')
# ax1.annotate(r'$80\degree S$',fontsize=10,xy=m(190,-80),xycoords='data',xytext = m(190,-80),textcoords='data')
# cs = m.contourf(x,y, 1-pval2, levels=[0.9, 1] ,colors='none',hatches=['.', None],alpha=0)
# ax1.set_title('Low Polynya',fontsize=18)
# m.drawcoastlines(color='blue')
#------------------------------公共部分-------------------------
cax = fig.add_axes([0.24, 0, 0.55, 0.1], aspect=0.05)
cb = fig.colorbar(im,orientation='horizontal',
                  ticks=tick_marks,
                  extend='both',cax = cax)
cb.ax.tick_params(labelsize=16)
cb.set_label("Sea level pressure (hPa)",fontsize=18)
# fig.suptitle('Composite SLP & UV10(ON)', fontsize=24)
# fig.tight_layout(rect=[0.02,0.1,0.98,0.95])
plt.savefig('/stu02/weizx24/figures/Figure8a_CMIP6_slp_uv.png' ,dpi=300,bbox_inches='tight')
# plt.show()

#endregion
'''

'''
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
plt.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=180,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-130,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
#-----------------------------第一张图------------------------------------------
plt.close()
fig=plt.figure(figsize=(11,7.5))
fig.subplots_adjust(top=0.99)
ax1 = fig.add_subplot(111)
x, y = m(*np.meshgrid(lons,lats))
im = m.contourf(x,y,corSON3,levels=levels1 ,extend='both',
                shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)
m.drawmeridians([90,120,150,180,210,240,270,300,330], linewidth=1.2,labels=[True,False,False,True],fmt='%g',latmax=80,fontsize=15)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,False],latmax=80,fontsize=15)
ax1.annotate(r'$70\degree S$',fontsize=15,xy=m(190,-70),xycoords='data',xytext = m(190,-70),textcoords='data')
ax1.annotate(r'$60\degree S$',fontsize=15,xy=m(190,-60),xycoords='data',xytext = m(190,-60),textcoords='data')
ax1.annotate(r'$80\degree S$',fontsize=15,xy=m(190,-80),xycoords='data',xytext = m(190,-80),textcoords='data')
cs = m.contourf(x,y, 1-pval3, levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
ax1.set_title('2m temperature (ON)', fontsize=24)
m.drawcoastlines(color='blue')

#------------------------------公共部分-------------------------
cax = fig.add_axes([0.24, 0, 0.55, 0.1], aspect=0.05)
cb = fig.colorbar(im,orientation='horizontal',ticks=tickmarks,extend='both',cax = cax)
cb.ax.tick_params(labelsize=16)
cb.set_label("temperature (K)",fontsize=18)
fig.tight_layout(rect=[0.02,0.1,0.98,0.95])
plt.savefig('/stu02/weizx24/figures/Figure8b_obs_t2m.png' ,dpi=300,bbox_inches='tight')
# plt.show()


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
plt.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=180,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-130,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
#-----------------------------第一张图------------------------------------------
plt.close()
fig=plt.figure(figsize=(11,7.5))
fig.subplots_adjust(top=0.99)
ax1 = fig.add_subplot(111)
x, y = m(*np.meshgrid(lons,lats))
im = m.contourf(x,y,tas_dif,
                levels=levels1 ,
                extend='both',
                shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)
m.drawmeridians([90,120,150,180,210,240,270,300,330], linewidth=1.2,labels=[True,False,False,True],fmt='%g',latmax=80,fontsize=15)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,False],latmax=80,fontsize=15)
ax1.annotate(r'$70\degree S$',fontsize=15,xy=m(190,-70),xycoords='data',xytext = m(190,-70),textcoords='data')
ax1.annotate(r'$60\degree S$',fontsize=15,xy=m(190,-60),xycoords='data',xytext = m(190,-60),textcoords='data')
ax1.annotate(r'$80\degree S$',fontsize=15,xy=m(190,-80),xycoords='data',xytext = m(190,-80),textcoords='data')
cs = m.contourf(x,y, 1-pval3, levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
ax1.set_title('Air temperature(ON)', fontsize=24)
m.drawcoastlines(color='blue')

#------------------------------公共部分-------------------------
cax = fig.add_axes([0.24, 0, 0.55, 0.1], aspect=0.05)
cb = fig.colorbar(im,orientation='horizontal',
                  ticks=tick_marks,
                  extend='both', extendfrac='auto',cax = cax)
cb.ax.tick_params(labelsize=16)
cb.set_label("temperature (K)",fontsize=18)
# fig.suptitle('Composite SLP & UV10(ON)', fontsize=24)
# fig.tight_layout(rect=[0.02,0.1,0.98,0.95])
plt.savefig('/stu02/weizx24/figures/Figure8b_CMIP6_t2m.png' ,dpi=600,bbox_inches='tight')
# plt.show()

#endregion
'''


#region sea ice advection
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
plt.close()
plt.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6
#x1, y1 = m(*np.meshgrid(lon,lat))
fig=plt.figure(figsize=(11,7.5))
fig.subplots_adjust(top=0.99)
ax = fig.add_subplot(111)
m = Basemap(projection='stere',resolution='h',
            lat_0=-90, lon_0=180,  lat_ts=(-90.+-55.)/2.,
            llcrnrlon=-130,urcrnrlon=20,llcrnrlat=-50,urcrnrlat=-72)
x1,y1 = m(*np.meshgrid(lon1,lat1))
im1 = m.contourf(x1,y1,adv_dif*100,levels=levels1,extend='both',
                shading='faceted', antialiased=True,cmap=cmaps.BlueWhiteOrangeRed)
cs = m.contourf(x1,y1, 1-pval3, levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
# Q1 = m.quiver(xx,yy,uproj*-1,vproj*-1,scale=11,color='purple',scale_units='inches')
# qk1 = plt.quiverkey(Q1, 0.2, -0.025, 2, '2 km/day', labelpos='W',fontproperties={'size':'large'})
m.drawmeridians([90,120,150,180,210,240,270,300,330], linewidth=1.2,labels=[True,False,False,True],fmt='%g',latmax=80,fontsize=15)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,False],latmax=80,fontsize=15)
ax.annotate(r'$70\degree S$',fontsize=15,xy=m(190,-70),xycoords='data',xytext = m(190,-70),textcoords='data')
ax.annotate(r'$60\degree S$',fontsize=15,xy=m(190,-60),xycoords='data',xytext = m(190,-60),textcoords='data')
ax.annotate(r'$80\degree S$',fontsize=15,xy=m(190,-80),xycoords='data',xytext = m(190,-80),textcoords='data')

#ax.annotate(r'$70\degree S$',fontsize=16,xy=m(10,-70),xycoords='data',xytext = m(10,-70),textcoords='data')
#ax.annotate(r'$60\degree S$',fontsize=16,xy=m(10,-60),xycoords='data',xytext = m(10,-60),textcoords='data')
#ax.annotate(r'$80\degree S$',fontsize=16,xy=m(10,-80),xycoords='data',xytext = m(10,-80),textcoords='data')
ax.set_title('Ice advection (ON)', fontsize=24)

# m.drawmeridians([160,230,300,20,90], linewidth=0.5, dashes=[1, 5],labels=[True,True,True, True])
m.drawcoastlines(color='blue')
m.fillcontinents(color='gray')
#cf = m.contour(x,y, fclm4,levels=np.array([980,986]),colors=['black'],alpha=0.5)
#m.drawparallels([-60],linewidth=1.)
#cs = m.contourf(x,y, p4, levels=[0.9, 0.95] ,colors='none',hatches=['++', None],alpha=0)
# m.drawcoastlines(color='blue')
# cax = fig.add_axes([0.2, 0.02, 0.6, 0.1], aspect=0.05)
# draw_latlon_polygon(m, [180,230], [-70,-70], '-.',color='black',linewidth=2,alpha=1)
cax = fig.add_axes([0.24, 0, 0.55, 0.1], aspect=0.05)
cb = fig.colorbar(im1,orientation='horizontal',ticks=tickmarks,extend='both', cax = cax)
cb.ax.tick_params(labelsize=16)
#cb.dividers.set_color('k')
#cb.dividers.set_linewidth(1)
#cb.dividers.set_linestyle('dotted')
cb.set_label('ice advection ' + '(% km/day)',fontsize=18)
#fig.suptitle('temperature [K]', fontsize=14)
#fig.tight_layout(rect=[0,0.15,1,0.95])
# fig.tight_layout(rect=[0,0.08,1,0.95])
#plt.savefig('fig/Reg_asl_t2m.png' ,dpi=300,bbox_inches='tight')
plt.savefig('/stu02/weizx24/figures/Figure8d_obs_ice_advection.png' ,dpi=300,bbox_inches='tight')


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
plt.close()
plt.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6
#x1, y1 = m(*np.meshgrid(lon,lat))
fig=plt.figure(figsize=(11,7.5))
fig.subplots_adjust(top=0.99)
ax = fig.add_subplot(111)
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
m.drawmeridians([90,120,150,180,210,240,270,300,330], linewidth=1.2,labels=[True,False,False,True],fmt='%g',latmax=80,fontsize=15)
m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,False],latmax=80,fontsize=15)
ax.annotate(r'$70\degree S$',fontsize=15,xy=m(190,-70),xycoords='data',xytext = m(190,-70),textcoords='data')
ax.annotate(r'$60\degree S$',fontsize=15,xy=m(190,-60),xycoords='data',xytext = m(190,-60),textcoords='data')
ax.annotate(r'$80\degree S$',fontsize=15,xy=m(190,-80),xycoords='data',xytext = m(190,-80),textcoords='data')

#ax.annotate(r'$70\degree S$',fontsize=16,xy=m(10,-70),xycoords='data',xytext = m(10,-70),textcoords='data')
#ax.annotate(r'$60\degree S$',fontsize=16,xy=m(10,-60),xycoords='data',xytext = m(10,-60),textcoords='data')
#ax.annotate(r'$80\degree S$',fontsize=16,xy=m(10,-80),xycoords='data',xytext = m(10,-80),textcoords='data')
ax.set_title('Ice advection (ON)', fontsize=24)

# m.drawmeridians([160,230,300,20,90], linewidth=0.5, dashes=[1, 5],labels=[True,True,True, True])
m.drawcoastlines(color='blue')
m.fillcontinents(color='gray')
#cf = m.contour(x,y, fclm4,levels=np.array([980,986]),colors=['black'],alpha=0.5)
#m.drawparallels([-60],linewidth=1.)
#cs = m.contourf(x,y, p4, levels=[0.9, 0.95] ,colors='none',hatches=['++', None],alpha=0)
# m.drawcoastlines(color='blue')
# cax = fig.add_axes([0.2, 0.02, 0.6, 0.1], aspect=0.05)
# draw_latlon_polygon(m, [180,230], [-70,-70], '-.',color='black',linewidth=2,alpha=1)
cax = fig.add_axes([0.24, 0, 0.55, 0.1], aspect=0.05)
cb = fig.colorbar(im1,orientation='horizontal',
                  ticks=tickmarks,
                  extend='both',cax = cax)
cb.ax.tick_params(labelsize=16)
#cb.dividers.set_color('k')
#cb.dividers.set_linewidth(1)
#cb.dividers.set_linestyle('dotted')
cb.set_label('ice advection ' + '(% km/day)',fontsize=18)
#fig.suptitle('temperature [K]', fontsize=14)
#fig.tight_layout(rect=[0,0.15,1,0.95])
# fig.tight_layout(rect=[0,0.08,1,0.95])
#plt.savefig('fig/Reg_asl_t2m.png' ,dpi=300,bbox_inches='tight') 
plt.savefig('/stu02/weizx24/figures/Figure8d_CMIP6_iceadvection.png' ,dpi=300,bbox_inches='tight')
#endregion
