'''sea ice concentration longterm evolution
first start with CESM2
'''
import numpy as np
import matplotlib.pyplot as plt
import glob
#import matplotlib.dates as mdates
from pylab import *
#from matplotlib import dates, ticker
#from scipy import stats
#import cmaps
import xarray as xr
#import matplotlib.path as mpath
#from datetime import datetime, timedelta
#from dateutil.relativedelta import relativedelta
import pandas as pd
#import scipy.signal as signal
# import cartopy.crs as ccrs

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
# from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.gridspec as gridspec
from scipy.stats import t
# from pyproj import Proj, transform
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import AxesGrid

start_time = time.time()

path_area_weight ='/stu02/weizx24/data/areacello_Ofx_CESM2-WACCM_piControl_r1i1p1f1_gn.nc'
areacella= xr.open_dataset(path_area_weight)
aw_xr = areacella['areacello']

thetao = xr.open_dataset('/stu02/weizx24/data/thetao_Omon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all_sellonlatbox_remapdis_zonmean.nc')['thetao'].mean('lon').sel(lev=slice(0,30000)).sel(lat=slice(-80.,-50.))
# .sel(lev=slice(0,30000)).sel(nlat=slice(-80,-50))

# 筛选最大最小年份
CMIP6_opwa = np.load('/stu02/weizx24/data/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']*1e-12
ross_sie_cmip6 = np.nanmean(CMIP6_opwa,axis=1)
Q1 = np.quantile(ross_sie_cmip6,0.05)
Q3 = np.quantile(ross_sie_cmip6,0.95)
sic_Feb = xr.open_dataset('/stu02/weizx24/data/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_Feb.nc')['siconc']
year_min = sic_Feb.time[ross_sie_cmip6<Q1].dt.year
year_max = sic_Feb.time[ross_sie_cmip6>Q3].dt.year

# year_max = thetao.time[ross_sie_cmip6<Q1].dt.year
thetao_max_lst = []
for i in range(len(year_max)):
    y = str(year_max[i].values).zfill(4)
#     print(y)
    thetao_max_lst.append(thetao.sel(time=y))
thetao_max_mean = np.array(thetao_max_lst).squeeze().mean(axis=0)
thetao_max_std = np.array(thetao_max_lst).squeeze().std(axis=0)

# year_min = thetao.time[ross_sie_cmip6<Q1].dt.year
thetao_min_lst = []
for i in range(len(year_min)):
    y = str(year_min[i].values).zfill(4)
#     print(y)
    thetao_min_lst.append(thetao.sel(time=y))
thetao_min_mean = np.array(thetao_min_lst).squeeze().mean(axis=0)
thetao_min_std = np.array(thetao_min_lst).squeeze().std(axis=0)

year_min_1 = sic_Feb.time[ross_sie_cmip6<Q1].dt.year+1
year_max_1 = sic_Feb.time[ross_sie_cmip6>Q3].dt.year+1

# year_min_1 = thetao.time[ross_sie_cmip6<Q1].dt.year
thetao_min_1_lst = []
for i in range(len(year_min_1)):
    y = str(year_min_1[i].values).zfill(4)
#     print(y)
    thetao_min_1_lst.append(thetao.sel(time=y))
thetao_min_1_mean = np.array(thetao_min_1_lst).squeeze().mean(axis=0)
thetao_min_1_std = np.array(thetao_min_1_lst).squeeze().std(axis=0)

# year_max_1 = thetao.time[ross_sie_cmip6<Q1].dt.year
thetao_max_1_lst = []
for i in range(len(year_max_1)):
    y = str(year_max_1[i].values).zfill(4)
#     print(y)
    thetao_max_1_lst.append(thetao.sel(time=y))
thetao_max_1_mean = np.array(thetao_max_1_lst).squeeze().mean(axis=0)
thetao_max_1_std = np.array(thetao_max_1_lst).squeeze().std(axis=0)

thetao_min_mean_all = np.concatenate([thetao_min_mean,thetao_min_1_mean],axis=0)[9:15]
thetao_max_mean_all = np.concatenate([thetao_max_mean,thetao_max_1_mean],axis=0)[9:15]
thetao_min_std_all = np.concatenate([thetao_min_std, thetao_min_1_std],axis=0)[9:15]
thetao_max_std_all = np.concatenate([thetao_max_std, thetao_max_1_std],axis=0)[9:15]

thetao_dif = thetao_max_mean_all - thetao_min_mean_all
# print(thetao_dif.shape)

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

num_y = 25

pre_clm_son = thetao_max_mean_all
pre_std_son = thetao_max_std_all

# thetao_clm_son = thetao_clm
# thetao_std_son = thetao_std

post_clm_son = thetao_min_mean_all
post_std_son = thetao_min_std_all


# pval1=np.zeros((6,332,316),'float')
# pval2=np.zeros((6,332,316),'float')
pval3=np.zeros((6,26,61),'float')

for k in range(6):
    for i in range(26):
        for j in range(61):
#             pval1[k,i,j] = independent_ttest(pre_clm_son[k,i,j], thetao_clm_son[k,i,j], pre_std_son[k,i,j],thetao_std_son[k,i,j], num_y, 30, 0.1)
#             pval2[k,i,j] = independent_ttest(post_clm_son[k,i,j], thetao_clm_son[k,i,j], post_std_son[k,i,j],thetao_std_son[k,i,j], num_y, 30, 0.1)
            pval3[k,i,j] = independent_ttest(pre_clm_son[k,i,j], post_clm_son[k,i,j], pre_std_son[k,i,j],post_std_son[k,i,j], num_y, num_y, 0.1)


#——————混合层深度——————————
mlotst = xr.open_dataset('/stu02/weizx24/data/mlotst_Omon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_all_sellonlatbox_remapdis.nc')['mlotst'].sel(lat=slice(-80,-50)).mean('lon')
# .sel(lev=slice(0,30000)).mean('lon')
lats = mlotst.lat
# lons = mlotst.lon

# year_max = mlotst.time[ross_sie_cmip6<Q1].dt.year
mlotst_max_lst = []
for i in range(len(year_max)):
    y = str(year_max[i].values).zfill(4)
#     print(y)
    mlotst_max_lst.append(mlotst.sel(time=y))
mlotst_max_mean = np.array(mlotst_max_lst).squeeze().mean(axis=0)
mlotst_max_std = np.array(mlotst_max_lst).squeeze().std(axis=0)

# year_min = mlotst.time[ross_sie_cmip6<Q1].dt.year
mlotst_min_lst = []
for i in range(len(year_min)):
    y = str(year_min[i].values).zfill(4)
#     print(y)
    mlotst_min_lst.append(mlotst.sel(time=y))
mlotst_min_mean = np.array(mlotst_min_lst).squeeze().mean(axis=0)
mlotst_min_std = np.array(mlotst_min_lst).squeeze().std(axis=0)

year_min_1 = sic_Feb.time[ross_sie_cmip6<Q1].dt.year+1
year_max_1 = sic_Feb.time[ross_sie_cmip6>Q3].dt.year+1

# year_min_1 = mlotst.time[ross_sie_cmip6<Q1].dt.year
mlotst_min_1_lst = []
for i in range(len(year_min_1)):
    y = str(year_min_1[i].values).zfill(4)
#     print(y)
    mlotst_min_1_lst.append(mlotst.sel(time=y))
mlotst_min_1_mean = np.array(mlotst_min_1_lst).squeeze().mean(axis=0)
mlotst_min_1_std = np.array(mlotst_min_1_lst).squeeze().std(axis=0)

# year_max_1 = mlotst.time[ross_sie_cmip6<Q1].dt.year
mlotst_max_1_lst = []
for i in range(len(year_max_1)):
    y = str(year_max_1[i].values).zfill(4)
#     print(y)
    mlotst_max_1_lst.append(mlotst.sel(time=y))
mlotst_max_1_mean = np.array(mlotst_max_1_lst).squeeze().mean(axis=0)
mlotst_max_1_std = np.array(mlotst_max_1_lst).squeeze().std(axis=0)

mlotst_min_mean_all = np.concatenate([mlotst_min_mean,mlotst_min_1_mean],axis=0)[9:15]
mlotst_max_mean_all = np.concatenate([mlotst_max_mean,mlotst_max_1_mean],axis=0)[9:15]
mlotst_min_std_all = np.concatenate([mlotst_min_std, mlotst_min_1_std],axis=0)[9:15]
mlotst_max_std_all = np.concatenate([mlotst_max_std, mlotst_max_1_std],axis=0)[9:15]

mlotst_dif = mlotst_max_mean_all - mlotst_min_mean_all
lat = thetao.lat
lev = thetao.lev*1e-2
np.savez('/stu02/weizx24/data/npz/Figure6_data.npz',lats=lats,lat=lat,lev=lev,thetao=thetao,mlotst_min_mean_all=mlotst_min_mean_all,mlotst_max_mean_all=mlotst_max_mean_all,thetao_dif=thetao_dif,pval3=pval3)
'''
#——————作图——————————


labelfont=20
tickfont=16
colorbarfont=18
levels1=np.linspace(-1.5,1.5,50)
tick_marks  = np.linspace(-1.5,1.5,5)

mpl.rcParams['hatch.color'] ='gray'
mpl.rcParams['hatch.linewidth'] = 0.8
fig = plt.figure(figsize=(20,8))
ax = fig.add_subplot(161)
figname = 'Oct '
ax.set_title(figname,fontsize=20)
im1 = ax.contourf(lat,lev,thetao_dif[0],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[0], levels=[0.9, 1] ,colors='none',hatches=['.', None],alpha=0)
ct = ax.plot(lats,mlotst_max_mean_all[0],color='k',linewidth=2,linestyle='--',label='max_years')
ct = ax.plot(lats,mlotst_min_mean_all[0],color='r',linewidth=2,linestyle='--',label='min_years')
ax.set_ylabel('Depth (meter)',fontsize=20)
ax.set_xlabel('Latitude',fontsize=20)
ax.set_xticks([-50,-55,-60,-65,-70,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()

#---------------------------------------
ax = fig.add_subplot(162)
figname = 'Nov '
ax.set_title(figname,fontsize=20)
im1 = ax.contourf(lat,lev,thetao_dif[1],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[1], levels=[0.9, 1] ,colors='none',hatches=['.', None],alpha=0)
ct = ax.plot(lats,mlotst_max_mean_all[1],color='k',linewidth=2,linestyle='--',label='max_years')
ct = ax.plot(lats,mlotst_min_mean_all[1],color='r',linewidth=2,linestyle='--',label='min_years')
# ax.set_ylabel('Depth (meter)',fontsize=20)
ax.set_xlabel('Latitude',fontsize=20)
ax.set_xticks([-50,-55,-60,-65,-70,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()
ax.get_yaxis().set_visible(False)
#---------------------------------------
ax = fig.add_subplot(163)
figname = 'Dec '
ax.set_title(figname,fontsize=20)
im1 = ax.contourf(lat,lev,thetao_dif[2],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[2], levels=[0.9, 1] ,colors='none',hatches=['.', None],alpha=0)
ct = ax.plot(lats,mlotst_max_mean_all[2],color='k',linewidth=2,linestyle='--',label='max_years')
ct = ax.plot(lats,mlotst_min_mean_all[2],color='r',linewidth=2,linestyle='--',label='min_years')
# ax.set_ylabel('Depth (meter)',fontsize=20)
ax.set_xlabel('Latitude',fontsize=20)
ax.set_xticks([-50,-55,-60,-65,-70,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()
ax.get_yaxis().set_visible(False)
#---------------------------------------
ax = fig.add_subplot(164)
figname = 'Jan '
ax.set_title(figname,fontsize=20)
im1 = ax.contourf(lat,lev,thetao_dif[3],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[3], levels=[0.9, 1] ,colors='none',hatches=['.', None],alpha=0)
ct = ax.plot(lats,mlotst_max_mean_all[3],color='k',linewidth=2,linestyle='--',label='max_years')
ct = ax.plot(lats,mlotst_min_mean_all[3],color='r',linewidth=2,linestyle='--',label='min_years')
ax.set_xlabel('Latitude',fontsize=20)
ax.set_xticks([-50,-55,-60,-65,-70,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()
ax.get_yaxis().set_visible(False)
#--------------------------------
ax = fig.add_subplot(165)
figname = 'Feb '
ax.set_title(figname,fontsize=20)
im1 = ax.contourf(lat,lev,thetao_dif[4],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[4], levels=[0.9, 1] ,colors='none',hatches=['.', None],alpha=0)
ct = ax.plot(lats,mlotst_max_mean_all[4],color='k',linewidth=2,linestyle='--',label='max_years')
ct = ax.plot(lats,mlotst_min_mean_all[4],color='r',linewidth=2,linestyle='--',label='min_years')
ax.set_xlabel('Latitude',fontsize=20)
ax.set_xticks([-50,-55,-60,-65,-70,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()
ax.get_yaxis().set_visible(False)
#----------------------------------
ax = fig.add_subplot(166)
figname = 'Mar '
ax.set_title(figname,fontsize=20)
im1 = ax.contourf(lat,lev,thetao_dif[5],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[5], levels=[0.9, 1] ,colors='none',hatches=['.', None],alpha=0)
ct = ax.plot(lats,mlotst_max_mean_all[5],color='k',linewidth=2,linestyle='--',label='largest_years')
ct = ax.plot(lats,mlotst_min_mean_all[5],color='r',linewidth=2,linestyle='--',label='smallest_years')
ax.set_xlabel('Latitude',fontsize=20)
ax.set_xticks([-50,-55,-60,-65,-70,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()
ax.get_yaxis().set_visible(False)

#————添加colorbar————
plt.legend(fontsize=16)
cax = fig.add_axes([1, 0.1, 0.015, 0.9], aspect=18)
cb = fig.colorbar(im1,orientation='vertical',
                  ticks =tick_marks,
                  extend='both', extendfrac=0.05,
                  cax = cax
                 )
cb.ax.tick_params(labelsize=colorbarfont, width=0)
cb.set_label('Ocean temperature (℃)',fontsize=20)
# cb.dividers.set_color('k')
#cb.outline.set_edgecolor('white')
# cb.dividers.set_linewidth(1)
fig.tight_layout(rect=[0,0.1,1,0.95])
plt.savefig('/stu02/weizx24/figures/CMIP6_ocean.png' ,dpi=600,bbox_inches='tight')
plt.show()

end_time = time.time()
elapsed_time = end_time - start_time
print(f"程序运行时间: {elapsed_time:.2f} 秒")
'''