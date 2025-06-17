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
mpl.rcParams['hatch.color'] ='gray'
mpl.rcParams['hatch.linewidth'] = 0.8

#设置全局变量
plt.rcParams['font.family'] = 'Arial'
# plt.rcParams['font.size'] = 9
# plt.rcParams['lines.linewidth'] = 1
# plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
# plt.rcParams['lines.markersize'] = 1 


file_obs = np.load('/stu02/weizx24/data/npz/Figure6_obs.npz')
pval3 = file_obs['pval3']
lat = file_obs['lat']
lev = file_obs['lev']
votemp_dif = file_obs['votemp_dif']
somx_max_mean = file_obs['somx_max_mean']
somx_min_mean = file_obs['somx_min_mean']

#region 出图方式二
labelfont=40
tickfont=40
colorbarfont=40
levels1=np.linspace(-1.5,1.5,50)
tick_marks  = np.linspace(-1.5,1.5,5)

mpl.rcParams['hatch.color'] ='gray'
mpl.rcParams['hatch.linewidth'] = 0.8

fig = plt.figure(figsize=(20,7))


#--------------------------------------
ax = fig.add_subplot(161)
figname = 'Oct '
ax.set_title(figname,fontsize=labelfont)

im1 = ax.contourf(lat,lev,votemp_dif[0],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[0], levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
ax.plot(lat,somx_max_mean[0],linewidth=2,color='k',linestyle='--',label='Large')
ax.plot(lat,somx_min_mean[0],linewidth=2,color='r',linestyle='--',label='Small')
ax.set_ylabel('Depth (meter)',fontsize=labelfont)
ax.set_xlabel('Latitude',fontsize=labelfont)
ax.set_xticks([-50,-55,-60,-65,-70,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()
ax.xaxis.set_visible(False)
# ax.yaxis.set_visible(False)
#--------------------------------------
ax = fig.add_subplot(162)
figname = 'Nov '
ax.set_title(figname,fontsize=labelfont)

im1 = ax.contourf(lat,lev,votemp_dif[1],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[1], levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
ax.plot(lat,somx_max_mean[1],linewidth=2,color='k',linestyle='--',label='Large')
ax.plot(lat,somx_min_mean[1],linewidth=2,color='r',linestyle='--',label='Small')
ax.set_ylabel('Depth (meter)',fontsize=labelfont)
ax.set_xlabel('Latitude',fontsize=labelfont)
ax.set_xticks([-50,-55,-60,-65,-70,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()
ax.yaxis.set_visible(False)
ax.xaxis.set_visible(False)

#--------------------------------------
ax = fig.add_subplot(163)
figname = 'Dec '
ax.set_title(figname,fontsize=labelfont)

im1 = ax.contourf(lat,lev,votemp_dif[2],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[2], levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
ax.plot(lat,somx_max_mean[2],linewidth=2,color='k',linestyle='--',label='Large')
ax.plot(lat,somx_min_mean[2],linewidth=2,color='r',linestyle='--',label='Small')
ax.set_ylabel('Depth (meter)',fontsize=labelfont)
ax.set_xlabel('Latitude',fontsize=labelfont)
ax.set_xticks([-50,-55,-60,-65,-70,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()
ax.yaxis.set_visible(False)
ax.xaxis.set_visible(False)
#---------------------------------------
ax = fig.add_subplot(164)
figname = 'Jan '
ax.set_title(figname,fontsize=labelfont)
im1 = ax.contourf(lat,lev,votemp_dif[3],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[3], levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
ax.plot(lat,somx_max_mean[3],linewidth=2,color='k',linestyle='--',label='Large')
ax.plot(lat,somx_min_mean[3],linewidth=2,color='r',linestyle='--',label='Small')
ax.set_xlabel('Latitude',fontsize=labelfont)
ax.set_xticks([-50,-55,-60,-65,-70,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()
ax.yaxis.set_visible(False)
ax.xaxis.set_visible(False)
#--------------------------------
ax = fig.add_subplot(165)
figname = 'Feb '
ax.set_title(figname,fontsize=labelfont)
im1 = ax.contourf(lat,lev,votemp_dif[4],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[4], levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
ax.plot(lat,somx_max_mean[4],linewidth=2,color='k',linestyle='--',label='Large')
ax.plot(lat,somx_min_mean[4],linewidth=2,color='r',linestyle='--',label='Small')
ax.set_xlabel('Latitude',fontsize=labelfont)
ax.set_xticks([-50,-55,-60,-65,-70,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()
ax.yaxis.set_visible(False)
ax.xaxis.set_visible(False)

#----------------------------------
ax = fig.add_subplot(166)
figname = 'Mar '
ax.set_title(figname,fontsize=labelfont)
im1 = ax.contourf(lat,lev,votemp_dif[5],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[5], levels=[0.95, 1] ,colors='none',hatches=['.', None],alpha=0)
ax.plot(lat,somx_max_mean[5],linewidth=2,color='k',linestyle='--',label='Large')
ax.plot(lat,somx_min_mean[5],linewidth=2,color='r',linestyle='--',label='Small')
ax.set_xlabel('Latitude',fontsize=labelfont)
ax.set_xticks([-50,-55,-60,-65,-70,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()
ax.yaxis.set_visible(False)
ax.xaxis.set_visible(False)
#————添加colorbar————
plt.legend(fontsize=27,loc='lower right')
cax = fig.add_axes([1, 0.1, 0.015, 0.9], aspect=18)
cb = fig.colorbar(im1,orientation='vertical',
                  ticks =tick_marks,
                  extend='both', extendfrac=0.05,
#                   drawedges=True,
                  cax = cax
                 )
cb.ax.tick_params(labelsize=colorbarfont, width=0)
cb.set_label('Ocean temperature (°C)',fontsize=labelfont)
# cb.dividers.set_color('k')
#cb.outline.set_edgecolor('white')
# cb.dividers.set_linewidth(1)
fig.tight_layout(rect=[0,0.1,1,0.95])
plt.savefig('/stu02/weizx24/figures/0924/Fig6/Figure6_obs.png' ,dpi=300,bbox_inches='tight')
plt.show()


#endregion
