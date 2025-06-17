'''sea ice concentration longterm evolution
first start with CESM2
'''
import numpy as np
import matplotlib.pyplot as plt
import glob
from pylab import *
import xarray as xr
import pandas as pd
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

file = np.load('/stu02/weizx24/data/npz/Figure6_data.npz')
thetao = file['thetao']
thetao_dif =file['thetao_dif']
lats = file['lats']
mlotst_max_mean_all = file['mlotst_max_mean_all']
mlotst_min_mean_all = file['mlotst_min_mean_all']
pval3 = file['pval3']
lat = file['lat']
lev = file['lev']


labelfont=40
tickfont=40
colorbarfont=40
levels1=np.linspace(-1.5,1.5,50)
tick_marks  = np.linspace(-1.5,1.5,5)

mpl.rcParams['hatch.color'] ='gray'
mpl.rcParams['hatch.linewidth'] = 0.8
plt.rcParams['font.family'] = 'Arial'
fig = plt.figure(figsize=(20,8))
ax = fig.add_subplot(161)
# figname = 'Oct '
# ax.set_title(figname,fontsize=24)
im1 = ax.contourf(lat,lev,thetao_dif[0],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[0], levels=[0.9, 1] ,colors='none',hatches=['.', None],alpha=0)
ct = ax.plot(lats,mlotst_max_mean_all[0],color='k',linewidth=2,linestyle='--',label='max_years')
ct = ax.plot(lats,mlotst_min_mean_all[0],color='r',linewidth=2,linestyle='--',label='min_years')
ax.set_ylabel('Depth (meter)',fontsize=labelfont)
# ax.set_xlabel('Latitude',fontsize=24)
# ax.set_xticks([-50,-55,-60,-65,-70,-75])
ax.set_xticks([-55,-65,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()

#---------------------------------------
ax = fig.add_subplot(162)
figname = 'Nov '
# ax.set_title(figname,fontsize=24)
im1 = ax.contourf(lat,lev,thetao_dif[1],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[1], levels=[0.9, 1] ,colors='none',hatches=['.', None],alpha=0)
ct = ax.plot(lats,mlotst_max_mean_all[1],color='k',linewidth=2,linestyle='--',label='max_years')
ct = ax.plot(lats,mlotst_min_mean_all[1],color='r',linewidth=2,linestyle='--',label='min_years')
# ax.set_ylabel('Depth (meter)',fontsize=24)
# ax.set_xlabel('Latitude',fontsize=24)
# ax.set_xticks([-50,-55,-60,-65,-70,-75])
ax.set_xticks([-55,-65,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()
ax.get_yaxis().set_visible(False)
#---------------------------------------
ax = fig.add_subplot(163)
figname = 'Dec '
# ax.set_title(figname,fontsize=24)
im1 = ax.contourf(lat,lev,thetao_dif[2],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[2], levels=[0.9, 1] ,colors='none',hatches=['.', None],alpha=0)
ct = ax.plot(lats,mlotst_max_mean_all[2],color='k',linewidth=2,linestyle='--',label='max_years')
ct = ax.plot(lats,mlotst_min_mean_all[2],color='r',linewidth=2,linestyle='--',label='min_years')
# ax.set_ylabel('Depth (meter)',fontsize=24)
ax.set_xlabel('Latitude',fontsize=labelfont)
# ax.set_xticks([-50,-55,-60,-65,-70,-75])
ax.set_xticks([-55,-65,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()
ax.get_yaxis().set_visible(False)
#---------------------------------------
ax = fig.add_subplot(164)
figname = 'Jan '
# ax.set_title(figname,fontsize=24)
im1 = ax.contourf(lat,lev,thetao_dif[3],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[3], levels=[0.9, 1] ,colors='none',hatches=['.', None],alpha=0)
ct = ax.plot(lats,mlotst_max_mean_all[3],color='k',linewidth=2,linestyle='--',label='max_years')
ct = ax.plot(lats,mlotst_min_mean_all[3],color='r',linewidth=2,linestyle='--',label='min_years')
# ax.set_xlabel('Latitude',fontsize=24)
# ax.set_xticks([-50,-55,-60,-65,-70,-75])
ax.set_xticks([-55,-65,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()
ax.get_yaxis().set_visible(False)
#--------------------------------
ax = fig.add_subplot(165)
figname = 'Feb '
# ax.set_title(figname,fontsize=24)
im1 = ax.contourf(lat,lev,thetao_dif[4],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[4], levels=[0.9, 1] ,colors='none',hatches=['.', None],alpha=0)
ct = ax.plot(lats,mlotst_max_mean_all[4],color='k',linewidth=2,linestyle='--',label='max_years')
ct = ax.plot(lats,mlotst_min_mean_all[4],color='r',linewidth=2,linestyle='--',label='min_years')
# ax.set_xlabel('Latitude',fontsize=24)
# ax.set_xticks([-50,-55,-60,-65,-70,-75])
ax.set_xticks([-55,-65,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()
ax.get_yaxis().set_visible(False)
#----------------------------------
ax = fig.add_subplot(166)
figname = 'Mar '
# ax.set_title(figname,fontsize=24)
im1 = ax.contourf(lat,lev,thetao_dif[5],
                 levels=levels1,
                 extend='both',shading='faceted', antialiased=True,cmap='RdBu_r')
cs = ax.contourf(lat,lev, 1-pval3[5], levels=[0.9, 1] ,colors='none',hatches=['.', None],alpha=0)
ct = ax.plot(lats,mlotst_max_mean_all[5],color='k',linewidth=2,linestyle='--',label='Large')
ct = ax.plot(lats,mlotst_min_mean_all[5],color='r',linewidth=2,linestyle='--',label='Small')
# ax.set_xlabel('Latitude',fontsize=24)
ax.set_xticks([-55,-65,-75])
ax.set_yticks(np.linspace(0,250,11)[::2])
ax.yaxis.set_tick_params(labelsize=tickfont)
ax.xaxis.set_tick_params(labelsize=tickfont)
plt.gca().invert_yaxis()
ax.get_yaxis().set_visible(False)

#————添加colorbar————
plt.legend(fontsize=27,loc='lower right')
cax = fig.add_axes([1, 0.1, 0.015, 0.9], aspect=18)
cb = fig.colorbar(im1,orientation='vertical',
                  ticks =tick_marks,
                  extend='both', extendfrac=0.05,
                  cax = cax
                 )
cb.ax.tick_params(labelsize=colorbarfont, width=0)
cb.set_label('Ocean temperature (°C)',fontsize=labelfont)
# cb.dividers.set_color('k')
#cb.outline.set_edgecolor('white')
# cb.dividers.set_linewidth(1)
fig.tight_layout(rect=[0,0.1,1,0.95])
plt.savefig('/stu02/weizx24/figures/0924/Fig6/Fig6_CMIP6_ocean.png' ,dpi=300,bbox_inches='tight')
plt.show()
