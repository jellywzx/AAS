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
import cmocean as cm
from scipy.stats import t

#设置全局变量
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 15  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6

#观测
file_obs = np.load('D:/npz/Figure7_obs_allmon.npz')
ross_Qnet_obs = file_obs['ross_Qnet_all']
ross_votemp_obs = file_obs['ross_votemp_all']
corre_obs = np.corrcoef(ross_Qnet_obs[:-1], ross_votemp_obs[1:])[0, 1]
print(f'观测数据的平均滞后相关系数是{corre_obs}')
plt.close()
fig = plt.figure(figsize=(12, 6))
ax1 = fig.add_subplot(121)
ax1.scatter(ross_Qnet_obs[:-1],ross_votemp_obs[1:])
plt.title('Observation')
plt.xlabel('Qnet (leading variable)')
plt.ylabel('UOT (lagging variable)')
plt.grid(True)
ax1.text(-2.8, 2.8, f'Lag Coefficient: {corre_obs}', fontsize=12, color='red', 
         bbox=dict(facecolor='white', alpha=0.5))
# plt.savefig('/stu02/weizx24/figures/0924/Qnet_votemp_leadlag.png')

#模式
file_mdl = np.load('D:/npz/Figure7_mdl_allmon.npz')
ross_Qnet1_mdl = file_mdl['ross_Qnet1_all']
ross_thetao_mdl = file_mdl['ross_thetao_all']
corre_mdl = np.corrcoef(ross_Qnet1_mdl[:-1], ross_thetao_mdl[1:])[0, 1]
print(f'模式数据的平均滞后相关系数是{corre_mdl}')

ax2 = fig.add_subplot(122)
ax2.scatter(ross_Qnet1_mdl[:-1],ross_thetao_mdl[1:])
plt.title('CESM2-WACCM-FV2')
plt.xlabel('Qnet (leading variable)')
plt.ylabel('UOT (lagging variable)')
plt.grid(True)
# ax.text(0.6, 0.2,'pearson r: '+str(round(p[0],2))+'\n'+'p value < 0.05',
#         transform=ax.transAxes, fontsize=20,
#         verticalalignment='top', bbox=props,)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax2.text(-2.8, 2.8, f'Observation lead-Lag Coefficient: {corre_mdl} ', fontsize=12, color='red', 
         bbox=dict(facecolor='white', alpha=0.5))
plt.subplots_adjust(wspace=0.3)
plt.savefig('C:/Users/fzjxw/python/Figures/majorrevision_Figure7.png',dpi=300)