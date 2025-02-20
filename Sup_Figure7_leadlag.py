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
from scipy.stats import pearsonr

#设置全局变量
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 15
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 16  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6

#观测
file_obs = np.load('D:/npz/Figure7_obs_allmon.npz')
ross_Qnet_obs = file_obs['ross_Qnet_all']
ross_votemp_obs = file_obs['ross_votemp_all']
corre_obs = np.corrcoef(ross_Qnet_obs[:-1], ross_votemp_obs[1:])[0, 1]
corre_obs,p_obs = pearsonr(ross_Qnet_obs[:-1], ross_votemp_obs[1:])
print(f'观测数据的平均滞后两个月的相关系数是{corre_obs}')
print(f'观测数据的p值是{p_obs}')
# plt.close()
fig = plt.figure(figsize=(8, 6))
ax1 = fig.add_subplot(111)
ax1.scatter(ross_Qnet_obs[:-1],ross_votemp_obs[1:],color='blue',label='Observation')
# plt.title('Observation')
# plt.xlabel('Qnet (leading)')
# plt.ylabel('UOT (lagging)')
# plt.grid(True)
# ax1.text(0.1, 0.9,'pearson r: '+str(round(corre_obs,2))+'\n'+'p value < 0.05',
#         transform=ax1.transAxes, verticalalignment='top', )
# ax1.text(-2.8, 2.8, f'Lag Coefficient: {corre_obs}', fontsize=12, color='red', 
#          bbox=dict(facecolor='white', alpha=0.5))
# plt.savefig('/stu02/weizx24/figures/0924/Qnet_votemp_leadlag.png')

#模式
file_mdl = np.load('D:/npz/Figure7_mdl_allmon.npz')
ross_Qnet1_mdl = file_mdl['ross_Qnet1_all']
ross_thetao_mdl = file_mdl['ross_thetao_all']
# corre_mdl = np.corrcoef(ross_Qnet1_mdl[:-1], ross_thetao_mdl[1:])[0, 1]
corre_mdl,p_mdl = pearsonr(ross_Qnet1_mdl[:-1], ross_thetao_mdl[1:])
print(f'模式数据的平均滞后相关系数是{corre_mdl}')
print(f'模式数据的p值是{p_mdl}')
# ax2 = fig.add_subplot(122)
ax1.scatter(ross_Qnet1_mdl[:-1],ross_thetao_mdl[1:],color='black',label='CESM2-WACCM-FV2')
# plt.title('CESM2-WACCM-FV2')
plt.xlabel('Qnet (leading)')
plt.ylabel('UOT (lagging)')
plt.grid(True,linestyle='--',color='lightgray')
ax1.text(0.04, 0.8,'Observation pearson r: '+str(round(corre_obs,2))+'(p < 0.05)'+'\n'+'Model pearson r: '+str(round(corre_mdl,2))+'(p < 0.05)',
        transform=ax1.transAxes, 
        verticalalignment='top',)
plt.legend(markerscale=10)
# props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# ax2.text(-2.8, 2.8, f'Observation lead-Lag Coefficient: {corre_mdl} ', fontsize=12, color='red', 
#          bbox=dict(facecolor='white', alpha=0.5))
plt.subplots_adjust(wspace=0.3)
plt.savefig('C:/Users/fzjxw/python/code/Figures/Qnet_votemp_leadlag.png',dpi=300)