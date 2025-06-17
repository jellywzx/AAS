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

#设置画图的全局变量
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 9
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 4


#region 观测
date = ['Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr']
file = np.load('/stu02/weizx24/data/npz/Figure7_obs_0924.npz')
ross_lw_dif = file['ross_lw_dif']
ross_sw_dif = file['ross_sw_dif']
ross_lhf_dif = file['ross_lhf_dif']
ross_shf_dif = file['ross_shf_dif']
ross_Qnet_dif = file['ross_Qnet_dif']
ross_votemp_dif = file['ross_votemp_dif']

#单独出图
plt.close()
fig = plt.figure(1, figsize=(150/25.4,4))
ax = fig.add_subplot(111)
ax.plot(date,ross_lw_dif,color='C2',label='LW',alpha=0.5)
ax.plot(date,ross_sw_dif,color='y',label='SW',alpha=0.5)
ax.plot(date,ross_lhf_dif,color='purple',label='LHF',alpha=0.5)
ax.plot(date,ross_shf_dif,color='k',label='SHF',alpha=0.5)
ax.plot(date,ross_Qnet_dif,color='b',marker='^',label='Qnet')
ax1 = ax.twinx()
# ax1.plot(date,ross_sst_dif,color='m',marker='8',markersize=13,linewidth=4,label='SST_dif')
ax1.plot(date,ross_votemp_dif,color='red',marker='8',label='UOT')
ax1.set_ylabel('Ocean temperature(°C)',color='red')
ax1.set_yticks(np.linspace(-0.8,0.8,13)[::3])
# ax1.yaxis.set_tick_params(labelsize=20)
ax1.tick_params('both', length=8, width=1.2, which='major',colors='red')
ax1.tick_params('both', length=4, width=0.5, which='minor',colors='red')
ax.set_yticks(np.linspace(-20,20,13)[::3])
ax.set_xlabel('Month')
ax.set_ylabel('Heat flux ' + r'(W m$^{\text{-2}}$)')
ax.set_title('(a) Observation')
plt.tight_layout(rect=[0,0,0.9,1])
ax.grid(color='lightgray',linestyle='--')
#给两个y轴的线条添加Legend
lines = ax.get_lines()+ax1.get_lines()
plt.legend(lines, [l.get_label() for l in lines],edgecolor='k',loc='lower left')
plt.savefig('/stu02/weizx24/figures/0924/Figure7_obs.png' ,dpi=300,bbox_inches='tight')

#region 模式
# date = ['Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr']
file_mdl = np.load('/stu02/weizx24/data/npz/Figure7_model_0924.npz')
ross_lw_dif = file_mdl['ross_lw_dif']
ross_sw_dif = file_mdl['ross_sw_dif']
ross_hfls_dif = file_mdl['ross_hfls_dif']
ross_hfss_dif = file_mdl['ross_hfss_dif']
ross_Qnet_dif = file_mdl['ross_Qnet_dif']
ross_thetao_dif = file_mdl['ross_thetao_dif']

#单独出图
plt.close()
fig = plt.figure(1, figsize=(150/25.4,4))
ax3 = fig.add_subplot(111)
ax3.plot(date,ross_lw_dif,color='C2',label='LW',alpha=0.5)
ax3.plot(date,ross_sw_dif,color='y',label='SW',alpha=0.5)
ax3.plot(date,ross_hfls_dif,color='purple', label='LHF',alpha=0.5)
ax3.plot(date,ross_hfss_dif,color='k',label='SHF',alpha=0.5)
ax3.plot(date,ross_Qnet_dif,color='b',marker='^',label='Qnet')
# ax3.plot(date,ross_Qnet1_dif,color='purple',marker='^',markersize=13, label='Qnet1')
ax4 = ax3.twinx()
# ax4.plot(date,ross_sst_anom,color='m',marker='8',markersize=13, label='SST_anom')
ax4.plot(date,ross_thetao_dif,color='red',marker='8',label='UOT')
ax4.set_yticks(np.linspace(-0.8,0.8,13)[::3])
ax4.set_ylabel('Ocean temperature (°C)',color='red')
# ax4.yax3is.set_tick_params(labelsize=20)
ax4.tick_params('both', length=8, width=1.2, which='major',colors='xkcd:red')
ax4.tick_params('both', length=4, width=0.5, which='minor',colors='xkcd:red')
ax3.set_yticks(np.linspace(-20,20,13)[::3])
ax3.set_xlabel('Month')
ax3.set_ylabel('Heat flux ' +r'(W m$^{\text{-2}}$)')

ax3.set_title('(b) CESM2-WACCM-FV2')
ax3.grid(color='lightgray',linestyle='--')
plt.tight_layout(rect=[0,0,0.9,1])
lines = ax3.get_lines()+ax4.get_lines()
plt.legend(lines, [l.get_label() for l in lines],edgecolor='k',loc='lower right')
plt.savefig('/stu02/weizx24/figures/0924/Figure7_model.png' ,dpi=300,bbox_inches='tight')
#endregion
