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
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.gridspec as gridspec
from scipy.stats import t
from pyproj import Proj, transform
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import AxesGrid


plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 9
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 

#观测
#仅罗斯海冰间湖面积
import os
folder_path = '/stu02/weizx24/data/opw_rec/nsidc/new/'
files = os.listdir(folder_path)
files.sort()
# files = np.delete(files,[-1])

a = list()
for i in files:
#     prini(i)
    a.append(pd.read_csv(folder_path+i)['eastross']+pd.read_csv(folder_path+i)['westross'])
a = np.array(a)
a_7days_mean = a[:,27:34].mean(axis=1)

years = np.linspace(1992,2021,30)
year_min = [1994., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2016., 2018., 2021.]
bool_max = np.isclose(years[:, None], year_max).any(axis=1)
bool_min = np.isclose(years[:, None], year_min).any(axis=1)

a_max = a_7days_mean[bool_max].mean(axis=0)
a_min = a_7days_mean[bool_min].mean(axis=0)
a_clm = a_7days_mean.mean(axis=0)

#因为模式的数据量太大了，所以只计算了这七天平均的冰间湖面积情况。
#如果要说明模式与观测的冰间湖面积大小区别，选取额12-01的进行对比。加上12-15的看一下，12-26的看一下。

#模式冰间湖数据
CMIP6_opwa = np.load('/stu02/weizx24/data/npz/CMIP6_opwa_1128-1204_20240528.npz')['CMIP6_opwa']
CMIP6_opwa_7day_mean = np.nanmean(np.array(CMIP6_opwa),axis=1)
#挑选最大最小年份
Q1 = np.quantile(CMIP6_opwa_7day_mean,0.05)
Q3 = np.quantile(CMIP6_opwa_7day_mean,0.95)
sic_Feb = xr.open_dataset('/stu02/weizx24/data/siconc_SImon_CESM2-WACCM-FV2_piControl_r1i1p1f1_gn_000101-049912_Feb.nc')['siconc']
year_min = sic_Feb.time[CMIP6_opwa_7day_mean<Q1].dt.year
year_max = sic_Feb.time[CMIP6_opwa_7day_mean>Q3].dt.year
CMIP6_opwa_max = CMIP6_opwa_7day_mean[CMIP6_opwa_7day_mean>Q3].mean(axis=0)*1e-12
CMIP6_opwa_min = CMIP6_opwa_7day_mean[CMIP6_opwa_7day_mean<Q1].mean(axis=0)*1e-12
CMIP6_opwa_clm = CMIP6_opwa_7day_mean.mean(axis=0)*1e-12

#timeseries变化图

# date_name = pd.date_range('1979-11-01','1979-12-31')
# date_list = list()
# for i in range(len(date_name)):
#     date_list.append(str(date_name[i])[5:10])

# fig = plt.figure(1, figsize=(6,6))
# ax4 = fig.add_subplot(111)
# ax4.plot(date_list[:-5],a_max[:-5],color='r',label='Largest')
# ax4.plot(date_list[:-5],a_clm[:-5],color='k',label='Climatology ')
# ax4.plot(date_list[:-5],a_min[:-5],color='b',label='Smallest')
# ax4.set_xlabel('Date')
# ax4.set_ylabel('Areas'+ r' (10$^{6}$ km$^{2}$)',)
# ax4.set_xticks(date_list[:-5][::5])
# ax4.set_yticks([0,0.2,0.4,0.6])
# ax4.grid(color='lightgray',linestyle='--',alpha=0.4)
# ax4.legend(edgecolor='k',loc='upper left')
# ax4.text(0, 1.05, '(d)', fontsize=10, transform=ax4.transAxes, va='top', ha='right')
# plt.subplots_adjust(hspace=0.35)
# plt.savefig('/stu02/weizx24/figures/0924/Figure3_all.png',dpi=300,bbox_inches='tight')

#bar plot
# obs_file = np.load('/stu02/weizx24/data/npz/Figure3_polynyaarea_obs.npz')
group1 = [a_clm,a_max,a_min]
# mdl_file = np.load('/stu02/weizx24/data/npz/Figure3_polynyaarea_mdl.npz')
group2 = [CMIP6_opwa_clm,CMIP6_opwa_max,CMIP6_opwa_min]

labels = ['Climatology', 'Large', 'Small',]
x = np.arange(len(labels))  # x轴的位置
width = 0.15  # 柱子的宽度
# 创建图表
fig, ax = plt.subplots()
# 画两组数据的柱状图
rects1 = ax.bar(x - width/2, group1, width, label='Observation')
rects2 = ax.bar(x + width/2, group2, width, label='CESM2-WACCM-FV2')
# ax.set_xlabel('Categories')
ax.set_ylabel('Areas'+ r' (10$^{6}$ km$^{2}$)')
# ax.set_title('Comparison ')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()
# plt.tight_layout()
plt.savefig('/stu02/weizx24/figures/0924/SupFig_polynya_barplot.png',dpi=300)
# plt.show()

