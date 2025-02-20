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
import os

folder_path = 'D:/opw_rec/nsidc/new/'
files = os.listdir(folder_path)
files.sort()
files
a = list()
for i in files:
#     prini(i)
    a.append(pd.read_csv(folder_path+i)['eastross']+pd.read_csv(folder_path+i)['westross'])
a = np.array(a)

#春季平均
spring_mean_polynya = np.nanmean(a[:,0:30],axis=1)
a_7daysmean = a[:,27:34].mean(axis=1)
correlation_matrix = np.corrcoef(spring_mean_polynya, a_7daysmean)
correlation_coefficient = correlation_matrix[0, 1]
print("相关系数为:", correlation_coefficient)

#设置全局变量
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 9
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['axes.titlesize'] = 10  # 轴标题字体大小
plt.rcParams['lines.markersize'] = 1 
mpl.rcParams['hatch.color'] ='springgreen'
mpl.rcParams['hatch.linewidth'] = 0.6
#来个默认设置的出图
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111)
ax.scatter(spring_mean_polynya,a_7daysmean,marker='o',color='k',s=20)
ax.set_ylabel('Spring average area'+ r' (10$^{6}$ km$^{2}$)',fontsize=20)
ax.set_xlabel('Average area from Nov 28 to Dec 4'+ r' (10$^{6}$ km$^{2}$)',fontsize=20)
ax.grid(linestyle='--',alpha=0.5)
ax.yaxis.set_tick_params(labelsize=20)
ax.xaxis.set_tick_params(labelsize=20)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.1, 0.9,'correlation coefficient: '+str(round(correlation_coefficient,2))+'\n'+'p value < 0.05',
        transform=ax.transAxes, fontsize=20,
        verticalalignment='top', bbox=props,)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
plt.tight_layout(rect=[0,0,0.9,1])
plt.savefig('C:/Users/fzjxw/python/code/Figures/majorrevision_areacorre.png' ,dpi=300,bbox_inches='tight')
plt.show()
