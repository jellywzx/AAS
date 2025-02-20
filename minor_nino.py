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
import numpy as np

#region 读取Nino数据
file_path = 'D:/ENSO_index/nino34_1979_2021.txt'
with open(file_path, 'r') as file:
    lines = file.readlines()
years = []
data = []
for line in lines:
    values = line.split()
    years.append(int(values[0]))  # 第一列是年份
    data.append([float(v) for v in values[1:]])  # 其余列是数据
# years数组是从1979到2022,总共44年)
years = np.array(years)
data = np.array(data)
# 输出验证
print("Years array:", years)
print("Data array shape:", data.shape)
print("Data array (first row):", data[0])
#endregion

#冰间湖面积较大和较小的年份
year_min = [1994., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2016., 2018., 2021.]
nino_max_lst = []
nino_min_lst = []
for i in range(44):
    if years[i] in year_max:
        nino_max_lst.append(data[i])
    if years[i] in year_min:
        nino_min_lst.append(data[i])
nino_max_arr = np.array(nino_max_lst)    
nino_min_arr = np.array(nino_min_lst)

nino_max_mean = np.nanmean(nino_max_arr,axis=0)

#作图
#October-November平均nino指数
nino_max_octnov = np.nanmean(nino_max_arr[:,9:11],axis=1)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
categories = ['2004', '2005', '2016', '2018','2021']
ax.bar(categories,nino_max_octnov)
ax.axhline(y=0,color='red',linestyle='--',)
ax.axhline(y=-0.5,color='red',linestyle='--',)
plt.show()

#November的nino指数
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
year_max_str = ['2004', '2005', '2016', '2018','2021']
ax.bar(year_max_str,nino_max_arr[:,10])
ax.axhline(y=0,color='red',linestyle='--',)
ax.axhline(y=-0.5,color='red',linestyle='--',)
plt.show()

# October的nino指数
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
year_max_str = ['2004', '2005', '2016', '2018','2021']
ax.bar(year_max_str,nino_max_arr[:,9])
ax.axhline(y=0,color='red',linestyle='--',)
ax.axhline(y=-0.5,color='red',linestyle='--',)
plt.show()