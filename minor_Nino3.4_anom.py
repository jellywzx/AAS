import numpy as np
import matplotlib.pyplot as plt
import glob
import matplotlib.dates as mdates
from pylab import *
from matplotlib import dates, ticker
from scipy import stats
import xarray as xr
import matplotlib.path as mpath
from datetime import datetime, timedelta, date
from dateutil.relativedelta import relativedelta
import pandas as pd
import scipy.signal as signal
import cmaps
from mpl_toolkits.basemap import Basemap
from scipy.stats import t
import os 
from scipy.stats import pearsonr

#region 读取Nino数据
file_path = 'D:/ENSO_index/Nino34_anom.txt'
with open(file_path, 'r') as file:
    lines = file.readlines()
years = []
data = []
for line in lines:
    values = line.split()
    years.append(int(values[0]))  # 第一列是年份
    data.append([float(v) for v in values[1:]])  # 其余列是数据

#只读取1979年-2024年的数据
years = np.array(years)[109:]
data = np.array(data)[109:]

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.bar(years,data)
ax.axhline(y=0,color='red',linestyle='--',)
ax.axhline(y=-0.5,color='red',linestyle='--',)
plt.show()