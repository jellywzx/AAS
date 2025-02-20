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

years = np.linspace(1992,2021,30)
year_min = [1994., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2016., 2018., 2021.]
bool_max = np.isclose(years[:, None], year_max).any(axis=1)
bool_min = np.isclose(years[:, None], year_min).any(axis=1)

a_max = a[bool_max].mean(axis=0)
a_min = a[bool_min].mean(axis=0)
a_clm = a.mean(axis=0)

# print('开始画图')
# plt.close()
date_name = pd.date_range('1979-11-01','1979-12-31')
date_list = list()
for i in range(len(date_name)):
    date_list.append(str(date_name[i])[5:10])

# fig = plt.figure(1, figsize=(160/25.4,3))
ax4 = fig.add_subplot(212)
ax4.plot(date_list[:-5],a_max[:-5],color='r',label='Large ')
ax4.plot(date_list[:-5],a_clm[:-5],color='k',label='Climatology ')
ax4.plot(date_list[:-5],a_min[:-5],color='b',label='Small ')
ax4.set_xlabel('Date')
ax4.set_ylabel('Areas'+ r' (10$^{6}$ km$^{2}$)',)
ax4.set_xticks(date_list[:-5][::5])
ax4.set_yticks([0,0.2,0.4,0.6])
ax4.grid(color='lightgray',linestyle='--',alpha=0.4)
ax4.legend(edgecolor='k',loc='upper left')
ax4.text(0, 1.05, '(d)', fontsize=10, transform=ax4.transAxes, va='top', ha='right')
plt.subplots_adjust(hspace=0.35)
plt.savefig('/stu02/weizx24/figures/0924/Figure3_all.png',dpi=300,bbox_inches='tight')
