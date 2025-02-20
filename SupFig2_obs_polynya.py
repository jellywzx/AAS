#本文的目的是：计算正文中Figure3的clm, max, min 12-01的平均冰间湖面积大小



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

sic = xr.open_dataset('/stu02/weizx24/data/Daily_SIC/process/nt_seaice_78-22.nc')['nsidc_nt_seaice_conc'].sel(time=slice('1992','2021'))
ftools='/stu02/weizx24/data/tools/'
with open(ftools+'pss25lats_v3.dat','rb') as flat:
    lats = np.fromfile(flat,dtype='<i4').reshape(332, 316)/100000.

with open(ftools+'pss25lons_v3.dat','rb') as flon:
    lons = np.fromfile(flon,dtype='<i4').reshape(332, 316)/100000.

with open(ftools+'pss25area_v3.dat','rb') as flon:
    area = np.fromfile(flon,dtype='<i4').reshape(332, 316)*1e-9

east_ross_mask = ((lons>-180.)&(lons<-130.)&(lats<-70.))
west_ross_mask = ((lons>160.)&(lats<-68.))
amun_mask = ((lons>-130.)&(lons<-90.)&(lats<-70.))
bell_mask = ((lons>-90.)&(lons<-60.)&(lats<-70.))

from global_land_mask import globe
globe_land_mask = globe.is_land(lats, lons)

year = np.linspace(1992,2021,30)
sic_list = []
for i in range(30):
    sel_time = str(int(year[i]))+'-12-01'
    sic_list.append(sic.sel(time=sel_time))
sic_1201_avg = np.nanmean(np.array(sic_list),axis=0)

year_min = [1994., 1997., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2010., 2016., 2018., 2021.]
num_y = 6

sic_list = []
for i in range(6):
    sel_time = str(int(year_max[i]))+'-12-01'
    sic_list.append(sic.sel(time=sel_time))
sic_max = np.nanmean(np.array(sic_list),axis=0)

year_min = [1994., 1997., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2010., 2016., 2018., 2021.]
num_y = 6
sic_list = []
for i in range(6):
    sel_time = str(int(year_min[i]))+'-12-01'
    sic_list.append(sic.sel(time=sel_time))
sic_min = np.nanmean(np.array(sic_list),axis=0)

#     mask_ross_nt = (((lons>-180.)&(lons<-90.)&(lats<-70.))|((lons>160.)&(lats<-71.5)))&(sic_1992[i]<0.15)
mask_ross_nt = (west_ross_mask|east_ross_mask)&(sic_max<0.15)
#     mask_ross_nt = ((lons>160.)&(lats<-69.)&(sic_1992[i]<0.15))|((lons>-180.)&(lons<-60.)&(lats<-70.)&(sic_1992[i]<0.15))
sic_ross_nt = np.where(mask_ross_nt,area,np.nan)
ross_sie_max = np.nansum(sic_ross_nt)