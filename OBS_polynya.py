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
np.set_printoptions(suppress=True)
np.set_printoptions(precision=9)
ftools='/stu02/weizx24/data/tools/'
with open(ftools+'pss25lats_v3.dat','rb') as flat:
    lats = np.fromfile(flat,dtype='<i4').reshape(332, 316)/100000.

with open(ftools+'pss25lons_v3.dat','rb') as flon:
    lons = np.fromfile(flon,dtype='<i4').reshape(332, 316)/100000.

with open(ftools+'pss25area_v3.dat','rb') as flon:
    area = np.fromfile(flon,dtype='<i4').reshape(332, 316)*1e-9

sic = xr.open_dataset('/stu02/weizx24/data/Daily_SIC/process/nt_seaice_78-22.nc')['nsidc_nt_seaice_conc']
sic_SO = sic[(sic.time.dt.month==9)|(sic.time.dt.month==10)].sel(time=slice('1992','2022'))
east_ross_mask = (((lons>-180.)&(lons<-170.)&(lats<-69.))|((lons>-170.)&(lons<-130.)&(lats<-70.)))
west_ross_mask = ((lons>160.)&(lats<-68.))
amun_mask = ((lons>-130.)&(lons<-90.)&(lats<-71.))
bell_mask = ((lons>-90.)&(lons<-60.)&(lats<-71.))
for i in range(len(sic_SO)):
#     mask_ross_nt = (((lons>-180.)&(lons<-90.)&(lats<-70.))|((lons>160.)&(lats<-71.5)))&(sic_1992[i]<0.15)
    mask_ross_nt = (west_ross_mask|east_ross_mask|amun_mask|bell_mask)&(sic_SO[i]<0.15)
#     mask_ross_nt = ((lons>160.)&(lats<-69.)&(sic_1992[i]<0.15))|((lons>-180.)&(lons<-60.)&(lats<-70.)&(sic_1992[i]<0.15))
    sic_ross_nt = np.where(mask_ross_nt,area,np.nan)
    
    fig = plt.figure(figsize=(32,6.5))
    levels1 = np.linspace(0,100,21)
    plt.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['hatch.color'] ='springgreen'
    mpl.rcParams['hatch.linewidth'] = 0.8

    pval1 = np.zeros((332,316),'float')
    #anomaly大于一倍标准差代表显著
    # pval1[np.where(abs(sic_2010_anom[i])>sic_std_conc[i]*1)]=0.99
    figname = str(sic_SO[i].time.values)[0:10]+' open water within ice pack'
    plt.title(figname,fontsize=16)
    m = Basemap(width=6000000,height=4000000,
                resolution='l',projection='stere',\
                lat_ts=-82,lat_0=-73.,lon_0=220.)
    x1,y1 = m(lons,lats)
    im1 = m.contourf(x1,y1,sic_ross_nt*100,
                         levels=levels1,
                         extend='both',
                        shading='faceted', antialiased=True,cmap='RdBu_r')
    m.contour(x1,y1,sic_SO[i],levels=[0.15],colors='g',linewidths=2)
    m.drawcoastlines()
    m.fillcontinents(color='lightgray')
    m.drawmeridians([90,120,150,180,210,240,270,300,330], linewidth=1.2,labels=[True,True,False,False],fmt='%g',latmax=80)
    m.drawparallels([-60,-70,-80], linewidth=1.2,labels=[False,False,False,True],latmax=80)

    fig.tight_layout(rect=[0,0.08,1,0.95])
    plt.savefig('/stu02/weizx24/figures/'+figname+'.png',dpi=600,bbox_inches='tight')
    plt.show()
