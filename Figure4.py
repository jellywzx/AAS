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
from global_land_mask import globe

def independent_ttest(mean1,mean2,std1,std2,n1,n2, alpha):
    # calculate means
    # calculate standard errors
    #se1, se2 = sem(data1), sem(data2)
    se1, se2 = std1/sqrt(n1), std2/sqrt(n2)
    # standard error on the difference between the samples
    sed = sqrt(se1**2.0 + se2**2.0)
    # calculate the t statistic
    t_stat = (mean1 - mean2) / sed
    # degrees of freedom
    df = n1 + n2 - 2
    # calculate the critical value
    cv = t.ppf(1.0 - alpha, df)
    # calculate the p-value
    p = (1.0 - t.cdf(abs(t_stat), df)) * 2.0
    # return everything
    return p

def draw_latlon_polygon(bmap, lons, lats, *args, **kwargs):
    '''Plot a polygon in lat/lon coordinates.

    bmap -- Basemap object.
    lons, lats -- Sequences of polygon vertices.
    *args, **kwargs -- Aditional arguments to pyplot.plot().

    You should use 'k-' in *args to draw the lines in color black.

    '''
    if len(lons) != len(lats):
        raise IndexError('lons and lats have different lenghts')
    if lons[-1] != lons[0] or lats[-1] != lats[0]:
        lons = np.concatenate((lons, lons[:1]))
        lats = np.concatenate((lats, lats[:1]))
    n = len(lons) - 1
    res = 10000
    for i in range(n):
        x = np.linspace(lons[i], lons[i + 1], res)
        y = np.linspace(lats[i], lats[i + 1], res)
        x, y = bmap(x, y)
        bmap.plot(x, y, *args, **kwargs)
        
ftools='/stu02/weizx24/data/tools/'
with open(ftools+'pss25lats_v3.dat','rb') as flat:
    lats = np.fromfile(flat,dtype='<i4').reshape(332, 316)/100000.

with open(ftools+'pss25lons_v3.dat','rb') as flon:
    lons = np.fromfile(flon,dtype='<i4').reshape(332, 316)/100000.

with open(ftools+'pss25area_v3.dat','rb') as flon:
    area = np.fromfile(flon,dtype='<i4').reshape(332, 316)*1e-9

year_min = [1994., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2016., 2018., 2021.]
num_y = 5

# year_min = [2002., 2003., 2006., 2007., 2014.]
# year_max = [2004., 2005., 2010., 2016., 2021.]

# year_min = [2002., 2003., 2007., 2014.]
# year_max = [2004., 2005.,2016., 2021.]

#region SIC
sic = xr.open_dataset('/stu02/weizx24/data/monthly_sic/seaice_conc_monthly_sh_197811_202212_sub_latlon.nc')['nsidc_nt_seaice_conc_monthly'].sel(time=slice('1992','2022'))

year = np.linspace(1992,2021,30)
sic_list = []
for i in range(30):
    start_time = str(int(year[i]))+'-10-01'
    end_time = str(int(year[i]+1))+'-03-31'
    sic_list.append(sic.sel(time=slice(start_time,end_time)))
#     print(sic_sel)

sic_clm = np.nanmean(np.array(sic_list),axis=0)
sic_std = np.nanstd(np.array(sic_list),axis=0)
sic_list = []
for i in range(5):
    start_time = str(int(year_max[i]))+'-10-01'
    end_time = str(int(year_max[i]+1))+'-03-31'
    sic_list.append(sic.sel(time=slice(start_time,end_time)))
#     print(sic_sel)
sic_max = np.nanmean(np.array(sic_list),axis=0)
sic_std_max = np.nanstd(np.array(sic_list),axis=0)
sic_list = []
for i in range(5):
    start_time = str(int(year_min[i]))+'-10-01'
    end_time = str(int(year_min[i]+1))+'-03-31'
    sic_list.append(sic.sel(time=slice(start_time,end_time)))
#     print(sic_sel)
sic_min = np.nanmean(np.array(sic_list),axis=0)
sic_std_min = np.nanstd(np.array(sic_list),axis=0)
pre_clm_son = sic_max
pre_std_son = sic_std_max

sic_clm_son = sic_clm
sic_std_son = sic_std

post_clm_son = sic_min
post_std_son = sic_std_min


#6个月份
pval1=np.zeros((6,332,316),'float')
pval2=np.zeros((6,332,316),'float')
pval3=np.zeros((6,332,316),'float')

for k in range(6):
    for i in range(332):
        for j in range(316):
            pval1[k,i,j] = independent_ttest(pre_clm_son[k,i,j], sic_clm_son[k,i,j], pre_std_son[k,i,j],sic_std_son[k,i,j], num_y, 30, 0.1)
            pval2[k,i,j] = independent_ttest(post_clm_son[k,i,j], sic_clm_son[k,i,j], post_std_son[k,i,j],sic_std_son[k,i,j], num_y, 30, 0.1)
            pval3[k,i,j] = independent_ttest(pre_clm_son[k,i,j], post_clm_son[k,i,j], pre_std_son[k,i,j],post_std_son[k,i,j], num_y, num_y, 0.1)

sic_dif = (sic_max-sic_clm)-(sic_min-sic_clm)
np.savez('/stu02/weizx24/data/npz/Figure4_SIC.npz',pval3=pval3,sic_dif=sic_dif,sic_clm=sic_clm,lons=lons,lats=lats)
#endregion

#region SST
sst = xr.open_dataset('/stu02/weizx24/data/sst.mnmean_0412.nc')['sst'].sel(time=slice('1992','2022')).sel(lat=slice(-50,-90))
# sst_Feb = sic[sic.time.dt.month==2]
# lons = sst_Feb.lon
# lats = sst_Feb.lat
# sic_clm = np.nanmean(sst_Feb,axis=0)
# sic_std = np.nanstd(sst_Feb,axis=0)
lons = sst.lon
lats = sst.lat
year = np.linspace(1992,2021,30)
sst_list = []
for i in range(30):
    start_time = str(int(year[i]))+'-10-01'
    end_time = str(int(year[i]+1))+'-03-31'
    sst_list.append(sst.sel(time=slice(start_time,end_time)))
#     print(sst_sel)

sst_clm = np.nanmean(np.array(sst_list),axis=0)
sst_std = np.nanstd(np.array(sst_list),axis=0)
sst_list = []
for i in range(5):
    start_time = str(int(year_max[i]))+'-10-01'
    end_time = str(int(year_max[i]+1))+'-03-31'
    sst_list.append(sst.sel(time=slice(start_time,end_time)))
#     print(sst_sel)
sst_max = np.nanmean(np.array(sst_list),axis=0)
sst_std_max = np.nanstd(np.array(sst_list),axis=0)
sst_list = []
for i in range(5):
    start_time = str(int(year_min[i]))+'-10-01'
    end_time = str(int(year_min[i]+1))+'-03-31'
    sst_list.append(sst.sel(time=slice(start_time,end_time)))
#     print(sst_sel)
sst_min = np.nanmean(np.array(sst_list),axis=0)
sst_std_min = np.nanstd(np.array(sst_list),axis=0)
pre_clm_son = sst_max
pre_std_son = sst_std_max

sst_clm_son = sst_clm
sst_std_son = sst_std

post_clm_son = sst_min
post_std_son = sst_std_min


pval1=np.zeros((6,40,360),'float')
pval2=np.zeros((6,40,360),'float')
pval3=np.zeros((6,40,360),'float')

for k in range(6):
    for i in range(40):
        for j in range(360):
            pval1[k,i,j] = independent_ttest(pre_clm_son[k,i,j], sst_clm_son[k,i,j], pre_std_son[k,i,j],sst_std_son[k,i,j], num_y, 30, 0.1)
            pval2[k,i,j] = independent_ttest(post_clm_son[k,i,j], sst_clm_son[k,i,j], post_std_son[k,i,j],sst_std_son[k,i,j], num_y, 30, 0.1)
            pval3[k,i,j] = independent_ttest(pre_clm_son[k,i,j], post_clm_son[k,i,j], pre_std_son[k,i,j],post_std_son[k,i,j], num_y, num_y, 0.1)

sst_dif = (sst_max-sst_clm)-(sst_min-sst_clm)
np.savez('/stu02/weizx24/data/npz/Figure4_SST.npz',sst_dif=sst_dif,sst_clm=sst_clm,pval3=pval3,lons=lons,lats=lats)
#endregion

#region votemper
year_min = [1994., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2016., 2018., 2021.]
num_y = 5

# year_min = [2002., 2003., 2006., 2007., 2014.]
# year_max = [2004., 2005., 2010., 2016., 2021.]

# year_min = [2002., 2003., 2007., 2014.]
# year_max = [2004., 2005.,2016., 2021.]
# votemper = xr.open_dataset('G:/00-18 votemper.nc')['votemperer'].sel(LEV=slice(0,50)).sel(latitude=slice(-90,-50)).mean(dim='LEV')
votemper = xr.open_dataset('/stu02/weizx24/data/ORAS5/79-22_votemp_1000m_30S.nc')['votemper'].sel(LEV=slice(0,50)).mean(dim = 'LEV').sel(time=slice('1992','2022'))
# votemper_clm = np.nanmean(votemper.sel(time=slice('1992','2021')).values.reshape((30,12,60,360)),axis=0)
# votemper_std = np.nanmean(votemper.sel(time=slice('1992','2021')).values.reshape((30,12,60,360)),axis=0)
lons = votemper.lon
lats = votemper.lat

# #夏季
# votemper_clm_sel = np.concatenate([votemper_clm[11:],votemper_clm[0:2]])
# votemper_std_sel = np.concatenate([votemper_std[11:],votemper_std[0:2]])

# votemp_22 = xr.open_dataset('G:ORAS5/votemper_19-22.nc')['votemper'].sel(deptht=slice(0,50)).mean('deptht').rename({'time_counter':'time'})
year = np.linspace(1992,2021,30)
votemper_list = []
for i in range(30):
    start_time = str(int(year[i]))+'-10-01'
    end_time = str(int(year[i]+1))+'-03-31'
    votemper_list.append(votemper.sel(time=slice(start_time,end_time)))
#     print(votemper_sel)

votemper_clm = np.nanmean(np.array(votemper_list),axis=0)
votemper_std = np.nanstd(np.array(votemper_list),axis=0)
votemper_list = []
for i in range(5):
    start_time = str(int(year_max[i]))+'-10-01'
    end_time = str(int(year_max[i]+1))+'-03-31'
    votemper_list.append(votemper.sel(time=slice(start_time,end_time)))
#     print(votemper_sel)
votemper_max = np.nanmean(np.array(votemper_list),axis=0)
votemper_std_max = np.nanstd(np.array(votemper_list),axis=0)
votemper_list = []
for i in range(5):
    start_time = str(int(year_min[i]))+'-10-01'
    end_time = str(int(year_min[i]+1))+'-03-31'
    votemper_list.append(votemper.sel(time=slice(start_time,end_time)))
#     print(votemper_sel)
votemper_min = np.nanmean(np.array(votemper_list),axis=0)
votemper_std_min = np.nanstd(np.array(votemper_list),axis=0)
pre_clm_son = votemper_max
pre_std_son = votemper_std  _max

votemper_clm_son = votemper_clm
votemper_std_son = votemper_std

post_clm_son = votemper_min
post_std_son = votemper_std_min


pval1=np.zeros((6,60,360),'float')
pval2=np.zeros((6,60,360),'float')
pval3=np.zeros((6,60,360),'float')

for k in range(6):
    for i in range(60):
        for j in range(360):
            pval1[k,i,j] = independent_ttest(pre_clm_son[k,i,j], votemper_clm_son[k,i,j], pre_std_son[k,i,j],votemper_std_son[k,i,j], num_y, 30, 0.1)
            pval2[k,i,j] = independent_ttest(post_clm_son[k,i,j], votemper_clm_son[k,i,j], post_std_son[k,i,j],votemper_std_son[k,i,j], num_y, 30, 0.1)
            pval3[k,i,j] = independent_ttest(pre_clm_son[k,i,j], post_clm_son[k,i,j], pre_std_son[k,i,j],post_std_son[k,i,j], num_y, num_y, 0.1)

votemper_dif = (votemper_max-votemper_clm)-(votemper_min-votemper_clm)
np.savez('/stu02/weizx24/data/npz/Figure4_votemper.npz',votemper_dif=votemper_dif,votemper_clm=votemper_clm,lons=lons,lats=lats,pval3=pval3)
#endregion

print('画图完成')