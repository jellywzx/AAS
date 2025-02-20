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
from scipy.stats import pearsonr
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
# import seaborn as sns


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


#region 冰间湖面积线形图
import os
folder_path = '/stu02/weizx24/data/opw/nsidc/new/'
files = os.listdir(folder_path)
files.sort()

# 仅罗斯海冰间湖面积
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

year_post = np.linspace(1992,2021,30)
date_name = pd.date_range('1979-11-01','1979-12-31')
date_list = list()
for i in range(len(date_name)):
    date_list.append(str(date_name[i])[5:10])

#设置全局变量
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 9
plt.rcParams['lines.linewidth'] = 1
plt.rcParams['lines.markersize'] =2
plt.rcParams['axes.linewidth'] = 0.6



colors = ['black','gray', 'silver', 'lightcoral' ,'red','chocolate','saddlebrown','darkorange','goldenrod','gold','olive','olivedrab', 
          'forestgreen','turquoise', 'teal', 'deepskyblue', 'steelblue' ,'cornflowerblue','blue','darkviolet','deeppink','palevioletred',
          'black','gray', 'silver', 'lightcoral' ,'red','chocolate','saddlebrown','darkorange','goldenrod','gold','olive','olivedrab', 
          'forestgreen','turquoise', 'teal', 'deepskyblue', 'steelblue' ,'cornflowerblue','blue','darkviolet','deeppink','palevioletred',
         ]
marker_list = ['o','^','s','p','2','P','*','h','3','H','+','D','d',
              'o','^','s','p','2','P','*','h','3','H','+','D','d',
              'o','^','s','p','2','P','*','h','3','H','+','D','d',
              'o','^','s','p','2','P','*','h','3','H','+','D','d',]

fig = plt.figure(figsize=(160/25.4,12))
ax1 = fig.add_subplot(311)
for i in range(len(a)):
    ax1.plot(date_list,a[i],
            marker=marker_list[i],
            color=colors[i],
            label=str(int(year_post[i])))
ax1.set_xlabel('Date')
ax1.set_ylabel('Areas'+ r' (10$^{6}$ km$^{2}$)')
# ax1.tick_params('both', length=8, width=1.2, which='major')
# ax1.tick_params('both', length=4, width=0.5, which='minor')
ax1.set_xticks(date_list[::5])
# ax.yaxis.set_tick_params(labelsize=20)
# ax.xaxis.set_tick_params(labelsize=20)
props = dict(boxstyl='round', facecolor='wheat', alpha=0.5)
ax1.grid(color='lightgray',linestyle='--')
plt.legend(bbox_to_anchor=(0,1),edgecolor='k',loc='upper left',ncol=3,fontsize=8)
ax1.text(0, 1.05, '(a)', fontsize=12, transform=ax1.transAxes, va='top', ha='right')


#——————————————图2——————————————————————————
#region 箱线图
dtime = pd.date_range('2001-11-01','2001-12-31',freq='D')
date_list = []
for i in range(len(dtime)):
    date_list.append(str(dtime[i])[5:10])
date_num = np.linspace(1,61,61)
# 仅罗斯海冰间湖面积
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

# # mpl.rcParams['boxplot.meanprops.color'] = 'g'
mpl.rcParams['boxplot.boxprops.color'] = 'black'
plt.rcParams['boxplot.boxprops.linewidth'] = 0.5
# mpl.rcParams['boxplot.meanprops.color'] = 'yellow'
mpl.rcParams['boxplot.whiskerprops.color'] = 'k'
plt.rcParams['boxplot.whiskerprops.linewidth'] = 0.5
mpl.rcParams['boxplot.medianprops.color'] = 'orange'
mpl.rcParams['boxplot.medianprops.linewidth'] = 0.5
plt.rcParams['boxplot.flierprops.markersize'] = 1.5
plt.rcParams['boxplot.flierprops.linewidth'] = 0.5
# fig,ax1 = plt.subplots(figsize=(16,8))

#图2
ax2 = fig.add_subplot(312)
bp = plt.boxplot(a,patch_artist=True)
plt.xticks(date_num[::5],date_list[::5])
# plt.yticks(fontsize=20)
ax2.set_ylabel('Areas'+ r' (10$^{6}$ km$^{2}$)')
ax2.set_xlabel('Date')
ax2.grid(axis='y',alpha=0.5)
plt.setp(bp['boxes'][27:34],color='red')
ax2.text(0, 1.05, '(b)', fontsize=12,transform=ax2.transAxes, va='top', ha='right')
#endregion

#region 趋势图
# 仅罗斯海冰间湖面积
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

data = a[:,27:34].mean(axis=1)
years = np.linspace(1992,2021,30)
aa, intercept, r_value, p1, std_err = stats.linregress(years, data)

color = list()
for i in range(30):
    color.append('k')
    
year_min = [1994., 2002., 2003., 2007., 2014.]
year_max = [2004., 2005., 2016., 2018., 2021.]

max_list = list()
min_list = list()
for i in range(len(years)):
    max_list.append(years[i] in year_max)
    min_list.append(years[i] in year_min)

max_mask = np.array(max_list)
min_mask = np.array(min_list)

color_arr = np.array(color)
color_arr[max_mask]='r'
color_arr[min_mask]='b'

df = pd.DataFrame({'years':years,'data':data,'color':color_arr})


ax3 = fig.add_subplot(313)
ax3.scatter(df['years'],df['data'],color=df['color'],s=10)
ax3.plot(years,data,color='k')
# plt.legend()
legend_elements = [Line2D([0], [0], marker='o', color='r', label='Large',
                          markerfacecolor='r'),
                    Line2D([0], [0], marker='o', color='b', label='Small',
                          markerfacecolor='b'),]
ax3.legend(handles=legend_elements, loc='lower right',)

#画出散点的拟合直线
aa, intercept, r_value, p1, std_err = stats.linregress(years, data)
y = aa*np.array(years)+intercept
print(aa,intercept)
ax3.plot(years,y,color='r',linestyle='--',linewidth=1)
# ax.set_ylim(0,0.4)
# ax.set_xlim(-5,5)
ax3.set_ylabel('Areas '+ r' (10$^{6}$ km$^{2}$)' )
ax3.set_xlabel('Year')
ax3.set_xticks(years[::4])
# ax.yaxis.set_tick_params(labelsize=20)
# ax.xaxis.set_tick_params(labelsize=20)
ax3.grid(linestyle='--',alpha=0.5)
# figname = 'Areas of Coastal Polynyas on 12-01'
# ax3.set_title(figname,)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#需要手动修改拟合函数的斜率和截距
ax3.text(0.05, 0.95,r'${y = 0.002x-3.65}$'+'\n'+'p value: '+str(round(p1,3)),
        transform=ax3.transAxes,
        verticalalignment='top', bbox=props,)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax3.text(0, 1.05, '(c)', fontsize=12,transform=ax3.transAxes, va='top', ha='right')
# fig.set_size_inches(150 / 25.4, 100 / 25.4)  # 转换为英寸，150mm 宽度示例
plt.subplots_adjust(hspace=0.2)
plt.savefig('/stu02/weizx24/figures/0924/Figure2_all.png' ,dpi=300,bbox_inches = 'tight')
#endregion