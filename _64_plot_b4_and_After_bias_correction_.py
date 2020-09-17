'''
Created on 17 Sep 2020

@author: hachem
'''
__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# =============================================================================

import os
import timeit
import time

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from matplotlib import rc
from matplotlib import rcParams
from pandas.plotting import register_matplotlib_converters

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears
from scipy.spatial import cKDTree


from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


from _00_additional_functions import (resample_intersect_2_dfs,
                                      select_convective_season,
                                      select_df_within_period,
                                      build_edf_fr_vals,
                                      get_cdf_part_abv_thr,
                                      plt_on_map_agreements,
                                      plt_correlation_with_distance)

from _01_2_read_hdf5 import HDF5


register_matplotlib_converters()

plt.ioff()

rc('font', size=16)
rc('font', family='arial')
rc('axes', labelsize=16)
rcParams['axes.labelpad'] = 10
plt.rcParams['axes.labelweight'] = 'bold'
#==============================================================================
def plot_histogram(x, density=True, nbins=10):
    hist, bins = np.histogram(x, bins=nbins, density=density)
    width = 0.5 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    return hist, width, center

#==============================================================================
main_dir = Path(r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\bias_correction_paper_revie2")
os.chdir(main_dir)


path_to_dwd_h5 = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\rain_bw_1hour"
    r"\netatmo_bw__60min_20152019.h5")

pws_hdf5 = HDF5(infile=path_to_dwd_h5)
pws_ids = pws_hdf5.get_all_names()

#df_stn = dwd_hdf5.get_pandas_dataframe(stn_name)
pws_plus = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
    r"\bias_correction_paper_revie2\netatmo_stn_70_ee_50_00_65_7e.csv")

pws_minus = (
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\bias_correction_paper_revie2"
    r"\netatmo_stn_70_ee_50_01_e5_a6.csv")

id_minus = '70_ee_50_01_e5_a6'
id_plus = '70_ee_50_00_65_7e'


netatmo_plus = pd.read_csv(
    pws_plus, index_col=0, sep=';',
    parse_dates=True, infer_datetime_format=True)

netatmo_minus = pd.read_csv(
    pws_minus, index_col=0, sep=';',
    parse_dates=True, infer_datetime_format=True)

netatmo_plus_b4 = pws_hdf5.get_pandas_dataframe_bet_dates(id_plus,
                                                          start_date=netatmo_plus.index[0],
                                                          end_date=netatmo_plus.index[-1])

netatmo_minus_b4 = pws_hdf5.get_pandas_dataframe_bet_dates(id_minus,
                                                          start_date=netatmo_minus.index[0],
                                                          end_date=netatmo_minus.index[-1])

netatmo_plus_b4_res = netatmo_plus_b4.resample('M').sum()
netatmo_plus_b4_res = netatmo_plus_b4_res[netatmo_plus_b4_res >1].dropna(how='all')
netatmo_plus_res = netatmo_plus.resample('M').sum()
netatmo_plus_res = netatmo_plus_res[netatmo_plus_res > 1].dropna(how='all')


netatmo_minus_b4_res = netatmo_minus_b4.resample('M').sum()
netatmo_minus_b4_res = netatmo_minus_b4_res[netatmo_minus_b4_res > 1].dropna(how='all')
netatmo_minus_res = netatmo_minus.resample('M').sum()

netatmo_minus_res = netatmo_minus_res[netatmo_minus_res >1].dropna(how='all')

import matplotlib.dates as mdates
# Set the locator
locator = mdates.MonthLocator()  # every month
# Specify the format - %b gives us Jan, Feb...
fmt = mdates.DateFormatter('%Y-%m')

plt.ioff()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 12),
                               sharex=True, sharey=True, dpi=100)

max_y = max(max(netatmo_minus_b4_res.values.ravel()), max(netatmo_plus_res.values.ravel()))
ax1.set_ylim([0, max_y + 3])
ax2.set_yticks(np.arange(0, max_y + 3, 25))

ms = 10
alpha_v = 0.25

ax1.plot(netatmo_minus_b4_res.index, netatmo_minus_b4_res.values.ravel(),
        alpha=0.95, color='blue', 
        linewidth=1, label='Before')
ax1.scatter(netatmo_minus_b4_res.index, netatmo_minus_b4_res.values.ravel(),s=ms,
        alpha=0.95, color='blue', marker='o')

ax1.plot(netatmo_minus_res.index, netatmo_minus_res.values.ravel(),
        alpha=0.95, color='r', 
        linewidth=1, label='After')
ax1.scatter(netatmo_minus_res.index, netatmo_minus_res.values.ravel(), s=ms,
        alpha=0.95, color='r', marker='+')
ax1.axvline(x='2015-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)
ax1.axvline(x='2016-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)
ax1.axvline(x='2017-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)
ax1.axvline(x='2018-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)
ax1.axvline(x='2019-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)

#ax1.axhline(y=netatmo_minus_b4_res.values.mean(), color='b', linestyle='--', alpha=0.25)
#ax1.axhline(y=netatmo_minus_res.values.mean(), color='r', linestyle='--', alpha=0.25)

ax1.set_xticks([
                   '2015-06-30 00:00:00',
                   '2016-06-30 00:00:00',
                   '2017-06-30 00:00:00',
                   '2018-06-30 00:00:00',
                   '2019-06-30 00:00:00'],  # 03 bor BW
                  minor=False)

#ax1.set_xlabel('Month', fontsize=18)
ax1.grid(alpha=0.35)
ax1.set_ylabel('[mm/month]', fontsize=18)
ax1.legend(loc='upper right')
ax1.title.set_text('Bias correction - PWS %s' % id_minus)

#ax1.xaxis.set_minor_locator(AutoMinorLocator())
#ax1.tick_params(which='both', width=2)
#ax1.tick_params(which='major', length=7)
#ax1.tick_params(which='minor', length=3, color='grey')
 
 
ax2.plot(netatmo_plus_b4_res.index, netatmo_plus_b4_res.values.ravel(),
        alpha=0.95, color='blue', 
        linewidth=1, label='Before')
ax2.plot(netatmo_plus_res.index, netatmo_plus_res.values.ravel(),
        alpha=0.95, color='r', 
        linewidth=1, label='After')

 
ax2.scatter(netatmo_plus_b4_res.index, netatmo_plus_b4_res.values.ravel(),
        alpha=0.95, color='blue', marker='o', s=ms)
ax2.scatter(netatmo_plus_res.index, netatmo_plus_res.values.ravel(), s=ms,
        alpha=0.95, color='r', marker='+')


ax2.axvline(x='2015-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)
ax2.axvline(x='2016-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)
ax2.axvline(x='2017-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)
ax2.axvline(x='2018-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)
ax2.axvline(x='2019-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)


#ax2.axhline(y=netatmo_plus_b4_res.values.mean(), color='b', linestyle='--', alpha=0.25)
#ax2.axhline(y=netatmo_plus_res.values.mean(), color='r', linestyle='--', alpha=0.25)

ax2.set_xticks([
                   '2015-06-30 00:00:00',
                   '2016-06-30 00:00:00',
                   '2017-06-30 00:00:00',
                   '2018-06-30 00:00:00',
                   '2019-06-30 00:00:00'],  # 03 bor BW
                  minor=False)

#ax1.set_xlabel('Month', fontsize=18)
ax2.grid(alpha=0.35)
#ax2.set_ylabel('[mm/month]', fontsize=18)
ax2.legend(loc='upper right')
ax2.title.set_text('Bias correction - PWS %s' % id_plus)
#'Increase in the sum [April-October]')

#ax2.xaxis.set_major_locator(AutoMinorLocator())
 
#ax2.tick_params(which='both', width=2)
#ax2.tick_params(which='major', length=7)
#ax2.tick_params(which='minor', length=3, color='grey')
 
 
fig.tight_layout()
fig.savefig(
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\bias_correction_paper_revie2\corr_.png",
    frameon=True, papertype='a4',
    bbox_inches='tight', pad_inches=.2)
plt.close()


#===============================================================================

netatmo_minus_b4_res_cm = netatmo_minus_b4_res.cumsum()
netatmo_minus_res_cm = netatmo_minus_res.cumsum()

netatmo_plus_b4_res_cm = netatmo_plus_b4_res.cumsum()

netatmo_plus_res_cmn = netatmo_plus_res.cumsum()


plt.ioff()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 12),
                               sharex=True, sharey=True, dpi=100)

max_y = max(max(netatmo_minus_b4_res_cm.values.ravel()), max(netatmo_plus_res_cmn.values.ravel()))
ax1.set_ylim([0, max_y + 3])
ax2.set_yticks(np.arange(0, max_y + 20, 100))

ms = 10
alpha_v = 0.25

ax1.plot(netatmo_minus_b4_res_cm.index, netatmo_minus_b4_res_cm.values.ravel(),
        alpha=0.95, color='blue', 
        linewidth=1, label='Before')
ax1.scatter(netatmo_minus_b4_res_cm.index, netatmo_minus_b4_res_cm.values.ravel(),s=ms,
        alpha=0.95, color='blue', marker='o')

ax1.plot(netatmo_minus_res_cm.index, netatmo_minus_res_cm.values.ravel(),
        alpha=0.95, color='r', 
        linewidth=1, label='After')
ax1.scatter(netatmo_minus_res_cm.index, netatmo_minus_res_cm.values.ravel(), s=ms,
        alpha=0.95, color='r', marker='+')
ax1.axvline(x='2015-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)
ax1.axvline(x='2016-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)
ax1.axvline(x='2017-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)
ax1.axvline(x='2018-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)
ax1.axvline(x='2019-08-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)

#ax1.axhline(y=netatmo_minus_b4_res.values.mean(), color='b', linestyle='--', alpha=0.25)
#ax1.axhline(y=netatmo_minus_res.values.mean(), color='r', linestyle='--', alpha=0.25)

ax1.set_xticks([
                   '2015-06-30 00:00:00',
                   '2016-06-30 00:00:00',
                   '2017-06-30 00:00:00',
                   '2018-06-30 00:00:00',
                   '2019-06-30 00:00:00'],  # 03 bor BW
                  minor=False)

#ax1.set_xlabel('Month', fontsize=18)
ax1.grid(alpha=0.35)
ax1.set_ylabel('[mm]', fontsize=18)
ax1.legend(loc='lower right')
ax1.title.set_text('Bias correction - PWS %s' % id_minus)

#ax1.xaxis.set_minor_locator(AutoMinorLocator())
#ax1.tick_params(which='both', width=2)
#ax1.tick_params(which='major', length=7)
#ax1.tick_params(which='minor', length=3, color='grey')
 
 
ax2.plot(netatmo_plus_b4_res_cm.index, netatmo_plus_b4_res_cm.values.ravel(),
        alpha=0.95, color='blue', 
        linewidth=1, label='Before')
ax2.plot(netatmo_plus_res_cmn.index, netatmo_plus_res_cmn.values.ravel(),
        alpha=0.95, color='r', 
        linewidth=1, label='After')

 
ax2.scatter(netatmo_plus_b4_res_cm.index, netatmo_plus_b4_res_cm.values.ravel(),
        alpha=0.95, color='blue', marker='o', s=ms)
ax2.scatter(netatmo_plus_res_cmn.index, netatmo_plus_res_cmn.values.ravel(), s=ms,
        alpha=0.95, color='r', marker='+')


ax2.axvline(x='2015-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)
ax2.axvline(x='2016-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)
ax2.axvline(x='2017-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)
ax2.axvline(x='2018-12-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)
ax2.axvline(x='2019-08-31 00:00:00', color='k', linestyle='--', alpha=alpha_v)


#ax2.axhline(y=netatmo_plus_b4_res.values.mean(), color='b', linestyle='--', alpha=0.25)
#ax2.axhline(y=netatmo_plus_res.values.mean(), color='r', linestyle='--', alpha=0.25)

ax2.set_xticks([
                   '2015-06-30 00:00:00',
                   '2016-06-30 00:00:00',
                   '2017-06-30 00:00:00',
                   '2018-06-30 00:00:00',
                   '2019-06-30 00:00:00'],  # 03 bor BW
                  minor=False)

#ax1.set_xlabel('Month', fontsize=18)
ax2.grid(alpha=0.35)
#ax2.set_ylabel('[mm/month]', fontsize=18)
ax2.legend(loc='lower right')
ax2.title.set_text('Bias correction - PWS %s' % id_plus)
#'Increase in the sum [April-October]')

#ax2.xaxis.set_major_locator(AutoMinorLocator())
 
#ax2.tick_params(which='both', width=2)
#ax2.tick_params(which='major', length=7)
#ax2.tick_params(which='minor', length=3, color='grey')
 
 
fig.tight_layout()
fig.savefig(
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\bias_correction_paper_revie2\corr_2.png",
    frameon=True, papertype='a4',
    bbox_inches='tight', pad_inches=.2)
plt.close()
#===============================================================================
#===============================================================================
# 
#===============================================================================
#===============================================================================
start = 2015
end = '2015'

netatmo_plus_b4_res = netatmo_plus_b4.resample('M').sum()
netatmo_plus_b4_res = netatmo_plus_b4_res[netatmo_plus_b4_res >1].dropna(how='all')
netatmo_plus_res = netatmo_plus.resample('M').sum()
netatmo_plus_res = netatmo_plus_res[netatmo_plus_res > 1].dropna(how='all')


netatmo_minus_b4_res = netatmo_minus_b4.resample('M').sum()
netatmo_minus_b4_res = netatmo_minus_b4_res[netatmo_minus_b4_res > 1].dropna(how='all')
netatmo_minus_res = netatmo_minus.resample('M').sum()

netatmo_minus_res = netatmo_minus_res[netatmo_minus_res >1].dropna(how='all')

one_year_minus = np.where(netatmo_minus_b4_res.index.year == start)[0]
one_year_plus = np.where(netatmo_plus_b4_res.index.year == start)[0]


netatmo_minus_b4_res = netatmo_minus_b4_res.iloc[one_year_minus]
netatmo_minus_res = netatmo_minus_res.iloc[one_year_minus]

netatmo_plus_b4_res = netatmo_plus_b4_res.iloc[one_year_plus]
netatmo_plus_res = netatmo_plus_res.iloc[one_year_plus]

labels =  {4:'April', 5:'Mai', 6:'June', 7:'July', 8:'August', 9:'September', 10:'October'}
labels_str = []

for _xx in netatmo_minus_b4_res.index.month:
    for ll in labels.keys():
        if _xx == ll:
            labels_str.append(labels[ll])
            

x = np.arange(len(netatmo_minus_b4_res.index))  # the label locations
x2 = np.arange(len(netatmo_plus_b4_res.index))  # the label locations
width = 0.4  # the width of the bars

plt.ioff()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 12),
                               sharex=True, sharey=True, dpi=100)

max_y = max(max(netatmo_minus_b4_res.values.ravel()),
            max(netatmo_plus_res.values.ravel()))
ax1.set_ylim([0, max_y + 3])
ax1.set_yticks(np.arange(0, max_y + 5, 10))


ms = 10
alpha_v = 0.25

ax1.bar(x-width/2,  netatmo_minus_b4_res.values.ravel(), 0.35,
        alpha=0.95, color='blue', 
        linewidth=1, label='Before')

ax1.bar(x+width/2 ,  netatmo_minus_res.values.ravel(),0.35,
        alpha=0.95, color='r',  
        linewidth=1, label='After')

#ax1.xaxis.set_minor_locator(AutoMinorLocator())

ax1.tick_params(which='both', width=2)
ax1.tick_params(which='major', length=7)
ax1.tick_params(which='minor', length=3, color='grey')

ax1.set_xticks(x)
ax1.set_xticklabels(labels_str)
#ax1.set_xlabel('Month', fontsize=18)
ax1.grid(alpha=0.35)
ax1.set_ylabel('[mm/month]', fontsize=18)
ax1.legend(loc='upper right')
ax1.title.set_text('Bias correction for %s - PWS ID1' % (end))

 
ax2.bar(x2 - width/2,  netatmo_plus_b4_res.values.ravel(), 0.35,
        alpha=0.95, color='b', 
        linewidth=1, label='Before')
ax2.bar(x2 + width/2,  netatmo_plus_res.values.ravel(), 0.35,
        alpha=0.95, color='r', 
        linewidth=1, label='After')
ax2.set_xticks(x)
ax2.set_xticklabels(labels_str)

#ax2.xaxis.set_minor_locator(AutoMinorLocator())

ax2.tick_params(which='both', width=2)
ax2.tick_params(which='major', length=7)
ax2.tick_params(which='minor', length=3, color='grey')

ax2.grid(alpha=0.35)
#ax2.set_ylabel('[mm/month]', fontsize=18)
ax2.legend(loc='upper right')
ax2.title.set_text('Bias correction for %s - PWS ID2' % (end))
 
 
fig.tight_layout()
fig.savefig(
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\bias_correction_paper_revie2\2015.png",
    frameon=True, papertype='a4',
    bbox_inches='tight', pad_inches=.2)
plt.close()




#===============================================================================