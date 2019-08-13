
# coding: utf-8

# In[6]:


import pandas as pd
import os

import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import curve_fit

# %%
# def func(x, a, b, c):
#    return a * np.exp(-b * x) + c

percent = '95'
# In[2]:


# df0 = r'/home/abbas/Desktop/home-office/year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_60min_dwd_netatmo_upper_99_percent_data_considered_neighbor_0_.csv'
#
# df1 = r'/home/abbas/Desktop/home-office/year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_60min_dwd_netatmo_upper_99_percent_data_considered_neighbor_1_.csv'
# df2 = r'/home/abbas/Desktop/home-office/year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_60min_dwd_netatmo_upper_99_percent_data_considered_neighbor_2_.csv'
# df3 = r'/home/abbas/Desktop/home-office/year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_60min_dwd_netatmo_upper_99_percent_data_considered_neighbor_3_.csv'
df0 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_60min_dwd_netatmo_upper_%s_percent_data_considered_neighbor_0_.csv" % percent
df1 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_60min_dwd_netatmo_upper_%s_percent_data_considered_neighbor_1_.csv" % percent
df2 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_60min_dwd_netatmo_upper_%s_percent_data_considered_neighbor_2_.csv" % percent
df3 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_60min_dwd_netatmo_upper_%s_percent_data_considered_neighbor_3_.csv" % percent
df4 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_60min_dwd_netatmo_upper_%s_percent_data_considered_neighbor_4_.csv" % percent

#df0 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_indicator_correlations_DWD_DWD_ppt_\year_allyears_df_dwd_correlationsfreq_60min_dwd_netatmo_upper_95_percent_data_considered_neighbor_1_.csv"
#df1 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_indicator_correlations_DWD_DWD_ppt_\year_allyears_df_dwd_correlationsfreq_60min_dwd_netatmo_upper_95_percent_data_considered_neighbor_2_.csv"


# In[3]:


in_df0 = pd.read_csv(df0, index_col=0, sep=';')
in_df1 = pd.read_csv(df1, index_col=0, sep=';')
in_df2 = pd.read_csv(df2, index_col=0, sep=';')
in_df3 = pd.read_csv(df3, index_col=0, sep=';')
in_df4 = pd.read_csv(df4, index_col=0, sep=';')

in_df0.dropna(how='any', inplace=True)
in_df1.dropna(how='any', inplace=True)
in_df2.dropna(how='any', inplace=True)
in_df3.dropna(how='any', inplace=True)
in_df4.dropna(how='any', inplace=True)

in_df0 = in_df0[in_df0['Bool_Spearman_Correlation'] < 1]
in_df1 = in_df1[in_df1['Bool_Spearman_Correlation'] < 1]
in_df2 = in_df2[in_df2['Bool_Spearman_Correlation'] < 1]
in_df3 = in_df3[in_df3['Bool_Spearman_Correlation'] < 1]
in_df4 = in_df4[in_df4['Bool_Spearman_Correlation'] < 1]
# In[5]:

s0 = in_df0.index.values
s1 = in_df1.index.values
s2 = in_df2.index.values
s3 = in_df3.index.values
s4 = in_df4.index.values

x0 = in_df0.loc[:, 'Distance to neighbor'].values
x1 = in_df1.loc[:, 'Distance to neighbor'].values
x2 = in_df2.loc[:, 'Distance to neighbor'].values
x3 = in_df3.loc[:, 'Distance to neighbor'].values
x4 = in_df4.loc[:, 'Distance to neighbor'].values

y0 = in_df0.loc[:, 'Bool_Spearman_Correlation'].values
y1 = in_df1.loc[:, 'Bool_Spearman_Correlation'].values
y2 = in_df2.loc[:, 'Bool_Spearman_Correlation'].values
y3 = in_df3.loc[:, 'Bool_Spearman_Correlation'].values
y4 = in_df4.loc[:, 'Bool_Spearman_Correlation'].values

# %% create funtion to fit and optimization scheme


def func(x, a, b, c, d):
    return a * x**3 + b * x**2 + c * x + d


def fit_curve_get_vals_below_curve(x, y, func, stns):
    popt, pcov = curve_fit(func, x, y)  # , bounds=(0, np.inf))
    y_fitted = func(x, *popt)

    lower_bound = y_fitted.max() * 0.10
    y_fitted_shifted = y_fitted - lower_bound

    xvals_below_curve = x[np.where(y <= y_fitted_shifted)]
    yvals_below_curve = y[np.where(y <= y_fitted_shifted)]

    stns_below_curve = stns[np.where(y <= y_fitted_shifted)]

    xvals_above_curve = x[np.where(y > y_fitted_shifted)]
    yvals_above_curve = y[np.where(y > y_fitted_shifted)]

    stns_above_curve = stns[np.where(y > y_fitted_shifted)]

    return (y_fitted_shifted, xvals_below_curve, yvals_below_curve,
            xvals_above_curve, yvals_above_curve, stns_below_curve,
            stns_above_curve)


# %% apply filter for every neighbor seperatly

(y_fitted_shifted0, xvals_below_curve0,
 yvals_below_curve0, xvals_above_curve0,
 yvals_above_curve0, stnbelow0, stnabove0) = fit_curve_get_vals_below_curve(
    x0, y0, func, s0)

(y_fitted_shifted1, xvals_below_curve1,
 yvals_below_curve1, xvals_above_curve1,
 yvals_above_curve1, stnbelow1, stnabove1) = fit_curve_get_vals_below_curve(
    x1, y1, func, s1)
(y_fitted_shifted2, xvals_below_curve2,
 yvals_below_curve2, xvals_above_curve2,
 yvals_above_curve2, stnbelow2, stnabove2) = fit_curve_get_vals_below_curve(
    x2, y2, func, s2)

(y_fitted_shifted3, xvals_below_curve3,
 yvals_below_curve3, xvals_above_curve3,
 yvals_above_curve3, stnbelow3, stnabove3) = fit_curve_get_vals_below_curve(
    x3, y3, func, s3)

(y_fitted_shifted4, xvals_below_curve4,
 yvals_below_curve4, xvals_above_curve4,
 yvals_above_curve4, stnbelow4, stnabove4) = fit_curve_get_vals_below_curve(
    x4, y4, func, s4)

# np.min(np.concatenate((y_fitted_shifted0, y_fitted_shifted1,
#  y_fitted_shifted2, y_fitted_shifted3)))
# np.intersect1d(s2, s3)
# %%

# intersect all stations

stns_keep_all = np.intersect1d(np.intersect1d(np.intersect1d(np.intersect1d(
    stnabove0, stnabove1),
    stnabove2), stnabove3), stnabove4)

stns_keep_al_sr = pd.Series(stns_keep_all)
stns_keep_al_sr.to_csv(
    (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
        r'\filter_Netamo_data_basedon_indicator_correlation'
        r'\keep_stns_all_neighbors_combined_%s_per_.csv' % percent),
    sep=';')
# %%
x0_abv = in_df0.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
y0_abv = in_df0.loc[stns_keep_all, 'Bool_Spearman_Correlation'].dropna().values

x1_abv = in_df1.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
y1_abv = in_df1.loc[stns_keep_all, 'Bool_Spearman_Correlation'].dropna().values

x2_abv = in_df2.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
y2_abv = in_df2.loc[stns_keep_all, 'Bool_Spearman_Correlation'].dropna().values

x3_abv = in_df3.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
y3_abv = in_df3.loc[stns_keep_all, 'Bool_Spearman_Correlation'].dropna().values

x4_abv = in_df4.loc[stns_keep_all, 'Distance to neighbor'].dropna().values
y4_abv = in_df4.loc[stns_keep_all, 'Bool_Spearman_Correlation'].dropna().values
#==============================================================================
#
#==============================================================================
plt.ioff()
plt.figure(figsize=(12, 12), dpi=100)


plt.scatter(x0_abv, y0_abv, c='r', alpha=0.5,
            marker='x', label='First Neighbor', s=28)
plt.scatter(x1_abv, y1_abv, c='orange', alpha=0.5,
            marker='.', label='Second Neighbor', s=28)
plt.scatter(x2_abv, y2_abv, c='g', alpha=0.5,
            marker='d', label='Third Neighbor', s=28)
plt.scatter(x3_abv, y3_abv, c='b', alpha=0.5,
            marker='1', label='Fourth Neighbor', s=28)
plt.scatter(x4_abv, y4_abv, c='m', alpha=0.5,
            marker='1', label='Fifth Neighbor', s=28)

plt.scatter(x0, y_fitted_shifted0, c='darkred', alpha=0.5,
            marker='x', label='Fitted curve 1', s=20)

plt.scatter(x1, y_fitted_shifted1, c='darkorange', alpha=0.5,
            marker='x', label='Fitted curve 2', s=20)

plt.scatter(x2, y_fitted_shifted2, c='darkgreen', alpha=0.5,
            marker='x', label='Fitted curve 3', s=20)

plt.scatter(x3, y_fitted_shifted3, c='darkblue', alpha=0.5,
            marker='x', label='Fitted curve 4', s=20)
plt.scatter(x4, y_fitted_shifted4, c='purple', alpha=0.5,
            marker='x', label='Fitted curve 5', s=20)
#
# plt.scatter(x0_abv, y0_abv, c='r', alpha=0.5,
#            marker='x', label='First Neighbor', s=28)
#
#
# plt.scatter(xvals_above_curve0, yvals_above_curve0, c='r',
#             alpha=0.5, marker='*', label='First Neighbor', s=28)
# plt.scatter(xvals_above_curve1, yvals_above_curve1, c='orange',
#             alpha=0.5, marker='.', label='Second Neighbor', s=28)
# plt.scatter(xvals_above_curve2, yvals_above_curve2, c='green',
#             alpha=0.5, marker='d', label='Third Neighbor', s=28)
# plt.scatter(xvals_above_curve3, yvals_above_curve3, c='b',
#             alpha=0.5, marker='1', label='Fourth Neighbor', s=28)
# #
#
# plt.scatter(xvals_below_curve0, yvals_below_curve0, c='pink',
#             alpha=0.5, marker='*')  # , label='First Neighbor', s=28)
# plt.scatter(xvals_below_curve1, yvals_below_curve1, c='y',
#             alpha=0.5, marker='.')  # , label='Second Neighbor', s=28)
# plt.scatter(xvals_below_curve2, yvals_below_curve2, c='lightgreen',
#             alpha=0.5, marker='d')  # , label='Third Neighbor', s=28)
# plt.scatter(xvals_below_curve3, yvals_below_curve3, c='lightblue',
#             alpha=0.5, marker='1')  # , label='Fourth Neighbor', s=28)
#

plt.xlim([0, max([x0.max(), x1.max(), x2.max(), x3.max()]) + 1000])
plt.ylim([-0.1, 1.1])
plt.xlabel('Distance (m)')
plt.ylabel('Indicator Spearman Correlation')
plt.legend(loc=0)
plt.grid(alpha=.25)
# plt.tight_layout()
plt.title('Keeping Netatmo %d stations: Indicator correlation'
          ' with distance for upper %s percent of data values' %
          (stns_keep_all.shape[0], percent))
# plt.savefig(os.path.join(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature', 'Netatmo_98percent_indic_corr.png'),
# frameon=True, papertype='a4',
# bbox_inches='tight', pad_inches=.2)
plt.show()
#
