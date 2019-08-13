
# coding: utf-8

# In[6]:

# TODO: NEED FIXING
import pandas as pd
import os

import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import curve_fit

# %%
# def func(x, a, b, c):
#    return a * np.exp(-b * x) + c


def func(x, a, b, c, d):
    return a * x**3 + b * x**2 + c * x + d


def fit_curve_get_vals_below_curve(x, y, func, stns):
    popt, pcov = curve_fit(func, x, y)
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

# df0 = r'/home/abbas/Desktop/home-office/year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_60min_dwd_netatmo_upper_99_percent_data_considered_neighbor_0_.csv'


# df1 = r'/home/abbas/Desktop/home-office/year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_60min_dwd_netatmo_upper_99_percent_data_considered_neighbor_1_.csv'
# df2 = r'/home/abbas/Desktop/home-office/year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_60min_dwd_netatmo_upper_99_percent_data_considered_neighbor_2_.csv'
# df3 = r'/home/abbas/Desktop/home-office/year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_60min_dwd_netatmo_upper_99_percent_data_considered_neighbor_3_.csv'
df0 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_60min_dwd_netatmo_upper_90_percent_data_considered_neighbor_0_.csv"
df1 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_60min_dwd_netatmo_upper_90_percent_data_considered_neighbor_1_.csv"
df2 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_60min_dwd_netatmo_upper_90_percent_data_considered_neighbor_2_.csv"
df3 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_60min_dwd_netatmo_upper_90_percent_data_considered_neighbor_3_.csv"
#df4 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_p5_freq_60min_dwd_netatmo_upper_98_percent_data_considered_neighbor_4_.csv"

#df0 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_indicator_correlations_DWD_DWD_ppt_\year_allyears_df_dwd_correlationsfreq_60min_dwd_netatmo_upper_95_percent_data_considered_neighbor_1_.csv"
#df1 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_indicator_correlations_DWD_DWD_ppt_\year_allyears_df_dwd_correlationsfreq_60min_dwd_netatmo_upper_95_percent_data_considered_neighbor_2_.csv"


# In[3]:


in_df0 = pd.read_csv(df0, index_col=0, sep=';')
in_df1 = pd.read_csv(df1, index_col=0, sep=';')
in_df2 = pd.read_csv(df2, index_col=0, sep=';')
in_df3 = pd.read_csv(df3, index_col=0, sep=';')
#in_df4 = pd.read_csv(df4, index_col=0, sep=';')

#in_df0.dropna(how='all', inplace=True)
#in_df1.dropna(how='all', inplace=True)
#in_df2.dropna(how='all', inplace=True)
#in_df3.dropna(how='all', inplace=True)

# In[4]:


dfs_all = pd.DataFrame(index=in_df0.index)


# In[5]:


# in_df0.head(5)
# In[6]:

dfs_all.loc[in_df0.index, 'stns0'] = in_df0.index
dfs_all.loc[in_df1.index, 'stns1'] = in_df1.index
dfs_all.loc[in_df2.index, 'stns2'] = in_df2.index
dfs_all.loc[in_df3.index, 'stns3'] = in_df3.index

dfs_all.loc[in_df0.index, 'distance0'] = in_df0.loc[:, 'Distance to neighbor']
dfs_all.loc[in_df1.index, 'distance1'] = in_df1.loc[:, 'Distance to neighbor']
dfs_all.loc[in_df2.index, 'distance2'] = in_df2.loc[:, 'Distance to neighbor']
dfs_all.loc[in_df3.index, 'distance3'] = in_df3.loc[:, 'Distance to neighbor']
#dfs_all.loc[in_df4.index, 'distance4'] = in_df4.loc[:, 'Distance to neighbor']

dfs_all.loc[in_df0.index, 'Indic_corr0'] = in_df0.loc[:,
                                                      'Bool_Spearman_Correlation']
dfs_all.loc[in_df1.index, 'Indic_corr1'] = in_df1.loc[:,
                                                      'Bool_Spearman_Correlation']
dfs_all.loc[in_df2.index, 'Indic_corr2'] = in_df2.loc[:,
                                                      'Bool_Spearman_Correlation']
dfs_all.loc[in_df3.index, 'Indic_corr3'] = in_df3.loc[:,
                                                      'Bool_Spearman_Correlation']
#dfs_all.loc[in_df4.index, 'Indic_corr4'] = in_df4.loc[:, 'Bool_Spearman_Correlation']

dfs_all.dropna(how='all', inplace=True)
# %%


# %%


s0 = dfs_all.stns0.values

xvals0 = dfs_all.distance0.values
y0 = dfs_all.Indic_corr0.values

s1, x1, y1 = dfs_all.stns1.values, dfs_all.distance1.values, dfs_all.Indic_corr1.values
s2, x2, y2 = dfs_all.stns2.values, dfs_all.distance2.values, dfs_all.Indic_corr2.values
s3, x3, y3 = dfs_all.stns3.values, dfs_all.distance3.values, dfs_all.Indic_corr3.values
# %%

#popt, pcov = curve_fit(func, xvals0, y0)
#print ("a = %s , b = %s, c = %s, d=%s" % (popt[0], popt[1], popt[2], popt[3] ))

stns = np.concatenate((s0, s1, s2, s3))
xs = np.concatenate((xvals0, x1, x2, x3))
ys = np.concatenate((y0, y1, y2, y3))

# remove nans
stns = stns[np.where(ys >= 0)]
xs = xs[np.where(ys >= 0)]
ys = ys[ys >= 0]


yfs, xvals_below_curves, yvals_below_curves, xa, ya, sb, sa = fit_curve_get_vals_below_curve(
    xs, ys, func, stns)

stns_to_keep = np.unique(sa)
stns_to_remove = np.unique(sb)

# TODO: NEED FIXING ASK PROF
s1_keep = np.intersect1d(sb, s2)
#%%
x0_abv = in_df0.loc[s1_keep, 'Distance to neighbor'].dropna().values
y0_abv = in_df0.loc[s1_keep, 'Bool_Spearman_Correlation'].dropna().values

x1_abv = in_df1.loc[s1_keep, 'Distance to neighbor']  # .dropna().values
y1_abv = in_df1.loc[s1_keep, 'Bool_Spearman_Correlation']  # .dropna().values

x2_abv = in_df2.loc[stns_to_keep, 'Distance to neighbor'].dropna().values
y2_abv = in_df2.loc[stns_to_keep, 'Bool_Spearman_Correlation'].dropna().values

x3_abv = in_df3.loc[stns_to_keep, 'Distance to neighbor'].dropna().values
y3_abv = in_df3.loc[stns_to_keep, 'Bool_Spearman_Correlation'].dropna().values

# In[11]: Additional filter on stations


def get_vals_abv_thr_based_on_2nd_arr(var1, var2, var2_thr):
    var1_abv = var1[np.where(var2 > var2_thr)]
    return var1_abv


# stns_to_keep = stns_to_keep[y0_abv > yfs[:y0_abv.shape[0]]]
# x0_abv = x0_abv[y0_abv > yfs]
# y0_abv = y0_abv[y0_abv > yfs]
#

plt.ioff()
plt.figure(figsize=(12, 8), dpi=100)


# plt.scatter(x0_abv, y0_abv, c='r', alpha=0.5,
#             marker='x', label='First Neighbor', s=28)
# plt.scatter(x1_abv, y1_abv, c='orange', alpha=0.5,
#             marker='.', label='Second Neighbor', s=28)
# plt.scatter(x2_abv, y2_abv, c='g', alpha=0.5,
#             marker='d', label='Third Neighbor', s=28)
# plt.scatter(x3_abv, y3_abv, c='b', alpha=0.5,
#             marker='1', label='Fourth Neighbor', s=28)

plt.scatter(xs, yfs, label="Polynomial function", alpha=.5, c='m', marker=',')
plt.scatter(xa, ya, label="Stns above Curve", alpha=.5, c='g', marker='d')


plt.scatter(xvals_below_curves, yvals_below_curves,
            label="Below Curve %d" % stns_to_remove.shape[0],
            alpha=.5, c='k', marker='X')

#plt.scatter(dfs_all.distance0, dfs_all.Indic_corr0, c='r', alpha=0.5, marker='x', label='First Neighbor', s=28)
#plt.scatter(dfs_all.distance1, dfs_all.Indic_corr1, c='orange', alpha=0.5, marker='.', label='Second Neighbor', s=28)
#plt.scatter(dfs_all.distance2, dfs_all.Indic_corr2, c='g', alpha=0.5, marker='d', label='Third Neighbor', s=28)
#plt.scatter(dfs_all.distance3, dfs_all.Indic_corr3, c='b', alpha=0.5, marker='1', label='Fourth Neighbor', s=28)
#plt.scatter(dfs_all.distance4, dfs_all.Indic_corr4, c='m', alpha=0.5, marker='2', label='Fifth Neighbor', s=28)

#plt.scatter([],[], label='Number of stations %d' % dfs_all.shape[0], marker='.')
plt.xlim([0, (dfs_all.distance3.max()) + 1000])
plt.ylim([-0.1, 1.1])
plt.xlabel('Distance (m)')
plt.ylabel('Indicator Spearman Correlation')
plt.legend(loc=0)
plt.grid(alpha=.25)
plt.tight_layout()
# plt.axis('equal')
plt.title('Keeping Netatmo %d stations: Indicator correlation with distance for upper 90percent of data values' %
          stns_to_keep.shape[0])
# plt.savefig(os.path.join(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature', 'Netatmo_98percent_indic_corr.png'),
#                frameon=True, papertype='a4',
#                bbox_inches='tight', pad_inches=.2)
plt.show()
