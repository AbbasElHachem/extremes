
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears

from _00_additional_functions import (build_edf_fr_vals, find_nearest)
plt.ioff()
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'axes.labelsize': 12})


#==============================================================================
# # In[2]:
#==============================================================================

min_orig_qnt_thr = 0.

temp_freq = '1440min'

kriging_with_uncertainty = False

if kriging_with_uncertainty:
    _acc_ = 'kriging_with_uncert_'
if not kriging_with_uncertainty:
    _acc_ = ''

add_to_path = ''

if temp_freq == '360min':
    add_to_path = '_Temporal_filter_used_'

if temp_freq == '720min':
    add_to_path = '_Temporal_filter_used__Temporal_filter_used_'

if temp_freq == '1440min':
    add_to_path = '_Temporal_filter_used__Temporal_filter_used__Temporal_filter_used_'


#==============================================================================
#
#==============================================================================


path_dwd_edf_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\edf_ppt_all_dwd_%s_.csv" % temp_freq

df_dwd_edf = pd.read_csv(path_dwd_edf_data, sep=';', index_col=0,
                         parse_dates=True, infer_datetime_format=True)

path_dwd_ppt_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\ppt_all_dwd_%s_.csv" % temp_freq

df_dwd_ppt = pd.read_csv(path_dwd_ppt_data, sep=';', index_col=0,
                         parse_dates=True, infer_datetime_format=True)

df_improvements = pd.DataFrame(index=df_dwd_edf.columns,
                               columns=['pearson_corr_dwd_',
                                        'spearman_corr_dwd_',
                                        'pearson_corr_dwd_netatmo',
                                        'spearman_corr_dwd_netatmo'])

for i in range(12):
    i = int(i)
    # i=0
    try:
        # path_interpolated_using_netatmo = (
        #    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo"
        #    r"\not_filtered_interpolated_quantiles_%s_data_basedon_qunatiles_%s_season_using_netatmo_only_grp_%d_.csv" %(temp_freq,path_acc, i))
        # _Temporal_filter_used_

        path_interpolated_using_dwd = (
            r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo"
            r"\%sinterpolated_quantiles_dwd_%s_data_Quantiles_Temporal_filter_used_%s_using_dwd_only_grp_%d_.csv" % (
                _acc_, temp_freq, add_to_path, i))

        path_interpolated_using_netatmo_dwd = (
            r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo"
            r"\%sinterpolated_quantiles_dwd_%s_data_Quantiles_Temporal_filter_used_%s_using_dwd_netamo_grp_%d_.csv" % (
                _acc_, temp_freq, add_to_path, i))

        # df_netatmo = pd.read_csv(path_interpolated_using_netatmo,
        # sep=';', index_col=0, parse_dates=True, infer_datetime_format=True)

        df_dwd = pd.read_csv(path_interpolated_using_dwd,
                             sep=';', index_col=0, parse_dates=True,
                             infer_datetime_format=True)

        df_netatmo_dwd = pd.read_csv(path_interpolated_using_netatmo_dwd,
                                     sep=';', index_col=0, parse_dates=True,
                                     infer_datetime_format=True)

        df_compare = pd.DataFrame(index=df_dwd.index)

        cmn_interpolated_events = df_netatmo_dwd.index.intersection(
            df_dwd.index)
        print('Total number of events is', cmn_interpolated_events.shape[0])
        for stn_ in df_dwd.columns:
            print(stn_)
            #stn_ = df_dwd.columns[0]
            for event_date in cmn_interpolated_events:
                # interpolated_quantile_netatmo = df_netatmo.loc[event_date, stn_]
                #event_date = cmn_interpolated_events[0]
                interpolated_quantile_dwd = df_dwd.loc[event_date, stn_]

                interpolated_quantile_netatmo_dwd = df_netatmo_dwd.loc[event_date, stn_]

                # transform back to precipitation
                original_quantile = df_dwd_edf.loc[event_date, stn_]

                ppt_stn, edf_stn_dist = build_edf_fr_vals(
                    df_dwd_ppt.loc[:, stn_].dropna().values.ravel())

                # find ppt corresponding to quantiles
                ppt_orig = df_dwd_ppt.loc[event_date, stn_]
                ppt_orig_fr_edf = find_nearest(ppt_stn, ppt_orig)

                edf_orig_fr_ppt = edf_stn_dist[
                    np.where(ppt_stn == ppt_orig_fr_edf)][0]

                # interpolated from DWD
                ppt_interp_fr_dwd = find_nearest(
                    ppt_stn, interpolated_quantile_dwd)

                # interpolated from DWD-Netatmo
                ppt_interp_fr_dwd_netatmo = find_nearest(
                    ppt_stn, interpolated_quantile_netatmo_dwd)

                if ppt_orig >= min_orig_qnt_thr:
                    # print(original_quantile)
                    df_compare.loc[event_date,
                                   'original_ppt'] = ppt_orig
                    df_compare.loc[event_date,
                                   'interpolated_ppt_dwd'] = ppt_interp_fr_dwd
                    # df_compare.loc[event_date,
                    #               'interpolated_quantile_netatmo'] = interpolated_quantile_netatmo
                    df_compare.loc[event_date,
                                   'interpolated_ppt_netatmo_dwd'] = ppt_interp_fr_dwd_netatmo

                else:
                    df_compare.loc[event_date,
                                   'original_ppt'] = np.nan
                    df_compare.loc[event_date,
                                   'interpolated_ppt_dwd'] = np.nan
                    # df_compare.loc[event_date,
                    #               'interpolated_quantile_netatmo'] = np.nan
                    df_compare.loc[event_date,
                                   'interpolated_ppt_netatmo_dwd'] = np.nan

            df_compare = df_compare[df_compare > 0]
            df_compare.dropna(how='any', inplace=True)

            values_x = df_compare['original_ppt'].values
            values_dwd = df_compare['interpolated_ppt_dwd'].values
            # values_netatmo =df_compare['interpolated_quantile_netatmo'].values
            values_netatmo_dwd = df_compare['interpolated_ppt_netatmo_dwd'].values

            # plot the stations in shapefile, look at the results of agreements

            # calculate correlations (pearson and spearman)
            corr_dwd = pears(values_x, values_dwd)[0]
            rho_dwd = spr(values_x, values_dwd)[0]

            #corr_netatmo = pears(values_x, values_netatmo)[0]
            #rho_netatmo = spr(values_x, values_netatmo)[0]

            corr_netatmo_dwd = pears(values_x, values_netatmo_dwd)[0]
            rho_netatmo_dwd = spr(values_x, values_netatmo_dwd)[0]

            df_improvements.loc[stn_, 'pearson_corr_dwd_'] = corr_dwd
            df_improvements.loc[stn_, 'spearman_corr_dwd_'] = rho_dwd
            df_improvements.loc[stn_,
                                'pearson_corr_dwd_netatmo'] = corr_netatmo_dwd
            df_improvements.loc[stn_,
                                'spearman_corr_dwd_netatmo'] = rho_netatmo_dwd

    except Exception as msg:
        print(msg)
        continue

df_improvements.dropna(how='all', inplace=True)


# In[41]:


df_improvements.head()


# In[42]:


# sum(i > j for (i, j) in zip(df_improvements.pearson_corr_dwd_netatmo.values,
#                            df_improvements.pearson_corr_dwd_.values))

dwd_ids_no_improvement_pearson = pd.DataFrame(data=df_improvements.index[np.where(
    df_improvements.pearson_corr_dwd_netatmo.values < df_improvements.pearson_corr_dwd_.values)],
    columns=['stns_no_impv_pearson_corr'])
dwd_ids_no_improvement_spr = pd.DataFrame(data=df_improvements.index[np.where(
    df_improvements.spearman_corr_dwd_netatmo.values < df_improvements.spearman_corr_dwd_.values)],
    columns=['stns_no_impv_spearman_corr'])

# dwd_ids_no_improvement_pearson.to_csv(
#    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\stns_no_impv_pearson_corr_daily_events_all.csv',sep=';')

# dwd_ids_no_improvement_spr.to_csv(
#    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\stns_no_impv_spearman_corr_daily_events_all.csv',sep=';')

df_improvements


# In[43]:


stations_with_improvements = sum(i >= j for (i, j) in zip(df_improvements.pearson_corr_dwd_netatmo.values,
                                                          df_improvements.pearson_corr_dwd_.values))

stations_without_improvements = sum(i < j for (i, j) in zip(df_improvements.pearson_corr_dwd_netatmo.values,
                                                            df_improvements.pearson_corr_dwd_.values))

percent_of_improvment = 100 * (stations_with_improvements /
                               df_improvements.pearson_corr_dwd_netatmo.shape[0])
#########################################################
plt.ioff()
fig = plt.figure(figsize=(24, 12), dpi=150)

ax = fig.add_subplot(111)

# ax.scatter(df_improvements.index,
#           df_improvements.pearson_corr_dwd_,
#           alpha=.8,
#           c='r',  # colors_arr,
#           s=15,
#           marker='d',
# cmap=plt.get_cmap('viridis'),
#          label='DWD Interpolation')

# ax.scatter(df_improvements.index,
#           df_improvements.pearson_corr_dwd_netatmo,
#           alpha=.8,
#           c='b',  # colors_arr,
#           s=15,
#           marker='*',
# cmap=plt.get_cmap('viridis'),
#           label='DWD-Netatmo Interpolation')

ax.plot(df_improvements.index,
        df_improvements.pearson_corr_dwd_,
        alpha=.8,
        c='b',  # colors_arr,
        marker='d',
        label='DWD Interpolation')
ax.plot(df_improvements.index,
        df_improvements.pearson_corr_dwd_netatmo,
        alpha=.8,
        c='r',  # colors_arr,
        marker='*',
        label='DWD-Netatmo Interpolation')

ax.set_title('Pearson Correlation Interpolated from DWD or DWD-Netatmo %s \n '
             'Rainfall of %s Extreme Events \n Stations with Improvemnts %d / %d, Percentage %0.0f'
             % (_acc_, temp_freq, stations_with_improvements,
                df_improvements.pearson_corr_dwd_netatmo.shape[0], percent_of_improvment))
ax.grid(alpha=0.25)
plt.setp(ax.get_xticklabels(), rotation=45)
ax.set_yticks(np.arange(0, 1.05, .10))
ax.set_xlabel('DWD Stations')
ax.legend(loc='upper right')
ax.set_ylabel('Pearson Correlation')

plt.savefig((r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo'
             r'\%soberserved_vs_interpolated_pears_correlationrainfall_extreme_%s_events_dwd_all_.png'
             % (_acc_, temp_freq)),
            frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
plt.close()


# In[44]:


stations_with_improvements = sum(i >= j for (i, j) in zip(df_improvements.spearman_corr_dwd_netatmo.values,
                                                          df_improvements.spearman_corr_dwd_.values))

percent_of_improvment = 100 * (stations_with_improvements /
                               df_improvements.spearman_corr_dwd_netatmo.shape[0])
#########################################################
plt.ioff()
fig = plt.figure(figsize=(24, 12), dpi=150)

ax = fig.add_subplot(111)

# ax.scatter(df_improvements.index,
#           df_improvements.spearman_corr_dwd_,
#           alpha=.8,
#           c='r',  # colors_arr,
#           s=15,
#           marker='d',
#           # cmap=plt.get_cmap('viridis'),
#           label='DWD Interpolation')

# ax.scatter(df_improvements.index,
#           df_improvements.spearman_corr_dwd_netatmo,
#           alpha=.8,
#           c='b',  # colors_arr,
#           s=15,
#           marker='*',
# cmap=plt.get_cmap('viridis'),
#           label='DWD-Netatmo Interpolation')

ax.plot(df_improvements.index,
        df_improvements.spearman_corr_dwd_,
        alpha=.8,
        c='b',  # colors_arr,
        marker='d',
        label='DWD Interpolation')
ax.plot(df_improvements.index,
        df_improvements.spearman_corr_dwd_netatmo,
        alpha=.8,
        c='r',  # colors_arr,
        marker='*',
        label='DWD-Netatmo Interpolation')

ax.set_title('Spearman Correlation Interpolated from DWD or DWD-Netatmo %s\n '
             'Rainfall of %s Extreme Events \n Stations with Improvemnts %d / %d, Percentage %0.0f'
             % (_acc_, temp_freq, stations_with_improvements,
                df_improvements.spearman_corr_dwd_netatmo.shape[0], percent_of_improvment))

plt.setp(ax.get_xticklabels(), rotation=45)
ax.grid(alpha=0.25)
ax.set_yticks(np.arange(0, 1.05, .10))
ax.set_xlabel('DWD Stations')
ax.legend(loc='upper right')
ax.set_ylabel('Spearman Correlation')

plt.savefig((r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo'
             r'\%soberserved_vs_interpolated_spearman_correlationrainfall_extreme_%s_events_dwd_all_.png' % (
                 _acc_, temp_freq)),
            frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
plt.close()


# In[43]:
