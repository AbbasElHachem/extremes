
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import fnmatch
from pathlib import Path

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears

from _00_additional_functions import (calculate_probab_ppt_below_thr,
                                      build_edf_fr_vals)
plt.ioff()
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'axes.labelsize': 12})


# In[2]:


# path_dwd_edf_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\edf_ppt_all_dwd_daily_.csv"
main_dir = Path(
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo')

#main_dir = Path(r'/home/IWS/hachem/Extremes')

min_orig_qnt_thr = 0.
_acc_ = ''

plot_not_filtered = False
plot_filtered = True


#==============================================================================
#
#==============================================================================
def list_all_full_path(ext, file_dir):
    """
    Purpose: To return full path of files in all dirs of a given folder with a
    -------  given extension in ascending order.

    Keyword arguments:
    ------------------
        ext (string) = Extension of the files to list
            e.g. '.txt', '.tif'.
        file_dir (string) = Full path of the folder in which the files
            reside.
    """
    new_list = []
    patt = '*' + ext
    for root, _, files in os.walk(file_dir):
        for elm in files:
            if fnmatch.fnmatch(elm, patt):
                full_path = os.path.join(root, elm)
                new_list.append(full_path)
    return(sorted(new_list))


#==============================================================================
# # ONLY PLOT NO FILTERTING
#==============================================================================
if plot_not_filtered:
    for temp_freq in ['60min', '360min', '720min', '1440min']:
        print(temp_freq)

        path_to_Quantiles_netatmo_no_flt___ = main_dir / (
            r'Qt_ok_ok_un_netatmo_no_flt___%s' % temp_freq)
        Quantiles_netatmo_no_flt___ = list_all_full_path(
            '.csv', path_to_Quantiles_netatmo_no_flt___)

        #########################################################
        path_dwd_edf_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\edf_ppt_all_dwd_%s_.csv" % temp_freq

        path_to_dwd_ppt = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\ppt_all_dwd_%s_.csv"
                           % temp_freq)

        df_dwd_edf = pd.read_csv(path_dwd_edf_data, sep=';', index_col=0,
                                 parse_dates=True, infer_datetime_format=True)

        # DWD ppt
        dwd_in_ppt_vals_df = pd.read_csv(
            path_to_dwd_ppt, sep=';', index_col=0, encoding='utf-8',
            parse_dates=True, infer_datetime_format=True)

        dwd_in_ppt_vals_df.dropna(how='all', axis=0, inplace=True)

        #########################################################
        df_improvements = pd.DataFrame(index=df_dwd_edf.columns,
                                       columns=['pearson_corr_dwd_', 'spearman_corr_dwd_',
                                                'pearson_corr_dwd_netatmo', 'spearman_corr_dwd_netatmo'])

        #########################################################
        path_interpolated_using_dwd_list = []
        path_interpolated_using_netatmo_dwd_list = []

        path_to_use = path_to_Quantiles_netatmo_no_flt___
        data_to_use = Quantiles_netatmo_no_flt___

        _interp_acc_ = str(r'%s' % (str(path_to_use).split('\\')[-1]))
        # for i in range(12):
        #   i = int(i)
        try:
            for df_file in data_to_use:

                if ('using_dwd_only_grp_') in df_file:
                    print(df_file)
                    path_interpolated_using_dwd_list.append(df_file)
                if ('using_dwd_netamo_grp_') in df_file:

                    print(df_file)
                    path_interpolated_using_netatmo_dwd_list.append(df_file)

        except Exception as msg:
            print(msg)
            continue
        #########################################################

        for (path_interpolated_using_dwd,
             path_interpolated_using_netatmo_dwd
             ) in zip(
                path_interpolated_using_dwd_list,
                path_interpolated_using_netatmo_dwd_list):

            df_dwd = pd.read_csv(path_interpolated_using_dwd, skiprows=[1],
                                 sep=';', index_col=0, parse_dates=True, infer_datetime_format=True)

            df_netatmo_dwd = pd.read_csv(path_interpolated_using_netatmo_dwd, skiprows=[1],
                                         sep=';', index_col=0, parse_dates=True, infer_datetime_format=True)

            df_compare = pd.DataFrame(index=df_dwd.index)

            cmn_interpolated_events = df_netatmo_dwd.index.intersection(
                df_dwd.index).intersection(df_dwd_edf.index)
            print('Total number of events is',
                  cmn_interpolated_events.shape[0])
            for stn_ in df_dwd.columns:
                # print(stn_)
                # print(df_dwd.columns.shape[0])
                #                 stn_ppt = dwd_in_ppt_vals_df.loc[:, stn_].dropna()
                #
                #                 p0 = calculate_probab_ppt_below_thr(stn_ppt.values,
                #                     0.1) / 2
                #                 print('P0', p0, p0 / 2)
                #
                #                 x0, y0 = build_edf_fr_vals(stn_df_no_nans.values)

                for event_date in cmn_interpolated_events:
                    # interpolated_quantile_netatmo = df_netatmo.loc[event_date, stn_]

                    interpolated_quantile_dwd = df_dwd.loc[event_date, stn_]

                    interpolated_quantile_netatmo_dwd = df_netatmo_dwd.loc[event_date, stn_]

                    original_quantile = df_dwd_edf.loc[event_date, stn_]

                    if original_quantile >= min_orig_qnt_thr:
                        # print(original_quantile)
                        df_compare.loc[event_date,
                                       'original_quantile'] = original_quantile
                        df_compare.loc[event_date,
                                       'interpolated_quantile_dwd'] = interpolated_quantile_dwd
                        # df_compare.loc[event_date,
                        #               'interpolated_quantile_netatmo'] = interpolated_quantile_netatmo
                        df_compare.loc[event_date,
                                       'interpolated_quantile_netatmo_dwd'] = interpolated_quantile_netatmo_dwd

                    else:
                        df_compare.loc[event_date,
                                       'original_quantile'] = np.nan
                        df_compare.loc[event_date,
                                       'interpolated_quantile_dwd'] = np.nan
                        # df_compare.loc[event_date,
                        #               'interpolated_quantile_netatmo'] = np.nan
                        df_compare.loc[event_date,
                                       'interpolated_quantile_netatmo_dwd'] = np.nan

                df_compare = df_compare[df_compare > 0]
                df_compare.dropna(how='any', inplace=True)

                values_x = df_compare['original_quantile'].values
                values_dwd = df_compare['interpolated_quantile_dwd'].values
                # values_netatmo =df_compare['interpolated_quantile_netatmo'].values
                values_netatmo_dwd = df_compare['interpolated_quantile_netatmo_dwd'].values

                # plot the stations in shapefile, look at the results of
                # agreements

                # calculate correlations (pearson and spearman) dwd-dwd
                corr_dwd = pears(values_x, values_dwd)[0]
                rho_dwd = spr(values_x, values_dwd)[0]

                # netatmo dwd
                corr_netatmo_dwd = pears(values_x, values_netatmo_dwd)[0]
                rho_netatmo_dwd = spr(values_x, values_netatmo_dwd)[0]

                df_improvements.loc[stn_, 'pearson_corr_dwd_'] = corr_dwd
                df_improvements.loc[stn_, 'spearman_corr_dwd_'] = rho_dwd
                df_improvements.loc[stn_,
                                    'pearson_corr_dwd_netatmo'] = corr_netatmo_dwd
                df_improvements.loc[stn_,
                                    'spearman_corr_dwd_netatmo'] = rho_netatmo_dwd

            df_improvements.dropna(how='all', inplace=True)

        #######################################################################
        dwd_ids_no_improvement_pearson = pd.DataFrame(
            data=df_improvements.index[np.where(
                df_improvements.pearson_corr_dwd_netatmo.values <
                df_improvements.pearson_corr_dwd_.values)],
            columns=['stns_no_impv_pearson_corr'])

        dwd_ids_no_improvement_spr = pd.DataFrame(
            data=df_improvements.index[np.where(
                df_improvements.spearman_corr_dwd_netatmo.values <
                df_improvements.spearman_corr_dwd_.values)],
            columns=['stns_no_impv_spearman_corr'])

        #########################################################
        stations_with_improvements = sum(i >= j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo.values,
            df_improvements.pearson_corr_dwd_.values))

        stations_without_improvements = sum(i < j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo.values,
            df_improvements.pearson_corr_dwd_.values))

        percent_of_improvment = 100 * (stations_with_improvements /
                                       df_improvements.pearson_corr_dwd_netatmo.shape[0])

        #########################################################

        stations_with_improvements_spr = sum(i >= j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo.values,
            df_improvements.spearman_corr_dwd_.values))

        percent_of_improvment_spr = 100 * (stations_with_improvements_spr /
                                           df_improvements.spearman_corr_dwd_netatmo.shape[0])

        stations_with_improvements_spr = sum(i >= j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo.values,
            df_improvements.spearman_corr_dwd_.values))

        percent_of_improvment_spr = 100 * (stations_with_improvements_spr /
                                           df_improvements.spearman_corr_dwd_netatmo.shape[0])
        #########################################################
        #########################################################

        mean_pearson_correlation_dwd_only = df_improvements.pearson_corr_dwd_.mean()
        mean_pearson_correlation_dwd_netatmo = df_improvements.pearson_corr_dwd_netatmo.mean()

        #########################################################
        mean_spr_correlation_dwd_only = df_improvements.spearman_corr_dwd_.mean()
        mean_spr_correlation_dwd_netatmo = df_improvements.spearman_corr_dwd_netatmo.mean()

        ########################
        plt.ioff()

        fig = plt.figure(figsize=(24, 12), dpi=150)

        ax = fig.add_subplot(111)

        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation %.2f' % np.round(
                    mean_pearson_correlation_dwd_only, 2))
        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_netatmo,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %.2f' % np.round(
                    mean_pearson_correlation_dwd_netatmo, 2))

        ax.set_title('Pearson Correlation Interpolated from DWD or DWD-Netatmo %s \n '
                     'Quantiles of %s Extreme Events %s \n Stations with Improvemnts %d / %d, Percentage %0.0f'
                     % (_acc_, temp_freq, _interp_acc_, stations_with_improvements,
                        df_improvements.pearson_corr_dwd_netatmo.shape[0], percent_of_improvment))
        ax.grid(alpha=0.25)
        plt.setp(ax.get_xticklabels(), rotation=45)
        ax.set_yticks(np.arange(0, 1.05, .10))
        ax.set_xlabel('DWD Stations')
        ax.legend(loc='lower right')
        ax.set_ylabel('Pearson Correlation')

        plt.savefig((path_to_use / (
                    r'pears_corr_%s_events_dwd_%s.png' % (temp_freq, _interp_acc_))),
                    frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

        ########################

        plt.ioff()
        fig = plt.figure(figsize=(24, 12), dpi=150)

        ax = fig.add_subplot(111)

        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation %.2f' % np.round(
                       mean_spr_correlation_dwd_only, 2))
        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_netatmo,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %.2f' % np.round(
                    mean_spr_correlation_dwd_netatmo, 2))

        ax.set_title('Spearman Correlation Interpolated from DWD or DWD-Netatmo %s\n '
                     'Quantiles of %s Extreme Events %s \n Stations with Improvemnts %d / %d, Percentage %0.0f'
                     % (_acc_, temp_freq, _interp_acc_, stations_with_improvements_spr,
                        df_improvements.spearman_corr_dwd_netatmo.shape[0], percent_of_improvment_spr))

        plt.setp(ax.get_xticklabels(), rotation=45)
        ax.grid(alpha=0.25)
        ax.set_yticks(np.arange(0, 1.05, .10))
        ax.set_xlabel('DWD Stations')
        ax.legend(loc='lower right')
        ax.set_ylabel('Spearman Correlation')

        plt.savefig((path_to_use / (
                    r'spr_corr_%s_events_dwd_%s.png' % (temp_freq, _interp_acc_))),
                    frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

        ########################

#==============================================================================
# # PLOT FILTERING
#==============================================================================
if plot_filtered:
    for temp_freq in ['60min', '360min', '720min', '1440min']:
        # '60min', '360min',  '1440min'
        print(temp_freq)

        path_to_Qt_ok_un_first_flt__temp_flt_1st_ = main_dir / (
            r'Qt_ok_ok_un_2_first_flt__temp_flt__1st_%s' % temp_freq)
        Qt_ok_un_first_flt__temp_flt_1st_ = list_all_full_path(
            '.csv', path_to_Qt_ok_un_first_flt__temp_flt_1st_)

        path_to_Qt_ok_un_first_flt__temp_flt_comb_ = main_dir / (
            r'Qt_ok_ok_un_2_first_flt__temp_flt__comb_%s' % temp_freq)
        Qt_ok_un_first_flt__temp_flt_comb_ = list_all_full_path(
            '.csv', path_to_Qt_ok_un_first_flt__temp_flt_comb_)

        path_to_Qt_ok_un_first_flt_1st_ = main_dir / (
            r'Qt_ok_ok_un_2_first_flt__1st_%s' % temp_freq)
        Qt_ok_un_first_flt_1st_ = list_all_full_path(
            '.csv', path_to_Qt_ok_un_first_flt_1st_)

        path_to_Qt_ok_un_first_flt_comb_ = main_dir / (
            r'Qt_ok_ok_un_2_first_flt__comb_%s' % temp_freq)
        Qt_ok_un_first_flt_comb_ = list_all_full_path(
            '.csv', path_to_Qt_ok_un_first_flt_comb_)

        #########################################################
        path_dwd_edf_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\edf_ppt_all_dwd_%s_.csv" % temp_freq

        df_dwd_edf = pd.read_csv(path_dwd_edf_data, sep=';', index_col=0,
                                 parse_dates=True, infer_datetime_format=True)

        #########################################################
        df_improvements = pd.DataFrame(
            index=df_dwd_edf.columns,
            columns=['pearson_corr_dwd_',
                     'spearman_corr_dwd_',
                     'pearson_corr_dwd_netatmo',
                     'spearman_corr_dwd_netatmo',
                     'pearson_corr_dwd_netatmo_unc',
                     'spearman_corr_dwd_netatmo_unc'])

        #########################################################
        path_interpolated_using_dwd_list = []
        path_interpolated_using_netatmo_dwd_list = []
        path_interpolated_using_netatmo_dwd_list_un = []

        path_to_use = path_to_Qt_ok_un_first_flt__temp_flt_1st_
        data_to_use = Qt_ok_un_first_flt__temp_flt_1st_

        _interp_acc_ = str(r'%s' % (str(path_to_use).split('\\')[-1]))
        # for i in range(12):
        #   i = int(i)
        try:
            for df_file in data_to_use:

                if ('using_dwd_only_grp_') in df_file:
                    print(df_file)
                    path_interpolated_using_dwd_list.append(df_file)
                if ('using_dwd_netamo_grp_') in df_file:

                    print(df_file)
                    path_interpolated_using_netatmo_dwd_list.append(df_file)

                if ('interpolated_quantiles_un_dwd_') in df_file:
                    print(df_file)
                    path_interpolated_using_netatmo_dwd_list_un.append(df_file)

        except Exception as msg:
            print(msg)
            continue
        #########################################################

        for (path_interpolated_using_dwd,
             path_interpolated_using_netatmo_dwd,
             path_interpolated_using_netatmo_dwd_unc) in zip(
                path_interpolated_using_dwd_list,
                path_interpolated_using_netatmo_dwd_list,
                path_interpolated_using_netatmo_dwd_list_un):

            df_dwd = pd.read_csv(path_interpolated_using_dwd, skiprows=[1],
                                 sep=';', index_col=0,
                                 parse_dates=True, infer_datetime_format=True)

            df_netatmo_dwd = pd.read_csv(path_interpolated_using_netatmo_dwd, skiprows=[1],
                                         sep=';', index_col=0,
                                         parse_dates=True, infer_datetime_format=True)

            df_netatmo_dwd_unc = pd.read_csv(path_interpolated_using_netatmo_dwd_unc, skiprows=[1],
                                             sep=';', index_col=0,
                                             parse_dates=True, infer_datetime_format=True)

            df_compare = pd.DataFrame(index=df_dwd.index)

            cmn_interpolated_events = df_netatmo_dwd.index.intersection(
                df_dwd.index).intersection(df_dwd_edf.index)

            for stn_ in df_dwd.columns:
                # print(stn_)
                # print(df_dwd.columns.shape[0])
                for event_date in cmn_interpolated_events:
                    # interpolated_quantile_netatmo = df_netatmo.loc[event_date, stn_]

                    interpolated_quantile_dwd = df_dwd.loc[event_date, stn_]

                    interpolated_quantile_netatmo_dwd = df_netatmo_dwd.loc[event_date, stn_]

                    interpolated_quantile_netatmo_dwd_unc = df_netatmo_dwd_unc.loc[event_date, stn_]

                    original_quantile = df_dwd_edf.loc[event_date, stn_]

                    if original_quantile >= min_orig_qnt_thr:
                        # print(original_quantile)
                        df_compare.loc[event_date,
                                       'original_quantile'] = original_quantile
                        df_compare.loc[event_date,
                                       'interpolated_quantile_dwd'] = interpolated_quantile_dwd
                        # df_compare.loc[event_date,
                        #               'interpolated_quantile_netatmo'] = interpolated_quantile_netatmo
                        df_compare.loc[event_date,
                                       'interpolated_quantile_netatmo_dwd'] = interpolated_quantile_netatmo_dwd
                        df_compare.loc[event_date,
                                       'interpolated_quantile_netatmo_dwd_unc'] = interpolated_quantile_netatmo_dwd_unc

                    else:
                        df_compare.loc[event_date,
                                       'original_quantile'] = np.nan
                        df_compare.loc[event_date,
                                       'interpolated_quantile_dwd'] = np.nan
                        # df_compare.loc[event_date,
                        #               'interpolated_quantile_netatmo'] = np.nan
                        df_compare.loc[event_date,
                                       'interpolated_quantile_netatmo_dwd'] = np.nan

                        df_compare.loc[event_date,
                                       'interpolated_quantile_netatmo_dwd_unc'] = np.nan

                df_compare = df_compare[df_compare > 0]
                df_compare.dropna(how='any', inplace=True)
                print('Total number of events is',
                      df_compare.shape[0])
                values_x = df_compare['original_quantile'].values
                values_dwd = df_compare['interpolated_quantile_dwd'].values
                # values_netatmo =df_compare['interpolated_quantile_netatmo'].values
                values_netatmo_dwd = df_compare['interpolated_quantile_netatmo_dwd'].values
                values_netatmo_dwd_unc = df_compare['interpolated_quantile_netatmo_dwd_unc'].values
                # plot the stations in shapefile, look at the results of
                # agreements

                # calculate correlations (pearson and spearman) dwd-dwd
                corr_dwd = pears(values_x, values_dwd)[0]
                rho_dwd = spr(values_x, values_dwd)[0]

                # netatmo dwd
                corr_netatmo_dwd = pears(values_x, values_netatmo_dwd)[0]
                rho_netatmo_dwd = spr(values_x, values_netatmo_dwd)[0]

                # netatmo dwd unc
                corr_netatmo_dwd_unc = pears(
                    values_x, values_netatmo_dwd_unc)[0]
                rho_netatmo_dwd_unc = spr(values_x, values_netatmo_dwd_unc)[0]

                df_improvements.loc[stn_, 'pearson_corr_dwd_'] = corr_dwd
                df_improvements.loc[stn_, 'spearman_corr_dwd_'] = rho_dwd
                df_improvements.loc[stn_,
                                    'pearson_corr_dwd_netatmo'] = corr_netatmo_dwd
                df_improvements.loc[stn_,
                                    'spearman_corr_dwd_netatmo'] = rho_netatmo_dwd

                df_improvements.loc[stn_,
                                    'pearson_corr_dwd_netatmo_unc'] = corr_netatmo_dwd_unc
                df_improvements.loc[stn_,
                                    'spearman_corr_dwd_netatmo_unc'] = rho_netatmo_dwd_unc

            df_improvements.dropna(how='all', inplace=True)

        #######################################################################

        dwd_ids_no_improvement_pearson = pd.DataFrame(
            data=df_improvements.index[np.where(
                df_improvements.pearson_corr_dwd_netatmo.values <
                df_improvements.pearson_corr_dwd_.values)],
            columns=['stns_no_impv_pearson_corr'])
        dwd_ids_no_improvement_spr = pd.DataFrame(
            data=df_improvements.index[np.where(
                df_improvements.spearman_corr_dwd_netatmo.values <
                df_improvements.spearman_corr_dwd_.values)],
            columns=['stns_no_impv_spearman_corr'])
        # netatmo dwd unc
        dwd_ids_no_improvement_pearson_unc = pd.DataFrame(
            data=df_improvements.index[np.where(
                df_improvements.pearson_corr_dwd_netatmo_unc.values <
                df_improvements.pearson_corr_dwd_.values)],
            columns=['stns_no_impv_pearson_corr_unc'])

        dwd_ids_no_improvement_spr_unc = pd.DataFrame(
            data=df_improvements.index[np.where(
                df_improvements.spearman_corr_dwd_netatmo_unc.values <
                df_improvements.spearman_corr_dwd_.values)],
            columns=['stns_no_impv_spearman_corr_unc'])

        # dwd_ids_no_improvement_pearson.to_csv((path_to_use / (
        # r'pears_corr_%s_events_dwd_%s_no_improv.csv' % (temp_freq,
        # _interp_acc_))),sep=';')

        # dwd_ids_no_improvement_spr.to_csv((path_to_use / (
        #            r'spr_corr_%s_events_dwd_%s_no_improv.csv' % (temp_freq, _interp_acc_))),
        #                                  sep=';')

        # dwd_ids_no_improvement_pearson.to_csv((path_to_use / (
        # r'pears_corr_%s_events_dwd_%s_no_improv.csv' % (temp_freq,
        # _interp_acc_))),sep=';')

        # dwd_ids_no_improvement_spr.to_csv((path_to_use / (
        #            r'spr_corr_%s_events_dwd_%s_no_improv.csv' % (temp_freq, _interp_acc_))),
        #                                  sep=';')
        #########################################################
        stations_with_improvements = sum(i >= j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo.values,
            df_improvements.pearson_corr_dwd_.values))

        stations_without_improvements = sum(i < j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo.values,
            df_improvements.pearson_corr_dwd_.values))

        percent_of_improvment = 100 * (stations_with_improvements /
                                       df_improvements.pearson_corr_dwd_netatmo.shape[0])

        #########################################################

        stations_with_improvements_spr = sum(i >= j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo.values,
            df_improvements.spearman_corr_dwd_.values))

        percent_of_improvment_spr = 100 * (stations_with_improvements_spr /
                                           df_improvements.spearman_corr_dwd_netatmo.shape[0])

        stations_with_improvements_spr = sum(i >= j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo.values,
            df_improvements.spearman_corr_dwd_.values))

        percent_of_improvment_spr = 100 * (stations_with_improvements_spr /
                                           df_improvements.spearman_corr_dwd_netatmo.shape[0])
        #########################################################

        # netatmo dwd unc
        stations_with_improvements_unc = sum(i >= j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo_unc.values,
            df_improvements.pearson_corr_dwd_.values))

        stations_without_improvements_unc = sum(i < j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo_unc.values,
            df_improvements.pearson_corr_dwd_.values))

        percent_of_improvment_unc = 100 * (stations_with_improvements_unc /
                                           df_improvements.pearson_corr_dwd_netatmo_unc.shape[0])
        #########################################################
        stations_with_improvements_spr_unc = sum(i >= j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo_unc.values,
            df_improvements.spearman_corr_dwd_.values))

        percent_of_improvment_spr_unc = 100 * (stations_with_improvements_spr_unc /
                                               df_improvements.spearman_corr_dwd_netatmo.shape[0])

        stations_with_improvements_spr_unc = sum(i >= j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo_unc.values,
            df_improvements.spearman_corr_dwd_.values))

        percent_of_improvment_spr_unc = 100 * (stations_with_improvements_spr_unc /
                                               df_improvements.spearman_corr_dwd_netatmo.shape[0])
        #########################################################

        mean_pearson_correlation_dwd_only = df_improvements.pearson_corr_dwd_.mean()
        mean_pearson_correlation_dwd_netatmo = df_improvements.pearson_corr_dwd_netatmo.mean()
        mean_pearson_correlation_dwd_netatmo_unc = df_improvements.pearson_corr_dwd_netatmo_unc.mean()

        #########################################################
        mean_spr_correlation_dwd_only = df_improvements.spearman_corr_dwd_.mean()
        mean_spr_correlation_dwd_netatmo = df_improvements.spearman_corr_dwd_netatmo.mean()
        mean_spr_correlation_dwd_netatmo_unc = df_improvements.spearman_corr_dwd_netatmo_unc.mean()

        ########################
        plt.ioff()

        fig = plt.figure(figsize=(24, 12), dpi=150)

        ax = fig.add_subplot(111)

        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation %.2f' % np.round(
                    mean_pearson_correlation_dwd_only, 2))
        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_netatmo,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %.2f' % np.round(
                    mean_pearson_correlation_dwd_netatmo, 2))

        ax.set_title('Pearson Correlation Interpolated from DWD or DWD-Netatmo %s \n '
                     'Quantiles of %s Extreme Events %s \n Stations with Improvemnts %d / %d, Percentage %0.0f'
                     % (_acc_, temp_freq, _interp_acc_, stations_with_improvements,
                        df_improvements.pearson_corr_dwd_netatmo.shape[0], percent_of_improvment))
        ax.grid(alpha=0.25)
        plt.setp(ax.get_xticklabels(), rotation=45)
        ax.set_yticks(np.arange(0, 1.05, .10))
        ax.set_xlabel('DWD Stations')
        ax.legend(loc='lower right')
        ax.set_ylabel('Pearson Correlation')

        plt.savefig((path_to_use / (
                    r'pears_corr_%s_events_dwd_%s.png' % (temp_freq, _interp_acc_))),
                    frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

        ########################

        plt.ioff()
        fig = plt.figure(figsize=(24, 12), dpi=150)

        ax = fig.add_subplot(111)

        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation %.2f' % np.round(
                       mean_spr_correlation_dwd_only, 2))
        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_netatmo,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %.2f' % np.round(
                    mean_spr_correlation_dwd_netatmo, 2))

        ax.set_title('Spearman Correlation Interpolated from DWD or DWD-Netatmo %s\n '
                     'Quantiles of %s Extreme Events %s \n Stations with Improvemnts %d / %d, Percentage %0.0f'
                     % (_acc_, temp_freq, _interp_acc_, stations_with_improvements_spr,
                        df_improvements.spearman_corr_dwd_netatmo.shape[0], percent_of_improvment_spr))

        plt.setp(ax.get_xticklabels(), rotation=45)
        ax.grid(alpha=0.25)
        ax.set_yticks(np.arange(0, 1.05, .10))
        ax.set_xlabel('DWD Stations')
        ax.legend(loc='lower right')
        ax.set_ylabel('Spearman Correlation')

        plt.savefig((path_to_use / (
                    r'spr_corr_%s_events_dwd_%s.png' % (temp_freq, _interp_acc_))),
                    frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

        ########################
        fig = plt.figure(figsize=(24, 12), dpi=150)

        ax = fig.add_subplot(111)

        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation %.2f' % np.round(
                    mean_pearson_correlation_dwd_only, 2))
        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_netatmo_unc,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %.2f' % np.round(
                    mean_pearson_correlation_dwd_netatmo_unc, 2))

        ax.set_title('Pearson Correlation Interpolated from DWD or DWD-Netatmo with Unc%s \n '
                     'Quantiles of %s Extreme Events %s \n Stations with Improvemnts %d / %d, Percentage %0.0f'
                     % (_acc_, temp_freq, _interp_acc_, stations_with_improvements_unc,
                        df_improvements.pearson_corr_dwd_netatmo.shape[0], percent_of_improvment_unc))
        ax.grid(alpha=0.25)
        plt.setp(ax.get_xticklabels(), rotation=45)
        ax.set_yticks(np.arange(0, 1.05, .10))
        ax.set_xlabel('DWD Stations')
        ax.legend(loc='lower right')
        ax.set_ylabel('Pearson Correlation')

        plt.savefig((path_to_use / (
                    r'pears_corr_%s_events_dwd_%s_unc.png' % (temp_freq, _interp_acc_))),
                    frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

        ########################

        fig = plt.figure(figsize=(24, 12), dpi=150)

        ax = fig.add_subplot(111)

        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation %.2f' % np.round(
                       mean_spr_correlation_dwd_only, 2))
        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_netatmo_unc,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation Unc %.2f' % np.round(
                    mean_spr_correlation_dwd_netatmo, 2))

        ax.set_title('Spearman Correlation Interpolated from DWD or DWD-Netatmo with Unc %s\n '
                     'Quantiles of %s Extreme Events %s \n Stations with Improvemnts %d / %d, Percentage %0.0f'
                     % (_acc_, temp_freq, _interp_acc_, stations_with_improvements_spr_unc,
                        df_improvements.spearman_corr_dwd_netatmo_unc.shape[0], percent_of_improvment_spr_unc))

        plt.setp(ax.get_xticklabels(), rotation=45)
        ax.grid(alpha=0.25)
        ax.set_yticks(np.arange(0, 1.05, .10))
        ax.set_xlabel('DWD Stations')
        ax.legend(loc='lower right')
        ax.set_ylabel('Spearman Correlation')

        plt.savefig((path_to_use / (
                    r'spr_corr_%s_events_dwd_%s_unc.png' % (temp_freq, _interp_acc_))),
                    frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()
