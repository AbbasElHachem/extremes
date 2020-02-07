
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import fnmatch
from pathlib import Path

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears
from scipy.stats import kendalltau as kdl
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

plt.ioff()
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'axes.labelsize': 16})


# In[12]:


main_dir = Path(
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo')


plot_filtered = True
plot_not_filtered = False


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


# In[14]:


main_dir


# In[26]:


if plot_filtered:
    # , '180min', '360min', '720min', '1440min']:
    for temp_freq in ['60min', '180min', '360min', '720min', '1440min']:
        print(temp_freq)

        path_to_Qt_ok_un_first_flt__temp_flt_1st_ = main_dir / (
            r'Final_results4/Ppt_ok_ok_un_new3_first_flt__temp_flt__1st_%s' % temp_freq)
        Qt_ok_un_first_flt__temp_flt_1st_ = list_all_full_path(
            '.csv', path_to_Qt_ok_un_first_flt__temp_flt_1st_)

        path_to_Qt_ok_un_first_flt__temp_flt_comb_ = main_dir / (
            r'Final_results/Ppt_ok_ok_un_new2_first_flt__temp_flt__comb_%s' % temp_freq)
    #         r'Qt_ok_ok_un_3_test_first_flt__temp_flt__comb_%s' % temp_freq)
        Qt_ok_un_first_flt__temp_flt_comb_ = list_all_full_path(
            '.csv', path_to_Qt_ok_un_first_flt__temp_flt_comb_)

        path_to_Qt_ok_un_first_flt_1st_ = main_dir / (
            r'Final_results/Ppt_ok_ok_un_new2_first_flt__1st_%s' % temp_freq)
        Qt_ok_un_first_flt_1st_ = list_all_full_path(
            '.csv', path_to_Qt_ok_un_first_flt_1st_)

        path_to_Qt_ok_un_first_flt_comb_ = main_dir / (
            r'Final_results/Ppt_ok_ok_un_new2_first_flt__comb_%s' % temp_freq)
        Qt_ok_un_first_flt_comb_ = list_all_full_path(
            '.csv', path_to_Qt_ok_un_first_flt_comb_)

    #     path_to_Qt_ok_un_first_flt_comb_ = main_dir / (
    #         r'Qt_ok_ok_un_3_first_flt__comb_%s' % temp_freq)
    #     Qt_ok_un_first_flt_comb_ = list_all_full_path(
    #         '.csv', path_to_Qt_ok_un_first_flt_comb_)

        path_to_Quantiles_netatmo_no_flt___ = main_dir / (
            r'Final_results/Ppt_ok_ok_un_new_netatmo_no_flt___%s' % temp_freq)

        Quantiles_netatmo_no_flt___ = list_all_full_path(
            '.csv', path_to_Quantiles_netatmo_no_flt___)

        path_dwd_ppt_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\ppt_all_dwd_%s_.csv" % temp_freq
        df_dwd_edf = pd.read_csv(path_dwd_ppt_data, sep=';', index_col=0,
                                 parse_dates=True, infer_datetime_format=True)
        # df_dwd_edf.dropna(inplace=True)

        print(df_dwd_edf.isna().sum().max())

        df_compare = pd.DataFrame(
            index=df_dwd_edf.columns,
            columns=['pearson_corr_dwd_',
                     'spearman_corr_dwd_',
                     'pearson_corr_dwd_netatmo',
                     'spearman_corr_dwd_netatmo',
                     'pearson_corr_dwd_netatmo_unc2perc',
                     'spearman_corr_dwd_netatmo_unc2perc',
                     'pearson_corr_dwd_netatmo_unc5perc',
                     'spearman_corr_dwd_netatmo_unc5perc',
                     'pearson_corr_dwd_netatmo_unc10perc',
                     'spearman_corr_dwd_netatmo_unc10perc',
                     'pearson_corr_dwd_netatmo_unc20perc',
                     'spearman_corr_dwd_netatmo_unc20perc'])

        min_orig_qnt_thr = 0.

        #########################################################

        path_to_use = path_to_Qt_ok_un_first_flt__temp_flt_1st_
        data_to_use = Qt_ok_un_first_flt__temp_flt_1st_

        _interp_acc_ = str(r'%s' % (str(path_to_use).split('\\')[-1]))
        # for i in range(12):
        #   i = int(i)
        try:
            for df_file in data_to_use:

                if ('dwd_only') in df_file:
                    print(df_file)
                    path_interpolated_using_dwd = df_file

                if ('dwd_netamo') in df_file and 'unc' not in df_file:
                    print(df_file)
                    path_interpolated_using_netatmo_dwd = df_file

                if ('dwd_netamo') in df_file and 'unc2perc' in df_file:
                    print(df_file)
                    path_interpolated_using_netatmo_dwd_unc2perc = df_file

                if ('dwd_netamo') in df_file and 'unc5perc' in df_file:
                    print(df_file)
                    path_interpolated_using_netatmo_dwd_unc5perc = df_file

                if ('dwd_netamo') in df_file and 'unc10perc' in df_file:
                    print(df_file)
                    path_interpolated_using_netatmo_dwd_unc10perc = df_file

                if ('dwd_netamo') in df_file and 'unc20perc' in df_file:
                    print(df_file)
                    path_interpolated_using_netatmo_dwd_unc20perc = df_file

    #             if ('combined_interp_dwd_%s' % temp_freq) in df_file:
    #                 print(df_file)
    #                 path_interpolated_using_dwd = df_file
    #
    #             if 'combined_interp_dwd_netatmo' in df_file and 'unc' not in df_file:
    #
    #                 print(df_file)
    #                 path_interpolated_using_netatmo_dwd = df_file
    #
    #             if 'combined_interp_dwd_netatmo_unc' in df_file:
    #                 print(df_file)
    #                 path_interpolated_using_netatmo_dwd_unc = df_file
        except Exception as msg:
            print(msg)
            continue
        #########################################################

        df_dwd = pd.read_csv(path_interpolated_using_dwd,
                             sep=';', index_col=0, parse_dates=True,
                             infer_datetime_format=True)

        # f_dwd.dropna(inplace=True)

        print(df_dwd.isna().sum().max())

        df_netatmo_dwd = pd.read_csv(path_interpolated_using_netatmo_dwd,
                                     sep=';', index_col=0,
                                     parse_dates=True,
                                     infer_datetime_format=True)
        # df_netatmo_dwd.dropna(inplace=True)

        print(df_netatmo_dwd.isna().sum().max())

        df_netatmo_dwd_unc2perc = pd.read_csv(path_interpolated_using_netatmo_dwd_unc2perc,
                                              sep=';', index_col=0,
                                              parse_dates=True, infer_datetime_format=True)

        df_netatmo_dwd_unc5perc = pd.read_csv(path_interpolated_using_netatmo_dwd_unc5perc,
                                              sep=';', index_col=0,
                                              parse_dates=True, infer_datetime_format=True)

        df_netatmo_dwd_unc10perc = pd.read_csv(path_interpolated_using_netatmo_dwd_unc10perc,
                                               sep=';', index_col=0,
                                               parse_dates=True, infer_datetime_format=True)

        df_netatmo_dwd_unc20perc = pd.read_csv(path_interpolated_using_netatmo_dwd_unc20perc,
                                               sep=';', index_col=0,
                                               parse_dates=True, infer_datetime_format=True)
    #     df_netatmo_dwd_unc.dropna(inplace=True)

        #########################################################

        df_compare = pd.DataFrame(index=df_dwd.index)

    #     df_dwd_edf = df_dwd_edf.loc[df_dwd_edf.index.intersection(
    #         df_netatmo_dwd.index), :]
    #
    #     df_dwd_edf = df_dwd_edf[df_dwd_edf > 0]
        # df_dwd_edf.dropna(how='all')
        # stns to remove from orig edf, because no data for 2015-2019
        df_dwd_edf.drop(columns=['00384', '13672', '05155'], inplace=True)

        try:

            for event_date in df_dwd.index.intersection(df_dwd_edf.index):
                (orig_edf_vals, dwd_interp_vals,
                 netatmo_dwd_interp_vals, netatmo_dwd_interp_vals_unc2perc,
                 netatmo_dwd_interp_vals_unc5perc,
                 netatmo_dwd_interp_vals_unc10perc,
                 netatmo_dwd_interp_vals_unc20perc) = [], [], [], [], [], [], []
                print(event_date)
                original_quantile = df_dwd_edf.loc[event_date, :]
                interpolated_quantile_dwd = df_dwd.loc[event_date, :]
                interpolated_quantile_netatmo_dwd = df_netatmo_dwd.loc[event_date, :]
                interpolated_quantile_netatmo_dwd_unc2perc = df_netatmo_dwd_unc2perc.loc[
                    event_date, :]
                interpolated_quantile_netatmo_dwd_unc5perc = df_netatmo_dwd_unc5perc.loc[
                    event_date, :]
                interpolated_quantile_netatmo_dwd_unc10perc = df_netatmo_dwd_unc10perc.loc[
                    event_date, :]
                interpolated_quantile_netatmo_dwd_unc20perc = df_netatmo_dwd_unc20perc.loc[
                    event_date, :]
                for stn_ in original_quantile.index:
                    # print(stn_)
                    edf_stn_orig = original_quantile.loc[stn_]
                    edf_stn_interp_dwd = interpolated_quantile_dwd.loc[stn_]
                    edf_stn_interp_netatmo_dwd = interpolated_quantile_netatmo_dwd.loc[stn_]
                    edf_stn_interp_netatmo_dwd_unc2perc = interpolated_quantile_netatmo_dwd_unc2perc.loc[
                        stn_]
                    edf_stn_interp_netatmo_dwd_unc5perc = interpolated_quantile_netatmo_dwd_unc5perc.loc[
                        stn_]
                    edf_stn_interp_netatmo_dwd_unc10perc = interpolated_quantile_netatmo_dwd_unc10perc.loc[
                        stn_]
                    edf_stn_interp_netatmo_dwd_unc20perc = interpolated_quantile_netatmo_dwd_unc20perc.loc[
                        stn_]
                    if ((edf_stn_orig >= 0) and
                        (edf_stn_interp_dwd >= 0) and
                            (edf_stn_interp_netatmo_dwd >= 0) and

                            (edf_stn_interp_netatmo_dwd_unc2perc >= 0)):

                        orig_edf_vals.append(edf_stn_orig)
                        dwd_interp_vals.append(edf_stn_interp_dwd)
                        netatmo_dwd_interp_vals.append(
                            edf_stn_interp_netatmo_dwd)

                        netatmo_dwd_interp_vals_unc2perc.append(
                            edf_stn_interp_netatmo_dwd_unc2perc)
                        netatmo_dwd_interp_vals_unc5perc.append(
                            edf_stn_interp_netatmo_dwd_unc5perc)
                        netatmo_dwd_interp_vals_unc10perc.append(
                            edf_stn_interp_netatmo_dwd_unc10perc)
                        netatmo_dwd_interp_vals_unc20perc.append(
                            edf_stn_interp_netatmo_dwd_unc20perc)

                corr_dwd = pears(orig_edf_vals, dwd_interp_vals)[0]
                rho_dwd = spr(orig_edf_vals, dwd_interp_vals)[0]
                # print(corr_dwd, rho_dwd)

                corr_netatmo_dwd = pears(
                    orig_edf_vals, netatmo_dwd_interp_vals)[0]
                rho_netatmo_dwd = spr(
                    orig_edf_vals, netatmo_dwd_interp_vals)[0]
                # print(corr_netatmo_dwd, rho_netatmo_dwd)
                corr_netatmo_dwd_unc2perc = pears(
                    orig_edf_vals, netatmo_dwd_interp_vals_unc2perc)[0]
                rho_netatmo_dwd_unc2perc = spr(
                    orig_edf_vals, netatmo_dwd_interp_vals_unc2perc)[0]

                corr_netatmo_dwd_unc5perc = pears(
                    orig_edf_vals, netatmo_dwd_interp_vals_unc5perc)[0]
                rho_netatmo_dwd_unc5perc = spr(
                    orig_edf_vals, netatmo_dwd_interp_vals_unc5perc)[0]

                corr_netatmo_dwd_unc10perc = pears(
                    orig_edf_vals, netatmo_dwd_interp_vals_unc10perc)[0]
                rho_netatmo_dwd_unc10perc = spr(
                    orig_edf_vals, netatmo_dwd_interp_vals_unc10perc)[0]

                corr_netatmo_dwd_unc20perc = pears(
                    orig_edf_vals, netatmo_dwd_interp_vals_unc20perc)[0]
                rho_netatmo_dwd_unc20perc = spr(
                    orig_edf_vals, netatmo_dwd_interp_vals_unc20perc)[0]

                df_compare.loc[event_date, 'pearson_corr_dwd_'] = corr_dwd
                df_compare.loc[event_date, 'spearman_corr_dwd_'] = rho_dwd
                df_compare.loc[event_date,
                               'pearson_corr_dwd_netatmo'] = corr_netatmo_dwd
                df_compare.loc[event_date,
                               'spearman_corr_dwd_netatmo'] = rho_netatmo_dwd

                df_compare.loc[event_date,
                               'pearson_corr_dwd_netatmo_unc2perc'] = corr_netatmo_dwd_unc2perc
                df_compare.loc[event_date,
                               'spearman_corr_dwd_netatmo_unc2perc'] = rho_netatmo_dwd_unc2perc

                df_compare.loc[event_date,
                               'pearson_corr_dwd_netatmo_unc5perc'] = corr_netatmo_dwd_unc5perc
                df_compare.loc[event_date,
                               'spearman_corr_dwd_netatmo_unc5perc'] = rho_netatmo_dwd_unc5perc

                df_compare.loc[event_date,
                               'pearson_corr_dwd_netatmo_unc10perc'] = corr_netatmo_dwd_unc10perc
                df_compare.loc[event_date,
                               'spearman_corr_dwd_netatmo_unc10perc'] = rho_netatmo_dwd_unc10perc

                df_compare.loc[event_date,
                               'pearson_corr_dwd_netatmo_unc20perc'] = corr_netatmo_dwd_unc20perc
                df_compare.loc[event_date,
                               'spearman_corr_dwd_netatmo_unc20perc'] = rho_netatmo_dwd_unc20perc

                # df_compare = df_compare[df_compare > 0]
                #df_compare.dropna(how='all', inplace=True)

        except Exception as msg:
            print(msg)

#         df_compare.loc['2016-06-24 22:00:00', :]
# ['2016-06-24 22:00:00'

        # difference dwd- dwdnetatmo
#         max_diff_pears = df_compare['pearson_corr_dwd_'] - \
#             df_compare['pearson_corr_dwd_netatmo']
#         max_diff_pears = pd.DataFrame(
#             index=max_diff_pears.index, data=max_diff_pears.values)
# #         sorted_vals = max_diff_pears.sort_values()
# #         sorted_vals.to_csv(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\correlations_events\corr_dwd_minus_dwdnetatmo_pearson.csv', sep=';')
#         max_diff_spr = df_compare['spearman_corr_dwd_'] - \
#             df_compare['spearman_corr_dwd_netatmo']
# #         sorted_vals = max_diff_spr.sort_values()
#         # orig_edf_vals, dwd_interp_vals, netatmo_dwd_interp_vals
#         max_diff_pears['spr'] = max_diff_spr.values.ravel()
#         sorted_vals = max_diff_pears.sort_values(['spr'])
#         sorted_vals.to_csv(
# r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\correlations_events\corr_dwd_minus_dwdnetatmo_spearman_daily.csv',
# sep=';')

        # hourly_events = ['2016-06-25 00:00:00',
#                  '2018-06-11 17:00:00',
#                  '2018-09-23 17:00:00',
#                  '2018-09-23 18:00:00',
#                  '2018-09-23 19:00:00']
#         hourly_events = ['2019-07-27 19:00:00',
#                          '2019-07-27 20:00:00',
#                          '2018-05-13 16:00:00',
#                          '2018-05-13 22:00:00',
#                          '2018-06-12 18:00:00',
#                          '2018-07-05 05:00:00',
#                          '2018-08-02 01:00:00',
#                          '2018-08-23 15:00:00',
#                          '2018-09-06 18:00:00',
#                          '2016-06-24 22:00:00',
#                          '2016-06-25 00:00:00']

# daily_events = ['2018-12-23 00:00:00',
#                 '2019-05-22 00:00:00']
        df_compare.sort_index(inplace=True)

#         df_compare.loc['2018-09-06 18:00:00', :]
#         for evt in hourly_events:
#             print(df_compare.loc[evt, :])
#         df_compare.sort_values(by='pearson_corr_dwd_', inplace=True)
#         df_compare.reset_index()

        # OK
        stations_with_improvements = sum(i >= j for (i, j) in zip(
            df_compare.pearson_corr_dwd_netatmo.values,
            df_compare.pearson_corr_dwd_.values))

        stations_without_improvements = sum(i < j for (i, j) in zip(
            df_compare.pearson_corr_dwd_netatmo.values,
            df_compare.pearson_corr_dwd_.values))

        percent_of_improvment = 100 * (stations_with_improvements /
                                       df_compare.pearson_corr_dwd_netatmo.shape[0])

        # OK with Unc

        stations_with_improvements_unc2perc = sum(i >= j for (i, j) in zip(
            df_compare.pearson_corr_dwd_netatmo_unc2perc.values,
            df_compare.pearson_corr_dwd_.values))

        stations_without_improvements_unc2perc = sum(i < j for (i, j) in zip(
            df_compare.pearson_corr_dwd_netatmo_unc2perc.values,
            df_compare.pearson_corr_dwd_.values))

        percent_of_improvment_unc2perc = 100 * (
            stations_with_improvements_unc2perc /
            df_compare.pearson_corr_dwd_netatmo_unc2perc.shape[0])

        stations_with_improvements_unc5perc = sum(i >= j for (i, j) in zip(
            df_compare.pearson_corr_dwd_netatmo_unc5perc.values,
            df_compare.pearson_corr_dwd_.values))

        stations_without_improvements_unc5perc = sum(i < j for (i, j) in zip(
            df_compare.pearson_corr_dwd_netatmo_unc5perc.values,
            df_compare.pearson_corr_dwd_.values))

        percent_of_improvment_unc5perc = 100 * (
            stations_with_improvements_unc5perc /
            df_compare.pearson_corr_dwd_netatmo_unc2perc.shape[0])

        stations_with_improvements_unc10perc = sum(i >= j for (i, j) in zip(
            df_compare.pearson_corr_dwd_netatmo_unc10perc.values,
            df_compare.pearson_corr_dwd_.values))

        stations_without_improvements_unc10perc = sum(i < j for (i, j) in zip(
            df_compare.pearson_corr_dwd_netatmo_unc10perc.values,
            df_compare.pearson_corr_dwd_.values))

        percent_of_improvment_unc10perc = 100 * (
            stations_with_improvements_unc10perc /
            df_compare.pearson_corr_dwd_netatmo_unc2perc.shape[0])

        stations_with_improvements_unc20perc = sum(i >= j for (i, j) in zip(
            df_compare.pearson_corr_dwd_netatmo_unc20perc.values,
            df_compare.pearson_corr_dwd_.values))

        stations_without_improvements_unc20perc = sum(i < j for (i, j) in zip(
            df_compare.pearson_corr_dwd_netatmo_unc20perc.values,
            df_compare.pearson_corr_dwd_.values))

        percent_of_improvment_unc20perc = 100 * (
            stations_with_improvements_unc20perc /
            df_compare.pearson_corr_dwd_netatmo_unc2perc.shape[0])

        ####
        mean_pearson_correlation_dwd_only = df_compare.pearson_corr_dwd_.mean()
        mean_pearson_correlation_dwd_netatmo = df_compare.pearson_corr_dwd_netatmo.mean()

        mean_pearson_correlation_dwd_netatmo_unc2perc = df_compare.pearson_corr_dwd_netatmo_unc2perc.mean()

        mean_pearson_correlation_dwd_netatmo_unc5perc = df_compare.pearson_corr_dwd_netatmo_unc5perc.mean()
        mean_pearson_correlation_dwd_netatmo_unc10perc = df_compare.pearson_corr_dwd_netatmo_unc10perc.mean()
        mean_pearson_correlation_dwd_netatmo_unc20perc = df_compare.pearson_corr_dwd_netatmo_unc20perc.mean()

        #########################################################
        mean_spr_correlation_dwd_only = df_compare.spearman_corr_dwd_.mean()
        mean_spr_correlation_dwd_netatmo = df_compare.spearman_corr_dwd_netatmo.mean()

        mean_spr_correlation_dwd_netatmo_unc2perc = df_compare.spearman_corr_dwd_netatmo_unc2perc.mean()
        mean_spr_correlation_dwd_netatmo_unc5perc = df_compare.spearman_corr_dwd_netatmo_unc5perc.mean()
        mean_spr_correlation_dwd_netatmo_unc10perc = df_compare.spearman_corr_dwd_netatmo_unc10perc.mean()
        mean_spr_correlation_dwd_netatmo_unc20perc = df_compare.spearman_corr_dwd_netatmo_unc20perc.mean()
        #########################################################
        plt.ioff()
        fig = plt.figure(figsize=(24, 12), dpi=150)

        ax = fig.add_subplot(111)

        ax.plot(df_compare.index,
                df_compare.pearson_corr_dwd_,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation %0.2f'
                % mean_pearson_correlation_dwd_only)
        ax.plot(df_compare.index,
                df_compare.pearson_corr_dwd_netatmo,
                alpha=.5,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %0.2f'
                % mean_pearson_correlation_dwd_netatmo)

        ax.plot(df_compare.index,
                df_compare.pearson_corr_dwd_netatmo_unc2perc,
                alpha=.5,
                c='g',  # colors_arr,
                marker='+',
                label='DWD-Netatmo Interpolation 2percUnc %0.2f'
                % mean_pearson_correlation_dwd_netatmo_unc2perc)

        ax.plot(df_compare.index,
                df_compare.pearson_corr_dwd_netatmo_unc5perc,
                alpha=.5,
                c='m',  # colors_arr,
                marker='1',
                label='DWD-Netatmo Interpolation 5percUnc %0.2f'
                % mean_pearson_correlation_dwd_netatmo_unc5perc)

        ax.plot(df_compare.index,
                df_compare.pearson_corr_dwd_netatmo_unc10perc,
                alpha=.5,
                c='c',  # colors_arr,
                marker='3',
                label='DWD-Netatmo Interpolation 10percUnc %0.2f'
                % mean_pearson_correlation_dwd_netatmo_unc10perc)

        ax.plot(df_compare.index,
                df_compare.pearson_corr_dwd_netatmo_unc20perc,
                alpha=.5,
                c='orange',  # colors_arr,
                marker=',',
                label='DWD-Netatmo Interpolation 20percUnc %0.2f'
                % mean_pearson_correlation_dwd_netatmo_unc20perc)

        ax.set_title('Pearson Correlation Interpolated Quantiles from DWD or DWD-Netatmo\n '
                     'Precipitation of %s Extreme Events %s\n Events with Improvemnts %d / %d, Percentage %0.0f\n'
                     'Events with Improvemnts with OK 2percUnc %d / %d, Percentage %0.0f'
                     '\nEvents with Improvemnts with OK 5percUnc %d / %d, Percentage %0.0f'
                     '\nEventswith Improvemnts with OK 10percUnc %d / %d, Percentage %0.0f'
                     '\nEvents with Improvemnts with OK 20percUnc %d / %d, Percentage %0.0f'
                     % (temp_freq, _interp_acc_,
                        stations_with_improvements,
                         df_compare.pearson_corr_dwd_netatmo.shape[0],
                        percent_of_improvment,
                        stations_with_improvements_unc2perc,
                        df_compare.pearson_corr_dwd_netatmo_unc2perc.shape[0],
                        percent_of_improvment_unc2perc,
                        stations_with_improvements_unc5perc,
                        df_compare.pearson_corr_dwd_netatmo_unc2perc.shape[0],
                        percent_of_improvment_unc5perc,
                        stations_with_improvements_unc10perc,
                        df_compare.pearson_corr_dwd_netatmo_unc2perc.shape[0],
                        percent_of_improvment_unc10perc,
                        stations_with_improvements_unc20perc,
                        df_compare.pearson_corr_dwd_netatmo_unc2perc.shape[0],
                        percent_of_improvment_unc20perc))
        ax.grid(alpha=0.25)
        plt.setp(ax.get_xticklabels(), rotation=45)

        # years = mdates.YearLocator(1)   # every year
        months = mdates.MonthLocator()  # every month

        #hrs = mdates.HourLocator()
        # years_fmt = mdates.DateFormatter('%Y-%m-%d')

        years_fmt = mdates.DateFormatter('%Y-%m-%d %H')
        ax.xaxis.set_major_locator(months)

        ax.xaxis.set_major_formatter(years_fmt)
        # ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))

        #ax.set_yticks(np.arange(0., 1.05, .10))
        ax.set_xlabel('Date of Event')
        ax.legend(loc='lower right')
        ax.set_ylabel('Pearson Correlation')

        plt.savefig(os.path.join(path_to_use,
                                 r'ppt_spatial_pears_corr_%s_events_dwd2.png' % temp_freq,
                                 ),
                    frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

        ####################################################
        stations_with_improvements = sum(i >= j for (i, j) in zip(
            df_compare.spearman_corr_dwd_netatmo.values,
            df_compare.spearman_corr_dwd_.values))

        stations_without_improvements = sum(i < j for (i, j) in zip(
            df_compare.spearman_corr_dwd_netatmo.values,
            df_compare.spearman_corr_dwd_.values))

        percent_of_improvment = 100 * (stations_with_improvements /
                                       df_compare.spearman_corr_dwd_netatmo.shape[0])

        # ok with Un

        stations_with_improvements_unc2perc = sum(i >= j for (i, j) in zip(
            df_compare.spearman_corr_dwd_netatmo_unc2perc.values,
            df_compare.spearman_corr_dwd_.values))

        stations_without_improvements_unc2perc = sum(i < j for (i, j) in zip(
            df_compare.spearman_corr_dwd_netatmo_unc2perc.values,
            df_compare.spearman_corr_dwd_.values))

        percent_of_improvment_unc2perc = 100 * (stations_with_improvements_unc2perc /
                                                df_compare.spearman_corr_dwd_netatmo_unc2perc.shape[0])

        stations_with_improvements_unc5perc = sum(i >= j for (i, j) in zip(
            df_compare.spearman_corr_dwd_netatmo_unc5perc.values,
            df_compare.spearman_corr_dwd_.values))

        stations_without_improvements_unc5perc = sum(i < j for (i, j) in zip(
            df_compare.spearman_corr_dwd_netatmo_unc5perc.values,
            df_compare.spearman_corr_dwd_.values))

        percent_of_improvment_unc5perc = 100 * (
            stations_with_improvements_unc5perc /
            df_compare.spearman_corr_dwd_netatmo_unc2perc.shape[0])

        stations_with_improvements_unc10perc = sum(i >= j for (i, j) in zip(
            df_compare.spearman_corr_dwd_netatmo_unc10perc.values,
            df_compare.spearman_corr_dwd_.values))

        stations_without_improvements_unc10perc = sum(i < j for (i, j) in zip(
            df_compare.spearman_corr_dwd_netatmo_unc10perc.values,
            df_compare.spearman_corr_dwd_.values))

        percent_of_improvment_unc10perc = 100 * (
            stations_with_improvements_unc10perc /
            df_compare.spearman_corr_dwd_netatmo_unc2perc.shape[0])

        stations_with_improvements_unc20perc = sum(i >= j for (i, j) in zip(
            df_compare.spearman_corr_dwd_netatmo_unc20perc.values,
            df_compare.spearman_corr_dwd_.values))

        stations_without_improvements_unc20perc = sum(i < j for (i, j) in zip(
            df_compare.spearman_corr_dwd_netatmo_unc20perc.values,
            df_compare.spearman_corr_dwd_.values))

        percent_of_improvment_unc20perc = 100 * (
            stations_with_improvements_unc20perc /
            df_compare.spearman_corr_dwd_netatmo_unc2perc.shape[0])
        #########################################################
        plt.ioff()
        fig = plt.figure(figsize=(24, 12), dpi=150)

        ax = fig.add_subplot(111)

        ax.plot(df_compare.index,
                df_compare.spearman_corr_dwd_,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation %0.2f'
                % mean_spr_correlation_dwd_only)
        ax.plot(df_compare.index,
                df_compare.spearman_corr_dwd_netatmo,
                alpha=.5,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %0.2f'
                % mean_spr_correlation_dwd_netatmo)

        ax.plot(df_compare.index,
                df_compare.spearman_corr_dwd_netatmo_unc2perc,
                alpha=.5,
                c='g',  # colors_arr,
                marker='+',
                label='DWD-Netatmo Interpolation 2percUnc %0.2f'
                % mean_spr_correlation_dwd_netatmo_unc2perc)

        ax.plot(df_compare.index,
                df_compare.spearman_corr_dwd_netatmo_unc5perc,
                alpha=.5,
                c='m',  # colors_arr,
                marker='1',
                label='DWD-Netatmo Interpolation 5percUnc %0.2f'
                % mean_spr_correlation_dwd_netatmo_unc5perc)

        ax.plot(df_compare.index,
                df_compare.spearman_corr_dwd_netatmo_unc10perc,
                alpha=.5,
                c='c',  # colors_arr,
                marker='3',
                label='DWD-Netatmo Interpolation 10percUnc %0.2f'
                % mean_spr_correlation_dwd_netatmo_unc10perc)

        ax.plot(df_compare.index,
                df_compare.spearman_corr_dwd_netatmo_unc20perc,
                alpha=.5,
                c='orange',  # colors_arr,
                marker=',',
                label='DWD-Netatmo Interpolation 20percUnc %0.2f'
                % mean_spr_correlation_dwd_netatmo_unc20perc)

        ax.set_title('Spearman Correlation Interpolated Quantiles from DWD or DWD-Netatmo \n '
                     'Rainfall of %s Extreme Events %s \n Events with Improvemnts %d / %d, Percentage %0.0f'
                     '\nEvents with Improvemnts with OK 2percUnc %d / %d, Percentage %0.0f'
                     '\nEvents with Improvemnts with OK 5percUnc %d / %d, Percentage %0.0f'
                     '\nEvents with Improvemnts with OK 10percUnc %d / %d, Percentage %0.0f'
                     '\nEvents with Improvemnts with OK 20percUnc %d / %d, Percentage %0.0f'
                     % (temp_freq, _interp_acc_,
                        stations_with_improvements,
                         df_compare.spearman_corr_dwd_netatmo.shape[0],
                        percent_of_improvment,
                        stations_with_improvements_unc2perc,
                        df_compare.spearman_corr_dwd_netatmo_unc2perc.shape[0],
                        percent_of_improvment_unc2perc,
                        stations_with_improvements_unc5perc,
                        df_compare.spearman_corr_dwd_netatmo_unc2perc.shape[0],
                        percent_of_improvment_unc5perc,
                        stations_with_improvements_unc10perc,
                        df_compare.spearman_corr_dwd_netatmo_unc2perc.shape[0],
                        percent_of_improvment_unc10perc,
                        stations_with_improvements_unc20perc,
                        df_compare.spearman_corr_dwd_netatmo_unc2perc.shape[0],
                        percent_of_improvment_unc20perc))
        ax.grid(alpha=0.25)
        plt.setp(ax.get_xticklabels(), rotation=60)

        #time_vals = df_compare.index.to_pydatetime()
        #dates_formatted = [pd.to_datetime(d, format='%Y-%m') for d in time_vals]

        # years = mdates.YearLocator(1)   # every year
        months = mdates.MonthLocator()  # every month
        #days = mdates.DayLocator(10)
        years_fmt = mdates.DateFormatter('%Y-%m-%d %H:%M:%S')
        # set ticks every week
        ax.xaxis.set_major_locator(mdates.MonthLocator())

        # set major ticks format

        # ax.set_xticks(dates_formatted[::10])
        # ax.xaxis.set_major_locator(days)
        ax.xaxis.set_major_formatter(years_fmt)
        # ax.xaxis.set_minor_formatter(years_fmt)

        #ax.set_yticks(np.arange(0.5, 1.05, .10))
        ax.set_xlabel('Date of Event')
        ax.legend(loc='lower left')
        ax.set_ylabel('Spearman Correlation')

        plt.savefig(os.path.join(path_to_use,
                                 r'ppt_spatial_spearm_corr_%s_events_dwd2.png' % (temp_freq)),
                    frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()


#==============================================================================
#
#==============================================================================

if plot_not_filtered:
        # , '180min', '360min', '720min', '1440min']:
    for temp_freq in ['60min', '180min', '360min', '720min', '1440min']:
        print(temp_freq)

        path_to_Quantiles_netatmo_no_flt___ = main_dir / (
            r'Final_results/Ppt_ok_ok_un_new_netatmo_no_flt___%s' % temp_freq)

        Quantiles_netatmo_no_flt___ = list_all_full_path(
            '.csv', path_to_Quantiles_netatmo_no_flt___)

        path_dwd_ppt_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\ppt_all_dwd_%s_.csv" % temp_freq
        df_dwd_edf = pd.read_csv(path_dwd_ppt_data, sep=';', index_col=0,
                                 parse_dates=True, infer_datetime_format=True)
        # df_dwd_edf.dropna(inplace=True)

        print(df_dwd_edf.isna().sum().max())

        df_compare = pd.DataFrame(
            index=df_dwd_edf.columns,
            columns=['pearson_corr_dwd_',
                     'spearman_corr_dwd_',
                     'pearson_corr_dwd_netatmo',
                     'spearman_corr_dwd_netatmo'])

        min_orig_qnt_thr = 0.

        #########################################################

        path_to_use = path_to_Quantiles_netatmo_no_flt___
        data_to_use = Quantiles_netatmo_no_flt___

        _interp_acc_ = str(r'%s' % (str(path_to_use).split('\\')[-1]))
        # for i in range(12):
        #   i = int(i)
        try:
            for df_file in data_to_use:

                if ('dwd_only') in df_file:
                    print(df_file)
                    path_interpolated_using_dwd = df_file

                if ('dwd_netamo') in df_file and 'unc' not in df_file:
                    print(df_file)
                    path_interpolated_using_netatmo_dwd = df_file

    #             if ('combined_interp_dwd_%s' % temp_freq) in df_file:
    #                 print(df_file)
    #                 path_interpolated_using_dwd = df_file
    #
    #             if 'combined_interp_dwd_netatmo' in df_file and 'unc' not in df_file:
    #
    #                 print(df_file)
    #                 path_interpolated_using_netatmo_dwd = df_file
    #
    #             if 'combined_interp_dwd_netatmo_unc' in df_file:
    #                 print(df_file)
    #                 path_interpolated_using_netatmo_dwd_unc = df_file
        except Exception as msg:
            print(msg)
            continue
        #########################################################

        df_dwd = pd.read_csv(path_interpolated_using_dwd,
                             sep=';', index_col=0, parse_dates=True,
                             infer_datetime_format=True)

        # f_dwd.dropna(inplace=True)

        print(df_dwd.isna().sum().max())

        df_netatmo_dwd = pd.read_csv(path_interpolated_using_netatmo_dwd,
                                     sep=';', index_col=0,
                                     parse_dates=True,
                                     infer_datetime_format=True)
        # df_netatmo_dwd.dropna(inplace=True)

        #########################################################

        df_compare = pd.DataFrame(index=df_dwd.index)

    #     df_dwd_edf = df_dwd_edf.loc[df_dwd_edf.index.intersection(
    #         df_netatmo_dwd.index), :]
    #
    #     df_dwd_edf = df_dwd_edf[df_dwd_edf > 0]
        # df_dwd_edf.dropna(how='all')

        try:

            for event_date in df_dwd.index.intersection(df_dwd_edf.index):
                (orig_edf_vals, dwd_interp_vals,
                 netatmo_dwd_interp_vals) = [], [], []
                print(event_date)
                original_quantile = df_dwd_edf.loc[event_date, :]
                interpolated_quantile_dwd = df_dwd.loc[event_date, :]
                interpolated_quantile_netatmo_dwd = df_netatmo_dwd.loc[event_date, :]

                for stn_ in original_quantile.index:
                    # print(stn_)
                    edf_stn_orig = original_quantile.loc[stn_]
                    edf_stn_interp_dwd = interpolated_quantile_dwd.loc[stn_]
                    edf_stn_interp_netatmo_dwd = interpolated_quantile_netatmo_dwd.loc[stn_]

                    if ((edf_stn_orig >= 0) and
                        (edf_stn_interp_dwd >= 0) and
                            (edf_stn_interp_netatmo_dwd >= 0)):
                        orig_edf_vals.append(edf_stn_orig)
                        dwd_interp_vals.append(edf_stn_interp_dwd)
                        netatmo_dwd_interp_vals.append(
                            edf_stn_interp_netatmo_dwd)

                corr_dwd = pears(orig_edf_vals, dwd_interp_vals)[0]
                rho_dwd = spr(orig_edf_vals, dwd_interp_vals)[0]
                # print(corr_dwd, rho_dwd)

                corr_netatmo_dwd = pears(
                    orig_edf_vals, netatmo_dwd_interp_vals)[0]
                rho_netatmo_dwd = spr(
                    orig_edf_vals, netatmo_dwd_interp_vals)[0]
                # print(corr_netatmo_dwd, rho_netatmo_dwd)

                df_compare.loc[event_date, 'pearson_corr_dwd_'] = corr_dwd
                df_compare.loc[event_date, 'spearman_corr_dwd_'] = rho_dwd
                df_compare.loc[event_date,
                               'pearson_corr_dwd_netatmo'] = corr_netatmo_dwd
                df_compare.loc[event_date,
                               'spearman_corr_dwd_netatmo'] = rho_netatmo_dwd

                # df_compare = df_compare[df_compare > 0]
                #df_compare.dropna(how='all', inplace=True)

        except Exception as msg:
            print(msg)

        # orig_edf_vals, dwd_interp_vals, netatmo_dwd_interp_vals

        df_compare.sort_index(inplace=True)
#         df_compare.loc['2018-06-11 16:00:00', :]
        # OK
        stations_with_improvements = sum(i >= j for (i, j) in zip(
            df_compare.pearson_corr_dwd_netatmo.values,
            df_compare.pearson_corr_dwd_.values))

        stations_without_improvements = sum(i < j for (i, j) in zip(
            df_compare.pearson_corr_dwd_netatmo.values,
            df_compare.pearson_corr_dwd_.values))

        percent_of_improvment = 100 * (stations_with_improvements /
                                       df_compare.pearson_corr_dwd_netatmo.shape[0])

        # OK with Unc
        # hourly_events = ['2016-06-25 00:00:00',
#                  '2018-06-11 17:00:00',
#                  '2018-09-23 17:00:00',
#                  '2018-09-23 18:00:00',
#                  '2018-09-23 19:00:00']

# daily_events = ['2018-12-23 00:00:00',
#                 '2019-05-22 00:00:00',
#                    '2018-05-14 00:00:00']
#         df_compare.loc['2019-07-29 19:00:00', :]
        ####
        mean_pearson_correlation_dwd_only = df_compare.pearson_corr_dwd_.mean()
        mean_pearson_correlation_dwd_netatmo = df_compare.pearson_corr_dwd_netatmo.mean()

        #########################################################
        mean_spr_correlation_dwd_only = df_compare.spearman_corr_dwd_.mean()
        mean_spr_correlation_dwd_netatmo = df_compare.spearman_corr_dwd_netatmo.mean()

        #########################################################

        df_compare.sort_values(by='pearson_corr_dwd_', inplace=True)
        df_compare.reset_index()
        plt.ioff()
        fig = plt.figure(figsize=(24, 12), dpi=150)

        ax = fig.add_subplot(111)

        ax.plot(df_compare.index,
                df_compare.pearson_corr_dwd_,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation %0.2f'
                % mean_pearson_correlation_dwd_only)
        ax.plot(df_compare.index,
                df_compare.pearson_corr_dwd_netatmo,
                alpha=.5,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %0.2f'
                % mean_pearson_correlation_dwd_netatmo)

        ax.set_title('Pearson Correlation Interpolated Quantiles from DWD or DWD-Netatmo\n '
                     'Precipitation of %s Intense Events %s\n Events with Improvemnts %d / %d, Percentage %0.0f\n'

                     % (temp_freq, _interp_acc_,
                        stations_with_improvements,
                         df_compare.pearson_corr_dwd_netatmo.shape[0],
                        percent_of_improvment
                        ))
        ax.grid(alpha=0.25)
        plt.setp(ax.get_xticklabels(), rotation=45)

        # years = mdates.YearLocator(1)   # every year
        months = mdates.MonthLocator()  # every month

        #hrs = mdates.HourLocator()
        # years_fmt = mdates.DateFormatter('%Y-%m-%d')

        years_fmt = mdates.DateFormatter('%Y-%m-%d %H')
        ax.xaxis.set_major_locator(months)

        ax.xaxis.set_major_formatter(years_fmt)
        # ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))

        #ax.set_yticks(np.arange(0.5, 1.05, .10))
        ax.set_xlabel('Date of Event')
        ax.legend(loc='lower left')
        ax.set_ylabel('Pearson Correlation')

        plt.savefig(os.path.join(path_to_use,
                                 r'ppt_spatial_pears_corr_%s_events_dwd2.png' % temp_freq,
                                 ),
                    frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

        ####################################################
        stations_with_improvements = sum(i >= j for (i, j) in zip(
            df_compare.spearman_corr_dwd_netatmo.values,
            df_compare.spearman_corr_dwd_.values))

        stations_without_improvements = sum(i < j for (i, j) in zip(
            df_compare.spearman_corr_dwd_netatmo.values,
            df_compare.spearman_corr_dwd_.values))

        percent_of_improvment = 100 * (stations_with_improvements /
                                       df_compare.spearman_corr_dwd_netatmo.shape[0])

        #########################################################
        plt.ioff()
        fig = plt.figure(figsize=(24, 12), dpi=150)

        ax = fig.add_subplot(111)

        ax.plot(df_compare.index,
                df_compare.spearman_corr_dwd_,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation %0.2f'
                % mean_spr_correlation_dwd_only)
        ax.plot(df_compare.index,
                df_compare.spearman_corr_dwd_netatmo,
                alpha=.5,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %0.2f'
                % mean_spr_correlation_dwd_netatmo)

        ax.set_title('Spearman Correlation Interpolated Quantiles from DWD or DWD-Netatmo \n '
                     'Rainfall of %s Intense Events %s \n Events with Improvemnts %d / %d, Percentage %0.0f'

                     % (temp_freq, _interp_acc_,
                        stations_with_improvements,
                         df_compare.spearman_corr_dwd_netatmo.shape[0],
                        percent_of_improvment))
        ax.grid(alpha=0.25)
        plt.setp(ax.get_xticklabels(), rotation=60)

        #time_vals = df_compare.index.to_pydatetime()
        #dates_formatted = [pd.to_datetime(d, format='%Y-%m') for d in time_vals]

        # years = mdates.YearLocator(1)   # every year
        months = mdates.MonthLocator()  # every month
        #days = mdates.DayLocator(10)
        years_fmt = mdates.DateFormatter('%Y-%m-%d %H:%M:%S')
        # set ticks every week
        ax.xaxis.set_major_locator(mdates.MonthLocator())

        # set major ticks format

        # ax.set_xticks(dates_formatted[::10])
        # ax.xaxis.set_major_locator(days)
        ax.xaxis.set_major_formatter(years_fmt)
        # ax.xaxis.set_minor_formatter(years_fmt)

        #ax.set_yticks(np.arange(0.5, 1.05, .10))
        ax.set_xlabel('Date of Event')
        ax.legend(loc='lower left')
        ax.set_ylabel('Spearman Correlation')

        plt.savefig(os.path.join(path_to_use,
                                 r'ppt_spatial_spearm_corr_%s_events_dwd2.png' % (temp_freq)),
                    frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()
