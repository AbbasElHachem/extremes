
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears

import fnmatch
from pathlib import Path

from _00_additional_functions import (list_all_full_path,
                                      build_edf_fr_vals, find_nearest,
                                      calculate_probab_ppt_below_thr)
plt.ioff()
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'axes.labelsize': 12})


#==============================================================================
# # In[2]:
#==============================================================================

main_dir = Path(
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo')

#main_dir = Path(r'/home/IWS/hachem/Extremes')

min_orig_qnt_thr = 0.
_acc_ = ''

plot_not_filtered = False  # False
plot_filtered = True  # True
#==============================================================================
#
#==============================================================================

if plot_filtered:
    # '60min', '180min', '360min', '720min', '1440min']:
    for temp_freq in ['60min', '180min', '360min', '720min', '1440min']:
        # '60min', '360min',  '1440min'
        print(temp_freq)

        path_to_Qt_ok_un_first_flt__temp_flt_1st_ = main_dir / (
            r'Final_results4/Ppt_ok_ok_un_new3_first_flt__temp_flt__1st_%s' % temp_freq)
        Qt_ok_un_first_flt__temp_flt_1st_ = list_all_full_path(
            '.csv', path_to_Qt_ok_un_first_flt__temp_flt_1st_)

        path_to_Qt_ok_un_first_flt__temp_flt_comb_ = main_dir / (
            r'Final_results/Ppt_ok_ok_un_new2_first_flt__temp_flt__comb_%s' % temp_freq)
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

        #########################################################
        path_dwd_edf_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\edf_ppt_all_dwd_%s_.csv" % temp_freq

        df_dwd_edf = pd.read_csv(path_dwd_edf_data, sep=';', index_col=0,
                                 parse_dates=True, infer_datetime_format=True)

        path_dwd_ppt_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\ppt_all_dwd_%s_.csv" % temp_freq

        df_dwd_ppt = pd.read_csv(path_dwd_ppt_data, sep=';', index_col=0,
                                 parse_dates=True, infer_datetime_format=True)

        #########################################################
        df_improvements = pd.DataFrame(
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

        except Exception as msg:
            print(msg)
            continue
            #########################################################

        df_dwd = pd.read_csv(path_interpolated_using_dwd,
                             sep=';', index_col=0, parse_dates=True,
                             infer_datetime_format=True)

        print(df_dwd.isna().sum().max())

        df_netatmo_dwd = pd.read_csv(path_interpolated_using_netatmo_dwd,
                                     sep=';', index_col=0,
                                     parse_dates=True,
                                     infer_datetime_format=True)

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

        cmn_interpolated_events = df_netatmo_dwd.index.intersection(
            df_dwd.index).intersection(df_dwd_edf.index)

        print('Total number of events is',
              cmn_interpolated_events.shape[0])

        # stns to remove from orig edf, because no data for 2015-2019
        df_dwd_edf.drop(columns=['00384', '13672', '05155'], inplace=True)

        for stn_ in df_dwd.columns:
            df_compare = pd.DataFrame(index=df_dwd.index)

            for event_date in cmn_interpolated_events:
                # print(event_date)
                #                 if str(event_date) == '2015-06-07 22:00:00':
                #                     print(event_date)
                # interpolated_quantile_netatmo = df_netatmo.loc[event_date, stn_]
                #event_date = cmn_interpolated_events[0]
                interpolated_ppt_dwd = df_dwd.loc[event_date, stn_]

                interpolated_ppt_netatmo_dwd = df_netatmo_dwd.loc[event_date, stn_]

                interpolated_quantile_netatmo_dwd_unc2perc = df_netatmo_dwd_unc2perc.loc[
                    event_date, stn_]
                interpolated_quantile_netatmo_dwd_unc5perc = df_netatmo_dwd_unc5perc.loc[
                    event_date, stn_]
                interpolated_quantile_netatmo_dwd_unc10perc = df_netatmo_dwd_unc10perc.loc[
                    event_date, stn_]
                interpolated_quantile_netatmo_dwd_unc20perc = df_netatmo_dwd_unc20perc.loc[
                    event_date, stn_]

                # original quantile when transforming ppt to edf
                original_ppt = df_dwd_ppt.loc[event_date, stn_]

                df_compare.loc[event_date,
                               'original_ppt'] = original_ppt

                df_compare.loc[
                    event_date,
                    'interpolated_ppt_dwd'] = interpolated_ppt_dwd

                df_compare.loc[
                    event_date,
                    'interpolated_ppt_netatmo_dwd'] = interpolated_ppt_netatmo_dwd
                df_compare.loc[
                    event_date,
                    'interpolated_ppt_netatmo_dwd_unc2perc'] = interpolated_quantile_netatmo_dwd_unc2perc

                df_compare.loc[
                    event_date,
                    'interpolated_ppt_netatmo_dwd_unc5perc'] = interpolated_quantile_netatmo_dwd_unc5perc

                df_compare.loc[
                    event_date,
                    'interpolated_ppt_netatmo_dwd_unc10perc'] = interpolated_quantile_netatmo_dwd_unc10perc

                df_compare.loc[
                    event_date,
                    'interpolated_ppt_netatmo_dwd_unc20perc'] = interpolated_quantile_netatmo_dwd_unc20perc

            df_compare = df_compare[df_compare >= 0]
            df_compare.dropna(how='any', inplace=True)
            #==============================================================
            # # PPT
            #==============================================================
            values_x_ppt = df_compare['original_ppt'].values
            values_dwd_ppt = df_compare['interpolated_ppt_dwd'].values

            values_netatmo_dwd_ppt = df_compare['interpolated_ppt_netatmo_dwd'].values
            values_netatmo_dwd_unc2perc = df_compare['interpolated_ppt_netatmo_dwd_unc2perc'].values

            values_netatmo_dwd_unc5perc = df_compare['interpolated_ppt_netatmo_dwd_unc5perc'].values
            values_netatmo_dwd_unc10perc = df_compare['interpolated_ppt_netatmo_dwd_unc10perc'].values
            values_netatmo_dwd_unc20perc = df_compare['interpolated_ppt_netatmo_dwd_unc20perc'].values

            # calculate sqared error between obsv and interpolated
            # ppt
            try:
                mse_dwd_interp_ppt = np.sqrt(np.square(
                    np.subtract(values_x_ppt, values_dwd_ppt)).mean())
            except Exception as msg:
                print(msg)

            mse_dwd_netatmo_interp_ppt = np.sqrt(np.square(
                np.subtract(values_x_ppt, values_netatmo_dwd_ppt)).mean())

            mse_dwd_netatmo_interp_unc_ppt_2perc = np.sqrt(np.square(
                np.subtract(values_x_ppt, values_netatmo_dwd_unc2perc)).mean())

            mse_dwd_netatmo_interp_unc_ppt_5perc = np.sqrt(np.square(
                np.subtract(values_x_ppt, values_netatmo_dwd_unc5perc)).mean())

            mse_dwd_netatmo_interp_unc_ppt_10perc = np.sqrt(np.square(
                np.subtract(values_x_ppt, values_netatmo_dwd_unc10perc)).mean())

            mse_dwd_netatmo_interp_unc_ppt_20perc = np.sqrt(np.square(
                np.subtract(values_x_ppt, values_netatmo_dwd_unc20perc)).mean())

            # calculate correlations (pearson and spearman)

            df_improvements.loc[stn_,
                                'mse_dwd_interp_ppt'] = mse_dwd_interp_ppt
            df_improvements.loc[stn_,
                                'mse_dwd_netatmo_interp_ppt'] = mse_dwd_netatmo_interp_ppt
            df_improvements.loc[stn_,
                                'mse_dwd_netatmo_interp_unc_ppt_2perc'] = mse_dwd_netatmo_interp_unc_ppt_2perc

            df_improvements.loc[stn_,
                                'mse_dwd_netatmo_interp_unc_ppt_5perc'] = mse_dwd_netatmo_interp_unc_ppt_5perc

            df_improvements.loc[stn_,
                                'mse_dwd_netatmo_interp_unc_ppt_10perc'] = mse_dwd_netatmo_interp_unc_ppt_10perc

            df_improvements.loc[stn_,
                                'mse_dwd_netatmo_interp_unc_ppt_20perc'] = mse_dwd_netatmo_interp_unc_ppt_20perc

        df_improvements.sort_values(by=['mse_dwd_interp_ppt'],
                                    ascending=True, inplace=True)
        df_improvements.dropna(how='all', inplace=True)
        #======================================================================
        # PPT RMSE
        #======================================================================
        stations_with_mse_improvements_ppt = sum(i < j for (i, j) in zip(
            df_improvements.mse_dwd_netatmo_interp_ppt.values,
            df_improvements.mse_dwd_interp_ppt.values))

        stations_without_mse_improvements_ppt = sum(i >= j for (i, j) in zip(
            df_improvements.mse_dwd_netatmo_interp_ppt.values,
            df_improvements.mse_dwd_interp_ppt.values))

        percent_of_mse_improvment_ppt = 100 * (stations_with_mse_improvements_ppt /
                                               df_improvements.mse_dwd_interp_ppt.shape[0])

        ##########################
        stations_with_mse_improvements_unc_ppt_2perc = sum(i < j for (i, j) in zip(
            df_improvements.mse_dwd_netatmo_interp_unc_ppt_2perc.values,
            df_improvements.mse_dwd_interp_ppt.values))

        stations_without_mse_improvements_unc_ppt_2perc = sum(i >= j for (i, j) in zip(
            df_improvements.mse_dwd_netatmo_interp_unc_ppt_2perc.values,
            df_improvements.mse_dwd_interp_ppt.values))

        percent_of_mse_improvment_unc_ppt_2perc = 100 * (stations_with_mse_improvements_unc_ppt_2perc /
                                                         df_improvements.mse_dwd_interp_ppt.shape[0])
        # 5perc
        stations_with_mse_improvements_unc_ppt_5perc = sum(i < j for (i, j) in zip(
            df_improvements.mse_dwd_netatmo_interp_unc_ppt_5perc.values,
            df_improvements.mse_dwd_interp_ppt.values))

        stations_without_mse_improvements_unc_ppt_5perc = sum(i >= j for (i, j) in zip(
            df_improvements.mse_dwd_netatmo_interp_unc_ppt_5perc.values,
            df_improvements.mse_dwd_interp_ppt.values))

        percent_of_mse_improvment_unc_ppt_5perc = 100 * (stations_with_mse_improvements_unc_ppt_5perc /
                                                         df_improvements.mse_dwd_interp_ppt.shape[0])

        # 10perc
        stations_with_mse_improvements_unc_ppt_10perc = sum(i < j for (i, j) in zip(
            df_improvements.mse_dwd_netatmo_interp_unc_ppt_10perc.values,
            df_improvements.mse_dwd_interp_ppt.values))

        stations_without_mse_improvements_unc_ppt_10perc = sum(i >= j for (i, j) in zip(
            df_improvements.mse_dwd_netatmo_interp_unc_ppt_10perc.values,
            df_improvements.mse_dwd_interp_ppt.values))

        percent_of_mse_improvment_unc_ppt_10perc = 100 * (stations_with_mse_improvements_unc_ppt_10perc /
                                                          df_improvements.mse_dwd_interp_ppt.shape[0])

        # 20perc
        stations_with_mse_improvements_unc_ppt_20perc = sum(i < j for (i, j) in zip(
            df_improvements.mse_dwd_netatmo_interp_unc_ppt_20perc.values,
            df_improvements.mse_dwd_interp_ppt.values))

        stations_without_mse_improvements_unc_ppt_20perc = sum(i >= j for (i, j) in zip(
            df_improvements.mse_dwd_netatmo_interp_unc_ppt_20perc.values,
            df_improvements.mse_dwd_interp_ppt.values))

        percent_of_mse_improvment_unc_ppt_20perc = 100 * (stations_with_mse_improvements_unc_ppt_20perc /
                                                          df_improvements.mse_dwd_interp_ppt.shape[0])
        #########################################################
        plt.ioff()
        fig = plt.figure(figsize=(24, 12), dpi=150)

        ax = fig.add_subplot(111)

        ax.plot(df_improvements.index,
                df_improvements.mse_dwd_interp_ppt,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation PPT RMSE mean = %0.4f' % np.mean(
                    df_improvements.mse_dwd_interp_ppt.values))
        ax.plot(df_improvements.index,
                df_improvements.mse_dwd_netatmo_interp_ppt,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation PPT RMSE mean = %0.4f' % np.mean(
                    df_improvements.mse_dwd_netatmo_interp_ppt.values))
        ax.plot(df_improvements.index,
                df_improvements.mse_dwd_netatmo_interp_unc_ppt_2perc,
                alpha=.5,
                c='g',  # colors_arr,
                marker='+',
                label='DWD-Netatmo Interpolation Unc 2perc PPT RMSE mean = %0.4f' % np.mean(
                    df_improvements.mse_dwd_netatmo_interp_unc_ppt_2perc.values))

        ax.plot(df_improvements.index,
                df_improvements.mse_dwd_netatmo_interp_unc_ppt_5perc,
                alpha=.5,
                c='m',  # colors_arr,
                marker='1',
                label='DWD-Netatmo Interpolation Unc 5perc PPT RMSE mean = %0.4f' % np.mean(
                    df_improvements.mse_dwd_netatmo_interp_unc_ppt_5perc.values))

        ax.plot(df_improvements.index,
                df_improvements.mse_dwd_netatmo_interp_unc_ppt_10perc,
                alpha=.5,
                c='c',  # colors_arr,
                marker='2',
                label='DWD-Netatmo Interpolation Unc 10perc PPT RMSE mean = %0.4f' % np.mean(
                    df_improvements.mse_dwd_netatmo_interp_unc_ppt_10perc.values))

        ax.plot(df_improvements.index,
                df_improvements.mse_dwd_netatmo_interp_unc_ppt_20perc,
                alpha=.5,
                c='orange',  # colors_arr,
                marker=',',
                label='DWD-Netatmo Interpolation Unc 20perc PPT RMSE mean = %0.4f' % np.mean(
                    df_improvements.mse_dwd_netatmo_interp_unc_ppt_20perc.values))

        ax.set_title('Root mean squared error Interpolated from DWD or DWD-Netatmo %s \n '
                     'PPT of %s Extreme Events \n Stations with Improvemnts %d / %d, Percentage %0.0f'
                     '\n Stations with Improvemnts 2perc Unc %d / %d, Percentage %0.0f'
                     '\n Stations with Improvemnts 5perc Unc %d / %d, Percentage %0.0f'
                     '\n Stations with Improvemnts 10perc Unc %d / %d, Percentage %0.0f'
                     '\n Stations with Improvemnts 20perc Unc %d / %d, Percentage %0.0f'
                     % (_acc_, temp_freq, stations_with_mse_improvements_ppt,
                        df_improvements.mse_dwd_interp_ppt.shape[0], percent_of_mse_improvment_ppt,
                         stations_with_mse_improvements_unc_ppt_2perc,
                        df_improvements.mse_dwd_interp_ppt.shape[0],
                         percent_of_mse_improvment_unc_ppt_2perc,

                         stations_with_mse_improvements_unc_ppt_5perc,
                        df_improvements.mse_dwd_interp_ppt.shape[0],
                         percent_of_mse_improvment_unc_ppt_5perc,

                         stations_with_mse_improvements_unc_ppt_10perc,
                        df_improvements.mse_dwd_interp_ppt.shape[0],
                         percent_of_mse_improvment_unc_ppt_10perc,

                         stations_with_mse_improvements_unc_ppt_20perc,
                        df_improvements.mse_dwd_interp_ppt.shape[0],
                         percent_of_mse_improvment_unc_ppt_20perc))
        ax.grid(alpha=0.25)
        plt.setp(ax.get_xticklabels(), rotation=45)
        #ax.set_yticks(np.arange(0, 1.05, .10))
        ax.set_xlabel('DWD Stations')
        ax.legend(loc='upper right')
        ax.set_ylabel('RMSE')

        plt.savefig((path_to_use / (
            r'_pptrmse_%s_events_dwd_%s_2.png' % (temp_freq, _interp_acc_))),
            frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()


#==============================================================================
# plot_not_filtered
#==============================================================================
if plot_not_filtered:
    for temp_freq in ['60min', '180min', '360min', '720min', '1440min']:  #
        print(temp_freq)

        path_to_Quantiles_netatmo_no_flt___ = main_dir / (
            r'Final_results/Ppt_ok_ok_un_new_netatmo_no_flt___%s' % temp_freq)

        Quantiles_netatmo_no_flt___ = list_all_full_path(
            '.csv', path_to_Quantiles_netatmo_no_flt___)

        #########################################################
        path_dwd_edf_data = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
                             r"\NetAtmo_BW\edf_ppt_all_dwd_%s_.csv"
                             % temp_freq)

        path_to_dwd_ppt = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
                           r"\NetAtmo_BW\ppt_all_dwd_%s_.csv"
                           % temp_freq)

        df_dwd_edf = pd.read_csv(path_dwd_edf_data, sep=';', index_col=0,
                                 parse_dates=True, infer_datetime_format=True)

        # DWD ppt
        path_dwd_ppt_data = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
                             r"\NetAtmo_BW\ppt_all_dwd_%s_.csv" % temp_freq)

        df_dwd_ppt = pd.read_csv(path_dwd_ppt_data, sep=';', index_col=0,
                                 parse_dates=True, infer_datetime_format=True)

        #########################################################
        df_improvements = pd.DataFrame(
            index=df_dwd_edf.columns,
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
                                 sep=';', index_col=0, parse_dates=True,
                                 infer_datetime_format=True)

            df_netatmo_dwd = pd.read_csv(path_interpolated_using_netatmo_dwd,
                                         skiprows=[1],
                                         sep=';', index_col=0, parse_dates=True,
                                         infer_datetime_format=True)

            cmn_interpolated_events = df_netatmo_dwd.index.intersection(
                df_dwd.index).intersection(df_dwd_edf.index)
            print('Total number of events is',
                  cmn_interpolated_events.shape[0])
            for stn_ in df_dwd.columns:
                print(stn_)
                df_compare = pd.DataFrame(index=df_dwd.index)
                for event_date in cmn_interpolated_events:
                    print(event_date)
                    #                 if str(event_date) == '2015-06-07 22:00:00':
                    #                     print(event_date)
                    # interpolated_quantile_netatmo = df_netatmo.loc[event_date, stn_]
                    #event_date = cmn_interpolated_events[0]
                    interpolated_ppt_dwd = df_dwd.loc[event_date, stn_]

                    interpolated_ppt_netatmo_dwd = df_netatmo_dwd.loc[event_date, stn_]

                    # original quantile when transforming ppt to edf
                    original_ppt = df_dwd_ppt.loc[event_date, stn_]

                    df_compare.loc[
                        event_date,
                        'original_ppt'] = original_ppt

                    df_compare.loc[
                        event_date,
                        'interpolated_ppt_dwd'] = interpolated_ppt_dwd
                    # df_compare.loc[event_date,
                    #               'interpolated_quantile_netatmo'] = interpolated_quantile_netatmo
                    df_compare.loc[
                        event_date,
                        'interpolated_ppt_netatmo_dwd'] = interpolated_ppt_netatmo_dwd

                df_compare = df_compare[df_compare >= 0]
                df_compare.dropna(how='any', inplace=True)
                #==============================================================
                # # PPT
                #==============================================================
                values_x_ppt = df_compare['original_ppt'].values
                values_dwd_ppt = df_compare['interpolated_ppt_dwd'].values

                values_netatmo_dwd_ppt = df_compare['interpolated_ppt_netatmo_dwd'].values

                # calculate sqared error between obsv and interpolated
                # ppt
                mse_dwd_interp_ppt = np.sqrt(np.square(
                    np.subtract(values_x_ppt, values_dwd_ppt)).mean())
                mse_dwd_netatmo_interp_ppt = np.sqrt(np.square(
                    np.subtract(values_x_ppt, values_netatmo_dwd_ppt)).mean())

                # calculate correlations (pearson and spearman)

                df_improvements.loc[
                    stn_,
                    'mse_dwd_interp_ppt'] = mse_dwd_interp_ppt
                df_improvements.loc[
                    stn_,
                    'mse_dwd_netatmo_interp_ppt'] = mse_dwd_netatmo_interp_ppt

                values_x = df_compare['original_ppt'].values
                values_dwd = df_compare['interpolated_ppt_dwd'].values
                # values_netatmo =df_compare['interpolated_quantile_netatmo'].values
                values_netatmo_dwd = df_compare['interpolated_ppt_netatmo_dwd'].values

                # calculate sqared error between obsv and interpolated
                # quantiles
                mse_dwd_interp = np.sqrt(np.square(
                    np.subtract(values_x, values_dwd)).mean())
                mse_dwd_netatmo_interp = np.sqrt(np.square(
                    np.subtract(values_x, values_netatmo_dwd)).mean())

                # calculate correlations (pearson and spearman)
                corr_dwd = pears(values_x, values_dwd)[0]
                rho_dwd = spr(values_x, values_dwd)[0]

                #corr_netatmo = pears(values_x, values_netatmo)[0]
                #rho_netatmo = spr(values_x, values_netatmo)[0]

                corr_netatmo_dwd = pears(values_x, values_netatmo_dwd)[0]
                rho_netatmo_dwd = spr(values_x, values_netatmo_dwd)[0]

                df_improvements.loc[stn_, 'pearson_corr_dwd_'] = corr_dwd
                df_improvements.loc[stn_, 'spearman_corr_dwd_'] = rho_dwd
                df_improvements.loc[
                    stn_,
                    'pearson_corr_dwd_netatmo'] = corr_netatmo_dwd
                df_improvements.loc[
                    stn_,
                    'spearman_corr_dwd_netatmo'] = rho_netatmo_dwd

                df_improvements.loc[
                    stn_,
                    'mse_dwd_interp'] = mse_dwd_interp
                df_improvements.loc[
                    stn_,
                    'mse_dwd_netatmo_interp'] = mse_dwd_netatmo_interp

            df_improvements.sort_values(by=['mse_dwd_interp_ppt'],
                                        ascending=True, inplace=True)
            df_improvements.dropna(how='all', inplace=True)

            # In[42]:

            #==================================================================
            # PPT RMSE
            #==================================================================
            stations_with_mse_improvements_ppt = sum(i < j for (i, j) in zip(
                df_improvements.mse_dwd_netatmo_interp_ppt.values,
                df_improvements.mse_dwd_interp_ppt.values))

            stations_without_mse_improvements_ppt = sum(i >= j for (i, j) in zip(
                df_improvements.mse_dwd_netatmo_interp_ppt.values,
                df_improvements.mse_dwd_interp_ppt.values))

            percent_of_mse_improvment_ppt = 100 * (
                stations_with_mse_improvements_ppt /
                df_improvements.mse_dwd_interp_ppt.shape[0])

            ##########################
#             stations_with_mse_improvements_unc_ppt = sum(i < j for (i, j) in zip(
#                 df_improvements.mse_dwd_netatmo_interp_unc_ppt.values,
#                 df_improvements.mse_dwd_interp_ppt.values))
#
#             stations_without_mse_improvements_unc_ppt = sum(i >= j for (i, j) in zip(
#                 df_improvements.mse_dwd_netatmo_interp_unc_ppt.values,
#                 df_improvements.mse_dwd_interp_ppt.values))

#             percent_of_mse_improvment_unc_ppt = 100 * (stations_with_mse_improvements_unc_ppt /
#                                                        df_improvements.mse_dwd_interp_ppt.shape[0])
            #########################################################
            mean_dwd = np.mean(df_improvements.mse_dwd_interp_ppt.values)
            mean_dwd_netatmo = np.mean(
                df_improvements.mse_dwd_netatmo_interp_ppt.values)

            # for plotting
            df_improvements['mse_dwd_netatmo_interp_ppt'].values[
                df_improvements['mse_dwd_netatmo_interp_ppt'] > 700] = 700

            plt.ioff()
            fig = plt.figure(figsize=(24, 12), dpi=150)

            ax = fig.add_subplot(111)

            ax.plot(df_improvements.index,
                    df_improvements.mse_dwd_interp_ppt,
                    alpha=.8,
                    c='b',  # colors_arr,
                    marker='d',
                    label='DWD Interpolation PPT RMSE mean = %0.4f' % mean_dwd)
            ax.plot(df_improvements.index,
                    df_improvements.mse_dwd_netatmo_interp_ppt,
                    alpha=.8,
                    c='r',  # colors_arr,
                    marker='*',
                    label='DWD-Netatmo Interpolation PPT RMSE mean = %0.4f'
                    % mean_dwd_netatmo)

            ax.set_title(
                'Root mean squared error Interpolated from DWD or DWD-Netatmo %s \n '
                'PPT of %s Extreme Events \n Stations with Improvemnts %d / %d, Percentage %0.0f'
                #                          '\n Stations with Improvemnts Unc %d / %d, Percentage %0.0f'
                % (_acc_, temp_freq, stations_with_mse_improvements_ppt,
                   df_improvements.mse_dwd_interp_ppt.shape[0],
                    percent_of_mse_improvment_ppt))
            ax.grid(alpha=0.25)
            plt.setp(ax.get_xticklabels(), rotation=45)
            #ax.set_yticks(np.arange(0, 1.05, .10))
            ax.set_xlabel('DWD Stations')
            ax.legend(loc='upper right')
            ax.set_ylabel('RMSE')

            plt.savefig((path_to_use / (
                r'_pptrmse_%s_events_dwd_%s_2.png' % (temp_freq, _interp_acc_))),
                frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
            plt.close()

#             #==================================================================
#             # QUANTILES RMSE
#             #==================================================================
#
#             # sum(i > j for (i, j) in zip(df_improvements.pearson_corr_dwd_netatmo.values,
#             #                            df_improvements.pearson_corr_dwd_.values))
#             dwd_ids_mse_no_improvement_pearson = pd.DataFrame(
#                 data=df_improvements.index[np.where(
#                     df_improvements.mse_dwd_netatmo_interp.values > df_improvements.mse_dwd_interp.values)],
#                 columns=['stns_no_impv_mse'])
#
#             # dwd_ids_no_improvement_pearson.to_csv(
#             #    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\stns_no_impv_pearson_corr_daily_events_all.csv',sep=';')
#
#             # dwd_ids_no_improvement_spr.to_csv(
#             #    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\stns_no_impv_spearman_corr_daily_events_all.csv',sep=';')
#
#             # In[43]:
#
#             stations_with_mse_improvements = sum(i < j for (i, j) in zip(
#                 df_improvements.mse_dwd_netatmo_interp.values,
#                 df_improvements.mse_dwd_interp.values))
#
#             stations_without_mse_improvements = sum(i >= j for (i, j) in zip(
#                 df_improvements.mse_dwd_netatmo_interp.values,
#                 df_improvements.mse_dwd_interp.values))
#
#             percent_of_mse_improvment = 100 * (stations_with_mse_improvements /
#                                                df_improvements.mse_dwd_interp.shape[0])
#
#             ##########################
#
#             #########################################################
#             plt.ioff()
#             fig = plt.figure(figsize=(24, 12), dpi=150)
#
#             ax = fig.add_subplot(111)
#
#             # ax.scatter(df_improvements.index,
#             #           df_improvements.pearson_corr_dwd_,
#             #           alpha=.8,
#             #           c='r',  # colors_arr,
#             #           s=15,
#             #           marker='d',
#             # cmap=plt.get_cmap('viridis'),
#             #          label='DWD Interpolation')
#
#             # ax.scatter(df_improvements.index,
#             #           df_improvements.pearson_corr_dwd_netatmo,
#             #           alpha=.8,
#             #           c='b',  # colors_arr,
#             #           s=15,
#             #           marker='*',
#             # cmap=plt.get_cmap('viridis'),
#             #           label='DWD-Netatmo Interpolation')
#
#             ax.plot(df_improvements.index,
#                     df_improvements.mse_dwd_interp,
#                     alpha=.8,
#                     c='b',  # colors_arr,
#                     marker='d',
#                     label='DWD Interpolation RMSE mean = %0.4f' % np.mean(
#                         df_improvements.mse_dwd_interp.values))
#             ax.plot(df_improvements.index,
#                     df_improvements.mse_dwd_netatmo_interp,
#                     alpha=.8,
#                     c='r',  # colors_arr,
#                     marker='*',
#                     label='DWD-Netatmo Interpolation RMSE mean = %0.4f' % np.mean(
#                         df_improvements.mse_dwd_netatmo_interp.values))
#
#             ax.set_title('Root mean squared error Interpolated from DWD or DWD-Netatmo %s \n '
#                          'Rainfall of %s Extreme Events \n Stations with Improvemnts %d / %d, Percentage %0.0f'
#
#                          % (_acc_, temp_freq, stations_with_mse_improvements,
#                             df_improvements.mse_dwd_interp.shape[0], percent_of_mse_improvment))
#             ax.grid(alpha=0.25)
#             plt.setp(ax.get_xticklabels(), rotation=45)
#             #ax.set_yticks(np.arange(0, 1.05, .10))
#             ax.set_xlabel('DWD Stations')
#             ax.legend(loc='upper right')
#             ax.set_ylabel('RMSE')
#
#             plt.savefig((path_to_use / (
#                 r'rmse_%s_events_dwd_%s_.png' % (temp_freq, _interp_acc_))),
#                 frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
#             plt.close()
#
#             #==================================================================
#             #
#             #==================================================================
#             dwd_ids_no_improvement_pearson = pd.DataFrame(
#                 data=df_improvements.index[np.where(
#                     df_improvements.pearson_corr_dwd_netatmo.values < df_improvements.pearson_corr_dwd_.values)],
#                 columns=['stns_no_impv_pearson_corr'])
#             dwd_ids_no_improvement_spr = pd.DataFrame(data=df_improvements.index[np.where(
#                 df_improvements.spearman_corr_dwd_netatmo.values < df_improvements.spearman_corr_dwd_.values)],
#                 columns=['stns_no_impv_spearman_corr'])
#
#             # dwd_ids_no_improvement_pearson.to_csv(
#             #    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\stns_no_impv_pearson_corr_daily_events_all.csv',sep=';')
#
#             # dwd_ids_no_improvement_spr.to_csv(
#             #    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\stns_no_impv_spearman_corr_daily_events_all.csv',sep=';')
#
#             # In[43]:
#
#             stations_with_improvements = sum(i >= j for (i, j) in zip(
#                 df_improvements.pearson_corr_dwd_netatmo.values,
#                 df_improvements.pearson_corr_dwd_.values))
#
#             stations_without_improvements = sum(i < j for (i, j) in zip(
#                 df_improvements.pearson_corr_dwd_netatmo.values,
#                 df_improvements.pearson_corr_dwd_.values))
#
#             percent_of_improvment = 100 * (stations_with_improvements /
#                                            df_improvements.pearson_corr_dwd_netatmo.shape[0])
#
#             #########################################################
#             plt.ioff()
#             fig = plt.figure(figsize=(24, 12), dpi=150)
#
#             ax = fig.add_subplot(111)
#
#             ax.plot(df_improvements.index,
#                     df_improvements.pearson_corr_dwd_,
#                     alpha=.8,
#                     c='b',  # colors_arr,
#                     marker='d',
#                     label='DWD Interpolation %0.4f' %
#                     df_improvements.pearson_corr_dwd_.mean())
#
#             ax.plot(df_improvements.index,
#                     df_improvements.pearson_corr_dwd_netatmo,
#                     alpha=.8,
#                     c='r',  # colors_arr,
#                     marker='*',
#                     label='DWD-Netatmo Interpolation %0.4f' %
#                     df_improvements.pearson_corr_dwd_netatmo.mean())
#
#             ax.set_title('Pearson Correlation Interpolated from DWD or DWD-Netatmo %s \n '
#                          'Rainfall of %s Extreme Events \n Stations with Improvemnts %d / %d, Percentage %0.0f'
#                          % (_acc_, temp_freq, stations_with_improvements,
#                             df_improvements.pearson_corr_dwd_netatmo.shape[0],
#                              percent_of_improvment))
#             ax.grid(alpha=0.25)
#             plt.setp(ax.get_xticklabels(), rotation=45)
#             ax.set_yticks(np.arange(0, 1.05, .10))
#             ax.set_xlabel('DWD Stations')
#             ax.legend(loc='upper right')
#             ax.set_ylabel('Pearson Correlation')
#
#             plt.savefig((path_to_use / (
#                 r'rmse_pears_corr_%s_events_dwd_%s_.png' % (temp_freq, _interp_acc_))),
#                 frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
#             plt.close()
#
#             # In[44]:
#
#             stations_with_improvements = sum(i >= j for (i, j) in zip(
#                 df_improvements.spearman_corr_dwd_netatmo.values,
#                 df_improvements.spearman_corr_dwd_.values))
#
#             percent_of_improvment = 100 * (stations_with_improvements /
#                                            df_improvements.spearman_corr_dwd_netatmo.shape[0])
#             #########################################################
#             plt.ioff()
#             fig = plt.figure(figsize=(24, 12), dpi=150)
#
#             ax = fig.add_subplot(111)
#
#             ax.plot(df_improvements.index,
#                     df_improvements.spearman_corr_dwd_,
#                     alpha=.8,
#                     c='b',  # colors_arr,
#                     marker='d',
#                     label='DWD Interpolation %0.4f' %
#                     df_improvements.spearman_corr_dwd_.mean())
#             ax.plot(df_improvements.index,
#                     df_improvements.spearman_corr_dwd_netatmo,
#                     alpha=.8,
#                     c='r',  # colors_arr,
#                     marker='*',
#                     label='DWD-Netatmo Interpolation %0.4f' %
#                     df_improvements.spearman_corr_dwd_netatmo.mean())
#
#             ax.set_title('Spearman Correlation Interpolated from DWD or DWD-Netatmo %s\n '
#                          'Rainfall of %s Extreme Events \n Stations with Improvemnts %d / %d, Percentage %0.0f'
#
#                          % (_acc_, temp_freq, stations_with_improvements,
#                             df_improvements.spearman_corr_dwd_netatmo.shape[0],
#                              percent_of_improvment
#                             ))
#
#             plt.setp(ax.get_xticklabels(), rotation=45)
#             ax.grid(alpha=0.25)
#             ax.set_yticks(np.arange(0, 1.05, .10))
#             ax.set_xlabel('DWD Stations')
#             ax.legend(loc='upper right')
#             ax.set_ylabel('Spearman Correlation')
#
#             plt.savefig((path_to_use / (
#                 r'rmse_spr_corr_%s_events_dwd_%s_.png' % (temp_freq, _interp_acc_))),
#                 frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
#             plt.close()
