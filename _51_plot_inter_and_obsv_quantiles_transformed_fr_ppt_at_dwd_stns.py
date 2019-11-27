
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

plot_not_filtered = False
plot_filtered = True
#==============================================================================
#
#==============================================================================

if plot_filtered:
    for temp_freq in ['60min', '360min', '720min', '1440min']:
        # '60min', '360min',  '1440min'
        print(temp_freq)

        path_to_Qt_ok_un_first_flt__temp_flt_1st_ = main_dir / (
            r'Ppt_ok_ok_un__first_flt__temp_flt__1st_%s' % temp_freq)
        Qt_ok_un_first_flt__temp_flt_1st_ = list_all_full_path(
            '.csv', path_to_Qt_ok_un_first_flt__temp_flt_1st_)

        path_to_Qt_ok_un_first_flt__temp_flt_comb_ = main_dir / (
            r'Ppt_ok_ok_un__first_flt__temp_flt__comb_%s' % temp_freq)
        Qt_ok_un_first_flt__temp_flt_comb_ = list_all_full_path(
            '.csv', path_to_Qt_ok_un_first_flt__temp_flt_comb_)

        path_to_Qt_ok_un_first_flt_1st_ = main_dir / (
            r'Ppt_ok_ok_un__first_flt__1st_%s' % temp_freq)
        Qt_ok_un_first_flt_1st_ = list_all_full_path(
            '.csv', path_to_Qt_ok_un_first_flt_1st_)

        path_to_Qt_ok_un_first_flt_comb_ = main_dir / (
            r'Ppt_ok_ok_un__first_flt__comb_%s' % temp_freq)
        Qt_ok_un_first_flt_comb_ = list_all_full_path(
            '.csv', path_to_Qt_ok_un_first_flt_comb_)

        path_to_Quantiles_netatmo_no_flt___ = Path(
            r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
            r'\oridinary_kriging_compare_DWD_Netatmo\Ppt_ok_ok_un_netatmo_no_flt__%s'
            % temp_freq)

        Quantiles_netatmo_no_flt___ = list_all_full_path(
            '.csv', path_to_Quantiles_netatmo_no_flt___)

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
                     'pearson_corr_dwd_netatmo_unc',
                     'spearman_corr_dwd_netatmo_unc'])

        #########################################################
        path_interpolated_using_dwd_list = []
        path_interpolated_using_netatmo_dwd_list = []
        path_interpolated_using_netatmo_dwd_list_un = []

        path_to_use = path_to_Qt_ok_un_first_flt__temp_flt_comb_
        data_to_use = Qt_ok_un_first_flt__temp_flt_comb_

        _interp_acc_ = str(r'%s' % (str(path_to_use).split('\\')[-1]))
        # for i in range(12):
        #   i = int(i)
        try:
            for df_file in data_to_use:

                if ('using_dwd_only_grp_') in df_file:
                    print(df_file)
                    path_interpolated_using_dwd_list.append(df_file)
                if ('using_dwd_netamo_grp_') in df_file and 'unc' not in df_file:

                    print(df_file)
                    path_interpolated_using_netatmo_dwd_list.append(df_file)

                if ('using_dwd_netamo_grp_') in df_file and 'unc' in df_file:
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

            print('Total number of events is',
                  cmn_interpolated_events.shape[0])
            for stn_ in df_dwd.columns:
                print(stn_)
                #stn_ = df_dwd.columns[0]
                for event_date in cmn_interpolated_events:
                    # print(event_date)
                    #                 if str(event_date) == '2015-06-07 22:00:00':
                    #                     print(event_date)
                    # interpolated_quantile_netatmo = df_netatmo.loc[event_date, stn_]
                    #event_date = cmn_interpolated_events[0]
                    interpolated_ppt_dwd = df_dwd.loc[event_date, stn_]

                    interpolated_ppt_netatmo_dwd = df_netatmo_dwd.loc[event_date, stn_]

                    interpolated_ppt_netatmo_dwd_unc = df_netatmo_dwd_unc.loc[event_date, stn_]
                    # original quantile when transforming ppt to edf
                    original_ppt = df_dwd_ppt.loc[event_date, stn_]

                    df_compare.loc[event_date,
                                   'original_ppt'] = original_ppt

                    df_compare.loc[event_date,
                                   'interpolated_ppt_dwd'] = interpolated_ppt_dwd
                    # df_compare.loc[event_date,
                    #               'interpolated_quantile_netatmo'] = interpolated_quantile_netatmo
                    df_compare.loc[event_date,
                                   'interpolated_ppt_netatmo_dwd'] = interpolated_ppt_netatmo_dwd
                    df_compare.loc[event_date,
                                   'interpolated_ppt_netatmo_dwd_unc'] = interpolated_ppt_netatmo_dwd_unc
                    if original_ppt >= 0:

                        #
                        original_ppt_stn = df_dwd_ppt.loc[
                            :,
                            stn_].dropna().values.ravel()

                        ppt_stn_dist, edf_stn_dist = build_edf_fr_vals(
                            original_ppt_stn)

                        # plt.plot(ppt_stn, edf_stn_dist)
                        edf_orig_ = df_dwd_edf.loc[event_date, stn_]

                        # interpolated from DWD
                        edf_interp_fr_dwd = edf_stn_dist[
                            np.where(ppt_stn_dist == find_nearest(
                                ppt_stn_dist,
                                interpolated_ppt_dwd))][0]

                        # interpolated from DWD-Netatmo
                        edf_interp_fr_dwd_netatmo = edf_stn_dist[
                            np.where(ppt_stn_dist == find_nearest(
                                ppt_stn_dist,
                                interpolated_ppt_netatmo_dwd))][0]

                        # interpolated from DWD-Netatmo Unc
                        edf_interp_fr_dwd_netatmo_unc = edf_stn_dist[
                            np.where(ppt_stn_dist == find_nearest(
                                ppt_stn_dist,
                                interpolated_ppt_netatmo_dwd_unc))][0]

                        if edf_orig_ >= 0.5:
                            # print(original_quantile)
                            df_compare.loc[event_date,
                                           'original_edf'] = edf_orig_
                            df_compare.loc[event_date,
                                           'interpolated_edf_dwd'] = edf_interp_fr_dwd
                            # df_compare.loc[event_date,
                            #               'interpolated_quantile_netatmo'] = interpolated_quantile_netatmo
                            df_compare.loc[event_date,
                                           'interpolated_edf_netatmo_dwd'] = edf_interp_fr_dwd_netatmo
                            df_compare.loc[event_date,
                                           'interpolated_edf_netatmo_dwd_unc'] = edf_interp_fr_dwd_netatmo_unc

                        else:
                            df_compare.loc[event_date,
                                           'original_edf'] = np.nan
                            df_compare.loc[event_date,
                                           'interpolated_edf_dwd'] = np.nan
                            # df_compare.loc[event_date,
                            #               'interpolated_quantile_netatmo'] = np.nan
                            df_compare.loc[event_date,
                                           'interpolated_edf_netatmo_dwd'] = np.nan
                            df_compare.loc[event_date,
                                           'interpolated_edf_netatmo_dwd_unc'] = np.nan
                    else:
                        df_compare.loc[event_date,
                                       'original_edf'] = np.nan
                        df_compare.loc[event_date,
                                       'interpolated_edf_dwd'] = np.nan
                        # df_compare.loc[event_date,
                        #               'interpolated_quantile_netatmo'] = np.nan
                        df_compare.loc[event_date,
                                       'interpolated_edf_netatmo_dwd'] = np.nan
                        df_compare.loc[event_date,
                                       'interpolated_edf_netatmo_dwd_unc'] = np.nan

                df_compare = df_compare[df_compare >= 0]
                df_compare.dropna(how='any', inplace=True)
                #==============================================================
                # # PPT
                #==============================================================
                values_x_ppt = df_compare['original_ppt'].values
                values_dwd_ppt = df_compare['interpolated_ppt_dwd'].values

                values_netatmo_dwd_ppt = df_compare['interpolated_ppt_netatmo_dwd'].values
                values_netatmo_dwd_unc_ppt = df_compare['interpolated_ppt_netatmo_dwd_unc'].values

                # calculate sqared error between obsv and interpolated
                # ppt
                mse_dwd_interp_ppt = np.square(
                    np.subtract(values_x_ppt, values_dwd_ppt)).mean()
                mse_dwd_netatmo_interp_ppt = np.square(
                    np.subtract(values_x_ppt, values_netatmo_dwd_ppt)).mean()
                mse_dwd_netatmo_interp_unc_ppt = np.square(
                    np.subtract(values_x_ppt, values_netatmo_dwd_unc_ppt)).mean()

                # calculate correlations (pearson and spearman)

                df_improvements.loc[stn_,
                                    'mse_dwd_interp_ppt'] = mse_dwd_interp_ppt
                df_improvements.loc[stn_,
                                    'mse_dwd_netatmo_interp_ppt'] = mse_dwd_netatmo_interp_ppt
                df_improvements.loc[stn_,
                                    'mse_dwd_netatmo_interp_unc_ppt'] = mse_dwd_netatmo_interp_unc_ppt

                #==============================================================
                # # Quantiles
                #==============================================================
                values_x = df_compare['original_edf'].values
                values_dwd = df_compare['interpolated_edf_dwd'].values
                # values_netatmo =df_compare['interpolated_quantile_netatmo'].values
                values_netatmo_dwd = df_compare['interpolated_edf_netatmo_dwd'].values
                values_netatmo_dwd_unc = df_compare['interpolated_edf_netatmo_dwd_unc'].values

                # calculate sqared error between obsv and interpolated
                # quantiles
                mse_dwd_interp = np.square(
                    np.subtract(values_x, values_dwd)).mean()
                mse_dwd_netatmo_interp = np.square(
                    np.subtract(values_x, values_netatmo_dwd)).mean()
                mse_dwd_netatmo_interp_unc = np.square(
                    np.subtract(values_x, values_netatmo_dwd_unc)).mean()

                # calculate correlations (pearson and spearman)
                corr_dwd = pears(values_x, values_dwd)[0]
                rho_dwd = spr(values_x, values_dwd)[0]

                #corr_netatmo = pears(values_x, values_netatmo)[0]
                #rho_netatmo = spr(values_x, values_netatmo)[0]

                corr_netatmo_dwd = pears(values_x, values_netatmo_dwd)[0]
                rho_netatmo_dwd = spr(values_x, values_netatmo_dwd)[0]

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

                df_improvements.loc[stn_,
                                    'mse_dwd_interp'] = mse_dwd_interp
                df_improvements.loc[stn_,
                                    'mse_dwd_netatmo_interp'] = mse_dwd_netatmo_interp
                df_improvements.loc[stn_,
                                    'mse_dwd_netatmo_interp_unc'] = mse_dwd_netatmo_interp_unc

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
        stations_with_mse_improvements_unc_ppt = sum(i < j for (i, j) in zip(
            df_improvements.mse_dwd_netatmo_interp_unc_ppt.values,
            df_improvements.mse_dwd_interp_ppt.values))

        stations_without_mse_improvements_unc_ppt = sum(i >= j for (i, j) in zip(
            df_improvements.mse_dwd_netatmo_interp_unc_ppt.values,
            df_improvements.mse_dwd_interp_ppt.values))

        percent_of_mse_improvment_unc_ppt = 100 * (stations_with_mse_improvements_unc_ppt /
                                                   df_improvements.mse_dwd_interp_ppt.shape[0])

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
                df_improvements.mse_dwd_interp_ppt,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation PPT RMSE mean = %0.3f' % np.mean(
                    df_improvements.mse_dwd_interp_ppt.values))
        ax.plot(df_improvements.index,
                df_improvements.mse_dwd_netatmo_interp_ppt,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation PPT RMSE mean = %0.3f' % np.mean(
                    df_improvements.mse_dwd_netatmo_interp_ppt.values))
        ax.plot(df_improvements.index,
                df_improvements.mse_dwd_netatmo_interp_unc_ppt,
                alpha=.8,
                c='g',  # colors_arr,
                marker='+',
                label='DWD-Netatmo Interpolation Unc PPT RMSE mean = %0.3f' % np.mean(
                    df_improvements.mse_dwd_netatmo_interp_unc_ppt.values))

        ax.set_title('Root mean squared error Interpolated from DWD or DWD-Netatmo %s \n '
                     'PPT of %s Extreme Events \n Stations with Improvemnts %d / %d, Percentage %0.0f'
                     '\n Stations with Improvemnts Unc %d / %d, Percentage %0.0f'
                     % (_acc_, temp_freq, stations_with_mse_improvements_ppt,
                        df_improvements.mse_dwd_interp_ppt.shape[0], percent_of_mse_improvment_ppt,
                         stations_with_mse_improvements_unc_ppt,
                        df_improvements.mse_dwd_interp_ppt.shape[0],
                         percent_of_mse_improvment_unc_ppt))
        ax.grid(alpha=0.25)
        plt.setp(ax.get_xticklabels(), rotation=45)
        #ax.set_yticks(np.arange(0, 1.05, .10))
        ax.set_xlabel('DWD Stations')
        ax.legend(loc='upper right')
        ax.set_ylabel('RMSE')

        plt.savefig((path_to_use / (
            r'_pptrmse_%s_events_dwd_%s_.png' % (temp_freq, _interp_acc_))),
            frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

        #======================================================================
        # QUANTILES
        #======================================================================

        # sum(i > j for (i, j) in zip(df_improvements.pearson_corr_dwd_netatmo.values,
        #                            df_improvements.pearson_corr_dwd_.values))
        dwd_ids_mse_no_improvement_pearson = pd.DataFrame(
            data=df_improvements.index[np.where(
                df_improvements.mse_dwd_netatmo_interp.values > df_improvements.mse_dwd_interp.values)],
            columns=['stns_no_impv_mse'])

        # dwd_ids_no_improvement_pearson.to_csv(
        #    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\stns_no_impv_pearson_corr_daily_events_all.csv',sep=';')

        # dwd_ids_no_improvement_spr.to_csv(
        #    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\stns_no_impv_spearman_corr_daily_events_all.csv',sep=';')

        # In[43]:

        stations_with_mse_improvements = sum(i < j for (i, j) in zip(
            df_improvements.mse_dwd_netatmo_interp.values,
            df_improvements.mse_dwd_interp.values))

        stations_without_mse_improvements = sum(i >= j for (i, j) in zip(
            df_improvements.mse_dwd_netatmo_interp.values,
            df_improvements.mse_dwd_interp.values))

        percent_of_mse_improvment = 100 * (stations_with_mse_improvements /
                                           df_improvements.mse_dwd_interp.shape[0])

        ##########################
        stations_with_mse_improvements_unc = sum(i < j for (i, j) in zip(
            df_improvements.mse_dwd_netatmo_interp_unc.values,
            df_improvements.mse_dwd_interp.values))

        stations_without_mse_improvements_unc = sum(i >= j for (i, j) in zip(
            df_improvements.mse_dwd_netatmo_interp_unc.values,
            df_improvements.mse_dwd_interp.values))

        percent_of_mse_improvment_unc = 100 * (stations_with_mse_improvements_unc /
                                               df_improvements.mse_dwd_interp.shape[0])

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
                df_improvements.mse_dwd_interp,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation RMSE mean = %0.3f' % np.mean(
                    df_improvements.mse_dwd_interp.values))
        ax.plot(df_improvements.index,
                df_improvements.mse_dwd_netatmo_interp,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation RMSE mean = %0.3f' % np.mean(
                    df_improvements.mse_dwd_netatmo_interp.values))
        ax.plot(df_improvements.index,
                df_improvements.mse_dwd_netatmo_interp_unc,
                alpha=.8,
                c='g',  # colors_arr,
                marker='+',
                label='DWD-Netatmo Interpolation Unc RMSE mean = %0.3f' % np.mean(
                    df_improvements.mse_dwd_netatmo_interp_unc.values))

        ax.set_title('Root mean squared error Interpolated from DWD or DWD-Netatmo %s \n '
                     'Quantiles of %s Extreme Events \n Stations with Improvemnts %d / %d, Percentage %0.0f'
                     '\n Stations with Improvemnts Unc %d / %d, Percentage %0.0f'
                     % (_acc_, temp_freq, stations_with_mse_improvements,
                        df_improvements.mse_dwd_interp.shape[0], percent_of_mse_improvment,
                         stations_with_mse_improvements_unc,
                        df_improvements.mse_dwd_interp.shape[0],
                         percent_of_mse_improvment_unc))
        ax.grid(alpha=0.25)
        plt.setp(ax.get_xticklabels(), rotation=45)
        #ax.set_yticks(np.arange(0, 1.05, .10))
        ax.set_xlabel('DWD Stations')
        ax.legend(loc='upper right')
        ax.set_ylabel('RMSE')

        plt.savefig((path_to_use / (
            r'rmse_%s_events_dwd_%s_.png' % (temp_freq, _interp_acc_))),
            frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

        #======================================================================
        #
        #======================================================================
        dwd_ids_no_improvement_pearson = pd.DataFrame(
            data=df_improvements.index[np.where(
                df_improvements.pearson_corr_dwd_netatmo.values < df_improvements.pearson_corr_dwd_.values)],
            columns=['stns_no_impv_pearson_corr'])
        dwd_ids_no_improvement_spr = pd.DataFrame(data=df_improvements.index[np.where(
            df_improvements.spearman_corr_dwd_netatmo.values < df_improvements.spearman_corr_dwd_.values)],
            columns=['stns_no_impv_spearman_corr'])

        dwd_ids_no_improvement_pearson_unc = pd.DataFrame(
            data=df_improvements.index[np.where(
                df_improvements.pearson_corr_dwd_netatmo_unc.values < df_improvements.pearson_corr_dwd_.values)],
            columns=['stns_no_impv_pearson_corr_unc'])
        dwd_ids_no_improvement_spr_unc = pd.DataFrame(data=df_improvements.index[np.where(
            df_improvements.spearman_corr_dwd_netatmo_unc.values < df_improvements.spearman_corr_dwd_.values)],
            columns=['stns_no_impv_spearman_corr_unc'])

        # dwd_ids_no_improvement_pearson.to_csv(
        #    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\stns_no_impv_pearson_corr_daily_events_all.csv',sep=';')

        # dwd_ids_no_improvement_spr.to_csv(
        #    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\stns_no_impv_spearman_corr_daily_events_all.csv',sep=';')

        # In[43]:

        stations_with_improvements = sum(i >= j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo.values,
            df_improvements.pearson_corr_dwd_.values))

        stations_without_improvements = sum(i < j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo.values,
            df_improvements.pearson_corr_dwd_.values))

        percent_of_improvment = 100 * (stations_with_improvements /
                                       df_improvements.pearson_corr_dwd_netatmo.shape[0])

        stations_with_improvements_unc = sum(i >= j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo_unc.values,
            df_improvements.pearson_corr_dwd_.values))

        stations_without_improvements_unc = sum(i < j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo_unc.values,
            df_improvements.pearson_corr_dwd_.values))

        percent_of_improvment_unc = 100 * (stations_with_improvements_unc /
                                           df_improvements.pearson_corr_dwd_netatmo.shape[0])
        #########################################################
        plt.ioff()
        fig = plt.figure(figsize=(24, 12), dpi=150)

        ax = fig.add_subplot(111)

        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation %0.3f' %
                df_improvements.pearson_corr_dwd_.mean())

        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_netatmo,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %0.3f' %
                df_improvements.pearson_corr_dwd_netatmo.mean())

        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_netatmo_unc,
                alpha=.8,
                c='g',  # colors_arr,
                marker='+',
                label='DWD-Netatmo Interpolation Unc %0.3f' %
                df_improvements.pearson_corr_dwd_netatmo_unc.mean())

        ax.set_title('Pearson Correlation Interpolated from DWD or DWD-Netatmo %s \n '
                     'Quantiles of %s Extreme Events \n Stations with Improvemnts %d / %d, Percentage %0.0f'
                     '\n Stations with Improvemnts Unc %d / %d, Percentage Unc %0.0f'
                     % (_acc_, temp_freq, stations_with_improvements,
                        df_improvements.pearson_corr_dwd_netatmo.shape[0], percent_of_improvment,
                         stations_with_improvements_unc,
                        df_improvements.pearson_corr_dwd_netatmo_unc.shape[0],
                         percent_of_improvment_unc))
        ax.grid(alpha=0.25)
        plt.setp(ax.get_xticklabels(), rotation=45)
        ax.set_yticks(np.arange(0, 1.05, .10))
        ax.set_xlabel('DWD Stations')
        ax.legend(loc='upper right')
        ax.set_ylabel('Pearson Correlation')

        plt.savefig((path_to_use / (
            r'rmse_pears_corr_%s_events_dwd_%s_.png' % (temp_freq, _interp_acc_))),
            frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

        # In[44]:

        stations_with_improvements = sum(i >= j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo.values,
            df_improvements.spearman_corr_dwd_.values))

        percent_of_improvment = 100 * (stations_with_improvements /
                                       df_improvements.spearman_corr_dwd_netatmo.shape[0])

        stations_with_improvements_unc = sum(i >= j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo_unc.values,
            df_improvements.spearman_corr_dwd_.values))

        percent_of_improvment_unc = 100 * (stations_with_improvements_unc /
                                           df_improvements.spearman_corr_dwd_netatmo.shape[0])
        #########################################################
        plt.ioff()
        fig = plt.figure(figsize=(24, 12), dpi=150)

        ax = fig.add_subplot(111)

        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation %0.3f' %
                df_improvements.spearman_corr_dwd_.mean())
        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_netatmo,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %0.3f' %
                df_improvements.spearman_corr_dwd_netatmo.mean())

        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_netatmo_unc,
                alpha=.8,
                c='g',  # colors_arr,
                marker='+',
                label='DWD-Netatmo Interpolation Unc %0.3f' %
                df_improvements.spearman_corr_dwd_netatmo_unc.mean())

        ax.set_title('Spearman Correlation Interpolated from DWD or DWD-Netatmo %s\n '
                     'Quantiles of %s Extreme Events \n Stations with Improvemnts %d / %d, Percentage %0.0f'
                     '\n Stations with Improvemnts Unc %d / %d, Percentage Unc %0.0f'
                     % (_acc_, temp_freq, stations_with_improvements,
                        df_improvements.spearman_corr_dwd_netatmo.shape[0], percent_of_improvment,
                        stations_with_improvements_unc,
                        df_improvements.spearman_corr_dwd_netatmo_unc.shape[0],
                         percent_of_improvment_unc))

        plt.setp(ax.get_xticklabels(), rotation=45)
        ax.grid(alpha=0.25)
        ax.set_yticks(np.arange(0, 1.05, .10))
        ax.set_xlabel('DWD Stations')
        ax.legend(loc='upper right')
        ax.set_ylabel('Spearman Correlation')

        plt.savefig((path_to_use / (
            r'rmse_spr_corr_%s_events_dwd_%s_.png' % (temp_freq, _interp_acc_))),
            frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

        # In[43]:

#==============================================================================
# plot_not_filtered
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
        path_dwd_ppt_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\ppt_all_dwd_%s_.csv" % temp_freq

        df_dwd_ppt = pd.read_csv(path_dwd_ppt_data, sep=';', index_col=0,
                                 parse_dates=True, infer_datetime_format=True)

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
                                 sep=';', index_col=0, parse_dates=True,
                                 infer_datetime_format=True)

            df_netatmo_dwd = pd.read_csv(path_interpolated_using_netatmo_dwd, skiprows=[1],
                                         sep=';', index_col=0, parse_dates=True,
                                         infer_datetime_format=True)

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
                    # print(event_date)
                    #                 if str(event_date) == '2015-06-07 22:00:00':
                    #                     print(event_date)
                    # interpolated_quantile_netatmo = df_netatmo.loc[event_date, stn_]
                    #event_date = cmn_interpolated_events[0]
                    interpolated_quantile_dwd = df_dwd.loc[event_date, stn_]

                    interpolated_quantile_netatmo_dwd = df_netatmo_dwd.loc[event_date, stn_]

                    # original quantile when transforming ppt to edf
                    original_quantile = df_dwd_edf.loc[event_date, stn_]

                    if original_quantile >= 0:
                        # calculate p0, since used when generating original
                        # quantiles

                         # original data of ppt
                        original_ppt_stn = df_dwd_ppt.loc[:,
                                                          stn_].dropna().values.ravel()
                        p0_stn = calculate_probab_ppt_below_thr(
                            original_ppt_stn, 0)

                        if np.isclose(interpolated_quantile_dwd,
                                      p0_stn / 2, atol=0.1):
                            ppt_interp_fr_dwd = 0

                        if np.isclose(interpolated_quantile_netatmo_dwd,
                                      p0_stn / 2, atol=0.1):
                            ppt_interp_fr_dwd_netatmo = 0

                        if np.isclose(original_quantile,
                                      p0_stn / 2, atol=0.1):
                            ppt_orig_fr_edf = 0

                        else:

                            ppt_stn_dist, edf_stn_dist = build_edf_fr_vals(
                                original_ppt_stn)

                            edf_orig_fr_dist = find_nearest(
                                edf_stn_dist, original_quantile)

                            # plt.plot(ppt_stn, edf_stn_dist)
                            ppt_orig_ = df_dwd_ppt.loc[event_date, stn_]

                            # interpolated from DWD
                            ppt_interp_fr_dwd = ppt_stn_dist[
                                np.where(edf_stn_dist == find_nearest(
                                    edf_stn_dist,
                                    interpolated_quantile_dwd))][0]

                            # interpolated from DWD-Netatmo
                            ppt_interp_fr_dwd_netatmo = ppt_stn_dist[
                                np.where(edf_stn_dist == find_nearest(
                                    edf_stn_dist,
                                    interpolated_quantile_netatmo_dwd))][0]

                            if ppt_orig_ >= 0.2:
                                # print(original_quantile)
                                df_compare.loc[event_date,
                                               'original_ppt'] = ppt_orig_
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

                # calculate sqared error between obsv and interpolated rainfall
                mse_dwd_interp = np.square(
                    np.subtract(values_x, values_dwd)).mean()
                mse_dwd_netatmo_interp = np.square(
                    np.subtract(values_x, values_netatmo_dwd)).mean()

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

                df_improvements.loc[stn_,
                                    'mse_dwd_interp'] = mse_dwd_interp
                df_improvements.loc[stn_,
                                    'mse_dwd_netatmo_interp'] = mse_dwd_netatmo_interp

            df_improvements.dropna(how='all', inplace=True)

            # In[42]:

            # sum(i > j for (i, j) in zip(df_improvements.pearson_corr_dwd_netatmo.values,
            #                            df_improvements.pearson_corr_dwd_.values))
            dwd_ids_mse_no_improvement_pearson = pd.DataFrame(
                data=df_improvements.index[np.where(
                    df_improvements.mse_dwd_netatmo_interp.values > df_improvements.mse_dwd_interp.values)],
                columns=['stns_no_impv_mse'])

            # dwd_ids_no_improvement_pearson.to_csv(
            #    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\stns_no_impv_pearson_corr_daily_events_all.csv',sep=';')

            # dwd_ids_no_improvement_spr.to_csv(
            #    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\stns_no_impv_spearman_corr_daily_events_all.csv',sep=';')

            # In[43]:

            stations_with_mse_improvements = sum(i < j for (i, j) in zip(
                df_improvements.mse_dwd_netatmo_interp.values,
                df_improvements.mse_dwd_interp.values))

            stations_without_mse_improvements = sum(i >= j for (i, j) in zip(
                df_improvements.mse_dwd_netatmo_interp.values,
                df_improvements.mse_dwd_interp.values))

            percent_of_mse_improvment = 100 * (stations_with_mse_improvements /
                                               df_improvements.mse_dwd_interp.shape[0])

            ##########################

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
                    df_improvements.mse_dwd_interp,
                    alpha=.8,
                    c='b',  # colors_arr,
                    marker='d',
                    label='DWD Interpolation RMSE mean = %0.1f' % np.mean(
                        df_improvements.mse_dwd_interp.values))
            ax.plot(df_improvements.index,
                    df_improvements.mse_dwd_netatmo_interp,
                    alpha=.8,
                    c='r',  # colors_arr,
                    marker='*',
                    label='DWD-Netatmo Interpolation RMSE mean = %0.1f' % np.mean(
                        df_improvements.mse_dwd_netatmo_interp.values))

            ax.set_title('Root mean squared error Interpolated from DWD or DWD-Netatmo %s \n '
                         'Rainfall of %s Extreme Events \n Stations with Improvemnts %d / %d, Percentage %0.0f'

                         % (_acc_, temp_freq, stations_with_mse_improvements,
                            df_improvements.mse_dwd_interp.shape[0], percent_of_mse_improvment))
            ax.grid(alpha=0.25)
            plt.setp(ax.get_xticklabels(), rotation=45)
            #ax.set_yticks(np.arange(0, 1.05, .10))
            ax.set_xlabel('DWD Stations')
            ax.legend(loc='upper right')
            ax.set_ylabel('RMSE')

            plt.savefig((path_to_use / (
                r'rmse_%s_events_dwd_%s_.png' % (temp_freq, _interp_acc_))),
                frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
            plt.close()

            #==================================================================
            #
            #==================================================================
            dwd_ids_no_improvement_pearson = pd.DataFrame(
                data=df_improvements.index[np.where(
                    df_improvements.pearson_corr_dwd_netatmo.values < df_improvements.pearson_corr_dwd_.values)],
                columns=['stns_no_impv_pearson_corr'])
            dwd_ids_no_improvement_spr = pd.DataFrame(data=df_improvements.index[np.where(
                df_improvements.spearman_corr_dwd_netatmo.values < df_improvements.spearman_corr_dwd_.values)],
                columns=['stns_no_impv_spearman_corr'])

            # dwd_ids_no_improvement_pearson.to_csv(
            #    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\stns_no_impv_pearson_corr_daily_events_all.csv',sep=';')

            # dwd_ids_no_improvement_spr.to_csv(
            #    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\stns_no_impv_spearman_corr_daily_events_all.csv',sep=';')

            # In[43]:

            stations_with_improvements = sum(i >= j for (i, j) in zip(
                df_improvements.pearson_corr_dwd_netatmo.values,
                df_improvements.pearson_corr_dwd_.values))

            stations_without_improvements = sum(i < j for (i, j) in zip(
                df_improvements.pearson_corr_dwd_netatmo.values,
                df_improvements.pearson_corr_dwd_.values))

            percent_of_improvment = 100 * (stations_with_improvements /
                                           df_improvements.pearson_corr_dwd_netatmo.shape[0])

            #########################################################
            plt.ioff()
            fig = plt.figure(figsize=(24, 12), dpi=150)

            ax = fig.add_subplot(111)

            ax.plot(df_improvements.index,
                    df_improvements.pearson_corr_dwd_,
                    alpha=.8,
                    c='b',  # colors_arr,
                    marker='d',
                    label='DWD Interpolation %0.2f' %
                    df_improvements.pearson_corr_dwd_.mean())

            ax.plot(df_improvements.index,
                    df_improvements.pearson_corr_dwd_netatmo,
                    alpha=.8,
                    c='r',  # colors_arr,
                    marker='*',
                    label='DWD-Netatmo Interpolation %0.2f' %
                    df_improvements.pearson_corr_dwd_netatmo.mean())

            ax.set_title('Pearson Correlation Interpolated from DWD or DWD-Netatmo %s \n '
                         'Rainfall of %s Extreme Events \n Stations with Improvemnts %d / %d, Percentage %0.0f'
                         % (_acc_, temp_freq, stations_with_improvements,
                            df_improvements.pearson_corr_dwd_netatmo.shape[0],
                             percent_of_improvment))
            ax.grid(alpha=0.25)
            plt.setp(ax.get_xticklabels(), rotation=45)
            ax.set_yticks(np.arange(0, 1.05, .10))
            ax.set_xlabel('DWD Stations')
            ax.legend(loc='upper right')
            ax.set_ylabel('Pearson Correlation')

            plt.savefig((path_to_use / (
                r'rmse_pears_corr_%s_events_dwd_%s_.png' % (temp_freq, _interp_acc_))),
                frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
            plt.close()

            # In[44]:

            stations_with_improvements = sum(i >= j for (i, j) in zip(
                df_improvements.spearman_corr_dwd_netatmo.values,
                df_improvements.spearman_corr_dwd_.values))

            percent_of_improvment = 100 * (stations_with_improvements /
                                           df_improvements.spearman_corr_dwd_netatmo.shape[0])
            #########################################################
            plt.ioff()
            fig = plt.figure(figsize=(24, 12), dpi=150)

            ax = fig.add_subplot(111)

            ax.plot(df_improvements.index,
                    df_improvements.spearman_corr_dwd_,
                    alpha=.8,
                    c='b',  # colors_arr,
                    marker='d',
                    label='DWD Interpolation %0.2f' %
                    df_improvements.spearman_corr_dwd_.mean())
            ax.plot(df_improvements.index,
                    df_improvements.spearman_corr_dwd_netatmo,
                    alpha=.8,
                    c='r',  # colors_arr,
                    marker='*',
                    label='DWD-Netatmo Interpolation %0.2f' %
                    df_improvements.spearman_corr_dwd_netatmo.mean())

            ax.set_title('Spearman Correlation Interpolated from DWD or DWD-Netatmo %s\n '
                         'Rainfall of %s Extreme Events \n Stations with Improvemnts %d / %d, Percentage %0.0f'

                         % (_acc_, temp_freq, stations_with_improvements,
                            df_improvements.spearman_corr_dwd_netatmo.shape[0],
                             percent_of_improvment
                            ))

            plt.setp(ax.get_xticklabels(), rotation=45)
            ax.grid(alpha=0.25)
            ax.set_yticks(np.arange(0, 1.05, .10))
            ax.set_xlabel('DWD Stations')
            ax.legend(loc='upper right')
            ax.set_ylabel('Spearman Correlation')

            plt.savefig((path_to_use / (
                r'rmse_spr_corr_%s_events_dwd_%s_.png' % (temp_freq, _interp_acc_))),
                frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
            plt.close()
