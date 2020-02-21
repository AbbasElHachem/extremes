
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
# import matplotlib.dates as mdates
import fnmatch
from pathlib import Path

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

plt.ioff()
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'axes.labelsize': 16})


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


if plot_filtered:
    # , '60min' '180min', '360min', '720min', '1440min']:
    for temp_freq in ['60min', '180min', '360min', '720min', '1440min']:
        print(temp_freq)

        path_to_Qt_ok_un_first_flt__temp_flt_1st_ = main_dir / (
            r'Final_results8/Ppt_ok_ok_un_new6_first_flt__temp_flt__1st_%s' % temp_freq)
        Qt_ok_un_first_flt__temp_flt_1st_ = list_all_full_path(
            '.csv', path_to_Qt_ok_un_first_flt__temp_flt_1st_)

        path_to_Qt_ok_un_first_flt__temp_flt_1st_2 = main_dir / (
            r'Final_results5/Ppt_ok_ok_un_new4_first_flt__temp_flt__1st_%s' % temp_freq)
        Qt_ok_un_first_flt__temp_flt_1st_2 = list_all_full_path(
            '.csv', path_to_Qt_ok_un_first_flt__temp_flt_1st_2)

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

        df_improvements = pd.DataFrame(
            index=df_dwd_edf.columns,
            columns=['pearson_corr_dwd_',
                     'spearman_corr_dwd_',
                     'pearson_corr_dwd_netatmo',
                     'spearman_corr_dwd_netatmo',
                     'pearson_corr_dwd_netatmo_unc05perc',
                     'spearman_corr_dwd_netatmo_unc05perc',
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

        path_to_use2 = path_to_Qt_ok_un_first_flt__temp_flt_1st_2
        data_to_use2 = Qt_ok_un_first_flt__temp_flt_1st_2

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

                if ('dwd_netamo') in df_file and 'unc05perc' in df_file:
                    print(df_file)
                    path_interpolated_using_netatmo_dwd_unc05perc = df_file

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

        try:
            for df_file in data_to_use2:

                if ('dwd_netamo') in df_file and 'unc' not in df_file:
                    print(df_file)
                    path_interpolated_using_netatmo_dwd2 = df_file

        except Exception as msg:
            print(msg)
            continue
        df_netatmo_dwd2 = pd.read_csv(path_interpolated_using_netatmo_dwd2,
                                      sep=';', index_col=0,
                                      parse_dates=True,
                                      infer_datetime_format=True)

        #======================================================================
        #
        #======================================================================
        df_dwd = pd.read_csv(path_interpolated_using_dwd,
                             sep=';', index_col=0, parse_dates=True,
                             infer_datetime_format=True)

        events_to_compare = df_netatmo_dwd2.index.intersection(df_dwd.index)
        # print(df_dwd.isna().sum().max())

        df_netatmo_dwd = pd.read_csv(path_interpolated_using_netatmo_dwd,
                                     sep=';', index_col=0,
                                     parse_dates=True,
                                     infer_datetime_format=True)

        print(df_netatmo_dwd.isna().sum().max())

        df_netatmo_dwd_unc05perc = pd.read_csv(
            path_interpolated_using_netatmo_dwd_unc05perc,
            sep=';', index_col=0,
            parse_dates=True, infer_datetime_format=True)

        df_netatmo_dwd_unc2perc = pd.read_csv(
            path_interpolated_using_netatmo_dwd_unc2perc,
            sep=';', index_col=0,
            parse_dates=True, infer_datetime_format=True)

        df_netatmo_dwd_unc5perc = pd.read_csv(
            path_interpolated_using_netatmo_dwd_unc5perc,
            sep=';', index_col=0,
            parse_dates=True, infer_datetime_format=True)

        df_netatmo_dwd_unc10perc = pd.read_csv(
            path_interpolated_using_netatmo_dwd_unc10perc,
            sep=';', index_col=0,
            parse_dates=True, infer_datetime_format=True)

        df_netatmo_dwd_unc20perc = pd.read_csv(
            path_interpolated_using_netatmo_dwd_unc20perc,
            sep=';', index_col=0,
            parse_dates=True, infer_datetime_format=True)
        #########################################################

        # stns to remove from orig edf, because no data for 2015-2019
        df_dwd_edf.drop(columns=['00384', '13672', '05155'], inplace=True)

        try:

            for stn_ in df_dwd_edf.columns:
                #                 if stn_ == '00384':
                #                     raise Exception
                #                     pass
                #                     print(stn_)
                df_compare = pd.DataFrame(index=df_dwd.index)

                # :
                for event_date in df_dwd.index.intersection(df_dwd_edf.index):
                    # print(event_date)
                    edf_stn_orig = df_dwd_edf.loc[event_date, stn_]
                    edf_stn_interp_dwd = df_dwd.loc[event_date, stn_]
                    edf_stn_interp_netatmo_dwd = df_netatmo_dwd.loc[event_date, stn_]
                    edf_stn_interp_netatmo_dwd_unc05perc = df_netatmo_dwd_unc05perc.loc[
                        event_date, stn_]
                    edf_stn_interp_netatmo_dwd_unc2perc = df_netatmo_dwd_unc2perc.loc[
                        event_date, stn_]
                    edf_stn_interp_netatmo_dwd_unc5perc = df_netatmo_dwd_unc5perc.loc[
                        event_date, stn_]
                    edf_stn_interp_netatmo_dwd_unc10perc = df_netatmo_dwd_unc10perc.loc[
                        event_date, stn_]
                    edf_stn_interp_netatmo_dwd_unc20perc = df_netatmo_dwd_unc20perc.loc[
                        event_date, stn_]
                    if ((edf_stn_orig >= 0) and
                        (edf_stn_interp_dwd >= 0) and
                            (edf_stn_interp_netatmo_dwd >= 0 and
                             (edf_stn_interp_netatmo_dwd_unc2perc >= 0) and
                             (edf_stn_interp_netatmo_dwd_unc05perc >= 0))):

                        df_compare.loc[
                            event_date,
                            'original_quantile'] = edf_stn_orig
                        df_compare.loc[
                            event_date,
                            'interpolated_quantile_dwd'] = edf_stn_interp_dwd

                        df_compare.loc[
                            event_date,
                            'interpolated_quantile_netatmo_dwd'] = edf_stn_interp_netatmo_dwd
                        df_compare.loc[
                            event_date,
                            'interpolated_quantile_netatmo_dwd_unc05perc'] = edf_stn_interp_netatmo_dwd_unc05perc

                        df_compare.loc[
                            event_date,
                            'interpolated_quantile_netatmo_dwd_unc2perc'] = edf_stn_interp_netatmo_dwd_unc2perc

                        df_compare.loc[
                            event_date,
                            'interpolated_quantile_netatmo_dwd_unc5perc'] = edf_stn_interp_netatmo_dwd_unc5perc

                        df_compare.loc[
                            event_date,
                            'interpolated_quantile_netatmo_dwd_unc10perc'] = edf_stn_interp_netatmo_dwd_unc10perc

                        df_compare.loc[
                            event_date,
                            'interpolated_quantile_netatmo_dwd_unc20perc'] = edf_stn_interp_netatmo_dwd_unc20perc

                    else:
                        df_compare.loc[event_date,
                                       'original_quantile'] = np.nan
                        df_compare.loc[event_date,
                                       'interpolated_quantile_dwd'] = np.nan

                        df_compare.loc[event_date,
                                       'interpolated_quantile_netatmo_dwd'] = np.nan
                        df_compare.loc[event_date,
                                       'interpolated_quantile_netatmo_dwd_unc05perc'] = np.nan

                        df_compare.loc[event_date,
                                       'interpolated_quantile_netatmo_dwd_unc2perc'] = np.nan

                        df_compare.loc[event_date,
                                       'interpolated_quantile_netatmo_dwd_unc5perc'] = np.nan

                        df_compare.loc[event_date,
                                       'interpolated_quantile_netatmo_dwd_unc10perc'] = np.nan
                        df_compare.loc[event_date,
                                       'interpolated_quantile_netatmo_dwd_unc20perc'] = np.nan

                df_compare = df_compare[df_compare > 0]
                df_compare.dropna(how='any', inplace=True)
                if df_compare.values.shape[0] > 1:  # at least 20evnts
                    values_x = df_compare['original_quantile'].values
                    values_dwd = df_compare['interpolated_quantile_dwd'].values

                    values_netatmo_dwd = df_compare['interpolated_quantile_netatmo_dwd'].values
                    values_netatmo_dwd_unc05perc = df_compare['interpolated_quantile_netatmo_dwd_unc05perc'].values

                    values_netatmo_dwd_unc2perc = df_compare['interpolated_quantile_netatmo_dwd_unc2perc'].values

                    values_netatmo_dwd_unc5perc = df_compare['interpolated_quantile_netatmo_dwd_unc5perc'].values
                    values_netatmo_dwd_unc10perc = df_compare['interpolated_quantile_netatmo_dwd_unc10perc'].values
                    values_netatmo_dwd_unc20perc = df_compare['interpolated_quantile_netatmo_dwd_unc20perc'].values
                    # plot the stations in shapefile, look at the results of
                    # agreements

                    # calculate correlations (pearson and spearman) dwd-dwd
                    corr_dwd = pears(values_x, values_dwd)[0]
                    rho_dwd = spr(values_x, values_dwd)[0]
                    print(corr_dwd)
#                     if corr_dwd < 0:
#                         raise Exception
#                         pass
                    # netatmo dwd
                    corr_netatmo_dwd = pears(values_x, values_netatmo_dwd)[0]
                    rho_netatmo_dwd = spr(values_x, values_netatmo_dwd)[0]
                    # netatmo dwd unc
                    corr_netatmo_dwd_unc05perc = pears(
                        values_x, values_netatmo_dwd_unc05perc)[0]
                    rho_netatmo_dwd_unc05perc = spr(
                        values_x, values_netatmo_dwd_unc05perc)[0]

                    corr_netatmo_dwd_unc2perc = pears(
                        values_x, values_netatmo_dwd_unc2perc)[0]
                    rho_netatmo_dwd_unc2perc = spr(
                        values_x, values_netatmo_dwd_unc2perc)[0]

                    corr_netatmo_dwd_unc5perc = pears(
                        values_x, values_netatmo_dwd_unc5perc)[0]
                    rho_netatmo_dwd_unc5perc = spr(
                        values_x, values_netatmo_dwd_unc5perc)[0]

                    corr_netatmo_dwd_unc10perc = pears(
                        values_x, values_netatmo_dwd_unc10perc)[0]
                    rho_netatmo_dwd_unc10perc = spr(
                        values_x, values_netatmo_dwd_unc10perc)[0]

                    corr_netatmo_dwd_unc20perc = pears(
                        values_x, values_netatmo_dwd_unc20perc)[0]
                    rho_netatmo_dwd_unc20perc = spr(
                        values_x, values_netatmo_dwd_unc20perc)[0]

                    df_improvements.loc[stn_, 'pearson_corr_dwd_'] = corr_dwd
                    df_improvements.loc[stn_, 'spearman_corr_dwd_'] = rho_dwd
                    df_improvements.loc[stn_,
                                        'pearson_corr_dwd_netatmo'] = corr_netatmo_dwd
                    df_improvements.loc[stn_,
                                        'spearman_corr_dwd_netatmo'] = rho_netatmo_dwd
                    df_improvements.loc[stn_,
                                        'pearson_corr_dwd_netatmo_unc05perc'] = corr_netatmo_dwd_unc05perc
                    df_improvements.loc[stn_,
                                        'spearman_corr_dwd_netatmo_unc05perc'] = rho_netatmo_dwd_unc05perc

                    df_improvements.loc[stn_,
                                        'pearson_corr_dwd_netatmo_unc2perc'] = corr_netatmo_dwd_unc2perc
                    df_improvements.loc[stn_,
                                        'spearman_corr_dwd_netatmo_unc2perc'] = rho_netatmo_dwd_unc2perc

                    df_improvements.loc[stn_,
                                        'pearson_corr_dwd_netatmo_unc5perc'] = corr_netatmo_dwd_unc5perc
                    df_improvements.loc[stn_,
                                        'spearman_corr_dwd_netatmo_unc5perc'] = rho_netatmo_dwd_unc5perc

                    df_improvements.loc[stn_,
                                        'pearson_corr_dwd_netatmo_unc10perc'] = corr_netatmo_dwd_unc10perc
                    df_improvements.loc[stn_,
                                        'spearman_corr_dwd_netatmo_unc10perc'] = rho_netatmo_dwd_unc10perc

                    df_improvements.loc[stn_,
                                        'pearson_corr_dwd_netatmo_unc20perc'] = corr_netatmo_dwd_unc20perc
                    df_improvements.loc[stn_,
                                        'spearman_corr_dwd_netatmo_unc20perc'] = rho_netatmo_dwd_unc20perc

            df_improvements.dropna(how='all', inplace=True)

            # df_compare = df_compare[df_compare > 0]
            #df_compare.dropna(how='all', inplace=True)

        except Exception as msg:
            print(msg)
        # dwd stn '07135', with weird results

        # hourly_events = ['2016-06-25 00:00:00',
#                  '2018-06-11 17:00:00',
#                  '2018-09-23 17:00:00',
#                  '2018-09-23 18:00:00',
#                  '2018-09-23 19:00:00']

# daily_events = ['2018-12-23 00:00:00',
#                 '2019-05-22 00:00:00']
        # OK
        stations_with_improvements = sum(i >= j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo.values,
            df_improvements.pearson_corr_dwd_.values))

        stations_without_improvements = sum(i < j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo.values,
            df_improvements.pearson_corr_dwd_.values))

        percent_of_improvment = 100 * (stations_with_improvements /
                                       df_improvements.pearson_corr_dwd_netatmo.shape[0])

        id_stns_no_impv = [id_stn for (id_stn, i, j) in zip(
            df_improvements.index,
            df_improvements.pearson_corr_dwd_netatmo.values,
            df_improvements.pearson_corr_dwd_.values) if i < j]

    #     # OK with Unc
        stations_with_improvements_unc05perc = sum(i >= j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo_unc05perc.values,
            df_improvements.pearson_corr_dwd_.values))
        stations_without_improvements_unc05perc = sum(i < j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo_unc05perc.values,
            df_improvements.pearson_corr_dwd_.values))
        percent_of_improvment_unc05perc = 100 * (
            stations_with_improvements_unc05perc /
            df_improvements.pearson_corr_dwd_netatmo_unc2perc.shape[0])

        stations_with_improvements_unc2perc = sum(i >= j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo_unc2perc.values,
            df_improvements.pearson_corr_dwd_.values))

        stations_without_improvements_unc2perc = sum(i < j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo_unc2perc.values,
            df_improvements.pearson_corr_dwd_.values))

        percent_of_improvment_unc2perc = 100 * (
            stations_with_improvements_unc2perc /
            df_improvements.pearson_corr_dwd_netatmo_unc2perc.shape[0])

        stations_with_improvements_unc5perc = sum(i >= j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo_unc5perc.values,
            df_improvements.pearson_corr_dwd_.values))

        stations_without_improvements_unc5perc = sum(i < j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo_unc5perc.values,
            df_improvements.pearson_corr_dwd_.values))

        percent_of_improvment_unc5perc = 100 * (
            stations_with_improvements_unc5perc /
            df_improvements.pearson_corr_dwd_netatmo_unc2perc.shape[0])

        stations_with_improvements_unc10perc = sum(i >= j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo_unc10perc.values,
            df_improvements.pearson_corr_dwd_.values))

        stations_without_improvements_unc10perc = sum(i < j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo_unc10perc.values,
            df_improvements.pearson_corr_dwd_.values))

        percent_of_improvment_unc10perc = 100 * (
            stations_with_improvements_unc10perc /
            df_improvements.pearson_corr_dwd_netatmo_unc2perc.shape[0])

        stations_with_improvements_unc20perc = sum(i >= j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo_unc20perc.values,
            df_improvements.pearson_corr_dwd_.values))

        stations_without_improvements_unc20perc = sum(i < j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo_unc20perc.values,
            df_improvements.pearson_corr_dwd_.values))

        percent_of_improvment_unc20perc = 100 * (
            stations_with_improvements_unc20perc /
            df_improvements.pearson_corr_dwd_netatmo_unc2perc.shape[0])

        ####
        mean_pearson_correlation_dwd_only = df_improvements.pearson_corr_dwd_.mean()
        mean_pearson_correlation_dwd_netatmo = df_improvements.pearson_corr_dwd_netatmo.mean()

        mean_pearson_correlation_dwd_netatmo_unc05perc = df_improvements.pearson_corr_dwd_netatmo_unc05perc.mean()
        mean_pearson_correlation_dwd_netatmo_unc2perc = df_improvements.pearson_corr_dwd_netatmo_unc2perc.mean()

        mean_pearson_correlation_dwd_netatmo_unc5perc = df_improvements.pearson_corr_dwd_netatmo_unc5perc.mean()
        mean_pearson_correlation_dwd_netatmo_unc10perc = df_improvements.pearson_corr_dwd_netatmo_unc10perc.mean()
        mean_pearson_correlation_dwd_netatmo_unc20perc = df_improvements.pearson_corr_dwd_netatmo_unc20perc.mean()
        #########################################################
        mean_spr_correlation_dwd_only = df_improvements.spearman_corr_dwd_.mean()
        mean_spr_correlation_dwd_netatmo = df_improvements.spearman_corr_dwd_netatmo.mean()

        mean_spr_correlation_dwd_netatmo_unc05perc = df_improvements.spearman_corr_dwd_netatmo_unc05perc.mean()
        mean_spr_correlation_dwd_netatmo_unc2perc = df_improvements.spearman_corr_dwd_netatmo_unc2perc.mean()
        mean_spr_correlation_dwd_netatmo_unc5perc = df_improvements.spearman_corr_dwd_netatmo_unc5perc.mean()
        mean_spr_correlation_dwd_netatmo_unc10perc = df_improvements.spearman_corr_dwd_netatmo_unc10perc.mean()
        mean_spr_correlation_dwd_netatmo_unc20perc = df_improvements.spearman_corr_dwd_netatmo_unc20perc.mean()
        #########################################################
        df_improvements.iloc[np.where(
            df_improvements.pearson_corr_dwd_.isna())]
        plt.ioff()
        fig = plt.figure(figsize=(32, 12), dpi=150)

        ax = fig.add_subplot(111)

        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation %0.2f'
                % mean_pearson_correlation_dwd_only)
        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_netatmo,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %0.2f'
                % mean_pearson_correlation_dwd_netatmo)
        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_netatmo_unc05perc,
                alpha=.8,
                c='k',  # colors_arr,
                marker='2',
                label='DWD-Netatmo Interpolation 05percUnc %0.2f'
                % mean_pearson_correlation_dwd_netatmo_unc05perc)
        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_netatmo_unc2perc,
                alpha=.8,
                c='g',  # colors_arr,
                marker='+',
                label='DWD-Netatmo Interpolation 2percUnc %0.2f'
                % mean_pearson_correlation_dwd_netatmo_unc2perc)

        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_netatmo_unc5perc,
                alpha=.8,
                c='m',  # colors_arr,
                marker='1',
                label='DWD-Netatmo Interpolation 5percUnc %0.2f'
                % mean_pearson_correlation_dwd_netatmo_unc5perc)

        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_netatmo_unc10perc,
                alpha=.8,
                c='c',  # colors_arr,
                marker='3',
                label='DWD-Netatmo Interpolation 10percUnc %0.2f'
                % mean_pearson_correlation_dwd_netatmo_unc10perc)

        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_netatmo_unc20perc,
                alpha=.8,
                c='orange',  # colors_arr,
                marker=',',
                label='DWD-Netatmo Interpolation 20percUnc %0.2f'
                % mean_pearson_correlation_dwd_netatmo_unc20perc)

        ax.set_title('Pearson Correlation Interpolated Rainfall from DWD or DWD-Netatmo\n '
                     'Precipitation of %s Intense Events \n Stations with Improvemnts %d / %d,'
                     ' Percentage %0.0f\n'
                     'Stations with Improvemnts with OK 05percUnc %d / %d, Percentage %0.0f'
                     '\nStations with Improvemnts with OK 2percUnc %d / %d, Percentage %0.0f'
                     '\nStations with Improvemnts with OK 5percUnc %d / %d, Percentage %0.0f'
                     '\nStations with Improvemnts with OK 10percUnc %d / %d, Percentage %0.0f'
                     '\nStations with Improvemnts with OK 20percUnc %d / %d, Percentage %0.0f'
                     % (temp_freq,  # _interp_acc_,
                        stations_with_improvements,
                         df_improvements.pearson_corr_dwd_netatmo.shape[0],
                        percent_of_improvment,
                        stations_with_improvements_unc05perc,
                        df_improvements.pearson_corr_dwd_netatmo_unc05perc.shape[0],
                        percent_of_improvment_unc05perc,
                        stations_with_improvements_unc2perc,
                        df_improvements.pearson_corr_dwd_netatmo_unc2perc.shape[0],
                        percent_of_improvment_unc2perc,
                        stations_with_improvements_unc5perc,
                        df_improvements.pearson_corr_dwd_netatmo_unc2perc.shape[0],
                        percent_of_improvment_unc5perc,
                        stations_with_improvements_unc10perc,
                        df_improvements.pearson_corr_dwd_netatmo_unc2perc.shape[0],
                        percent_of_improvment_unc10perc,
                        stations_with_improvements_unc20perc,
                        df_improvements.pearson_corr_dwd_netatmo_unc2perc.shape[0],
                        percent_of_improvment_unc20perc))

        plt.setp(ax.get_xticklabels(), rotation=90)
        ax.grid(alpha=0.25)
        ax.set_yticks(np.arange(0, 1.05, .10))
        ax.set_xlabel('DWD Stations')
        ax.legend(loc='lower left')
        ax.set_ylabel('Pearson Correlation')

        plt.savefig(os.path.join(path_to_use,
                                 r'ppt_temporal_pears_corr_%s_events_dwd_cr.png' % temp_freq,
                                 ),
                    frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

        ####################################################
        stations_with_improvements = sum(i >= j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo.values,
            df_improvements.spearman_corr_dwd_.values))

        stations_without_improvements = sum(i < j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo.values,
            df_improvements.spearman_corr_dwd_.values))

        percent_of_improvment = 100 * (stations_with_improvements /
                                       df_improvements.spearman_corr_dwd_netatmo.shape[0])

        # ok with Un
        stations_with_improvements_unc05perc = sum(i >= j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo_unc05perc.values,
            df_improvements.spearman_corr_dwd_.values))

        stations_without_improvements_unc05perc = sum(i < j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo_unc05perc.values,
            df_improvements.spearman_corr_dwd_.values))

        percent_of_improvment_unc05perc = 100 * (
            stations_with_improvements_unc05perc /
            df_improvements.spearman_corr_dwd_netatmo_unc2perc.shape[0])

        stations_with_improvements_unc2perc = sum(i >= j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo_unc2perc.values,
            df_improvements.spearman_corr_dwd_.values))

        stations_without_improvements_unc2perc = sum(i < j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo_unc2perc.values,
            df_improvements.spearman_corr_dwd_.values))

        percent_of_improvment_unc2perc = 100 * (
            stations_with_improvements_unc2perc /
            df_improvements.spearman_corr_dwd_netatmo_unc2perc.shape[0])

        stations_with_improvements_unc5perc = sum(i >= j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo_unc5perc.values,
            df_improvements.spearman_corr_dwd_.values))

        stations_without_improvements_unc5perc = sum(i < j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo_unc5perc.values,
            df_improvements.spearman_corr_dwd_.values))

        percent_of_improvment_unc5perc = 100 * (
            stations_with_improvements_unc5perc /
            df_improvements.spearman_corr_dwd_netatmo_unc2perc.shape[0])

        stations_with_improvements_unc10perc = sum(i >= j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo_unc10perc.values,
            df_improvements.spearman_corr_dwd_.values))

        stations_without_improvements_unc10perc = sum(i < j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo_unc10perc.values,
            df_improvements.spearman_corr_dwd_.values))

        percent_of_improvment_unc10perc = 100 * (
            stations_with_improvements_unc10perc /
            df_improvements.spearman_corr_dwd_netatmo_unc2perc.shape[0])

        stations_with_improvements_unc20perc = sum(i >= j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo_unc20perc.values,
            df_improvements.spearman_corr_dwd_.values))

        stations_without_improvements_unc20perc = sum(i < j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo_unc20perc.values,
            df_improvements.spearman_corr_dwd_.values))

        percent_of_improvment_unc20perc = 100 * (
            stations_with_improvements_unc20perc /
            df_improvements.spearman_corr_dwd_netatmo_unc2perc.shape[0])
        #########################################################
        plt.ioff()
        fig = plt.figure(figsize=(32, 12), dpi=150)

        ax = fig.add_subplot(111)

        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation %0.2f'
                % mean_spr_correlation_dwd_only)
        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_netatmo,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %0.2f'
                % mean_spr_correlation_dwd_netatmo)
        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_netatmo_unc05perc,
                alpha=.8,
                c='k',  # colors_arr,
                marker='2',
                label='DWD-Netatmo Interpolation 05percUnc %0.2f'
                % mean_spr_correlation_dwd_netatmo_unc2perc)
        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_netatmo_unc2perc,
                alpha=.8,
                c='g',  # colors_arr,
                marker='+',
                label='DWD-Netatmo Interpolation 2percUnc %0.2f'
                % mean_spr_correlation_dwd_netatmo_unc2perc)

        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_netatmo_unc5perc,
                alpha=.8,
                c='m',  # colors_arr,
                marker='1',
                label='DWD-Netatmo Interpolation 5percUnc %0.2f'
                % mean_spr_correlation_dwd_netatmo_unc5perc)

        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_netatmo_unc10perc,
                alpha=.8,
                c='c',  # colors_arr,
                marker='3',
                label='DWD-Netatmo Interpolation 10percUnc %0.2f'
                % mean_spr_correlation_dwd_netatmo_unc10perc)

        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_netatmo_unc20perc,
                alpha=.8,
                c='orange',  # colors_arr,
                marker=',',
                label='DWD-Netatmo Interpolation 20percUnc %0.2f'
                % mean_spr_correlation_dwd_netatmo_unc20perc)

        ax.set_title('Spearman Correlation Interpolated Rainfall from DWD or DWD-Netatmo \n '
                     'Rainfall of %s Intense Events \n Stations with Improvemnts %d / %d, Percentage %0.0f'
                     '\nStations with Improvemnts with OK 05percUnc %d / %d, Percentage %0.0f'
                     '\nStations with Improvemnts with OK 2percUnc %d / %d, Percentage %0.0f'
                     '\nStations with Improvemnts with OK 5percUnc %d / %d, Percentage %0.0f'
                     '\nStations with Improvemnts with OK 10percUnc %d / %d, Percentage %0.0f'
                     '\nStations with Improvemnts with OK 20percUnc %d / %d, Percentage %0.0f'
                     % (temp_freq,  # _interp_acc_,
                        stations_with_improvements,
                         df_improvements.spearman_corr_dwd_netatmo.shape[0],
                        percent_of_improvment,
                        stations_with_improvements_unc05perc,
                        df_improvements.spearman_corr_dwd_netatmo_unc05perc.shape[0],
                        percent_of_improvment_unc05perc,
                        stations_with_improvements_unc2perc,
                        df_improvements.spearman_corr_dwd_netatmo_unc2perc.shape[0],
                        percent_of_improvment_unc2perc,
                        stations_with_improvements_unc5perc,
                        df_improvements.spearman_corr_dwd_netatmo_unc2perc.shape[0],
                        percent_of_improvment_unc5perc,
                        stations_with_improvements_unc10perc,
                        df_improvements.spearman_corr_dwd_netatmo_unc2perc.shape[0],
                        percent_of_improvment_unc10perc,
                        stations_with_improvements_unc20perc,
                        df_improvements.spearman_corr_dwd_netatmo_unc2perc.shape[0],
                        percent_of_improvment_unc20perc))

        plt.setp(ax.get_xticklabels(), rotation=90)
        ax.grid(alpha=0.25)
        ax.set_yticks(np.arange(0, 1.05, .10))
        ax.set_xlabel('DWD Stations')

        ax.legend(loc='lower left')
        ax.set_ylabel('Spearman Correlation')

        plt.savefig(os.path.join(path_to_use,
                                 r'ppt_temporal_spearm_corr_%s_events_dwd_cr.png' % (temp_freq)),
                    frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

#==============================================================================
#
#==============================================================================
if plot_not_filtered:
        # , '60min' '180min', '360min', '720min', '1440min']:
    for temp_freq in ['60min', '180min', '360min', '720min', '1440min']:
        print(temp_freq)

        path_to_Quantiles_netatmo_no_flt___ = main_dir / (
            r'Ppt_ok_ok_un_new_netatmo_no_flt___%s' % temp_freq)

        Quantiles_netatmo_no_flt___ = list_all_full_path(
            '.csv', path_to_Quantiles_netatmo_no_flt___)

        path_dwd_ppt_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\ppt_all_dwd_%s_.csv" % temp_freq
        df_dwd_edf = pd.read_csv(path_dwd_ppt_data, sep=';', index_col=0,
                                 parse_dates=True, infer_datetime_format=True)
        # df_dwd_edf.dropna(inplace=True)

        print(df_dwd_edf.isna().sum().max())

        df_improvements = pd.DataFrame(
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

        print(df_dwd.isna().sum().max())

        df_netatmo_dwd = pd.read_csv(path_interpolated_using_netatmo_dwd,
                                     sep=';', index_col=0,
                                     parse_dates=True,
                                     infer_datetime_format=True)

        print(df_netatmo_dwd.isna().sum().max())

        try:

            for stn_ in df_dwd_edf.columns:
                df_compare = pd.DataFrame(index=df_dwd.index)
                for event_date in df_dwd.index.intersection(df_dwd_edf.index):
                    #                     print(event_date)
                    edf_stn_orig = df_dwd_edf.loc[event_date, stn_]
                    edf_stn_interp_dwd = df_dwd.loc[event_date, stn_]
                    edf_stn_interp_netatmo_dwd = df_netatmo_dwd.loc[event_date, stn_]

                    if ((edf_stn_orig >= 0) and
                        (edf_stn_interp_dwd >= 0) and
                            (edf_stn_interp_netatmo_dwd >= 0)):

                        df_compare.loc[
                            event_date,
                            'original_quantile'] = edf_stn_orig
                        df_compare.loc[
                            event_date,
                            'interpolated_quantile_dwd'] = edf_stn_interp_dwd

                        df_compare.loc[
                            event_date,
                            'interpolated_quantile_netatmo_dwd'] = edf_stn_interp_netatmo_dwd

                    else:
                        df_compare.loc[event_date,
                                       'original_quantile'] = np.nan
                        df_compare.loc[event_date,
                                       'interpolated_quantile_dwd'] = np.nan

                        df_compare.loc[event_date,
                                       'interpolated_quantile_netatmo_dwd'] = np.nan

                df_compare = df_compare[df_compare > 0]
                df_compare.dropna(how='any', inplace=True)

                if df_compare.values.shape[0] > 20:
                    values_x = df_compare['original_quantile'].values
                    values_dwd = df_compare['interpolated_quantile_dwd'].values

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

            # df_compare = df_compare[df_compare > 0]
            #df_compare.dropna(how='all', inplace=True)

        except Exception as msg:
            print(msg)

        # OK
        stations_with_improvements = sum(i >= j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo.values,
            df_improvements.pearson_corr_dwd_.values))

        stations_without_improvements = sum(i < j for (i, j) in zip(
            df_improvements.pearson_corr_dwd_netatmo.values,
            df_improvements.pearson_corr_dwd_.values))

        percent_of_improvment = 100 * (stations_with_improvements /
                                       df_improvements.pearson_corr_dwd_netatmo.shape[0])

        ####
        mean_pearson_correlation_dwd_only = df_improvements.pearson_corr_dwd_.mean()
        mean_pearson_correlation_dwd_netatmo = df_improvements.pearson_corr_dwd_netatmo.mean()

        #########################################################
        mean_spr_correlation_dwd_only = df_improvements.spearman_corr_dwd_.mean()
        mean_spr_correlation_dwd_netatmo = df_improvements.spearman_corr_dwd_netatmo.mean()

        #########################################################
        plt.ioff()
        fig = plt.figure(figsize=(24, 12), dpi=150)

        ax = fig.add_subplot(111)

        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation %0.2f'
                % mean_pearson_correlation_dwd_only)
        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_netatmo,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %0.2f'
                % mean_pearson_correlation_dwd_netatmo)

        ax.set_title('Pearson Correlation Interpolated Quantiles from DWD or DWD-Netatmo\n '
                     'Precipitation of %s Intense Events \n Events with Improvemnts %d / %d, Percentage %0.0f\n'

                     % (temp_freq,  # _interp_acc_,
                        stations_with_improvements,
                         df_improvements.pearson_corr_dwd_netatmo.shape[0],
                        percent_of_improvment))
        ax.grid(alpha=0.25)
        plt.setp(ax.get_xticklabels(), rotation=60)
        ax.grid(alpha=0.25)
        ax.set_yticks(np.arange(0, 1.05, .10))
        ax.set_xlabel('DWD Stations')
        ax.legend(loc='lower right')
        ax.set_ylabel('Pearson Correlation')

        plt.savefig(os.path.join(path_to_use,
                                 r'ppt_temporal_pears_corr_%s_events_dwd.png' % temp_freq,
                                 ),
                    frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

        ####################################################
        stations_with_improvements = sum(i >= j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo.values,
            df_improvements.spearman_corr_dwd_.values))

        stations_without_improvements = sum(i < j for (i, j) in zip(
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
                label='DWD Interpolation %0.2f'
                % mean_spr_correlation_dwd_only)
        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_netatmo,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %0.2f'
                % mean_spr_correlation_dwd_netatmo)

        ax.set_title('Spearman Correlation Interpolated Quantiles from DWD or DWD-Netatmo \n '
                     'Rainfall of %s Intense Events  \n Events with Improvemnts %d / %d, Percentage %0.0f'

                     % (temp_freq,  # _interp_acc_,
                        stations_with_improvements,
                         df_improvements.spearman_corr_dwd_netatmo.shape[0],
                        percent_of_improvment,
                        ))
        ax.grid(alpha=0.25)
        plt.setp(ax.get_xticklabels(), rotation=60)
        ax.grid(alpha=0.25)
        ax.set_yticks(np.arange(0, 1.05, .10))
        ax.set_xlabel('DWD Stations')
        ax.legend(loc='lower right')
        ax.set_ylabel('Spearman Correlation')

        plt.savefig(os.path.join(path_to_use,
                                 r'ppt_temporal_spearm_corr_%s_events_dwd.png' % (temp_freq)),
                    frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()
