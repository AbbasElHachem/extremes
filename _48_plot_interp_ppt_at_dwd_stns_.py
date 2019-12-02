
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

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

plt.ioff()
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'axes.labelsize': 12})


# In[12]:


temp_freq = '60min'  # daily

main_dir = Path(
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo')


# In[13]:


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


for temp_freq in ['60min', '360min', '720min', '1440min']:
    print(temp_freq)

    path_to_Qt_ok_un_first_flt__temp_flt_1st_ = main_dir / (
        r'Ppt_ok_ok_un_new_first_flt__temp_flt__1st_%s' % temp_freq)
    Qt_ok_un_first_flt__temp_flt_1st_ = list_all_full_path(
        '.csv', path_to_Qt_ok_un_first_flt__temp_flt_1st_)

    path_to_Qt_ok_un_first_flt__temp_flt_comb_ = main_dir / (
        r'Ppt_ok_ok_un_new_first_flt__temp_flt__comb_%s' % temp_freq)
    Qt_ok_un_first_flt__temp_flt_comb_ = list_all_full_path(
        '.csv', path_to_Qt_ok_un_first_flt__temp_flt_comb_)

    path_to_Qt_ok_un_first_flt_1st_ = main_dir / (
        r'Ppt_ok_ok_un_new_first_flt__1st_%s' % temp_freq)
    Qt_ok_un_first_flt_1st_ = list_all_full_path(
        '.csv', path_to_Qt_ok_un_first_flt_1st_)

    path_to_Qt_ok_un_first_flt_comb_ = main_dir / (
        r'Ppt_ok_ok_un_new_first_flt__comb_%s' % temp_freq)
    Qt_ok_un_first_flt_comb_ = list_all_full_path(
        '.csv', path_to_Qt_ok_un_first_flt_comb_)

    path_to_Quantiles_netatmo_no_flt___ = main_dir / (
        r'Ppt_ok_ok_un_new_netatmo_no_flt___%s' % temp_freq)

    Quantiles_netatmo_no_flt___ = list_all_full_path(
        '.csv', path_to_Quantiles_netatmo_no_flt___)

    path_dwd_ppt_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\ppt_all_dwd_%s_.csv" % temp_freq
    df_dwd_edf = pd.read_csv(path_dwd_ppt_data, sep=';', index_col=0,
                             parse_dates=True, infer_datetime_format=True)
    # df_dwd_edf.dropna(inplace=True)

    print(df_dwd_edf.isna().sum().max())

    df_improvements = pd.DataFrame(index=df_dwd_edf.columns,
                                   columns=['pearson_corr_dwd_',
                                            'spearman_corr_dwd_',
                                            'pearson_corr_dwd_netatmo',
                                            'spearman_corr_dwd_netatmo',
                                            'pearson_corr_dwd_netatmo_unc',
                                            'spearman_corr_dwd_netatmo_unc'])

    min_orig_qnt_thr = 0.

    #########################################################

    path_to_use = path_to_Qt_ok_un_first_flt_comb_
    data_to_use = Qt_ok_un_first_flt_comb_

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

            if ('dwd_netamo') in df_file:
                print(df_file)
                path_interpolated_using_netatmo_dwd_unc = df_file

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

    df_netatmo_dwd_unc = pd.read_csv(path_interpolated_using_netatmo_dwd_unc,
                                     sep=';', index_col=0,
                                     parse_dates=True, infer_datetime_format=True)

    #########################################################

    df_compare = pd.DataFrame(index=df_dwd.index)

    try:

        for stn_ in df_dwd_edf.columns:

            for event_date in df_dwd.index.intersection(df_dwd_edf.index):
                print(event_date)
                edf_stn_orig = df_dwd_edf.loc[event_date, stn_]
                edf_stn_interp_dwd = df_dwd.loc[event_date, stn_]
                edf_stn_interp_netatmo_dwd = df_netatmo_dwd.loc[event_date, stn_]
                edf_stn_interp_netatmo_dwd_unc = df_netatmo_dwd_unc.loc[event_date, stn_]

                if ((edf_stn_orig >= 0) and
                    (edf_stn_interp_dwd >= 0) and
                        (edf_stn_interp_netatmo_dwd >= 0)
                        and (edf_stn_interp_netatmo_dwd_unc >= 0)):

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
                        'interpolated_quantile_netatmo_dwd_unc'] = edf_stn_interp_netatmo_dwd_unc

                else:
                    df_compare.loc[event_date,
                                   'original_quantile'] = np.nan
                    df_compare.loc[event_date,
                                   'interpolated_quantile_dwd'] = np.nan

                    df_compare.loc[event_date,
                                   'interpolated_quantile_netatmo_dwd'] = np.nan
                    df_compare.loc[event_date,
                                   'interpolated_quantile_netatmo_dwd_unc'] = np.nan

            df_compare = df_compare[df_compare > 0]
            df_compare.dropna(how='any', inplace=True)

            values_x = df_compare['original_quantile'].values
            values_dwd = df_compare['interpolated_quantile_dwd'].values

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
            corr_netatmo_dwd_unc = pears(values_x, values_netatmo_dwd_unc)[0]
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

#     # OK with Unc
    stations_with_improvements_unc = sum(i >= j for (i, j) in zip(
        df_improvements.pearson_corr_dwd_netatmo_unc.values,
        df_improvements.pearson_corr_dwd_.values))

    stations_without_improvements_unc = sum(i < j for (i, j) in zip(
        df_improvements.pearson_corr_dwd_netatmo_unc.values,
        df_improvements.pearson_corr_dwd_.values))

    percent_of_improvment_unc = 100 * (
        stations_with_improvements_unc /
        df_improvements.pearson_corr_dwd_netatmo_unc.shape[0])

    ####
    mean_pearson_correlation_dwd_only = df_improvements.pearson_corr_dwd_.mean()
    mean_pearson_correlation_dwd_netatmo = df_improvements.pearson_corr_dwd_netatmo.mean()
    mean_pearson_correlation_dwd_netatmo_unc = df_improvements.pearson_corr_dwd_netatmo_unc.mean()

    #########################################################
    mean_spr_correlation_dwd_only = df_improvements.spearman_corr_dwd_.mean()
    mean_spr_correlation_dwd_netatmo = df_improvements.spearman_corr_dwd_netatmo.mean()
    mean_spr_correlation_dwd_netatmo_unc = df_improvements.spearman_corr_dwd_netatmo_unc.mean()
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

    ax.plot(df_improvements.index,
            df_improvements.pearson_corr_dwd_netatmo_unc,
            alpha=.8,
            c='g',  # colors_arr,
            marker='+',
            label='DWD-Netatmo Interpolation Unc %0.2f'
            % mean_pearson_correlation_dwd_netatmo_unc)

    ax.set_title('Pearson Correlation Interpolated Quantiles from DWD or DWD-Netatmo\n '
                 'Precipitation of %s Extreme Events %s\n Events with Improvemnts %d / %d, Percentage %0.0f\n'
                 'Events with Improvemnts with OK Unc %d / %d, Percentage %0.0f'
                 % (temp_freq, _interp_acc_,
                    stations_with_improvements,
                     df_improvements.pearson_corr_dwd_netatmo.shape[0],
                    percent_of_improvment,
                    stations_with_improvements_unc,
                    df_improvements.pearson_corr_dwd_netatmo_unc.shape[
                        0],
                    percent_of_improvment_unc))
    ax.grid(alpha=0.25)
    plt.setp(ax.get_xticklabels(), rotation=45)
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

    # ok with Un
    stations_with_improvements_unc = sum(i >= j for (i, j) in zip(
        df_improvements.spearman_corr_dwd_netatmo_unc.values,
        df_improvements.spearman_corr_dwd_.values))

    stations_without_improvements_unc = sum(i < j for (i, j) in zip(
        df_improvements.spearman_corr_dwd_netatmo_unc.values,
        df_improvements.spearman_corr_dwd_.values))

    percent_of_improvment_unc = 100 * (stations_with_improvements_unc /
                                       df_improvements.spearman_corr_dwd_netatmo_unc.shape[0])
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
    ax.plot(df_improvements.index,
            df_improvements.spearman_corr_dwd_netatmo_unc,
            alpha=.8,
            c='g',  # colors_arr,
            marker='+',
            label='DWD-Netatmo Interpolation Unc %0.2f'
            % mean_spr_correlation_dwd_netatmo_unc)

    ax.set_title('Spearman Correlation Interpolated Quantiles from DWD or DWD-Netatmo \n '
                 'Rainfall of %s Extreme Events %s \n Events with Improvemnts %d / %d, Percentage %0.0f'
                 '\n Events with Improvemnts OK with Unc %d / %d, Percentage %0.0f'
                 % (temp_freq, _interp_acc_,
                    stations_with_improvements,
                     df_improvements.spearman_corr_dwd_netatmo.shape[0],
                    percent_of_improvment,
                    stations_with_improvements_unc,
                    df_improvements.spearman_corr_dwd_netatmo_unc.shape[
                        0],
                    percent_of_improvment_unc))
    ax.grid(alpha=0.25)
    plt.setp(ax.get_xticklabels(), rotation=45)
    ax.grid(alpha=0.25)
    ax.set_yticks(np.arange(0, 1.05, .10))
    ax.set_xlabel('DWD Stations')
    ax.legend(loc='lower right')
    ax.set_ylabel('Spearman Correlation')

    plt.savefig(os.path.join(path_to_use,
                             r'ppt_temporal_spearm_corr_%s_events_dwd.png' % (temp_freq)),
                frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
    plt.close()
