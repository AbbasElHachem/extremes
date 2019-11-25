
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

plt.ioff()
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'axes.labelsize': 12})


# In[2]:


# path_dwd_edf_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\edf_ppt_all_dwd_daily_.csv"
main_dir = Path(
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo')


min_orig_qnt_thr = 0.


# In[3]:


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


# In[13]:


for temp_freq in ['60min', '360min', '720min', '1440min']:
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

    path_to_Quantiles_netatmo_no_flt___ = main_dir / (
        r'Qt_ok_ok_un_netatmo_no_flt___%s' % temp_freq)
    Quantiles_netatmo_no_flt___ = list_all_full_path(
        '.csv', path_to_Quantiles_netatmo_no_flt___)

    #########################################################
    path_dwd_edf_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\edf_ppt_all_dwd_%s_.csv" % temp_freq

    df_dwd_edf = pd.read_csv(path_dwd_edf_data, sep=';', index_col=0,
                             parse_dates=True, infer_datetime_format=True)

    #########################################################
    df_improvements = pd.DataFrame(index=df_dwd_edf.columns,
                                   columns=['pearson_corr_dwd_', 'spearman_corr_dwd_',
                                            'pearson_corr_dwd_netatmo', 'spearman_corr_dwd_netatmo'])

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
    stns_dwd = df_dwd_edf.columns

#     path_interpolated_using_netatmo_dwd_list_un =
    df_combined_interp_dwd = pd.DataFrame(columns=stns_dwd)
    df_combined_interp_dwd_netatmo = pd.DataFrame(columns=stns_dwd)
    df_combined_interp_dwd_netatmo_unc = pd.DataFrame(columns=stns_dwd)
    for (path_interpolated_using_dwd,
         path_interpolated_using_netatmo_dwd,
         path_interpolated_using_netatmo_dwd_unc
         ) in zip(
            path_interpolated_using_dwd_list,
            path_interpolated_using_netatmo_dwd_list,
            path_interpolated_using_netatmo_dwd_list_un):

        df_dwd = pd.read_csv(path_interpolated_using_dwd,
                             sep=';', index_col=0,
                             parse_dates=True,
                             infer_datetime_format=True)

        df_netatmo_dwd = pd.read_csv(path_interpolated_using_netatmo_dwd,
                                     sep=';', index_col=0,
                                     parse_dates=True,
                                     infer_datetime_format=True)

        df_netatmo_dwd_unc = pd.read_csv(path_interpolated_using_netatmo_dwd_unc,
                                         sep=';', index_col=0,
                                         parse_dates=True,
                                         infer_datetime_format=True)
        for stn_ in df_dwd.columns:
            # print(stn_)
            for event_date in df_dwd.index[1:]:
                # interpolated_quantile_netatmo = df_netatmo.loc[event_date, stn_]

                interpolated_quantile_dwd = df_dwd.loc[event_date, stn_]

                df_combined_interp_dwd.loc[event_date,
                                           stn_] = interpolated_quantile_dwd

        for stn_ in df_netatmo_dwd.columns:
            print(stn_)
            for event_date in df_netatmo_dwd.index[1:]:
                interpolated_quantile_netatmo_dwd = df_netatmo_dwd.loc[event_date, stn_]
                df_combined_interp_dwd_netatmo.loc[event_date,
                                                   stn_] = interpolated_quantile_netatmo_dwd

        for stn_ in df_netatmo_dwd_unc.columns:
            print(stn_)
            for event_date in df_netatmo_dwd_unc.index[1:]:
                interpolated_quantile_netatmo_dwd_unc = df_netatmo_dwd_unc.loc[event_date, stn_]
                df_combined_interp_dwd_netatmo_unc.loc[event_date,
                                                       stn_] = interpolated_quantile_netatmo_dwd_unc
    df_combined_interp_dwd.to_csv(path_to_use / (
        r'df_combined_interp_dwd_%s_events_at_dwd_%s.csv' % (temp_freq, _interp_acc_)),
        sep=';')
    df_combined_interp_dwd_netatmo.to_csv(path_to_use / (
        r'df_combined_interp_dwd_netatmo_%s_events_at_dwd_%s.csv' % (temp_freq, _interp_acc_)),
        sep=';')
    df_combined_interp_dwd_netatmo_unc.to_csv(path_to_use / (
        r'df_combined_interp_dwd_netatmo_unc_%s_events_at_dwd_%s.csv' % (temp_freq, _interp_acc_)),
        sep=';')


# In[12]:


# path_interpolated_using_netatmo_dwd_unc
