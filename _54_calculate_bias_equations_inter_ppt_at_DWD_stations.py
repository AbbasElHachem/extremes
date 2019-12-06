
# coding: utf-8

# In[1]:

'''


'''

import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

import fnmatch
from pathlib import Path

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

plt.ioff()
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'axes.labelsize': 12})


# In[12]:


main_dir = Path(
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo')


out_dir_bias_dfs = main_dir / 'bias_calculation'
if not os.path.exists(out_dir_bias_dfs):
    os.mkdir(out_dir_bias_dfs)
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

# In[26]:


df_overall_bias = pd.DataFrame(index=['60min', '360min',
                                      '720min', '1440min'])

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

    #########################################################
    path_dwd_edf_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\edf_ppt_all_dwd_%s_.csv" % temp_freq

    df_dwd_edf = pd.read_csv(path_dwd_edf_data, sep=';', index_col=0,
                             parse_dates=True, infer_datetime_format=True)

    path_dwd_ppt_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\ppt_all_dwd_%s_.csv" % temp_freq

    df_dwd_ppt = pd.read_csv(path_dwd_ppt_data, sep=';', index_col=0,
                             parse_dates=True, infer_datetime_format=True)

    #########################################################

    path_to_use = path_to_Quantiles_netatmo_no_flt___
    data_to_use = Quantiles_netatmo_no_flt___

    _interp_acc_ = str(r'%s' % (str(path_to_use).split('\\')[-1]))
    # for i in range(12):
    #   i = int(i)
    try:
        for df_file in data_to_use:

            if ('dwd_only') in df_file:
                path_interpolated_using_dwd = df_file

            if ('dwd_netamo') in df_file and 'unc' not in df_file:
                path_interpolated_using_netatmo_dwd = df_file

            if ('dwd_netamo') in df_file:
                path_interpolated_using_netatmo_dwd_unc = df_file

    except Exception as msg:
        print(msg)
        raise Exception

    #########################################################

    df_dwd = pd.read_csv(path_interpolated_using_dwd,
                         sep=';', index_col=0, parse_dates=True,
                         infer_datetime_format=True)

    df_netatmo_dwd = pd.read_csv(path_interpolated_using_netatmo_dwd,
                                 sep=';', index_col=0,
                                 parse_dates=True,
                                 infer_datetime_format=True)

#     df_netatmo_dwd_unc = pd.read_csv(path_interpolated_using_netatmo_dwd_unc,
#                                      sep=';', index_col=0,
#                                      parse_dates=True,
#                                      infer_datetime_format=True)

    nbr_stns = 111

    #==========================================================================
    #
    #==========================================================================

    def calculate_overall_bias(df_obsv, df_interp):
        nbr_events = []
        sum_over_stns = []
        for stn in df_interp.columns:

            if df_interp.loc[:, stn].dropna(how='all').values.size > 0:
                sum_over_events = []

                for event_date in df_interp.index:

                    obsv = df_obsv.loc[event_date, stn]
                    interp = df_interp.loc[event_date, stn]

                    if ((obsv >= 0) and (interp >= 0)):
                        sum_over_events.append(obsv - interp)

            if len(sum_over_events) > 0:
                nbr_events.append(len(sum_over_events))
                sum_over_stns.append(np.sum(sum_over_events))

        mean_nbr_events = np.mean(nbr_events)
        # print('mean_nbr_events: ', mean_nbr_events)
        B1 = np.round((1 / (nbr_stns * mean_nbr_events)) *
                      np.square(np.sum(sum_over_stns)), 6)
        return B1

    #==========================================================================
    #
    #==========================================================================
    def calculate_temporal_bias(df_obsv, df_interp):
        nbr_events = []
        sum_over_stns = []
        for stn in df_interp.columns:
            if df_interp.loc[:, stn].dropna(how='all').values.size > 0:
                sum_over_events = []

                for event_date in df_interp.index:
                    obsv = df_obsv.loc[event_date, stn]
                    interp = df_interp.loc[event_date, stn]

                    if ((obsv >= 0) and (interp >= 0)):
                        sum_over_events.append(obsv - interp)

            if len(sum_over_events) > 0:
                nbr_events.append(len(sum_over_events))
                sum_over_stns.append(
                    np.square(np.sum(sum_over_events) / np.mean(nbr_events)))
        # print('mean_nbr_events: ', np.mean(nbr_events))
        B2 = np.round((1 / (nbr_stns)) * (np.sum(sum_over_stns)), 6)
        return B2
    #==========================================================================
    #
    #==========================================================================

    def calculate_spatial_bias(df_obsv, df_interp):
        nbr_events = 0
        sum_over_events = []
        for event_date in df_interp.index:
            sum_over_stns = []
            for stn in df_interp.columns:
                if df_interp.loc[:, stn].dropna(how='all').values.size > 0:

                    obsv = df_obsv.loc[event_date, stn]
                    interp = df_interp.loc[event_date, stn]

                    if ((obsv >= 0) and (interp >= 0)):
                        sum_over_stns.append(obsv - interp)

            if len(sum_over_stns) > 0:

                sum_over_events.append(
                    np.square(np.sum(sum_over_stns) / nbr_stns))

                nbr_events += 1
        # print('mean_nbr_events: ', nbr_events)
        B3 = np.round((1 / (nbr_events)) * (np.sum(sum_over_events)), 6)
        return B3

    #==========================================================================
    #
    #==========================================================================

    def calculate_squared_error(df_obsv, df_interp):
        nbr_events = []
        sum_over_stns = []
        for stn in df_interp.columns:

            if df_interp.loc[:, stn].dropna(how='all').values.size > 0:
                sum_over_events = []

                for event_date in df_interp.index:

                    obsv = df_obsv.loc[event_date, stn]
                    interp = df_interp.loc[event_date, stn]

                    if ((obsv >= 0) and (interp >= 0)):
                        sum_over_events.append(
                            np.square(obsv - interp))

            if len(sum_over_events) > 0:
                nbr_events.append(len(sum_over_events))
                sum_over_stns.append(np.sum(sum_over_events))

        mean_nbr_events = np.mean(nbr_events)
        # print('mean_nbr_events: ', mean_nbr_events)
        B4 = np.round((1 / (nbr_stns * mean_nbr_events)) *
                      np.sum(sum_over_stns), 6)
        return B4

    #==========================================================================
    #
    #==========================================================================

    B1_dwd = calculate_overall_bias(df_dwd_ppt, df_dwd)
    B1_dwd_netatmo = calculate_overall_bias(df_dwd_ppt, df_netatmo_dwd)
#     B1_dwd_netatmo_unc = calculate_overall_bias(df_dwd_edf, df_netatmo_dwd_unc)
    print('*\n+-B1 DWD Interp:', B1_dwd)
    print('*\n+-B1 DWD Netatmo Interp:', B1_dwd_netatmo)
#     print('*\n+-B1 DWD Netatmo Unc Interp:', B1_dwd_netatmo_unc)

    B2_dwd = calculate_temporal_bias(df_dwd_ppt, df_dwd)
    B2_dwd_netatmo = calculate_temporal_bias(df_dwd_ppt, df_netatmo_dwd)
#     B2_dwd_netatmo_unc = calculate_temporal_bias(
#         df_dwd_edf, df_netatmo_dwd_unc)
    print('*\n+-B2 DWD Interp:', B2_dwd)
    print('*\n+-B2 DWD Netatmo Interp:', B2_dwd_netatmo)
#     print('*\n+-B2 DWD Netatmo Unc Interp:', B2_dwd_netatmo_unc)

    B3_dwd = calculate_spatial_bias(df_dwd_ppt, df_dwd)
    B3_dwd_netatmo = calculate_spatial_bias(df_dwd_ppt, df_netatmo_dwd)
#     B3_dwd_netatmo_unc = calculate_spatial_bias(df_dwd_edf, df_netatmo_dwd_unc)
    print('*\n+-B3 DWD Interp:', B3_dwd)
    print('*\n+-B3 DWD Netatmo Interp:', B3_dwd_netatmo)
#     print('*\n+-B3 DWD Netatmo Unc Interp:', B3_dwd_netatmo_unc)

    B4_dwd = calculate_squared_error(df_dwd_ppt, df_dwd)
    B4_dwd_netatmo = calculate_squared_error(df_dwd_ppt, df_netatmo_dwd)
#     B4_dwd_netatmo_unc = calculate_squared_error(
#         df_dwd_edf, df_netatmo_dwd_unc)
    print('*\n+-B4 DWD Interp:', B4_dwd)
    print('*\n+-B4 DWD Netatmo Interp:', B4_dwd_netatmo)
#     print('*\n+-B4 DWD Netatmo Unc Interp:', B4_dwd_netatmo_unc)

    df_overall_bias.loc[temp_freq, 'B1_dwd'] = B1_dwd
    df_overall_bias.loc[temp_freq, 'B1_dwd_netatmo'] = B1_dwd_netatmo
#     df_overall_bias.loc[temp_freq, 'B1_dwd_netatmo_unc'] = B1_dwd_netatmo_unc

    df_overall_bias.loc[temp_freq, 'B2_dwd'] = B2_dwd
    df_overall_bias.loc[temp_freq, 'B2_dwd_netatmo'] = B2_dwd_netatmo
#     df_overall_bias.loc[temp_freq, 'B2_dwd_netatmo_unc'] = B2_dwd_netatmo_unc

    df_overall_bias.loc[temp_freq, 'B3_dwd'] = B3_dwd
    df_overall_bias.loc[temp_freq, 'B3_dwd_netatmo'] = B3_dwd_netatmo
#     df_overall_bias.loc[temp_freq, 'B3_dwd_netatmo_unc'] = B3_dwd_netatmo_unc

    df_overall_bias.loc[temp_freq, 'B4_dwd'] = B4_dwd
    df_overall_bias.loc[temp_freq, 'B4_dwd_netatmo'] = B4_dwd_netatmo
#     df_overall_bias.loc[temp_freq, 'B4_dwd_netatmo_unc'] = B4_dwd_netatmo_unc

#==============================================================================
#
#==============================================================================
print('Done')
df_overall_bias.to_csv(
    Path(out_dir_bias_dfs / ('%s.csv' % _interp_acc_)), sep=';')
