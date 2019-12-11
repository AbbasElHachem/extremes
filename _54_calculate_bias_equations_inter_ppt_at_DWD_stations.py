
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


df_overall_bias = pd.DataFrame(index=['60min', '180min',
                                      '360min', '720min', '1440min'])
# '60min', '180min', '360min', '720min', '1440min']:
for temp_freq in ['60min']:
    # '60min', '360min',  '1440min'
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
#     path_dwd_edf_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\edf_ppt_all_dwd_%s_.csv" % temp_freq

#     df_dwd_edf = pd.read_csv(path_dwd_edf_data, sep=';', index_col=0,
#                              parse_dates=True, infer_datetime_format=True)

    path_dwd_ppt_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\ppt_all_dwd_%s_.csv" % temp_freq

    df_dwd_ppt = pd.read_csv(path_dwd_ppt_data, sep=';', index_col=0,
                             parse_dates=True, infer_datetime_format=True)

    #########################################################

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

        cmn_interpolated_events = df_netatmo_dwd.index.intersection(
            df_dwd.index).intersection(df_dwd_ppt.index)

        df_improvements = pd.DataFrame(
            index=df_dwd_ppt.columns)

#         df_improvements_b3 = pd.DataFrame(
#             index=df_dwd.index)

#         df_dwd_ppt = df_dwd_ppt.loc[df_dwd_ppt.index.intersection(
#             cmn_interpolated_events), :].dropna(how='all')

        print('Total number of events is',
              cmn_interpolated_events.shape[0])

        for stn_ in df_dwd_ppt.columns:
            #                 print(stn_)
            print(stn_)
            df_dwd_ppt_stn = df_dwd_ppt.loc[df_dwd_ppt.index.intersection(
                cmn_interpolated_events), stn_].dropna(how='all')
            print(df_dwd_ppt_stn.size)
            if df_dwd_ppt_stn.size > 0:
                df_compare = pd.DataFrame(index=df_dwd_ppt_stn.index)
                #stn_ = df_dwd.columns[0]
                for event_date in df_dwd_ppt_stn.index:
                    print(event_date)
                    # original ppt
                    original_ppt = df_dwd_ppt_stn.loc[event_date]

                    #event_date = cmn_interpolated_events[0]
                    interpolated_ppt_dwd = df_dwd.loc[event_date, stn_]
                    #print(stn_, interpolated_ppt_dwd, original_ppt)

                    interpolated_ppt_netatmo_dwd = df_netatmo_dwd.loc[event_date, stn_]

                    interpolated_ppt_netatmo_dwd_unc = df_netatmo_dwd_unc.loc[event_date, stn_]

                    if (interpolated_ppt_dwd >= 0) and (
                            interpolated_ppt_netatmo_dwd >= 0) and (
                            original_ppt >= 0):

                        df_compare.loc[event_date,
                                       'diff_dwd_dwd'] = (
                                           original_ppt - interpolated_ppt_dwd)
                        print(df_compare.loc[event_date,
                                             'diff_dwd_dwd'])
#                         df_compare.loc[event_date, 'diff_dwd_netatmo'] = (
#                             original_ppt -
#                             interpolated_ppt_netatmo_dwd)
#                         df_compare.loc[event_date, 'diff_dwd_netatmo_unc'] = (
#                             original_ppt -
#                             interpolated_ppt_netatmo_dwd_unc)
#
                df_compare.dropna(how='all', inplace=True)
                if df_compare.size > 0:
                    #             values_x_ppt = df_compare['original_ppt'].values
                    #             values_dwd_ppt = df_compare['interpolated_ppt_dwd'].values

                    #             values_netatmo_dwd_ppt = df_compare['interpolated_ppt_netatmo_dwd'].values
                    #             values_netatmo_dwd_unc_ppt = df_compare['interpolated_ppt_netatmo_dwd_unc'].values

                    diff_dwd_obs_interp = df_compare.loc[
                        :, 'diff_dwd_dwd'].values.sum()
#                     diff_dwd_obs_dwd_netatmo_interp = df_compare.loc[
#                         :, 'diff_dwd_netatmo'].values.sum()
#                     diff_dwd_obs_dwd_netatmo_interp_unc = df_compare.loc[
#                         :, 'diff_dwd_netatmo_unc'].values.sum()

                    # calculate sum per stn

                    df_improvements.loc[stn_,
                                        'dwd_interp_ppt_B1'] = diff_dwd_obs_interp
#                     df_improvements.loc[stn_,
#                                         'dwd_netatmo_interp_ppt_B1'] = diff_dwd_obs_dwd_netatmo_interp
#                     df_improvements.loc[stn_,
#                                         'dwd_netatmo_interp_unc_ppt_B1'] = diff_dwd_obs_dwd_netatmo_interp_unc

    #==========================================================================
    # CALCULATE TERMS
    #==========================================================================
        break
        df_improvements.dropna(how='any', inplace=True)

        B1_dwd = (
            np.square(
                df_improvements.loc[:,
                                    'dwd_interp_ppt_B1'].values.sum()) /
            (cmn_interpolated_events.shape[0] * 111))

#     B1_netatmo = (
#         np.square(
#             df_improvements.loc[:,
#                                 'dwd_netatmo_interp_ppt_B1'].values.sum()) /
#         (cmn_interpolated_events.shape[0] * 111))
#
#     B1_netatmo_unc = (
#         np.square(
#             df_improvements.loc[:, 'dwd_netatmo_interp_unc_ppt_B1'].values.sum()) /
#         (cmn_interpolated_events.shape[0] * 111))

        df_overall_bias.loc[temp_freq, 'B1_dwd'] = B1_dwd

        print(temp_freq, ' B1_dwd ', B1_dwd)
#     df_overall_bias.loc[temp_freq, 'B1_dwd_netatmo'] = B1_netatmo
#     df_overall_bias.loc[temp_freq, 'B1_dwd_netatmo_unc'] = B1_netatmo_unc


# #==============================================================================
#
#==============================================================================
print('Done')
df_overall_bias.to_csv(
    Path(out_dir_bias_dfs / ('%s.csv' % _interp_acc_)), sep=';')
