
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
# '60min', '180min', '360min', '720min', '1440min']: '60min', '180min', '360min', '720min',
for temp_freq in ['60min', '180min', '360min', '720min', '1440min']:
    # '60min', '360min',  '1440min'
    print(temp_freq)

    path_to_Qt_ok_un_first_flt__temp_flt_1st_ = main_dir / (
        r'Ppt_ok_ok_un_new2_first_flt__temp_flt__1st_%s' % temp_freq)
    Qt_ok_un_first_flt__temp_flt_1st_ = list_all_full_path(
        '.csv', path_to_Qt_ok_un_first_flt__temp_flt_1st_)

    path_to_Qt_ok_un_first_flt__temp_flt_comb_ = main_dir / (
        r'Ppt_ok_ok_un_new2_first_flt__temp_flt__comb_%s' % temp_freq)
    Qt_ok_un_first_flt__temp_flt_comb_ = list_all_full_path(
        '.csv', path_to_Qt_ok_un_first_flt__temp_flt_comb_)

    path_to_Qt_ok_un_first_flt_1st_ = main_dir / (
        r'Ppt_ok_ok_un_new2_first_flt__1st_%s' % temp_freq)
    Qt_ok_un_first_flt_1st_ = list_all_full_path(
        '.csv', path_to_Qt_ok_un_first_flt_1st_)

    path_to_Qt_ok_un_first_flt_comb_ = main_dir / (
        r'Ppt_ok_ok_un_new2_first_flt__comb_%s' % temp_freq)
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

    path_to_use = path_to_Qt_ok_un_first_flt__temp_flt_1st_
    data_to_use = Qt_ok_un_first_flt__temp_flt_1st_

    _interp_acc_ = str(r'%s' % (str(path_to_use).split('\\')[-1]))
    # for i in range(12):
    #   i = int(i)
    try:
        for df_file in data_to_use:

            if ('using_dwd_only_grp_') in df_file:
                print(df_file)
                path_interpolated_using_dwd = df_file
            if ('using_dwd_netamo_grp_') in df_file and 'perc' not in df_file:

                print(df_file)
                path_interpolated_using_netatmo_dwd = df_file

            if ('using_dwd_netamo_grp_') in df_file and 'unc10perc' in df_file:
                print(df_file)
                path_interpolated_using_netatmo_dwd_unc = df_file

    except Exception as msg:
        print(msg)
        continue
        #########################################################

    df_dwd = pd.read_csv(path_interpolated_using_dwd, skiprows=[1],
                         sep=';', index_col=0,
                         parse_dates=True,
                         infer_datetime_format=True).dropna(how='all')

    df_dwd_ppt = df_dwd_ppt.loc[df_dwd_ppt.index.intersection(
        df_dwd.index), :].dropna(how='all')

#     print(df_dwd.describe())

    df_dwd_ppt = df_dwd_ppt[df_dwd.columns]

#     print(df_dwd_ppt.describe())
    pass

    # for stn for event calculate difference

    df_netatmo_dwd = pd.read_csv(path_interpolated_using_netatmo_dwd,
                                 skiprows=[1],
                                 sep=';', index_col=0,
                                 parse_dates=True,
                                 infer_datetime_format=True).dropna(how='all')

    try:
        df_netatmo_dwd_unc = pd.read_csv(path_interpolated_using_netatmo_dwd_unc,
                                         skiprows=[1],
                                         sep=';', index_col=0,
                                         parse_dates=True,
                                         infer_datetime_format=True).dropna(how='all')
    except Exception as msg:
        df_netatmo_dwd_unc = df_netatmo_dwd
#     cmn_interpolated_events = df_dwd.index.intersection(df_dwd_ppt.index)
#
    df_sum_per_stn = pd.DataFrame(index=df_dwd_ppt.columns)

    df_improvements_b3 = pd.DataFrame(
        index=df_dwd.index)

    cmn_events = df_dwd_ppt.index.intersection(
        df_dwd.index).intersection(df_netatmo_dwd.index).intersection(
            df_netatmo_dwd_unc.index)

    plt.ioff()

    fig = plt.figure(figsize=(24, 12), dpi=200)

    ax = fig.add_subplot(111)

    fig2 = plt.figure(figsize=(24, 12), dpi=200)

    ax2 = fig2.add_subplot(111)

#     fig3 = plt.figure(figsize=(24, 12), dpi=200)
#
#     ax3 = fig3.add_subplot(111)

    for stn_ in df_dwd_ppt.columns:

        df_dwd_ppt_stn_obs = df_dwd_ppt.loc[:, stn_].dropna(how='all')
        df_dwd_ppt_stn_interp = df_dwd.loc[:, stn_].dropna(how='all')
        # .dropna(how='all')
        df_netatmo_ppt_stn_interp = df_netatmo_dwd.loc[:, stn_].dropna(
            how='all')

        df_netatmo_unc_ppt_stn_interp = df_netatmo_dwd_unc.loc[:, stn_].dropna(
            how='all')
#         df_netatmo_ppt_stn_interp.dropna(how='all', inplace=True)
#         df_dwd_ppt_stn_obs = df_dwd_ppt_stn_obs.loc[
#             df_netatmo_ppt_stn_interp.index]
#         df_dwd_ppt_stn_interp = df_dwd_ppt_stn_interp.loc[
#             df_netatmo_ppt_stn_interp.index]
        if (df_dwd_ppt_stn_obs.size > 0) and (
                df_dwd_ppt_stn_interp.size > 0) and (
                df_netatmo_ppt_stn_interp.size > 0) and (
                df_netatmo_unc_ppt_stn_interp.size > 0):

            df_sum_events = pd.DataFrame(index=df_dwd_ppt_stn_obs.index)

            for event_date in df_netatmo_ppt_stn_interp.index:
                if (event_date in df_dwd_ppt_stn_obs.index) and (
                        event_date in df_dwd_ppt_stn_interp.index) and (
                        event_date in df_netatmo_unc_ppt_stn_interp.index):

                    obsv_ppt = df_dwd_ppt_stn_obs.loc[event_date]
                    interp_ppt = df_dwd_ppt_stn_interp.loc[event_date]
                    interp_ppt_netatmo = df_netatmo_ppt_stn_interp.loc[event_date]
                    interp_ppt_netatmo_unc = df_netatmo_unc_ppt_stn_interp.loc[event_date]

                    if (obsv_ppt >= 0) and (
                            interp_ppt >= 0) and (
                            interp_ppt_netatmo >= 0) and (
                            interp_ppt_netatmo_unc >= 0):

                        diff_obsv_interp = obsv_ppt - interp_ppt
                        diff_obsv_interp_netatmo = obsv_ppt - interp_ppt_netatmo
                        diff_obsv_interp_netatmo_unc = obsv_ppt - interp_ppt_netatmo_unc

                        df_sum_events.loc[event_date,
                                          'diff_dwd_dwd'] = diff_obsv_interp
                        df_sum_events.loc[event_date,
                                          'diff_dwd_netatmo'] = diff_obsv_interp_netatmo
                        df_sum_events.loc[event_date,
                                          'diff_dwd_netatmo_unc'] = diff_obsv_interp_netatmo_unc

                        diff_obsv_interp_abs = np.abs(obsv_ppt - interp_ppt)
                        diff_obsv_interp_netatmo_abs = np.abs(
                            obsv_ppt - interp_ppt_netatmo)
                        diff_obsv_interp_netatmo_unc_abs = np.abs(
                            obsv_ppt - interp_ppt_netatmo_unc)

                        df_sum_events.loc[event_date,
                                          'diff_dwd_dwd_abs'] = diff_obsv_interp_abs
                        df_sum_events.loc[event_date,
                                          'diff_dwd_netatmo_abs'] = diff_obsv_interp_netatmo_abs
                        df_sum_events.loc[event_date,
                                          'diff_dwd_netatmo_unc_abs'] = diff_obsv_interp_netatmo_unc_abs

                    else:
                        df_sum_events.loc[event_date,
                                          'diff_dwd_dwd'] = np.nan
                        df_sum_events.loc[event_date,
                                          'diff_dwd_netatmo'] = np.nan

                        df_sum_events.loc[event_date,
                                          'diff_dwd_netatmo_unc'] = np.nan
                        df_sum_events.loc[event_date,
                                          'diff_dwd_dwd_abs'] = np.nan
                        df_sum_events.loc[event_date,
                                          'diff_dwd_netatmo_abs'] = np.nan

                        df_sum_events.loc[event_date,
                                          'diff_dwd_netatmo_unc_abs'] = np.nan

            df_sum_events.dropna(how='all', inplace=True)

            if df_sum_events.size > 0:
                #                 ax.plot(df_sum_events.loc[:, 'diff_dwd_dwd'].index,
                #                         df_sum_events.loc[:, 'diff_dwd_dwd'].values,
                #                         # marker='+',
                #                         color='b', alpha=0.1)
                #                 ax.plot(df_sum_events.loc[:, 'diff_dwd_netatmo'].index,
                #                         df_sum_events.loc[:, 'diff_dwd_netatmo'].values,
                #                         # marker='d',
                #                         color='r', alpha=0.1)

                #                 df_sums_sorted = df_sum_events.sort_values(
                #                     by=['diff_dwd_netatmo'])
                #
                #                 idx_int = np.arange(0,
                #                                     df_sums_sorted.loc[:,
                #                                                        'diff_dwd_netatmo'].index.size, 1)
                #                 ax3.plot(idx_int,
                #                          df_sums_sorted.loc[:, 'diff_dwd_dwd'].values,
                #                          # marker='+',
                #                          color='b', alpha=0.25)
                #                 ax3.plot(idx_int,
                #                          df_sums_sorted.loc[:, 'diff_dwd_netatmo'].values,
                #                          # marker='d',
                #                          color='r', alpha=0.25)

                #                 ax.plot(df_sum_events.loc[:, 'diff_dwd_netatmo_unc'].index,
                #                         df_sum_events.loc[:, 'diff_dwd_netatmo_unc'].values,
                #                         marker='o',
                #                         color='r', alpha=0.1)

                # difference obs-inter square stn for B1
                diff_dwd_obs_interp = df_sum_events.loc[
                    :, 'diff_dwd_dwd'].values.sum()
                diff_netatmo_obs_interp = df_sum_events.loc[
                    :, 'diff_dwd_netatmo'].values.sum()
                diff_netatmo_unc_obs_interp = df_sum_events.loc[
                    :, 'diff_dwd_netatmo_unc'].values.sum()

                # for B%

                diff_dwd_obs_interp_abs = df_sum_events.loc[
                    :, 'diff_dwd_dwd_abs'].values.sum()
                diff_netatmo_obs_interp_abs = df_sum_events.loc[
                    :, 'diff_dwd_netatmo_abs'].values.sum()
                diff_netatmo_unc_obs_interp_abs = df_sum_events.loc[
                    :, 'diff_dwd_netatmo_unc_abs'].values.sum()

                # difference obs-inter square stn for B2
                diff_dwd_obs_interp_b2 = (
                    diff_dwd_obs_interp /
                    df_dwd_ppt.index.shape[0])**2

                diff_dwd_obs_dwd_netatmo_interp_b2 = (
                    diff_netatmo_obs_interp /
                    df_dwd_ppt.index.shape[0])**2

                diff_dwd_obs_dwd_netatmo_interp_unc_b2 = (
                    diff_netatmo_unc_obs_interp /
                    df_dwd_ppt.index.shape[0])**2

                # difference obs-inter square for B4
                diff_dwd_obs_interp_square = sum(
                    df_sum_events.loc[:, 'diff_dwd_dwd'].values**2)

                diff_dwd_obs_dwd_netatmo_interp_square = np.sum(
                    df_sum_events.loc[:, 'diff_dwd_netatmo'].values**2)

                diff_dwd_obs_dwd_netatmo_interp_unc_square = np.sum(
                    df_sum_events.loc[:, 'diff_dwd_netatmo_unc'].values**2)

                #==============================================================
                # B1
                #==============================================================
                df_sum_per_stn.loc[
                    stn_,
                    'dwd_interp_ppt_B1'] = diff_dwd_obs_interp

                df_sum_per_stn.loc[
                    stn_,
                    'dwd_netatmo_interp_ppt_B1'] = diff_netatmo_obs_interp

                df_sum_per_stn.loc[
                    stn_,
                    'dwd_netatmo_unc_interp_ppt_B1'] = diff_netatmo_unc_obs_interp

                #==============================================================
                # B% _abs
                #==============================================================
                df_sum_per_stn.loc[
                    stn_,
                    'dwd_interp_ppt_B1_abs'] = diff_dwd_obs_interp_abs

                df_sum_per_stn.loc[
                    stn_,
                    'dwd_netatmo_interp_ppt_B1_abs'] = diff_netatmo_obs_interp_abs

                df_sum_per_stn.loc[
                    stn_,
                    'dwd_netatmo_unc_interp_ppt_B1_abs'] = diff_netatmo_unc_obs_interp_abs

                #==============================================================
                # B2
                #==============================================================
                df_sum_per_stn.loc[
                    stn_,
                    'dwd_interp_ppt_B2'] = diff_dwd_obs_interp_b2
                df_sum_per_stn.loc[
                    stn_,
                    'dwd_netatmo_interp_ppt_B2'] = diff_dwd_obs_dwd_netatmo_interp_b2
                df_sum_per_stn.loc[
                    stn_,
                    'dwd_netatmo_interp_unc_ppt_B2'] = diff_dwd_obs_dwd_netatmo_interp_unc_b2
                #==============================================================
                # B4
                #==============================================================
                df_sum_per_stn.loc[
                    stn_,
                    'dwd_interp_ppt_square_B4'] = diff_dwd_obs_interp_square
                df_sum_per_stn.loc[
                    stn_,
                    'dwd_netatmo_interp_ppt_square_B4'] = diff_dwd_obs_dwd_netatmo_interp_square
                df_sum_per_stn.loc[
                    stn_,
                    'dwd_netatmo_interp_unc_ppt_square_B4'] = diff_dwd_obs_dwd_netatmo_interp_unc_square

#     ax.grid(alpha=0.25)
#     ax.set_ylim([-200, 200])
#     ax.set_title(
#         'Difference Obsv-Interp for every station;'
#         '\n%s'
#         '\nDWD:blue, DWD-Netatmo:red' % _interp_acc_)
#     ax.set_xlabel('Extreme Events')
#     #ax.legend(loc='upper right')
#     ax.set_ylabel('(Z-Z*) mm/%s' % temp_freq)
#
#     fig.savefig((out_dir_bias_dfs / (
#         r'diff_obs_int_%s_events_dwd_%s_.png' % (temp_freq, _interp_acc_))),
#         frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)

#     ax3.grid(alpha=0.25)
#     ax3.set_ylim([-200, 200])
#     ax3.set_title(
#         'Sorted Difference Obsv-Interp for every station;'
#         '\n%s'
#         '\nDWD:blue, DWD-Netatmo:red' % _interp_acc_)
#     ax3.set_xlabel('Event number')
#     #ax.legend(loc='upper right')
#     ax3.set_ylabel('(Z-Z*) mm/%s' % temp_freq)
#
#     fig3.savefig((out_dir_bias_dfs / (
#         r'sorted_diff_obs_int_%s_events_dwd_%s_.png' % (temp_freq, _interp_acc_))),
#         frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)

    df_sum_per_stn.dropna(how='all', inplace=True)

#     ax2.plot(df_sum_per_stn.index,
#              df_sum_per_stn.loc[:, 'dwd_interp_ppt_B1'].values,
#              color='b', marker='+', alpha=0.75,
#              label='DWD %0.1f' % df_sum_per_stn.loc[
#                  :,
#                  'dwd_interp_ppt_B1'].values.sum())
#     mean_dwd_val = df_sum_per_stn.loc[:, 'dwd_interp_ppt_B1'].mean()
#     mean_dwd_vals = [mean_dwd_val for i in range(len(df_sum_per_stn.index))]
#
#     ax2.plot(df_sum_per_stn.index,
#              mean_dwd_vals, marker='_',
#              color='darkblue', alpha=0.75)
#
#     ax2.plot(df_sum_per_stn.index,
#              df_sum_per_stn.loc[:, 'dwd_netatmo_interp_ppt_B1'].values,
#              color='g', marker='o', alpha=0.75,
#              label='DWD-Netatmo %0.1f' % df_sum_per_stn.loc[
#                  :,
#                  'dwd_netatmo_interp_ppt_B1'].values.sum())
#
#     mean_netatmo_val = df_sum_per_stn.loc[:,
#                                           'dwd_netatmo_interp_ppt_B1'].mean()
#     mean_netatmo_vals = [
#         mean_netatmo_val for i in range(len(df_sum_per_stn.index))]
#
#     ax2.plot(df_sum_per_stn.index,
#              mean_netatmo_vals, marker='_',
#              color='darkgreen', alpha=0.75)
#
#     ax2.plot(df_sum_per_stn.index,
#              df_sum_per_stn.loc[:, 'dwd_netatmo_unc_interp_ppt_B1'].values,
#              color='r', marker='d', alpha=0.75,
#              label='DWD-Netatmo Unc %0.1f' % df_sum_per_stn.loc[
#                  :,
#                  'dwd_netatmo_unc_interp_ppt_B1'].values.sum())
#
#     ax2.legend(loc='upper right')
#     ax2.grid(alpha=0.25)
#     #ax2.set_ylim([-400, 200])
#     ax2.set_title(
#         'Sum Difference Obsv-Interp for every station;'
#         '%s'
#         ' \nDWD:blue, DWD-Netatmo:green, DWD-NEtatmo Und:red' % _interp_acc_)
#     ax2.set_xlabel('Extreme Events')
#     ax2.set_xticklabels(labels=df_sum_per_stn.index, rotation=45)
#     #ax.legend(loc='upper right')
#     ax2.set_ylabel('sum(Z-Z*) mm/%s' % temp_freq)
#
#     fig2.savefig((out_dir_bias_dfs / (
#         r'sum_diff_obs_int_%s_events_dwd_%s_.png' % (temp_freq, _interp_acc_))),
#         frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
#     plt.close()

    #======================================================================
    # calculte B3
    #======================================================================

    for event_date in cmn_events:
        df_compare_b3 = pd.DataFrame(index=df_dwd.columns)
        for stn_ in df_dwd.columns:

            interpolated_ppt_dwd = df_dwd.loc[event_date, stn_]

            interpolated_ppt_netatmo_dwd = df_netatmo_dwd.loc[event_date, stn_]

            interpolated_ppt_netatmo_dwd_unc = df_netatmo_dwd_unc.loc[event_date, stn_]
            # original quantile when transforming ppt to edf
            original_ppt = df_dwd_ppt.loc[event_date, stn_]

            if (interpolated_ppt_dwd >= 0) and (
                    interpolated_ppt_netatmo_dwd >= 0) and (
                        original_ppt >= 0):

                df_compare_b3.loc[stn_,
                                  'original_ppt'] = original_ppt

                df_compare_b3.loc[stn_,
                                  'interpolated_ppt_dwd'] = interpolated_ppt_dwd

                df_compare_b3.loc[stn_,
                                  'interpolated_ppt_netatmo_dwd'] = interpolated_ppt_netatmo_dwd
                df_compare_b3.loc[stn_,
                                  'interpolated_ppt_netatmo_dwd_unc'] = interpolated_ppt_netatmo_dwd_unc
                df_compare_b3.loc[stn_,
                                  'diff_dwd_dwd'] = original_ppt - interpolated_ppt_dwd

                df_compare_b3.loc[stn_, 'diff_dwd_netatmo'] = original_ppt - \
                    interpolated_ppt_netatmo_dwd
                df_compare_b3.loc[stn_, 'diff_dwd_netatmo_unc'] = original_ppt - \
                    interpolated_ppt_netatmo_dwd_unc
        # df_compare_b3 = df_compare_b3[df_compare_b3 >= 0]
        df_compare_b3.dropna(how='all', inplace=True)

        # difference obs-inter for B3

        diff_dwd_obs_interp_B3 = np.square(
            df_compare_b3.loc[:,
                              'diff_dwd_dwd'].sum() / df_compare_b3.index.shape[0])
        diff_dwd_obs_dwd_netatmo_interp_B3 = np.square(
            df_compare_b3.loc[:, 'diff_dwd_netatmo'].sum() /
            df_compare_b3.index.shape[0])

        diff_dwd_obs_dwd_netatmo_interp_unc_B3 = np.square(
            df_compare_b3.loc[:, 'diff_dwd_netatmo_unc'].sum() /
            df_compare_b3.index.shape[0])

        df_improvements_b3.loc[
            event_date,
            'dwd_interp_ppt_B3'] = diff_dwd_obs_interp_B3
        df_improvements_b3.loc[
            event_date,
            'dwd_netatmo_interp_ppt_B3'] = diff_dwd_obs_dwd_netatmo_interp_B3
        df_improvements_b3.loc[
            event_date,
            'dwd_netatmo_interp_unc_ppt_B3'] = diff_dwd_obs_dwd_netatmo_interp_unc_B3

    #==========================================================================
    # calculate TERMS
    #==========================================================================

    B1_dwd = ((df_sum_per_stn.loc[:, 'dwd_interp_ppt_B1'].values.sum()**2) /
              (df_dwd_ppt.index.shape[0] * 111))

    B1_netatmo = ((df_sum_per_stn.loc[:, 'dwd_netatmo_interp_ppt_B1'].values.sum()**2) /
                  (df_dwd_ppt.index.shape[0] * 111))

    B1_netatmo_unc = ((df_sum_per_stn.loc[:, 'dwd_netatmo_unc_interp_ppt_B1'].values.sum()**2) /
                      (df_dwd_ppt.index.shape[0] * 111))

    B6_dwd = (np.abs(df_sum_per_stn.loc[:, 'dwd_interp_ppt_B1'].values.sum()) /
              (df_dwd_ppt.index.shape[0] * 111))

    B6_netatmo = np.abs((df_sum_per_stn.loc[:, 'dwd_netatmo_interp_ppt_B1'].values.sum()) /
                        (df_dwd_ppt.index.shape[0] * 111))

    B6_netatmo_unc = np.abs((df_sum_per_stn.loc[:, 'dwd_netatmo_unc_interp_ppt_B1'].values.sum()) /
                            (df_dwd_ppt.index.shape[0] * 111))

    # B!_abs
    B1_dwd_abs = (df_sum_per_stn.loc[:, 'dwd_interp_ppt_B1_abs'].values.sum()**2 /
                  (df_dwd_ppt.index.shape[0] * 111))

    B1_netatmo_abs = (df_sum_per_stn.loc[:, 'dwd_netatmo_interp_ppt_B1_abs'].values.sum()**2 /
                      (df_dwd_ppt.index.shape[0] * 111))

    B1_netatmo_unc_abs = (df_sum_per_stn.loc[:, 'dwd_netatmo_unc_interp_ppt_B1_abs'].values.sum()**2 /
                          (df_dwd_ppt.index.shape[0] * 111))

    # B2
    B2_dwd = (
        df_sum_per_stn.loc[
            :,
            'dwd_interp_ppt_B2'].dropna().values.sum() /
        (df_sum_per_stn.index.shape[0]))

    B2_netatmo = (
        df_sum_per_stn.loc[
            :,
            'dwd_netatmo_interp_ppt_B2'].dropna().values.sum() /
        (df_sum_per_stn.index.shape[0]))

    B2_netatmo_unc = (
        df_sum_per_stn.loc[
            :,
            'dwd_netatmo_interp_unc_ppt_B2'].dropna().values.sum() /
        (df_sum_per_stn.index.shape[0]))
#     B3
    B3_dwd = (
        df_improvements_b3.loc[
            :,
            'dwd_interp_ppt_B3'].dropna().values.sum() /
        (cmn_events.shape[0]))

    B3_netatmo = (
        df_improvements_b3.loc[
            :,
            'dwd_netatmo_interp_ppt_B3'].dropna().values.sum() /
        (cmn_events.shape[0]))
    B3_netatmo_unc = (
        df_improvements_b3.loc[
            :,
            'dwd_netatmo_interp_unc_ppt_B3'].dropna().values.sum() /
        (cmn_events.shape[0]))

    # B4
    B4_dwd = (
        df_sum_per_stn.loc[
            :,
            'dwd_interp_ppt_square_B4'].dropna().values.sum() /
        (cmn_events.shape[0] * df_sum_per_stn.index.shape[0]))

    B4_netatmo = (
        df_sum_per_stn.loc[
            :,
            'dwd_netatmo_interp_ppt_square_B4'].dropna().values.sum() /
        (cmn_events.shape[0] * df_sum_per_stn.index.shape[0]))

    B4_netatmo_unc = (
        df_sum_per_stn.loc[
            :,
            'dwd_netatmo_interp_unc_ppt_square_B4'].dropna().values.sum() /
        (cmn_events.shape[0] * df_sum_per_stn.index.shape[0]))
    pass
#
    df_overall_bias.loc[temp_freq, 'B1_dwd'] = B1_dwd
    df_overall_bias.loc[temp_freq, 'B1_dwd_netatmo'] = B1_netatmo
    df_overall_bias.loc[temp_freq, 'B1_dwd_netatmo_unc'] = B1_netatmo_unc

    df_overall_bias.loc[temp_freq, 'B2_dwd'] = B2_dwd
    df_overall_bias.loc[temp_freq, 'B2_dwd_netatmo'] = B2_netatmo
    df_overall_bias.loc[temp_freq, 'B2_dwd_netatmo_unc'] = B2_netatmo_unc

    df_overall_bias.loc[temp_freq, 'B3_dwd'] = B3_dwd
    df_overall_bias.loc[temp_freq, 'B3_dwd_netatmo'] = B3_netatmo
    df_overall_bias.loc[temp_freq, 'B3_dwd_netatmo_unc'] = B3_netatmo_unc

    df_overall_bias.loc[temp_freq, 'B4_dwd'] = B4_dwd
    df_overall_bias.loc[temp_freq, 'B4_dwd_netatmo'] = B4_netatmo
    df_overall_bias.loc[temp_freq, 'B4_dwd_netatmo_unc'] = B4_netatmo_unc

    df_overall_bias.loc[temp_freq, 'B6_dwd'] = B6_dwd
    df_overall_bias.loc[temp_freq, 'B6_dwd_netatmo'] = B6_netatmo
    df_overall_bias.loc[temp_freq, 'B6_dwd_netatmo_unc'] = B6_netatmo_unc

    df_overall_bias.loc[temp_freq, 'B1_dwd_abs'] = B1_dwd_abs
    df_overall_bias.loc[temp_freq, 'B1_dwd_netatmo_abs'] = B1_netatmo_abs
    df_overall_bias.loc[temp_freq,
                        'B1_dwd_netatmo_unc_abs'] = B1_netatmo_unc_abs


#=====================================================================

print('Done')
df_overall_bias.to_csv(
    Path(out_dir_bias_dfs / ('abs%s.csv' % _interp_acc_)), sep=';')
