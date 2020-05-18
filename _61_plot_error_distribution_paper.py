
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

path_to_use = Path(
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\180520_error_dist_paper')

plot_filtered = True


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

def plot_histogram(x):
    hist, bins = np.histogram(x, bins=10, density=True)
    width = 0.5 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    return hist, width, center
# In[26]:


if plot_filtered:
    # , '180min', '360min', '720min', '1440min']:
    for temp_freq in ['60min', '180min', '360min', '720min', '1440min']:
        print(temp_freq)

        path_dwd_ppt_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\ppt_all_dwd_%s_.csv" % temp_freq
        df_dwd_edf = pd.read_csv(path_dwd_ppt_data, sep=';', index_col=0,
                                 parse_dates=True, infer_datetime_format=True)

        #########################################################

        dwd_int_path = (
            r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
            r"\oridinary_kriging_compare_DWD_Netatmo\Final_results"
            r"\Ppt_ok_ok_un_new2_first_flt__temp_flt__1st_%s"
            r"\interpolated_ppt_dwd_%s_data_Ppt_ok_ok_un_new2_first_flt__temp_flt__using_dwd_only_grp_11_1st.csv"
            % (temp_freq, temp_freq))

        C1_path = (
            r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
            r"\oridinary_kriging_compare_DWD_Netatmo"
            r"\Final_results"
            r"\Ppt_ok_ok_un_new_netatmo_no_flt___%s"
            r"\interpolated_ppt_dwd_%s_data_Ppt_ok_ok_un_2_netatmo_no_flt__using_dwd_netamo_grp_11_.csv"
            % (temp_freq, temp_freq))

        C2_path = (
            r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
            r"\oridinary_kriging_compare_DWD_Netatmo"
            r"\Final_results"
            r"\Ppt_ok_ok_un_new2_first_flt__1st_%s"
            r"\interpolated_ppt_dwd_%s_data_Ppt_ok_ok_un_new2_first_flt__using_dwd_netamo_grp_11_1st.csv"
            % (temp_freq, temp_freq))

        C3_path = (
            r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
            r"\oridinary_kriging_compare_DWD_Netatmo"
            r"\Final_results"
            r"\Ppt_ok_ok_un_new2_first_flt__temp_flt__1st_%s"
            r"\interpolated_ppt_dwd_%s_data_Ppt_ok_ok_un_new2_first_flt__temp_flt__using_dwd_netamo_grp_11_1st.csv"
            % (temp_freq, temp_freq))

        C4_path = (
            r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
            r"\oridinary_kriging_compare_DWD_Netatmo"
            r"\Final_results"
            r"\Ppt_ok_ok_un_new2_first_flt__temp_flt__1st_%s"
            r"\interpolated_ppt_dwd_%s_data_Ppt_ok_ok_un_new2_first_flt__temp_flt__using_dwd_netamo_grp_11_1st_unc10perc.csv"
            % (temp_freq, temp_freq))

        df_dwd = pd.read_csv(dwd_int_path,
                             sep=';', index_col=0, parse_dates=True,
                             infer_datetime_format=True)

        c1_df = pd.read_csv(
            C1_path,
            sep=';', index_col=0,
            parse_dates=True,
            infer_datetime_format=True)

        c2_df = pd.read_csv(
            C2_path,
            sep=';', index_col=0,
            parse_dates=True,
            infer_datetime_format=True)

        c3_df = pd.read_csv(
            C3_path,
            sep=';', index_col=0,
            parse_dates=True,
            infer_datetime_format=True)

        c4_df = pd.read_csv(
            C4_path,
            sep=';', index_col=0,
            parse_dates=True,
            infer_datetime_format=True)

        #########################################################

        df_compare = pd.DataFrame(
            index=df_dwd.index,
            columns=['mean_error_dwd',
                     'min_error_dwd',
                     'max_error_dwd',
                     'mean_error_c1',
                     'min_error_c1',
                     'max_error_c1',
                     'mean_error_c2',
                     'min_error_c2',
                     'max_error_c2',
                     'mean_error_c3',
                     'min_error_c3',
                     'max_error_c3',
                     'mean_error_c4',
                     'min_error_c4',
                     'max_error_c4'])

    #     df_dwd_edf = df_dwd_edf.loc[df_dwd_edf.index.intersection(
    #         df_netatmo_dwd.index), :]
    #
    #     df_dwd_edf = df_dwd_edf[df_dwd_edf > 0]
        # df_dwd_edf.dropna(how='all')
        # stns to remove from orig edf, because no data for 2015-2019
        df_dwd_edf.drop(columns=['00384', '13672', '05155'], inplace=True)

        try:

            for event_date in df_dwd.index.intersection(df_dwd_edf.index):
                (orig_ppt, dwd_int_ppt, c1_int_ppt,
                    c2_int_ppt, c3_int_ppt, c4_int_ppt
                 ) = [], [], [], [], [], []

                print(event_date)

                original_ppt = df_dwd_edf.loc[event_date, :]
                dwd_interp = df_dwd.loc[event_date, :]
                c1_int = c1_df.loc[event_date, :]
                c2_int = c2_df.loc[event_date, :]
                c3_int = c3_df.loc[event_date, :]
                c4_int = c4_df.loc[event_date, :]

                for stn_ in original_ppt.index:
                    # print(stn_)
                    ppt_stn_orig = original_ppt.loc[stn_]
                    dwd_int_stn = dwd_interp.loc[stn_]
                    c1_int_stn = c1_int.loc[stn_]
                    c2_int_stn = c2_int.loc[stn_]
                    c3_int_stn = c3_int.loc[stn_]
                    c4_int_stn = c4_int.loc[stn_]

                    if ((ppt_stn_orig >= 0) and
                        (dwd_int_stn >= 0) and
                        (c1_int_stn >= 0) and
                            (c2_int_stn >= 0) and
                            (c3_int_stn >= 0) and
                            (c4_int_stn >= 0)):

                        orig_ppt.append(ppt_stn_orig)
                        dwd_int_ppt.append(dwd_int_stn)
                        c1_int_ppt.append(c1_int_stn)
                        c2_int_ppt.append(
                            c2_int_stn)

                        c3_int_ppt.append(
                            c3_int_stn)
                        c4_int_ppt.append(
                            c4_int_stn)

                # error all dwd stns
                mean_error_dwd = np.round(np.mean(
                    np.array(orig_ppt) -
                    np.array(dwd_int_ppt)), 2)
                max_error_dwd = np.round(np.max(
                    np.array(orig_ppt) -
                    np.array(dwd_int_ppt)), 2)
                min_error_dwd = np.round(np.min(
                    np.array(orig_ppt) -
                    np.array(dwd_int_ppt)), 2)

                mean_error_c1 = np.round(np.mean(
                    np.array(orig_ppt) -
                    np.array(c1_int_ppt)), 2)
                max_error_c1 = np.round(np.max(
                    np.array(orig_ppt) -
                    np.array(c1_int_ppt)), 2)
                min_error_c1 = np.round(np.min(
                    np.array(orig_ppt) -
                    np.array(c1_int_ppt)), 2)

                mean_error_c2 = np.round(np.mean(
                    np.array(orig_ppt) -
                    np.array(c2_int_ppt)), 2)
                max_error_c2 = np.round(np.max(
                    np.array(orig_ppt) -
                    np.array(c2_int_ppt)), 2)
                min_error_c2 = np.round(np.min(
                    np.array(orig_ppt) -
                    np.array(c2_int_ppt)), 2)

                mean_error_c3 = np.round(np.mean(
                    np.array(orig_ppt) -
                    np.array(c3_int_ppt)), 2)
                max_error_c3 = np.round(np.max(
                    np.array(orig_ppt) -
                    np.array(c3_int_ppt)), 2)
                min_error_c3 = np.round(np.min(
                    np.array(orig_ppt) -
                    np.array(c3_int_ppt)), 2)

                mean_error_c4 = np.round(np.mean(
                    np.array(orig_ppt) -
                    np.array(c4_int_ppt)), 2)
                max_error_c4 = np.round(np.max(
                    np.array(orig_ppt) -
                    np.array(c4_int_ppt)), 2)
                min_error_c4 = np.round(np.min(
                    np.array(orig_ppt) -
                    np.array(c4_int_ppt)), 2)

                df_compare.loc[event_date, 'mean_error_dwd'] = mean_error_dwd
                df_compare.loc[event_date, 'max_error_dwd'] = max_error_dwd
                df_compare.loc[event_date, 'min_error_dwd'] = min_error_dwd

                df_compare.loc[event_date,
                               'mean_error_c1'] = mean_error_c1

                df_compare.loc[event_date,
                               'max_error_c1'] = max_error_c1
                df_compare.loc[event_date,
                               'min_error_c1'] = min_error_c1

                df_compare.loc[event_date,
                               'mean_error_c2'] = mean_error_c2
                df_compare.loc[event_date,
                               'max_error_c2'] = max_error_c2

                df_compare.loc[event_date,
                               'min_error_c2'] = min_error_c2

                df_compare.loc[event_date,
                               'mean_error_c3'] = mean_error_c3
                df_compare.loc[event_date,
                               'max_error_c3'] = max_error_c3

                df_compare.loc[event_date,
                               'min_error_c3'] = min_error_c3

                df_compare.loc[event_date,
                               'mean_error_c4'] = mean_error_c4
                df_compare.loc[event_date,
                               'max_error_c4'] = max_error_c4

                df_compare.loc[event_date,
                               'min_error_c4'] = min_error_c4

                # df_compare = df_compare[df_compare > 0]
                #df_compare.dropna(how='all', inplace=True)

        except Exception as msg:
            print(msg)

        df_compare.sort_index(inplace=True)
        df_compare.dropna(how='all', inplace=True)

        df_compare[df_compare > 30] = 30
        df_compare[df_compare < -30] = -30
        #======================================================================
        #
        #======================================================================
        plt.ioff()
        fig = plt.figure(figsize=(24, 12), dpi=150)

        ax4 = fig.add_subplot(111)

        hist1, width1, center1 = plot_histogram(
            df_compare.mean_error_c1.values)

        hist2, width2, center2 = plot_histogram(
            df_compare.mean_error_c2.values)
        hist3, width3, center3 = plot_histogram(
            df_compare.mean_error_c3.dropna().values)
        hist4, width4, center4 = plot_histogram(
            df_compare.mean_error_c4.dropna().values)

        ax4.bar(center1, hist1, align='center', width=width1,
                alpha=0.5, color='b', edgecolor='darkblue',
                linewidth=1,
                label='C1')

        ax4.bar(center2, hist2, align='center', width=width2,
                alpha=0.5, color='r', edgecolor='darkred',
                linewidth=1,
                label='C2')

        ax4.bar(center3, hist3, align='center', width=width3,
                alpha=0.5, color='g', edgecolor='darkgreen',
                linewidth=1,
                label='C3')

        ax4.bar(center4, hist4, align='center', width=width4,
                alpha=0.5, color='orange', edgecolor='darkorange',
                linewidth=1,
                label='C4')

        ax4.set_xlabel('Error [mm/%s]' % (temp_freq))
        ax4.set_ylabel('Normed Frequency')
        ax4.grid(alpha=0.5)
        ax4.legend(loc='upper right')
        ax4.title.set_text('Distribution of mean spatial error'
                           ' (each event all stations)')

        plt.savefig(
            os.path.join(
                path_to_use,
                r'hist_error_spatial_dist_%s_2.png' % temp_freq,
            ),
            frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

        #======================================================================
        #
        #======================================================================

        df_compare_stns = pd.DataFrame(
            index=df_dwd.columns,
            columns=['mean_error_dwd',
                     'min_error_dwd',
                     'max_error_dwd',
                     'mean_error_c1',
                     'min_error_c1',
                     'max_error_c1',
                     'mean_error_c2',
                     'min_error_c2',
                     'max_error_c2',
                     'mean_error_c3',
                     'min_error_c3',
                     'max_error_c3',
                     'mean_error_c4',
                     'min_error_c4',
                     'max_error_c4'])

        # stns to remove from orig edf, because no data for 2015-2019

        try:

            for dwd_stn in df_dwd.columns.intersection(df_dwd_edf.columns):
                (orig_ppt, dwd_int_ppt, c1_int_ppt,
                    c2_int_ppt, c3_int_ppt, c4_int_ppt
                 ) = [], [], [], [], [], []

                print(dwd_stn)

                original_ppt = df_dwd_edf.loc[:, dwd_stn]
                dwd_interp = df_dwd.loc[:, dwd_stn]
                c1_int = c1_df.loc[:, dwd_stn]
                c2_int = c2_df.loc[:, dwd_stn]
                c3_int = c3_df.loc[:, dwd_stn]
                c4_int = c4_df.loc[:, dwd_stn]

                for event_ in original_ppt.index.intersection(c1_int.index):
                    # print(stn_)
                    ppt_stn_orig = original_ppt.loc[event_]
                    dwd_int_stn = dwd_interp.loc[event_]
                    c1_int_stn = c1_int.loc[event_]
                    c2_int_stn = c2_int.loc[event_]
                    c3_int_stn = c3_int.loc[event_]
                    c4_int_stn = c4_int.loc[event_]

                    if ((ppt_stn_orig >= 0) and
                        (dwd_int_stn >= 0) and
                        (c1_int_stn >= 0) and
                            (c2_int_stn >= 0) and
                            (c3_int_stn >= 0) and
                            (c4_int_stn >= 0)):

                        orig_ppt.append(ppt_stn_orig)
                        dwd_int_ppt.append(dwd_int_stn)
                        c1_int_ppt.append(c1_int_stn)
                        c2_int_ppt.append(
                            c2_int_stn)

                        c3_int_ppt.append(
                            c3_int_stn)
                        c4_int_ppt.append(
                            c4_int_stn)

                # error all dwd stns
                mean_error_dwd = np.round(np.mean(
                    np.array(orig_ppt) -
                    np.array(dwd_int_ppt)), 2)
                max_error_dwd = np.round(np.max(
                    np.array(orig_ppt) -
                    np.array(dwd_int_ppt)), 2)
                min_error_dwd = np.round(np.min(
                    np.array(orig_ppt) -
                    np.array(dwd_int_ppt)), 2)

                mean_error_c1 = np.round(np.mean(
                    np.array(orig_ppt) -
                    np.array(c1_int_ppt)), 2)
                max_error_c1 = np.round(np.max(
                    np.array(orig_ppt) -
                    np.array(c1_int_ppt)), 2)
                min_error_c1 = np.round(np.min(
                    np.array(orig_ppt) -
                    np.array(c1_int_ppt)), 2)

                mean_error_c2 = np.round(np.mean(
                    np.array(orig_ppt) -
                    np.array(c2_int_ppt)), 2)
                max_error_c2 = np.round(np.max(
                    np.array(orig_ppt) -
                    np.array(c2_int_ppt)), 2)
                min_error_c2 = np.round(np.min(
                    np.array(orig_ppt) -
                    np.array(c2_int_ppt)), 2)

                mean_error_c3 = np.round(np.mean(
                    np.array(orig_ppt) -
                    np.array(c3_int_ppt)), 2)
                max_error_c3 = np.round(np.max(
                    np.array(orig_ppt) -
                    np.array(c3_int_ppt)), 2)
                min_error_c3 = np.round(np.min(
                    np.array(orig_ppt) -
                    np.array(c3_int_ppt)), 2)

                mean_error_c4 = np.round(np.mean(
                    np.array(orig_ppt) -
                    np.array(c4_int_ppt)), 2)
                max_error_c4 = np.round(np.max(
                    np.array(orig_ppt) -
                    np.array(c4_int_ppt)), 2)
                min_error_c4 = np.round(np.min(
                    np.array(orig_ppt) -
                    np.array(c4_int_ppt)), 2)

                df_compare_stns.loc[dwd_stn, 'mean_error_dwd'] = mean_error_dwd
                df_compare_stns.loc[dwd_stn, 'max_error_dwd'] = max_error_dwd
                df_compare_stns.loc[dwd_stn, 'min_error_dwd'] = min_error_dwd

                df_compare_stns.loc[dwd_stn,
                                    'mean_error_c1'] = mean_error_c1

                df_compare_stns.loc[dwd_stn,
                                    'max_error_c1'] = max_error_c1
                df_compare_stns.loc[dwd_stn,
                                    'min_error_c1'] = min_error_c1

                df_compare_stns.loc[dwd_stn,
                                    'mean_error_c2'] = mean_error_c2
                df_compare_stns.loc[dwd_stn,
                                    'max_error_c2'] = max_error_c2

                df_compare_stns.loc[dwd_stn,
                                    'min_error_c2'] = min_error_c2

                df_compare_stns.loc[dwd_stn,
                                    'mean_error_c3'] = mean_error_c3
                df_compare_stns.loc[event_date,
                                    'max_error_c3'] = max_error_c3

                df_compare_stns.loc[dwd_stn,
                                    'min_error_c3'] = min_error_c3

                df_compare_stns.loc[dwd_stn,
                                    'mean_error_c4'] = mean_error_c4
                df_compare_stns.loc[event_date,
                                    'max_error_c4'] = max_error_c4

                df_compare_stns.loc[dwd_stn,
                                    'min_error_c4'] = min_error_c4

                # df_compare = df_compare[df_compare > 0]
                #df_compare.dropna(how='all', inplace=True)

        except Exception as msg:
            print(msg)

        # df_compare_stns.sort_index(inplace=True)
        df_compare_stns.dropna(how='all', inplace=True)
        df_compare_stns[df_compare_stns > 30] = 30
        df_compare_stns[df_compare_stns < -30] = -30

        #======================================================================
        #
        #======================================================================
        plt.ioff()
        fig = plt.figure(figsize=(24, 12), dpi=150)

        ax4 = fig.add_subplot(111)

        hist1, width1, center1 = plot_histogram(
            df_compare_stns.mean_error_c1.dropna().values)

        hist2, width2, center2 = plot_histogram(
            df_compare_stns.mean_error_c2.dropna().values)
        hist3, width3, center3 = plot_histogram(
            df_compare_stns.mean_error_c3.dropna().values)
        hist4, width4, center4 = plot_histogram(
            df_compare_stns.mean_error_c4.dropna().values)

        ax4.bar(center1, hist1, align='center', width=width1,
                alpha=0.5, color='b', edgecolor='darkblue',
                linewidth=1,
                label='C1')

        ax4.bar(center2, hist2, align='center', width=width2,
                alpha=0.5, color='r', edgecolor='darkred',
                linewidth=1,
                label='C2')

        ax4.bar(center3, hist3, align='center', width=width3,
                alpha=0.5, color='g', edgecolor='darkgreen',
                linewidth=1,
                label='C3')

        ax4.bar(center4, hist4, align='center', width=width4,
                alpha=0.5, color='orange', edgecolor='darkorange',
                linewidth=1,
                label='C4')

        ax4.set_xlabel('Error [mm/%s]' % (temp_freq))
        ax4.set_ylabel('Normed Frequency')
        ax4.grid(alpha=0.5)
        ax4.legend(loc='upper right')
        ax4.title.set_text('Distribution of mean temporal error'
                           ' (each station all events)')

        plt.savefig(
            os.path.join(
                path_to_use,
                r'hist_error_temp_dist_%s_2.png' % temp_freq,
            ),
            frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

        #======================================================================
#
#==============================================================================
#         # OK
#         plt.ioff()
#         fig = plt.figure(figsize=(24, 12), dpi=150)
#
#         ax = fig.add_subplot(111)
#
#         ax.plot(df_compare.index,
#                 df_compare.mean_error_dwd,
#                 alpha=.8,
#                 c='b',  # colors_arr,
#                 marker='d',
#                 label='DWD Mean Err %0.2f'
#                 % np.mean(df_compare.mean_error_dwd))
#
#         ax.plot(df_compare.index,
#                 df_compare.max_error_dwd,
#                 alpha=.8,
#                 c='darkblue',  # colors_arr,
#                 marker='d',
#                 label='DWD Max Err %0.2f'
#                 % np.mean(df_compare.max_error_dwd))
#         ax.plot(df_compare.index,
#                 df_compare.min_error_dwd,
#                 alpha=.8,
#                 c='lightblue',  # colors_arr,
#                 marker='d',
#                 label='DWD Min Err %0.2f'
#                 % np.mean(df_compare.min_error_dwd))
#         #======================================================================
#         #
#         #======================================================================
#         ax.plot(df_compare.index,
#                 df_compare.mean_error_dwd_netatmo,
#                 alpha=.5,
#                 c='r',  # colors_arr,
#                 marker='*',
#                 label='DWD-Netatmo Mean Err %0.2f'
#                 % np.mean(df_compare.mean_error_dwd_netatmo))
#
#         ax.plot(df_compare.index,
#                 df_compare.max_error_dwd_netatmo,
#                 alpha=.5,
#                 c='darkred',  # colors_arr,
#                 marker='*',
#                 label='DWD-Netatmo Max Err %0.2f'
#                 % np.mean(df_compare.max_error_dwd_netatmo))
#
#         ax.plot(df_compare.index,
#                 df_compare.min_error_dwd_netatmo,
#                 alpha=.5,
#                 c='salmon',  # colors_arr,
#                 marker='*',
#                 label='DWD-Netatmo Min Err %0.2f'
#                 % np.mean(df_compare.min_error_dwd_netatmo))
#
#         #======================================================================
#         #
#         #======================================================================
#         ax.plot(df_compare.index,
#                 df_compare.mean_error_dwd_netatmo_unc,
#                 alpha=.5,
#                 c='g',  # colors_arr,
#                 marker='x',
#                 label='DWD-Netatmo Unc Mean Err %0.2f'
#                 % np.mean(df_compare.mean_error_dwd_netatmo_unc))
#
#         ax.plot(df_compare.index,
#                 df_compare.max_error_dwd_netatmo_unc,
#                 alpha=.5,
#                 c='darkgreen',  # colors_arr,
#                 marker='x',
#                 label='DWD-Netatmo Unc Max Err %0.2f'
#                 % np.mean(df_compare.max_error_dwd_netatmo_unc))
#
#         ax.plot(df_compare.index,
#                 df_compare.min_error_dwd_netatmo_unc,
#                 alpha=.5,
#                 c='lightgreen',  # colors_arr,
#                 marker='x',
#                 label='DWD-Netatmo Unc Min Err %0.2f'
#                 % np.mean(df_compare.min_error_dwd_netatmo_unc))


#         ax.set_title('Pearson Correlation Interpolated Quantiles from DWD or DWD-Netatmo\n '
#                      'Precipitation of %s Extreme Events %s\n Events with Improvemnts %d / %d, Percentage %0.0f\n'
#                      'Events with Improvemnts with OK 2percUnc %d / %d, Percentage %0.0f'
#                      '\nEvents with Improvemnts with OK 5percUnc %d / %d, Percentage %0.0f'
#                      '\nEventswith Improvemnts with OK 10percUnc %d / %d, Percentage %0.0f'
#                      '\nEvents with Improvemnts with OK 20percUnc %d / %d, Percentage %0.0f'
#                      % (temp_freq, _interp_acc_,
#                         stations_with_improvements,
#                          df_compare.pearson_corr_dwd_netatmo.shape[0],
#                         percent_of_improvment,
#                         stations_with_improvements_unc2perc,
#                         df_compare.pearson_corr_dwd_netatmo_unc2perc.shape[0],
#                         percent_of_improvment_unc2perc,
#                         stations_with_improvements_unc5perc,
#                         df_compare.pearson_corr_dwd_netatmo_unc2perc.shape[0],
#                         percent_of_improvment_unc5perc,
#                         stations_with_improvements_unc10perc,
#                         df_compare.pearson_corr_dwd_netatmo_unc2perc.shape[0],
#                         percent_of_improvment_unc10perc,
#                         stations_with_improvements_unc20perc,
#                         df_compare.pearson_corr_dwd_netatmo_unc2perc.shape[0],
#                         percent_of_improvment_unc20perc))
#         ax.grid(alpha=0.25)
#         plt.setp(ax.get_xticklabels(), rotation=45)
#
#         # years = mdates.YearLocator(1)   # every year
#         months = mdates.MonthLocator()  # every month
#
#         #hrs = mdates.HourLocator()
#         # years_fmt = mdates.DateFormatter('%Y-%m-%d')
#
#         years_fmt = mdates.DateFormatter('%Y-%m-%d %H')
#         ax.xaxis.set_major_locator(months)
#
#         ax.xaxis.set_major_formatter(years_fmt)
#         # ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))
#
#         #ax.set_yticks(np.arange(0., 1.05, .10))
#         ax.set_xlabel('Date of Event')
#         ax.legend(loc='lower right')
#         ax.set_ylabel('Error compared to observed DWD values')
#
#         plt.savefig(os.path.join(path_to_use,
#                                  r'error_spatial_distribution_%s_events_dwd2.png' % temp_freq,
#                                  ),
#                     frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
#         plt.close()
