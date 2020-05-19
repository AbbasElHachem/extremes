
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np
import warnings

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import fnmatch
from pathlib import Path

import scipy.stats as stats

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears

from pandas.plotting import register_matplotlib_converters

from _00_additional_functions import build_edf_fr_vals
from _00_additional_functions import FitKernelDensity, KernelDensityEstimate

with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)

FitKernel = FitKernelDensity()

fit_kernel_fct = KernelDensityEstimate()


register_matplotlib_converters()

plt.ioff()
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'axes.labelsize': 16})


path_to_use = Path(
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\180520_error_dist_paper')

plot_filtered = True

temp_freq = '60min'  # '1440min'  # '60min'

events_60min = ['2018-06-11 16:00:00', '2018-09-06 18:00:00']
events_1440min = ['2018-05-14 00:00:00', '2019-07-28 00:00:00']

if temp_freq == '60min':
    event_list = events_60min

if temp_freq == '1440min':
    event_list = events_1440min
#==============================================================================
#
#==============================================================================


def plot_histogram(x):
    hist, bins = np.histogram(x, bins=10, density=False)
    width = 0.5 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    return hist, width, center
# In[26]:


if plot_filtered:
    # , '180min', '360min', '720min', '1440min']:
    # , '180min', '360min', '720min', '1440min']:
    for temp_freq in ['60min']:
        print(temp_freq)

        path_dwd_ppt_data = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
                             r"\NetAtmo_BW\ppt_all_dwd_%s_.csv" % temp_freq)
        df_dwd_edf = pd.read_csv(path_dwd_ppt_data, sep=';', index_col=0,
                                 parse_dates=True, infer_datetime_format=True)

        #########################################################

        dwd_int_path = (
            "X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
            r"\oridinary_kriging_compare_DWD_Netatmo\Final_results"
            r"\Ppt_ok_ok_un_new_netatmo_no_flt___%s"
            r"\interpolated_ppt_dwd_%s_data_Ppt_ok_ok_un_2_netatmo_no_flt__using_dwd_only_grp_11_.csv"
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
                     'mean_error_c1',
                     'mean_error_c2',
                     'mean_error_c3',
                     'mean_error_c4',
                     'pears_dwd',
                     'pears_c3',
                     'spr_dwd',
                     'spr_c3'])

    #     df_dwd_edf = df_dwd_edf.loc[df_dwd_edf.index.intersection(
    #         df_netatmo_dwd.index), :]
    #
    #     df_dwd_edf = df_dwd_edf[df_dwd_edf > 0]
        # df_dwd_edf.dropna(how='all')
        # stns to remove from orig edf, because no data for 2015-2019
        df_dwd_edf.drop(columns=['00384', '13672', '05155'], inplace=True)

        df_dwd_edf = df_dwd_edf.loc[df_dwd.index, :]
        try:

            # df_dwd.index.intersection(df_dwd_edf.index):
            for event_date in df_dwd.index.intersection(df_dwd_edf.index):
                # event_list:
                (orig_ppt, dwd_int_ppt, c1_int_ppt,
                    c2_int_ppt, c3_int_ppt, c4_int_ppt
                 ) = [], [], [], [], [], []

                print(event_date)

                original_ppt = df_dwd_edf.loc[event_date, :].dropna(how='all')
                dwd_interp = df_dwd.loc[event_date,
                                        original_ppt.index]
                c1_int = c1_df.loc[event_date,
                                   original_ppt.index]
                c2_int = c2_df.loc[event_date,
                                   original_ppt.index]
                c3_int = c3_df.loc[event_date,
                                   original_ppt.index]
                c4_int = c4_df.loc[event_date,
                                   original_ppt.index]

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
                        c2_int_ppt.append(c2_int_stn)

                        c3_int_ppt.append(c3_int_stn)
                        c4_int_ppt.append(c4_int_stn)

                # corr
                corr_dwd = np.round(pears(orig_ppt, dwd_int_ppt)[0], 2)
                corr_dwd_net = np.round(pears(orig_ppt, c3_int_ppt)[0], 2)

                spr_dwd = np.round(spr(orig_ppt, dwd_int_ppt)[0], 2)
                spr_dwd_net = np.round(spr(orig_ppt, c3_int_ppt)[0], 2)

                data_c0 = np.array(orig_ppt) - np.array(dwd_int_ppt)
                data_c1 = np.array(orig_ppt) - np.array(c1_int_ppt)
                data_c2 = np.array(orig_ppt) - np.array(c2_int_ppt)
                data_c3 = np.array(orig_ppt) - np.array(c3_int_ppt)
                data_c4 = np.array(orig_ppt) - np.array(c4_int_ppt)

                d_opt0, no0, pp0 = fit_kernel_fct.fit_kernel_to_data(data_c0)
                d_opt1, no1, pp1 = fit_kernel_fct.fit_kernel_to_data(data_c1)
                d_opt2, no2, pp2 = fit_kernel_fct.fit_kernel_to_data(data_c2)
                d_opt3, no3, pp3 = fit_kernel_fct.fit_kernel_to_data(data_c3)
                d_opt4, no4, pp4 = fit_kernel_fct.fit_kernel_to_data(data_c4)

                plt.ioff()
                fig = plt.figure(figsize=(24, 12), dpi=100)

                ax3 = fig.add_subplot(111)
                ax3.plot(pp0, no0, c='k', label='C0', marker='_', alpha=0.85)
                ax3.plot(pp1, no1, c='r', label='C1', marker='+', alpha=0.5)
                ax3.plot(pp2, no2, c='b', label='C2', marker='1', alpha=0.5)
                ax3.plot(pp3, no3, c='g', label='C3', marker='o', alpha=0.5)
                ax3.plot(pp4, no4, c='orange', label='C4',
                         marker='.', alpha=0.5)

                ax3.set_xlabel('Error [mm/%s]' % (temp_freq))
                ax3.set_ylabel('Density')
                ax3.grid(alpha=0.5)
                ax3.legend(loc='upper right')
                ax3.title.set_text('Error distribution event %s\n'
                                   'Pears Corr. C0=%0.2f--C3=%.2f'
                                   ' Spear Corr.C0=%0.2f--C3=%.2f'
                                   % (event_date, corr_dwd, corr_dwd_net,
                                      spr_dwd, spr_dwd_net))
                plt.savefig(
                    os.path.join(
                        path_to_use,
                        r'error_dist_%s_%s1.png'
                        % (temp_freq, str(event_date).replace(':', '_')),
                    ),
                    frameon=True, papertype='a4',
                    bbox_inches='tight', pad_inches=.2)
                plt.close()
#                 a = np.array(orig_ppt) - np.array(dwd_int_ppt)
#                 b = np.array(orig_ppt) - np.array(c3_int_ppt)

#                 plt.ioff()
#                 plt.plot(a, c='r')
#                 plt.plot(b, c='b')
#                 plt.show()

                # error all dwd stns
                mean_error_dwd = np.round(np.mean(
                    np.array(orig_ppt) - np.array(dwd_int_ppt)), 2)

                mean_error_c1 = np.round(np.mean(
                    np.array(orig_ppt) - np.array(c1_int_ppt)), 2)

                mean_error_c2 = np.round(np.mean(
                    np.array(orig_ppt) - np.array(c2_int_ppt)), 2)

                mean_error_c3 = np.round(np.mean(
                    np.array(orig_ppt) - np.array(c3_int_ppt)), 2)

                mean_error_c4 = np.round(np.mean(
                    np.array(orig_ppt) - np.array(c4_int_ppt)), 2)

                df_compare.loc[event_date, 'mean_error_dwd'] = mean_error_dwd

                df_compare.loc[event_date,
                               'mean_error_c1'] = mean_error_c1

                df_compare.loc[event_date,
                               'mean_error_c2'] = mean_error_c2

                df_compare.loc[event_date,
                               'mean_error_c3'] = mean_error_c3

                df_compare.loc[event_date,
                               'mean_error_c4'] = mean_error_c4

                df_compare.loc[event_date,
                               'pears_dwd'] = corr_dwd
                df_compare.loc[event_date,
                               'pears_c3'] = corr_dwd_net

                df_compare.loc[event_date,
                               'spr_dwd'] = corr_dwd
                df_compare.loc[event_date,
                               'spr_c3'] = corr_dwd_net

                # df_compare = df_compare[df_compare > 0]
                #df_compare.dropna(how='all', inplace=True)

        except Exception as msg:
            print(msg)

        df_compare.sort_index(inplace=True)
        df_compare.dropna(how='all', inplace=True)


#         df_compare[df_compare > 10] = 10
#         df_compare[df_compare < -10] = -10

#         np.where(df_compare.pears_dwd > df_compare.pears_dwd_netatmo)
#         plt.ioff()
#         df_compare.pears_dwd.plot(legend=False, c='r')
#         df_compare.pears_dwd_netatmo.plot(legend=False, c='b')
#         plt.show()
        pass
        #======================================================================
        #
#         #======================================================================
#         plt.ioff()
#         fig = plt.figure(figsize=(24, 12), dpi=150)
#
#         ax4 = fig.add_subplot(111)
#
#         hist1, width1, center1 = plot_histogram(
#             df_compare.mean_error_c1.values)
#
#         hist2, width2, center2 = plot_histogram(
#             df_compare.mean_error_c2.values)
#         hist3, width3, center3 = plot_histogram(
#             df_compare.mean_error_c3.dropna().values)
#         hist4, width4, center4 = plot_histogram(
#             df_compare.mean_error_c4.dropna().values)
#
#         ax4.bar(center1, hist1, align='center', width=width1,
#                 alpha=0.5, color='b', edgecolor='darkblue',
#                 linewidth=1,
#                 label='C1')
#
#         ax4.bar(center2, hist2, align='center', width=width2,
#                 alpha=0.5, color='r', edgecolor='darkred',
#                 linewidth=1,
#                 label='C2')
#
#         ax4.bar(center3, hist3, align='center', width=width3,
#                 alpha=0.5, color='g', edgecolor='darkgreen',
#                 linewidth=1,
#                 label='C3')
#
#         ax4.bar(center4, hist4, align='center', width=width4,
#                 alpha=0.5, color='orange', edgecolor='darkorange',
#                 linewidth=1,
#                 label='C4')
#
#         ax4.set_xlabel('Error [mm/%s]' % (temp_freq))
#         ax4.set_ylabel('Normed Frequency')
#         ax4.grid(alpha=0.5)
#         ax4.legend(loc='upper right')
#         ax4.title.set_text('Distribution of mean spatial error'
#                            ' (each event all stations)')
#
#         plt.savefig(
#             os.path.join(
#                 path_to_use,
#                 r'hist_error_spatial_dist_%s_2.png' % temp_freq,
#             ),
#             frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
#         plt.close()
        #======================================================================
        #
        #======================================================================
#         xc0, yc0 = build_edf_fr_vals(df_compare.mean_error_dwd.dropna().values)
#         xc1, yc1 = build_edf_fr_vals(df_compare.mean_error_c1.dropna().values)
#         xc2, yc2 = build_edf_fr_vals(df_compare.mean_error_c2.dropna().values)
#         xc3, yc3 = build_edf_fr_vals(df_compare.mean_error_c3.dropna().values)
#         xc4, yc4 = build_edf_fr_vals(df_compare.mean_error_c4.dropna().values)
#
#         plt.ioff()
#         fig = plt.figure(figsize=(24, 12), dpi=150)
#
#         ax3 = fig.add_subplot(111)
#         ax3.plot(xc0, yc0, c='k', label='DWD', marker='_', alpha=0.5)
#         ax3.plot(xc1, yc1, c='r', label='C1', marker='+', alpha=0.5)
#         ax3.plot(xc2, yc2, c='b', label='C2', marker='1', alpha=0.5)
#         ax3.plot(xc3, yc3, c='g', label='C3', marker='o', alpha=0.5)
#         ax3.plot(xc4, yc4, c='orange', label='C4',
#                  marker='.', alpha=0.5)
#
#         ax3.set_xlabel('Error [mm/%s]' % (temp_freq))
#         ax3.set_ylabel('Frequency')
#         ax3.grid(alpha=0.5)
#         ax3.legend(loc='lower right')
#         ax3.title.set_text('Distribution of mean spatial error'
#                            ' (each event all stations)')
#
#         plt.savefig(
#             os.path.join(
#                 path_to_use,
#                 r'cdf_hist_error_spatial_dist_%s_2.png' % temp_freq,
#             ),
#             frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
#         plt.close()
        #======================================================================

        #======================================================================
        #
        #======================================================================
#         df_compare_stns = pd.DataFrame(
#             index=df_dwd.columns,
#             columns=['mean_error_dwd',
#                      'mean_error_c1',
#                      'mean_error_c2',
#                      'mean_error_c3',
#                      'mean_error_c4'])
#
#         # stns to remove from orig edf, because no data for 2015-2019
#
#         try:
#
#             for dwd_stn in df_dwd.columns.intersection(df_dwd_edf.columns):
#                 (orig_ppt, dwd_int_ppt, c1_int_ppt,
#                     c2_int_ppt, c3_int_ppt, c4_int_ppt
#                  ) = [], [], [], [], [], []
#
#                 print(dwd_stn)
#
#                 original_ppt = df_dwd_edf.loc[:, dwd_stn]
#                 dwd_interp = df_dwd.loc[:, dwd_stn]
#                 c1_int = c1_df.loc[:, dwd_stn]
#                 c2_int = c2_df.loc[:, dwd_stn]
#                 c3_int = c3_df.loc[:, dwd_stn]
#                 c4_int = c4_df.loc[:, dwd_stn]
#
#                 for event_ in original_ppt.index.intersection(c1_int.index):
#                     # print(stn_)
#                     ppt_stn_orig = original_ppt.loc[event_]
#                     dwd_int_stn = dwd_interp.loc[event_]
#                     c1_int_stn = c1_int.loc[event_]
#                     c2_int_stn = c2_int.loc[event_]
#                     c3_int_stn = c3_int.loc[event_]
#                     c4_int_stn = c4_int.loc[event_]
#
#                     if ((ppt_stn_orig >= 0) and
#                         (dwd_int_stn >= 0) and
#                         (c1_int_stn >= 0) and
#                             (c2_int_stn >= 0) and
#                             (c3_int_stn >= 0) and
#                             (c4_int_stn >= 0)):
#
#                         orig_ppt.append(ppt_stn_orig)
#                         dwd_int_ppt.append(dwd_int_stn)
#                         c1_int_ppt.append(c1_int_stn)
#                         c2_int_ppt.append(
#                             c2_int_stn)
#
#                         c3_int_ppt.append(
#                             c3_int_stn)
#                         c4_int_ppt.append(
#                             c4_int_stn)
#
#                 # error all dwd stns
#                 mean_error_dwd = np.round(np.mean(
#                     np.array(orig_ppt) -
#                     np.array(dwd_int_ppt)), 2)
#
#                 mean_error_c1 = np.round(np.mean(
#                     np.array(orig_ppt) -
#                     np.array(c1_int_ppt)), 2)
#
#                 mean_error_c2 = np.round(np.mean(
#                     np.array(orig_ppt) -
#                     np.array(c2_int_ppt)), 2)
#
#                 mean_error_c3 = np.round(np.mean(
#                     np.array(orig_ppt) -
#                     np.array(c3_int_ppt)), 2)
#
#                 mean_error_c4 = np.round(np.mean(
#                     np.array(orig_ppt) -
#                     np.array(c4_int_ppt)), 2)
#
#                 df_compare_stns.loc[dwd_stn, 'mean_error_dwd'] = mean_error_dwd
#
#                 df_compare_stns.loc[dwd_stn,
#                                     'mean_error_c1'] = mean_error_c1
#
#                 df_compare_stns.loc[dwd_stn,
#                                     'mean_error_c2'] = mean_error_c2
#
#                 df_compare_stns.loc[dwd_stn,
#                                     'mean_error_c3'] = mean_error_c3
#
#                 df_compare_stns.loc[dwd_stn, 'mean_error_c4'] = mean_error_c4
#
#                 # df_compare = df_compare[df_compare > 0]
#                 #df_compare.dropna(how='all', inplace=True)
#
#         except Exception as msg:
#             print(msg)
#
#         # df_compare_stns.sort_index(inplace=True)
#         df_compare_stns.dropna(how='all', inplace=True)
#         df_compare_stns[df_compare_stns > 10] = 10
#         df_compare_stns[df_compare_stns < -10] = -10
#
#         #======================================================================
#         #
#         #======================================================================
#         plt.ioff()
# #         fig = plt.figure(figsize=(24, 12), dpi=150)
# #
# #         ax4 = fig.add_subplot(111)
# #
# #         hist1, width1, center1 = plot_histogram(
# #             df_compare_stns.mean_error_c1.dropna().values)
# #
# #         hist2, width2, center2 = plot_histogram(
# #             df_compare_stns.mean_error_c2.dropna().values)
# #         hist3, width3, center3 = plot_histogram(
# #             df_compare_stns.mean_error_c3.dropna().values)
# #         hist4, width4, center4 = plot_histogram(
# #             df_compare_stns.mean_error_c4.dropna().values)
# #
# #         ax4.bar(center1, hist1, align='center', width=width1,
# #                 alpha=0.5, color='b', edgecolor='darkblue',
# #                 linewidth=1,
# #                 label='C1')
# #
# #         ax4.bar(center2, hist2, align='center', width=width2,
# #                 alpha=0.5, color='r', edgecolor='darkred',
# #                 linewidth=1,
# #                 label='C2')
# #
# #         ax4.bar(center3, hist3, align='center', width=width3,
# #                 alpha=0.5, color='g', edgecolor='darkgreen',
# #                 linewidth=1,
# #                 label='C3')
# #
# #         ax4.bar(center4, hist4, align='center', width=width4,
# #                 alpha=0.5, color='orange', edgecolor='darkorange',
# #                 linewidth=1,
# #                 label='C4')
# #
# #         ax4.set_xlabel('Error [mm/%s]' % (temp_freq))
# #         ax4.set_ylabel('Normed Frequency')
# #         ax4.grid(alpha=0.5)
# #         ax4.legend(loc='upper right')
# #         ax4.title.set_text('Distribution of mean temporal error'
# #                            ' (each station all events)')
# #
# #         plt.savefig(
# #             os.path.join(
# #                 path_to_use,
# #                 r'hist_error_temp_dist_%s_2.png' % temp_freq,
# #             ),
# #             frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
# #         plt.close()
#
#         xc0, yc0 = build_edf_fr_vals(
#             df_compare_stns.mean_error_dwd.dropna().values)
#         xc1, yc1 = build_edf_fr_vals(
#             df_compare_stns.mean_error_c1.dropna().values)
#         xc2, yc2 = build_edf_fr_vals(
#             df_compare_stns.mean_error_c2.dropna().values)
#         xc3, yc3 = build_edf_fr_vals(
#             df_compare_stns.mean_error_c3.dropna().values)
#         xc4, yc4 = build_edf_fr_vals(
#             df_compare_stns.mean_error_c4.dropna().values)
#
#         plt.ioff()
#         fig = plt.figure(figsize=(24, 12), dpi=150)
#
#         ax3 = fig.add_subplot(111)
#         ax3.plot(xc0, yc0, c='k', label='DWD', marker='_', alpha=0.5)
#         ax3.plot(xc1, yc1, c='r', label='C1', marker='+', alpha=0.5)
#         ax3.plot(xc2, yc2, c='b', label='C2', marker='1', alpha=0.5)
#         ax3.plot(xc3, yc3, c='g', label='C3', marker='o', alpha=0.5)
#         ax3.plot(xc4, yc4, c='orange', label='C4',
#                  marker='.', alpha=0.5)
#
#         ax3.set_xlabel('Error [mm/%s]' % (temp_freq))
#         ax3.set_ylabel('Frequency')
#         ax3.grid(alpha=0.5)
#         ax3.legend(loc='lower right')
#         ax3.title.set_text('Distribution of mean temporal error'
#                            ' (each station all events)')
#
#         plt.savefig(
#             os.path.join(
#                 path_to_use,
#                 r'cdf_hist_error_temp_dist_%s_2.png' % temp_freq,
#             ),
#             frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
#         plt.close()
#
#         #======================================================================
# #
# #==============================================================================
# #         # OK
# #         plt.ioff()
# #         fig = plt.figure(figsize=(24, 12), dpi=150)
# #
# #         ax = fig.add_subplot(111)
# #
# #         ax.plot(df_compare.index,
# #                 df_compare.mean_error_dwd,
# #                 alpha=.8,
# #                 c='b',  # colors_arr,
# #                 marker='d',
# #                 label='DWD Mean Err %0.2f'
# #                 % np.mean(df_compare.mean_error_dwd))
# #
# #         ax.plot(df_compare.index,
# #                 df_compare.max_error_dwd,
# #                 alpha=.8,
# #                 c='darkblue',  # colors_arr,
# #                 marker='d',
# #                 label='DWD Max Err %0.2f'
# #                 % np.mean(df_compare.max_error_dwd))
# #         ax.plot(df_compare.index,
# #                 df_compare.min_error_dwd,
# #                 alpha=.8,
# #                 c='lightblue',  # colors_arr,
# #                 marker='d',
# #                 label='DWD Min Err %0.2f'
# #                 % np.mean(df_compare.min_error_dwd))
# #         #=============================================================
# #         #
# #         #=============================================================
# #         ax.plot(df_compare.index,
# #                 df_compare.mean_error_dwd_netatmo,
# #                 alpha=.5,
# #                 c='r',  # colors_arr,
# #                 marker='*',
# #                 label='DWD-Netatmo Mean Err %0.2f'
# #                 % np.mean(df_compare.mean_error_dwd_netatmo))
# #
# #         ax.plot(df_compare.index,
# #                 df_compare.max_error_dwd_netatmo,
# #                 alpha=.5,
# #                 c='darkred',  # colors_arr,
# #                 marker='*',
# #                 label='DWD-Netatmo Max Err %0.2f'
# #                 % np.mean(df_compare.max_error_dwd_netatmo))
# #
# #         ax.plot(df_compare.index,
# #                 df_compare.min_error_dwd_netatmo,
# #                 alpha=.5,
# #                 c='salmon',  # colors_arr,
# #                 marker='*',
# #                 label='DWD-Netatmo Min Err %0.2f'
# #                 % np.mean(df_compare.min_error_dwd_netatmo))
# #
# #         #=============================================================
# #         #
# #         #=============================================================
# #         ax.plot(df_compare.index,
# #                 df_compare.mean_error_dwd_netatmo_unc,
# #                 alpha=.5,
# #                 c='g',  # colors_arr,
# #                 marker='x',
# #                 label='DWD-Netatmo Unc Mean Err %0.2f'
# #                 % np.mean(df_compare.mean_error_dwd_netatmo_unc))
# #
# #         ax.plot(df_compare.index,
# #                 df_compare.max_error_dwd_netatmo_unc,
# #                 alpha=.5,
# #                 c='darkgreen',  # colors_arr,
# #                 marker='x',
# #                 label='DWD-Netatmo Unc Max Err %0.2f'
# #                 % np.mean(df_compare.max_error_dwd_netatmo_unc))
# #
# #         ax.plot(df_compare.index,
# #                 df_compare.min_error_dwd_netatmo_unc,
# #                 alpha=.5,
# #                 c='lightgreen',  # colors_arr,
# #                 marker='x',
# #                 label='DWD-Netatmo Unc Min Err %0.2f'
# #                 % np.mean(df_compare.min_error_dwd_netatmo_unc))
#
#
# #         ax.set_title('Pearson Correlation Interpolated Quantiles from DWD or DWD-Netatmo\n '
# #                      'Precipitation of %s Extreme Events %s\n Events with Improvemnts %d / %d, Percentage %0.0f\n'
# #                      'Events with Improvemnts with OK 2percUnc %d / %d, Percentage %0.0f'
# #                      '\nEvents with Improvemnts with OK 5percUnc %d / %d, Percentage %0.0f'
# #                      '\nEventswith Improvemnts with OK 10percUnc %d / %d, Percentage %0.0f'
# #                      '\nEvents with Improvemnts with OK 20percUnc %d / %d, Percentage %0.0f'
# #                      % (temp_freq, _interp_acc_,
# #                         stations_with_improvements,
# #                          df_compare.pearson_corr_dwd_netatmo.shape[0],
# #                         percent_of_improvment,
# #                         stations_with_improvements_unc2perc,
# #                         df_compare.pearson_corr_dwd_netatmo_unc2perc.shape[0],
# #                         percent_of_improvment_unc2perc,
# #                         stations_with_improvements_unc5perc,
# #                         df_compare.pearson_corr_dwd_netatmo_unc2perc.shape[0],
# #                         percent_of_improvment_unc5perc,
# #                         stations_with_improvements_unc10perc,
# #                         df_compare.pearson_corr_dwd_netatmo_unc2perc.shape[0],
# #                         percent_of_improvment_unc10perc,
# #                         stations_with_improvements_unc20perc,
# #                         df_compare.pearson_corr_dwd_netatmo_unc2perc.shape[0],
# #                         percent_of_improvment_unc20perc))
# #         ax.grid(alpha=0.25)
# #         plt.setp(ax.get_xticklabels(), rotation=45)
# #
# #         # years = mdates.YearLocator(1)   # every year
# #         months = mdates.MonthLocator()  # every month
# #
# #         #hrs = mdates.HourLocator()
# #         # years_fmt = mdates.DateFormatter('%Y-%m-%d')
# #
# #         years_fmt = mdates.DateFormatter('%Y-%m-%d %H')
# #         ax.xaxis.set_major_locator(months)
# #
# #         ax.xaxis.set_major_formatter(years_fmt)
# #         # ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))
# #
# #         #ax.set_yticks(np.arange(0., 1.05, .10))
# #         ax.set_xlabel('Date of Event')
# #         ax.legend(loc='lower right')
# #         ax.set_ylabel('Error compared to observed DWD values')
# #
# #         plt.savefig(os.path.join(path_to_use,
# #                                  r'error_spatial_distribution_%s_events_dwd2.png' % temp_freq,
# #                                  ),
# #                     frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
# #         plt.close()
