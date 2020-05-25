
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np


import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from pathlib import Path

import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears

from mpl_toolkits.axes_grid1 import make_axes_locatable

from pandas.plotting import register_matplotlib_converters

from _00_additional_functions import KernelDensityEstimate

fit_kernel_fct = KernelDensityEstimate()


register_matplotlib_converters()

plt.ioff()

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


path_to_use = Path(
    r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo\180520_error_dist_paper')


events_60min = ['2018-06-11 16:00:00', '2018-09-06 18:00:00']
events_1440min = ['2018-05-14 00:00:00', '2019-07-28 00:00:00']


#==============================================================================
#
#==============================================================================
# for saving results
no_60 = {evt: {'C0': [], 'C1': [], 'C2': [], 'C3': [], 'C4': []}
         for evt in events_60min}
pp_60 = {evt: {'C0': [], 'C1': [], 'C2': [], 'C3': [], 'C4': []}
         for evt in events_60min}
no_1440 = {evt: {'C0': [], 'C1': [], 'C2': [], 'C3': [], 'C4': []}
           for evt in events_1440min}
pp_1440 = {evt: {'C0': [], 'C1': [], 'C2': [], 'C3': [], 'C4': []}
           for evt in events_1440min}


for temp_freq in ['60min', '1440min']:
    print(temp_freq)

    path_dwd_ppt_data = (r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes"
                         r"\NetAtmo_BW\ppt_all_dwd_%s_.csv" % temp_freq)
    df_dwd_edf = pd.read_csv(path_dwd_ppt_data,
                             sep=';', index_col=0,
                             parse_dates=True,
                             infer_datetime_format=True)

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

    # stns to remove from orig edf, because no data for 2015-2019
    df_dwd_edf.drop(columns=['00384', '13672', '05155'], inplace=True)

    df_dwd_edf = df_dwd_edf.loc[df_dwd.index, :]

    if temp_freq == '60min':
        event_list = events_60min

    if temp_freq == '1440min':
        event_list = events_1440min

    try:

        # df_dwd.index.intersection(df_dwd_edf.index):
        for event_date in event_list:
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

            #==================================================================
            # PLOT EACH EVENT
            #==================================================================
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

            if temp_freq == '60min':
                no_60[event_date]['C0'].append(np.array(no0))
                no_60[event_date]['C1'].append(np.array(no1))
                no_60[event_date]['C2'].append(np.array(no2))
                no_60[event_date]['C3'].append(np.array(no3))
                no_60[event_date]['C4'].append(np.array(no4))

                pp_60[event_date]['C0'].append(pp0)
                pp_60[event_date]['C1'].append(pp1)
                pp_60[event_date]['C2'].append(pp2)
                pp_60[event_date]['C3'].append(pp3)
                pp_60[event_date]['C4'].append(pp4)

            if temp_freq == '1440min':
                no_1440[event_date]['C0'].append(np.array(no0))
                no_1440[event_date]['C1'].append(np.array(no1))
                no_1440[event_date]['C2'].append(np.array(no2))
                no_1440[event_date]['C3'].append(np.array(no3))
                no_1440[event_date]['C4'].append(np.array(no4))

                pp_1440[event_date]['C0'].append(pp0)
                pp_1440[event_date]['C1'].append(pp1)
                pp_1440[event_date]['C2'].append(pp2)
                pp_1440[event_date]['C3'].append(pp3)
                pp_1440[event_date]['C4'].append(pp4)

    except Exception as msg:
        print(msg)


# get start and end date of events, for labeling
#==============================================================================
# # make Subplot
#==============================================================================
strt_ev1_60 = str(pd.to_datetime(events_60min[0]) - pd.Timedelta(hours=1))
strt_ev2_60 = str(pd.to_datetime(events_60min[1]) - pd.Timedelta(hours=1))

strt_ev1_1440 = str(pd.to_datetime(events_1440min[0]) -
                    pd.Timedelta(days=1))  # .split(' ')[0]
strt_ev2_1440 = str(pd.to_datetime(events_1440min[1]) -
                    pd.Timedelta(days=1))  # .split(' ')[0]

plt.ioff()


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,
                                             figsize=(24, 14), dpi=100)


ax1.plot(pp_60[events_60min[0]]['C0'][0], no_60[events_60min[0]]['C0'][0],
         c='k',  marker='_', alpha=0.85)
ax1.plot(pp_60[events_60min[0]]['C1'][0], no_60[events_60min[0]]['C1'][0],
         c='r',  marker='+', alpha=0.65)
ax1.plot(pp_60[events_60min[0]]['C2'][0], no_60[events_60min[0]]['C2'][0],
         c='b',  marker='1', alpha=0.65)
ax1.plot(pp_60[events_60min[0]]['C3'][0], no_60[events_60min[0]]['C3'][0],
         c='g',  marker='o', alpha=0.65)
ax1.plot(pp_60[events_60min[0]]['C4'][0], no_60[events_60min[0]]['C4'][0],
         c='orange', marker='.', alpha=0.65)

ax1.axvline(x=0, color='k', linestyle='--', alpha=0.35)

ax2.plot(pp_60[events_60min[1]]['C0'][0], no_60[events_60min[1]]['C0'][0],
         c='k',  marker='_', alpha=0.85)
ax2.plot(pp_60[events_60min[1]]['C1'][0], no_60[events_60min[1]]['C1'][0],
         c='r', marker='+', alpha=0.65)
ax2.plot(pp_60[events_60min[1]]['C2'][0], no_60[events_60min[1]]['C2'][0],
         c='b',  marker='1', alpha=0.65)
ax2.plot(pp_60[events_60min[1]]['C3'][0], no_60[events_60min[1]]['C3'][0],
         c='g',  marker='o', alpha=0.65)
ax2.plot(pp_60[events_60min[1]]['C4'][0], no_60[events_60min[1]]['C4'][0],
         c='orange',  marker='.', alpha=0.65)

ax2.axvline(x=0, color='k', linestyle='--', alpha=0.35)

ax3.plot(pp_1440[events_1440min[0]]['C0'][0],
         no_1440[events_1440min[0]]['C0'][0],
         c='k', marker='_', alpha=0.85)
ax3.plot(pp_1440[events_1440min[0]]['C1'][0],
         no_1440[events_1440min[0]]['C1'][0],
         c='r',  marker='+', alpha=0.65)
ax3.plot(pp_1440[events_1440min[0]]['C2'][0],
         no_1440[events_1440min[0]]['C2'][0],
         c='b',  marker='1', alpha=0.65)
ax3.plot(pp_1440[events_1440min[0]]['C3'][0],
         no_1440[events_1440min[0]]['C3'][0],
         c='g',  marker='o', alpha=0.65)
ax3.plot(pp_1440[events_1440min[0]]['C4'][0],
         no_1440[events_1440min[0]]['C4'][0],
         c='orange',  marker='.', alpha=0.65)

ax3.axvline(x=0, color='k', linestyle='--', alpha=0.35)

ax4.plot(pp_1440[events_1440min[1]]['C0'][0],
         no_1440[events_1440min[1]]['C0'][0],
         c='k', label='C0', marker='_', alpha=0.85)
ax4.plot(pp_1440[events_1440min[1]]['C1'][0],
         no_1440[events_1440min[1]]['C1'][0],
         c='r', label='C1', marker='+', alpha=0.65)
ax4.plot(pp_1440[events_1440min[1]]['C2'][0],
         no_1440[events_1440min[1]]['C2'][0],
         c='b', label='C2', marker='1', alpha=0.65)
ax4.plot(pp_1440[events_1440min[1]]['C3'][0],
         no_1440[events_1440min[1]]['C3'][0],
         c='g', label='C3', marker='o', alpha=0.65)
ax4.plot(pp_1440[events_1440min[1]]['C4'][0],
         no_1440[events_1440min[1]]['C4'][0],
         c='orange', label='C4', marker='.', alpha=0.5)

ax4.axvline(x=0, color='k', linestyle='-', alpha=0.35)

ax1.title.set_text('Hourly event between %s and %s'
                   % (strt_ev1_60, events_60min[0]))
ax2.title.set_text('Hourly event between %s and %s'
                   % (strt_ev2_60, events_60min[1]))
ax3.title.set_text('Daily event between %s and %s'
                   % (strt_ev1_1440,
                      events_1440min[0]))  # .split(' ')[0]
ax4.title.set_text('Daily event between %s and %s'
                   % (strt_ev2_1440,
                      events_1440min[1])
                   )

#'Pears Corr. C0=%0.2f--C3=%.2f'
#' Spear Corr.C0=%0.2f--C3=%.2f'
#% (event_date, corr_dwd, corr_dwd_net,
#   spr_dwd, spr_dwd_net))

ax1.grid(color='k', linestyle='--', alpha=0.2)
ax2.grid(color='k', linestyle='--', alpha=0.2)
ax3.grid(color='k', linestyle='--', alpha=0.2)
ax4.grid(color='k', linestyle='--', alpha=0.2)

ax1.set_ylabel('Density')
ax3.set_ylabel('Density')
ax3.set_xlabel('Error [mm]')
ax4.set_xlabel('Error [mm]')


# axins = inset_axes(ax4,  # here using axis of the lowest plot
#                    width="4%",  # width = 5% of parent_bbox width
#                    height="150%",  # height : 340% good for a (4x4) Grid
#                    loc='lower left',
#                    bbox_to_anchor=(1.09, 0.35, 1, 1),
#                    bbox_transform=ax4.transAxes,
#                    borderpad=0,
#                    )

# for ax in fig.get_axes():
#    ax.label_outer()
# fig.tight_layout()

# handles, labels = ax4.get_legend_handles_labels()
# unique_labels = np.unique(labels)

# ax_legend = fig.add_axes([0.1725, 0.02525, 0.68, 0.0225], zorder=3)
ax4.legend(loc=0)
# fig.tight_layout()


plt.savefig(
    os.path.join(
        path_to_use,
        r'error_dist_%s_%s21.png'
        % (temp_freq, str(event_date).replace(':', '_')),
    ),
    frameon=True, papertype='a4',
    bbox_inches='tight', pad_inches=.2)
plt.close()
