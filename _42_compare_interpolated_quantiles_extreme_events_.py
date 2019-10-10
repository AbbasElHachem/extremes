import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from scipy.stats import spearmanr as spr
from scipy.stats import pearsonr as pears

plt.ioff()
plt.rcParams.update({'font.size': 12})
plt.rcParams.update({'axes.labelsize': 12})

path_dwd_edf_data = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\NetAtmo_BW\edf_ppt_all_dwd_daily_.csv"

df_dwd_edf = pd.read_csv(path_dwd_edf_data, sep=';', index_col=0,
                         parse_dates=True, infer_datetime_format=True)

for i in range(12):
    try:
        i = int(i)
        path_interpolated_using_netatmo = (
            r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo"
            r"\interpolated_quantiles_daily_data_basedon_qunatiles_2019_08_08_season_using_netatmo_only_grp_%d_.csv" % i)

        path_interpolated_using_dwd = (
            r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo"
            r"\interpolated_quantiles_dwd_daily_data_basedon_qunatiles_2019_08_08_season_using_dwd_only_grp_%d_.csv" % i)

        path_interpolated_using_netatmo_dwd = (
            r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo"
            r"\interpolated_quantiles_dwd_daily_data_basedon_quantiles_2019_08_08_season_using_dwd_netamo_grp_%d_.csv" % i)

        df_netatmo = pd.read_csv(path_interpolated_using_netatmo,
                                 sep=';', index_col=0, parse_dates=True,
                                 infer_datetime_format=True)

        df_dwd = pd.read_csv(path_interpolated_using_dwd,
                             sep=';', index_col=0, parse_dates=True,
                             infer_datetime_format=True)

        df_netatmo_dwd = pd.read_csv(path_interpolated_using_netatmo_dwd,
                                     sep=';', index_col=0, parse_dates=True,
                                     infer_datetime_format=True)

        df_compare = pd.DataFrame(index=df_netatmo.index)

        for stn_ in df_netatmo.columns:
            print(stn_)
            for event_date in df_netatmo.index[1:]:
                interpolated_quantile_netatmo = df_netatmo.loc[event_date, stn_]

                interpolated_quantile_dwd = df_dwd.loc[event_date, stn_]

                interpolated_quantile_netatmo_dwd = df_netatmo_dwd.loc[event_date, stn_]

                original_quantile = df_dwd_edf.loc[event_date, stn_]

                if original_quantile >= 0.6:
                    df_compare.loc[event_date,
                                   'original_quantile'] = original_quantile
                    df_compare.loc[event_date,
                                   'interpolated_quantile_dwd'] = interpolated_quantile_dwd
                    df_compare.loc[event_date,
                                   'interpolated_quantile_netatmo'] = interpolated_quantile_netatmo
                    df_compare.loc[event_date,
                                   'interpolated_quantile_netatmo_dwd'] = interpolated_quantile_netatmo_dwd

                else:
                    df_compare.loc[event_date,
                                   'original_quantile'] = np.nan
                    df_compare.loc[event_date,
                                   'interpolated_quantile_dwd'] = np.nan
                    df_compare.loc[event_date,
                                   'interpolated_quantile_netatmo'] = np.nan
                    df_compare.loc[event_date,
                                   'interpolated_quantile_netatmo_dwd'] = interpolated_quantile_netatmo_dwd

            df_compare = df_compare[df_compare > 0]
            df_compare.dropna(how='any', inplace=True)

            values_x = df_compare['original_quantile'].values
            values_dwd = df_compare['interpolated_quantile_dwd'].values
            values_netatmo = df_compare['interpolated_quantile_netatmo'].values
            values_netatmo_dwd = df_compare['interpolated_quantile_netatmo_dwd'].values

            # plot the stations in shapefile, look at the results of agreements

            # calculate correlations (pearson and spearman)
            corr_dwd = pears(values_x, values_dwd)[0]
            rho_dwd = spr(values_x, values_dwd)[0]

            corr_netatmo = pears(values_x, values_netatmo)[0]
            rho_netatmo = spr(values_x, values_netatmo)[0]

            corr_netatmo_dwd = pears(values_x, values_netatmo_dwd)[0]
            rho_netatmo_dwd = spr(values_x, values_netatmo_dwd)[0]

            #########################################################
            plt.ioff()
            fig = plt.figure(figsize=(24, 12), dpi=150)

            ax = fig.add_subplot(111)

            # plot 45 deg line
            _min = min(values_x.min(), values_dwd.min())
            _max = max(values_x.max(), values_dwd.max())

            ax.plot([_min, _max], [_min, _max],
                    c='k', linestyle='--', alpha=0.4)

            # set plot limit
            ax.set_xlim(_min - 0.01, 1.01)
            ax.set_ylim(_min - 0.01, 1.01)

            ax.scatter(values_x,
                       values_dwd,
                       alpha=.8,
                       c='r',  # colors_arr,
                       s=15,
                       marker='d',
                       # cmap=plt.get_cmap('viridis'),
                       label='DWD Interpolated %d Events' % values_dwd.shape[0])

            ax.set_title('Observed and Interpolated Quantiles for Daily Extreme Events \n DWD Station %s \n'
                         'Pearson Cor=%0.3f; Spearman Cor=%0.3f'
                         % (stn_, corr_dwd, rho_dwd))
            ax.grid(alpha=0.25)

            ax.set_xlabel('Original Quantiles')
            ax.legend(loc='lower right')
            ax.set_ylabel('Interpolated Quantiles')

            plt.savefig((r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo'
                         r'\oberserved_vs_interpolated_quantiles_extreme_daily_events_stn_%s_dwd.png' % (stn_)),
                        frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
            plt.close()

            ###################################################################

            fig = plt.figure(figsize=(24, 12), dpi=150)

            ax = fig.add_subplot(111)

            _min = min(values_x.min(), values_netatmo.min())
            _max = max(values_x.max(), values_netatmo.max())
            ax.plot([_min, _max], [_min, _max],
                    c='k', linestyle='--', alpha=0.4)

            # set plot limit
            ax.set_xlim(_min - 0.01, 1.01)
            ax.set_ylim(_min - 0.01, 1.01)

            ax.scatter(values_x,
                       values_netatmo,
                       alpha=.8,
                       c='b',  # colors_arr,
                       s=15,
                       marker='d',
                       # cmap=plt.get_cmap('viridis'),
                       label='Netatmo Interpolated %d Events' % values_netatmo.shape[0])

            ax.set_title('Observed and Interpolated Quantiles for Daily Extreme Events \n DWD Station %s \n'
                         'Pearson Cor=%0.3f; Spearman Cor=%0.3f'
                         % (stn_, corr_netatmo, rho_netatmo))
            ax.grid(alpha=0.25)

            # ax.set_xlim([-0.1, max(interpolated_vals.values.max(),
            #                       ppt_vals_season.max()) + 2])
            ax.set_xlabel('Original Quantiles')
            ax.legend(loc='lower right')
            ax.set_ylabel('Interpolated Quantiles')

            plt.savefig((r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo'
                         r'\oberserved_vs_interpolated_quantiles_extreme_daily_events_stn_%s_netatmo.png' % (stn_)),
                        frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
            plt.close()

            ###################################################################

            fig = plt.figure(figsize=(24, 12), dpi=150)

            ax = fig.add_subplot(111)

            _min = min(values_x.min(), values_netatmo_dwd.min())
            _max = max(values_x.max(), values_netatmo_dwd.max())
            ax.plot([_min, _max], [_min, _max],
                    c='k', linestyle='--', alpha=0.4)

            ax.scatter(values_x,
                       values_netatmo_dwd,
                       alpha=.8,
                       c='g',  # colors_arr,
                       s=15,
                       marker='d',
                       # cmap=plt.get_cmap('viridis'),
                       label='DWD-Netatmo Interpolated %d Events' % values_netatmo_dwd.shape[0])

            ax.set_title('Observed and Interpolated Quantiles for Daily Extreme Events \n DWD Station %s \n'
                         'Pearson Cor=%0.3f; Spearman Cor=%0.3f'
                         % (stn_, corr_netatmo_dwd, rho_netatmo_dwd))
            ax.grid(alpha=0.25)
            ax.set_xlim(_min - 0.01, 1.01)
            ax.set_ylim(_min - 0.01, 1.01)
            # ax.set_xlim([-0.1, max(interpolated_vals.values.max(),
            #                       ppt_vals_season.max()) + 2])
            ax.set_xlabel('Original Quantiles')
            ax.legend(loc='lower right')
            ax.set_ylabel('Interpolated Quantiles')

            plt.savefig((r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo'
                         r'\oberserved_vs_interpolated_quantiles_extreme_daily_events_stn_%s_netatmo_dwd.png' % (stn_)),
                        frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
            plt.close()

            plt.ioff()

            ###################################################################

#             fig = plt.figure(figsize=(24, 12), dpi=150)
#
#             ax = fig.add_subplot(111)
#
#             # plot 45 deg line
#             _min = min(values_x.min(), values_dwd.min())
#             _max = max(values_x.max(), values_dwd.max())
#
#             ax.plot([_min, _max], [_min, _max],
#                     c='k', linestyle='--', alpha=0.4)
#
#             # set plot limit
#             ax.set_xlim(_min - 0.01, 1.01)
#             ax.set_ylim(_min - 0.01, 1.01)
#
#             ax.scatter(values_x,
#                        values_dwd,
#                        alpha=.6,
#                        c='r',  # colors_arr,
#                        s=15,
#                        marker='d',
#                        # cmap=plt.get_cmap('viridis'),
#                        label='DWD Interpolated ')
#             ax.scatter(values_x,
#                        values_netatmo,
#                        alpha=.6,
#                        c='b',  # colors_arr,
#                        s=15,
#                        marker='x',
#                        # cmap=plt.get_cmap('viridis'),
#                        label='Netatmo Interpolated ')
#             ax.scatter(values_x,
#                        values_netatmo_dwd,
#                        alpha=.6,
#                        c='g',  # colors_arr,
#                        s=15,
#                        marker='o',
#                        # cmap=plt.get_cmap('viridis'),
#                        label='DWD-Netatmo Interpolated  ')
#             ax.set_title('Observed and Interpolated Quantiles for Daily Extreme Events \n DWD Station %s \n'
#                          #'Pearson Cor=%0.3f; Spearman Cor=%0.3f'
#                          % (stn_))  # , corr_dwd, rho_dwd))
#             ax.grid(alpha=0.25)
#
#             ax.set_xlabel('Original Quantiles')
#             ax.legend(loc='lower right')
#             ax.set_ylabel('Interpolated Quantiles')
#
#             plt.savefig((r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\oridinary_kriging_compare_DWD_Netatmo'
#                          r'\oberserved_vs_interpolated_quantiles_extreme_daily_events_stn_%s_dwdnetatmo_.png' % (stn_)),
#                         frameon=True, papertype='a4', bbox_inches='tight', pad_inches=.2)
#             plt.close()
#             break
    except Exception as msg:
        print(msg)
