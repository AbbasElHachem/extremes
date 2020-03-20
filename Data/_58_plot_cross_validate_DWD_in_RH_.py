
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

from _00_additional_functions import list_all_full_path
register_matplotlib_converters()

plt.ioff()
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'axes.labelsize': 16})


#main_dir = Path(r'X:\staff\elhachem\2020_10_03_Rheinland_Pfalz')

main_dir = Path(
    r'/run/media/abbas/EL Hachem 2019/home_office/2020_10_03_Rheinland_Pfalz')

plot_results_per_stn = True
plot_results_per_event = True

used_data_acc = r'98_20_5'



if plot_results_per_stn:
    # , '60min' '180min', '360min', '720min', '1440min']:
    for temp_freq in ['60min']:
        print(temp_freq)

        path_to_netatmo_interpolation = (main_dir / r'ppt_cross_valid_RH_60min' /
            (r'df_interpolated_netatmo_only_%s_data_%s.csv' % (
                temp_freq, used_data_acc)))

        path_to_dwd_interpolation = (main_dir / 
             r'ppt_cross_valid_RH_60min' /
            (r'df_interpolated_dwd_only_%s_data_%s.csv' %(
                 temp_freq, used_data_acc)))
        
        path_to_dwd_netatmo_interpolation = (main_dir /
             r'ppt_cross_valid_RH_60min' /
            (r'df_interpolated_dwd_netatmos_comb_%s_data_%s.csv' %(
                 temp_freq, used_data_acc)))
                
    #     path_to_Qt_ok_un_first_flt_comb_ = main_dir / (
    #         r'Qt_ok_ok_un_3_first_flt__comb_%s' % temp_freq)
    #     Qt_ok_un_first_flt_comb_ = list_all_full_path(
    #         '.csv', path_to_Qt_ok_un_first_flt_comb_)


        path_dwd_ppt_data = main_dir / (
            r"ppt_dwd_2014_2019_%s_no_freezing_5deg.csv" % temp_freq)
        df_dwd_ppt = pd.read_csv(path_dwd_ppt_data, sep=';', index_col=0,
                                 parse_dates=True, infer_datetime_format=True)
        # df_dwd_edf.dropna(inplace=True)

        df_improvements = pd.DataFrame(
            index=df_dwd_ppt.columns,
            columns=['pearson_corr_dwd_',
                     'spearman_corr_dwd_',
                     'pearson_corr_netatmo_',
                     'spearman_corr_netatmo_',
                     'pearson_corr_dwd_netatmo',
                     'spearman_corr_dwd_netatmo',
                     ])

        min_orig_qnt_thr = 0.


        #======================================================================
        #
        #======================================================================
        df_dwd = pd.read_csv(path_to_dwd_interpolation,
                             sep=';', index_col=0, parse_dates=True,
                             infer_datetime_format=True)

        df_netatmo = pd.read_csv(path_to_netatmo_interpolation,
                                     sep=';', index_col=0,
                                     parse_dates=True,
                                     infer_datetime_format=True)

        df_netatmo_dwd = pd.read_csv(path_to_dwd_netatmo_interpolation,
                                     sep=';', index_col=0,
                                     parse_dates=True,
                                     infer_datetime_format=True)

        print(df_netatmo_dwd.isna().sum().max())


        try:

            for stn_ in df_dwd_ppt.columns:
                #                 if stn_ == '00384':
                #                     raise Exception
                #                     pass
                #                     print(stn_)
                df_compare = pd.DataFrame(index=df_dwd.index)

                # :
                for event_date in df_dwd.index.intersection(df_dwd_ppt.index):
                    # print(event_date)

                    # event_date = '2019-07-27 20:00:00'
                    ppt_stn_orig = df_dwd_ppt.loc[event_date, stn_]
                    ppt_stn_interp_dwd = df_dwd.loc[event_date, stn_]
                    ppt_stn_interp_netatmo = df_netatmo.loc[event_date, stn_]
                    ppt_stn_interp_netatmo_dwd = df_netatmo_dwd.loc[event_date, stn_]
#                     edf_stn_interp_netatmo_dwd_unc05perc = df_netatmo_dwd_unc05perc.loc[
#                         event_date, stn_]
                    if ((ppt_stn_orig >= 0) and
                        (ppt_stn_interp_dwd >= 0) and
                            (ppt_stn_interp_netatmo >= 0 and
                             (ppt_stn_interp_netatmo_dwd >= 0)
                             )):
                        # and(edf_stn_interp_netatmo_dwd_unc05perc >= 0)
                        df_compare.loc[
                            event_date,
                            'original_quantile'] = ppt_stn_orig
                        df_compare.loc[
                            event_date,
                            'interpolated_quantile_dwd'] = ppt_stn_interp_dwd
                        
                        df_compare.loc[
                            event_date,
                            'interpolated_quantile_netatmo'] = ppt_stn_interp_netatmo
                            
                        df_compare.loc[
                            event_date,
                            'interpolated_quantile_netatmo_dwd'] = ppt_stn_interp_netatmo_dwd
                        

                    else:
                        df_compare.loc[event_date,
                                       'original_quantile'] = np.nan
                        df_compare.loc[event_date,
                                       'interpolated_quantile_dwd'] = np.nan

                        df_compare.loc[event_date,
                                       'interpolated_quantile_netatmo'] = np.nan
                        df_compare.loc[event_date,
                                       'interpolated_quantile_netatmo_dwd'] = np.nan

                df_compare = df_compare[df_compare > 0]
                df_compare.dropna(how='any', inplace=True)
                
                if df_compare.values.shape[0] > 1:  # at least 20evnts
                    values_x = df_compare['original_quantile'].values
                    values_dwd = df_compare['interpolated_quantile_dwd'].values
                    values_netatmo = df_compare['interpolated_quantile_netatmo'].values

                    values_netatmo_dwd = df_compare['interpolated_quantile_netatmo_dwd'].values



                    # calculate correlations (pearson and spearman) dwd-dwd
                    corr_dwd = pears(values_x, values_dwd)[0]
                    rho_dwd = spr(values_x, values_dwd)[0]
                    print(corr_dwd)
                    # netatmo 
                    corr_netatmo = pears(
                        values_x, values_netatmo)[0]
                    rho_netatmo = spr(
                        values_x, values_netatmo)[0]

                    # netatmo dwd
                    corr_netatmo_dwd = pears(values_x, values_netatmo_dwd)[0]
                    rho_netatmo_dwd = spr(values_x, values_netatmo_dwd)[0]
                    


                    df_improvements.loc[stn_, 'pearson_corr_dwd_'] = corr_dwd
                    df_improvements.loc[stn_, 'spearman_corr_dwd_'] = rho_dwd
                    
                    df_improvements.loc[stn_,
                                        'pearson_corr_netatmo_'] = corr_netatmo
                    df_improvements.loc[stn_,
                                        'spearman_corr_netatmo_'] = rho_netatmo
                                        
                    df_improvements.loc[stn_,
                                        'pearson_corr_dwd_netatmo'] = corr_netatmo_dwd
                    df_improvements.loc[stn_,
                                        'spearman_corr_dwd_netatmo'] = rho_netatmo_dwd
                    


            df_improvements.dropna(how='all', inplace=True)



        except Exception as msg:
            print(msg)

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

    #     # OK with netatmo
        stations_with_improvements_netatmo= sum(i >= j for (i, j) in zip(
            df_improvements.pearson_corr_netatmo_.values,
            df_improvements.pearson_corr_dwd_.values))
        stations_without_improvements_netatmo = sum(i < j for (i, j) in zip(
            df_improvements.pearson_corr_netatmo_.values,
            df_improvements.pearson_corr_dwd_.values))
        percent_of_improvment_netatmo= 100 * (
            stations_with_improvements_netatmo /
            df_improvements.pearson_corr_dwd_netatmo.shape[0])


        ######################################################
        mean_pearson_correlation_dwd_only = df_improvements.pearson_corr_dwd_.mean()
        mean_pearson_correlation_netatmo_only = df_improvements.pearson_corr_netatmo_.mean()
        mean_pearson_correlation_dwd_netatmo = df_improvements.pearson_corr_dwd_netatmo.mean()

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
                df_improvements.pearson_corr_netatmo_,
                alpha=.8,
                c='k',  # colors_arr,
                marker='2',
                label='Netatmo Interpolation %0.2f'
                % mean_pearson_correlation_netatmo_only)
        
        ax.plot(df_improvements.index,
                df_improvements.pearson_corr_dwd_netatmo,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %0.2f'
                % mean_pearson_correlation_dwd_netatmo)

        ax.set_title('Pearson Correlation Interpolated Rainfall'
                     ' from DWD or Netatmo or DWD-Netatmo\n '
                     
                     'Precipitation of %s Intense Events \n'
                     'Stations with Improvemnts with OK Netatmo no Filter %d / %d, Percentage %0.0f'
                     '\nStations with Improvemnts with OK DWD-Netatmo %d / %d,'
                     ' Percentage %0.0f\n'
                     
                     %(temp_freq, stations_with_improvements_netatmo,
                        df_improvements.pearson_corr_netatmo_.shape[0],
                        percent_of_improvment_netatmo,
                         stations_with_improvements,
                         df_improvements.pearson_corr_dwd_netatmo.shape[0],
                        percent_of_improvment ))

        plt.setp(ax.get_xticklabels(), rotation=90)
        ax.grid(alpha=0.25)
        ax.set_yticks(np.arange(0, 1.05, .10))
        ax.set_xlabel('DWD Stations')
        ax.legend(loc='lower left')
        ax.set_ylabel('Pearson Correlation')

        plt.savefig((main_dir / r'ppt_cross_valid_RH_60min' / (
             r'ppt_temporal_pears_corr_%s_events_%s.png'
                                  % (temp_freq, used_data_acc ))),
                    papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

#===============================================================================
# 
#===============================================================================

        stations_with_improvements = sum(i >= j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo.values,
            df_improvements.spearman_corr_dwd_.values))

        stations_without_improvements = sum(i < j for (i, j) in zip(
            df_improvements.spearman_corr_dwd_netatmo.values,
            df_improvements.spearman_corr_dwd_.values))

        percent_of_improvment = 100 * (
            stations_with_improvements /
           df_improvements.spearman_corr_dwd_netatmo.shape[0])

        id_stns_no_impv = [id_stn for (id_stn, i, j) in zip(
            df_improvements.index,
            df_improvements.spearman_corr_dwd_netatmo.values,
            df_improvements.spearman_corr_dwd_.values) if i < j]

    #     # OK with netatmo
        stations_with_improvements_netatmo= sum(i >= j for (i, j) in zip(
            df_improvements.spearman_corr_netatmo_.values,
            df_improvements.spearman_corr_dwd_.values))
        stations_without_improvements_netatmo = sum(i < j for (i, j) in zip(
            df_improvements.spearman_corr_netatmo_.values,
            df_improvements.spearman_corr_dwd_.values))
        percent_of_improvment_netatmo= 100 * (
            stations_with_improvements_netatmo /
            df_improvements.spearman_corr_dwd_netatmo.shape[0])

        ####
        mean_spearman_correlation_dwd_only = df_improvements.spearman_corr_dwd_.mean()
        mean_spearman_correlation_netatmo_only = df_improvements.spearman_corr_netatmo_.mean()

        mean_spearman_correlation_dwd_netatmo = df_improvements.spearman_corr_dwd_netatmo.mean()

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
                % mean_spearman_correlation_dwd_only)
        
        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_netatmo_,
                alpha=.8,
                c='k',  # colors_arr,
                marker='2',
                label='Netatmo Interpolation %0.2f'
                % mean_spearman_correlation_netatmo_only)
        
        ax.plot(df_improvements.index,
                df_improvements.spearman_corr_dwd_netatmo,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %0.2f'
                % mean_spearman_correlation_dwd_netatmo)

        ax.set_title('Spearman Correlation Interpolated Rainfall from DWD or Netatmo or DWD-Netatmo\n '
                     
                     'Precipitation of %s Intense Events \n'
                     'Stations with Improvemnts with OK Netatmo no Filter %d / %d, Percentage %0.0f'
                     '\nStations with Improvemnts with OK DWD-Netatmo %d / %d,'
                     ' Percentage %0.0f\n'
                     
                     % (temp_freq, stations_with_improvements_netatmo,
                        df_improvements.spearman_corr_netatmo_.shape[0],
                        percent_of_improvment_netatmo,
                         stations_with_improvements,
                         df_improvements.spearman_corr_dwd_netatmo.shape[0],
                        percent_of_improvment ))

        plt.setp(ax.get_xticklabels(), rotation=90)
        ax.grid(alpha=0.25)
        ax.set_yticks(np.arange(0, 1.05, .10))
        ax.set_xlabel('DWD Stations')
        ax.legend(loc='lower left')
        ax.set_ylabel('Spearman Correlation')

        plt.savefig((main_dir / r'ppt_cross_valid_RH_60min' / (
             r'ppt_temporal_spr_corr_%s_events_%s.png'
                                  % (temp_freq, used_data_acc ))),
                    papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()
        
        
#===============================================================================
# 
#===============================================================================

if plot_results_per_event:
    # , '60min' '180min', '360min', '720min', '1440min']:
    for temp_freq in ['60min']:
        print(temp_freq)

        path_to_netatmo_interpolation = (main_dir / r'ppt_cross_valid_RH_60min' /
            (r'df_interpolated_netatmo_only_%s_data_%s.csv' % (
                temp_freq, used_data_acc)))

        path_to_dwd_interpolation = (main_dir / 
             r'ppt_cross_valid_RH_60min' /
            (r'df_interpolated_dwd_only_%s_data_%s.csv' %(
                 temp_freq, used_data_acc)))
        
        path_to_dwd_netatmo_interpolation = (main_dir /
             r'ppt_cross_valid_RH_60min' /
            (r'df_interpolated_dwd_netatmos_comb_%s_data_%s.csv' %(
                 temp_freq, used_data_acc)))
                
    #     path_to_Qt_ok_un_first_flt_comb_ = main_dir / (
    #         r'Qt_ok_ok_un_3_first_flt__comb_%s' % temp_freq)
    #     Qt_ok_un_first_flt_comb_ = list_all_full_path(
    #         '.csv', path_to_Qt_ok_un_first_flt_comb_)


        path_dwd_ppt_data = main_dir / (
            r"ppt_dwd_2014_2019_%s_no_freezing_5deg.csv" % temp_freq)
        df_dwd_ppt = pd.read_csv(path_dwd_ppt_data, sep=';', index_col=0,
                                 parse_dates=True, infer_datetime_format=True)
        # df_dwd_edf.dropna(inplace=True)

        df_improvements = pd.DataFrame(
            index=df_dwd_ppt.index,
            columns=['pearson_corr_dwd_',
                     'spearman_corr_dwd_',
                     'pearson_corr_netatmo_',
                     'spearman_corr_netatmo_',
                     'pearson_corr_dwd_netatmo',
                     'spearman_corr_dwd_netatmo',
                     ])
        
        #======================================================================
        #
        #======================================================================
        df_dwd = pd.read_csv(path_to_dwd_interpolation,
                             sep=';', index_col=0, parse_dates=True,
                             infer_datetime_format=True)

        df_netatmo = pd.read_csv(path_to_netatmo_interpolation,
                                     sep=';', index_col=0,
                                     parse_dates=True,
                                     infer_datetime_format=True)

        df_netatmo_dwd = pd.read_csv(path_to_dwd_netatmo_interpolation,
                                     sep=';', index_col=0,
                                     parse_dates=True,
                                     infer_datetime_format=True)
        df_dwd_ppt = df_dwd_ppt.loc[df_dwd.index, df_dwd.columns]
        
        
        #=======================================================================
        # plt.scatter(df_dwd_ppt.values.ravel(), df_dwd_ppt.values.ravel())
        # plt.scatter(df_dwd_ppt.values.ravel(), df_dwd.values.ravel())
        # plt.scatter(df_dwd_ppt.values.ravel(), df_netatmo.values.ravel())
        # plt.scatter(df_dwd_ppt.values.ravel(), df_netatmo_dwd.values.ravel())
        #=======================================================================
        df_compare = pd.DataFrame(index=df_dwd.index)
        try:

            for event_date in df_dwd.index:
                (orig_edf_vals, dwd_interp_vals,
                 netatmo_interp_vals,
                 netatmo_dwd_interp_vals) = [], [], [], []
                print(event_date)
                original_quantile = df_dwd_ppt.loc[event_date, :]
                interpolated_quantile_dwd = df_dwd.loc[event_date, :]
                interpolated_quantile_netatmo_dwd = df_netatmo_dwd.loc[event_date, :]
                interpolated_quantile_netatmo = df_netatmo.loc[
                    event_date, :]

                for stn_ in original_quantile.index:
                    # print(stn_)
                    edf_stn_orig = original_quantile.loc[stn_]
                    edf_stn_interp_dwd = interpolated_quantile_dwd.loc[stn_]
                    edf_stn_interp_netatmo_dwd = interpolated_quantile_netatmo_dwd.loc[stn_]
                    edf_stn_interp_netatmo = interpolated_quantile_netatmo.loc[stn_]

                    if ((edf_stn_orig >= 0) and
                        (edf_stn_interp_dwd >= 0) and
                            (edf_stn_interp_netatmo_dwd >= 0) and

                            (edf_stn_interp_netatmo >= 0)):

                        orig_edf_vals.append(edf_stn_orig)
                        dwd_interp_vals.append(edf_stn_interp_dwd)
                        netatmo_dwd_interp_vals.append(
                            edf_stn_interp_netatmo_dwd)

                        netatmo_interp_vals.append(
                            edf_stn_interp_netatmo)

                corr_dwd = pears(orig_edf_vals, dwd_interp_vals)[0]
                rho_dwd = spr(orig_edf_vals, dwd_interp_vals)[0]
                # print(corr_dwd, rho_dwd)

                corr_netatmo_dwd = pears(
                    orig_edf_vals, netatmo_dwd_interp_vals)[0]
                rho_netatmo_dwd = spr(
                    orig_edf_vals, netatmo_dwd_interp_vals)[0]
                # print(corr_netatmo_dwd, rho_netatmo_dwd)
                corr_netatmo = pears(
                    orig_edf_vals, netatmo_interp_vals)[0]
                rho_netatmo= spr(
                    orig_edf_vals, netatmo_interp_vals)[0]



                df_compare.loc[event_date, 'pearson_corr_dwd_'] = corr_dwd
                df_compare.loc[event_date, 'spearman_corr_dwd_'] = rho_dwd
                df_compare.loc[event_date,
                               'pearson_corr_dwd_netatmo'] = corr_netatmo_dwd
                df_compare.loc[event_date,
                               'spearman_corr_dwd_netatmo'] = rho_netatmo_dwd

                df_compare.loc[event_date,
                               'pearson_corr_netatmo_'] = corr_netatmo
                df_compare.loc[event_date,
                               'spearman_corr_netatmo_'] = rho_netatmo


                # df_compare = df_compare[df_compare > 0]
                #df_compare.dropna(how='all', inplace=True)

        except Exception as msg:
            print(msg)
        stations_with_improvements = sum(i >= j for (i, j) in zip(
            df_compare.pearson_corr_dwd_netatmo.values,
            df_compare.pearson_corr_dwd_.values))

        stations_without_improvements = sum(i < j for (i, j) in zip(
            df_compare.pearson_corr_dwd_netatmo.values,
            df_compare.pearson_corr_dwd_.values))

        percent_of_improvment = 100 * (stations_with_improvements /
                                       df_improvements.pearson_corr_dwd_netatmo.shape[0])

        id_stns_no_impv = [id_stn for (id_stn, i, j) in zip(
            df_compare.index,
            df_compare.pearson_corr_dwd_netatmo.values,
            df_compare.pearson_corr_dwd_.values) if i < j]

    #     # OK with netatmo
        stations_with_improvements_netatmo= sum(i >= j for (i, j) in zip(
            df_compare.pearson_corr_netatmo_.values,
            df_compare.pearson_corr_dwd_.values))
        stations_without_improvements_netatmo = sum(i < j for (i, j) in zip(
            df_compare.pearson_corr_netatmo_.values,
            df_compare.pearson_corr_dwd_.values))
        percent_of_improvment_netatmo= 100 * (
            stations_with_improvements_netatmo /
            df_compare.pearson_corr_dwd_netatmo.shape[0])


        ######################################################
        mean_pearson_correlation_dwd_only = df_improvements.pearson_corr_dwd_.mean()
        mean_pearson_correlation_netatmo_only = df_improvements.pearson_corr_netatmo_.mean()
        mean_pearson_correlation_dwd_netatmo = df_improvements.pearson_corr_dwd_netatmo.mean()

        #########################################################

        plt.ioff()
        fig = plt.figure(figsize=(32, 12), dpi=150)

        ax = fig.add_subplot(111)

        ax.plot(df_compare.index,
                df_compare.pearson_corr_dwd_,
                alpha=.8,
                c='b',  # colors_arr,
                marker='d',
                label='DWD Interpolation %0.2f'
                % mean_pearson_correlation_dwd_only)
        
        ax.plot(df_compare.index,
                df_compare.pearson_corr_netatmo_,
                alpha=.8,
                c='k',  # colors_arr,
                marker='2',
                label='Netatmo Interpolation %0.2f'
                % mean_pearson_correlation_netatmo_only)
        
        ax.plot(df_compare.index,
                df_compare.pearson_corr_dwd_netatmo,
                alpha=.8,
                c='r',  # colors_arr,
                marker='*',
                label='DWD-Netatmo Interpolation %0.2f'
                % mean_pearson_correlation_dwd_netatmo)

        ax.set_title('Pearson Correlation Interpolated Rainfall'
                     ' from DWD or Netatmo or DWD-Netatmo\n '
                     
                     'Precipitation of %s Intense Events \n'
                     'Stations with Improvemnts with OK Netatmo no Filter %d / %d, Percentage %0.0f'
                     '\nStations with Improvemnts with OK DWD-Netatmo %d / %d,'
                     ' Percentage %0.0f\n'
                     
                     %(temp_freq, stations_with_improvements_netatmo,
                        df_compare.pearson_corr_netatmo_.shape[0],
                        percent_of_improvment_netatmo,
                         stations_with_improvements,
                         df_compare.pearson_corr_dwd_netatmo.shape[0],
                        percent_of_improvment ))
        import matplotlib.dates as mdates
        plt.setp(ax.get_xticklabels(), rotation=90)
        ax.grid(alpha=0.25)
        ax.set_yticks(np.arange(0, 1.05, .10))
        years_fmt = mdates.DateFormatter('%Y-%m-%d %H')
        months = mdates.MonthLocator() 
        ax.xaxis.set_major_locator(months)

        ax.xaxis.set_major_formatter(years_fmt)
        # ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H'))

        #ax.set_yticks(np.arange(0., 1.05, .10))
        ax.set_xlabel('Date of Event')
        ax.legend(loc='lower right')
        ax.legend(loc='lower left')
        ax.set_ylabel('Pearson Correlation')

        plt.savefig((main_dir / r'ppt_cross_valid_RH_60min' / (
             r'ppt_spatial_pears_corr_%s_events_%s.png'
                                  % (temp_freq, used_data_acc ))),
                    papertype='a4', bbox_inches='tight', pad_inches=.2)
        plt.close()

#===============================================================================
#  
#  
#         stations_with_improvements = sum(i >= j for (i, j) in zip(
#             df_improvements.spearman_corr_dwd_netatmo.values,
#             df_improvements.spearman_corr_dwd_.values))
#  
#         stations_without_improvements = sum(i < j for (i, j) in zip(
#             df_improvements.spearman_corr_dwd_netatmo.values,
#             df_improvements.spearman_corr_dwd_.values))
#  
#         percent_of_improvment = 100 * (
#             stations_with_improvements /
#            df_improvements.spearman_corr_dwd_netatmo.shape[0])
#  
#         id_stns_no_impv = [id_stn for (id_stn, i, j) in zip(
#             df_improvements.index,
#             df_improvements.spearman_corr_dwd_netatmo.values,
#             df_improvements.spearman_corr_dwd_.values) if i < j]
#  
#     #     # OK with netatmo
#         stations_with_improvements_netatmo= sum(i >= j for (i, j) in zip(
#             df_improvements.spearman_corr_netatmo_.values,
#             df_improvements.spearman_corr_dwd_.values))
#         stations_without_improvements_netatmo = sum(i < j for (i, j) in zip(
#             df_improvements.spearman_corr_netatmo_.values,
#             df_improvements.spearman_corr_dwd_.values))
#         percent_of_improvment_netatmo= 100 * (
#             stations_with_improvements_netatmo /
#             df_improvements.spearman_corr_dwd_netatmo.shape[0])
#  
#         ####
#         mean_spearman_correlation_dwd_only = df_improvements.spearman_corr_dwd_.mean()
#         mean_spearman_correlation_netatmo_only = df_improvements.spearman_corr_netatmo_.mean()
#  
#         mean_spearman_correlation_dwd_netatmo = df_improvements.spearman_corr_dwd_netatmo.mean()
#  
#         #########################################################
#  
#         plt.ioff()
#         fig = plt.figure(figsize=(32, 12), dpi=150)
#  
#         ax = fig.add_subplot(111)
#  
#         ax.plot(df_improvements.index,
#                 df_improvements.spearman_corr_dwd_,
#                 alpha=.8,
#                 c='b',  # colors_arr,
#                 marker='d',
#                 label='DWD Interpolation %0.2f'
#                 % mean_spearman_correlation_dwd_only)
#          
#         ax.plot(df_improvements.index,
#                 df_improvements.spearman_corr_netatmo_,
#                 alpha=.8,
#                 c='k',  # colors_arr,
#                 marker='2',
#                 label='Netatmo Interpolation %0.2f'
#                 % mean_spearman_correlation_netatmo_only)
#          
#         ax.plot(df_improvements.index,
#                 df_improvements.spearman_corr_dwd_netatmo,
#                 alpha=.8,
#                 c='r',  # colors_arr,
#                 marker='*',
#                 label='DWD-Netatmo Interpolation %0.2f'
#                 % mean_spearman_correlation_dwd_netatmo)
#  
#         ax.set_title('Spearman Correlation Interpolated Rainfall from DWD or Netatmo or DWD-Netatmo\n '
#                       
#                      'Precipitation of %s Intense Events \n'
#                      'Stations with Improvemnts with OK Netatmo no Filter %d / %d, Percentage %0.0f'
#                      '\nStations with Improvemnts with OK DWD-Netatmo %d / %d,'
#                      ' Percentage %0.0f\n'
#                       
#                      % (temp_freq, stations_with_improvements_netatmo,
#                         df_improvements.spearman_corr_netatmo_.shape[0],
#                         percent_of_improvment_netatmo,
#                          stations_with_improvements,
#                          df_improvements.spearman_corr_dwd_netatmo.shape[0],
#                         percent_of_improvment ))
#  
#         plt.setp(ax.get_xticklabels(), rotation=90)
#         ax.grid(alpha=0.25)
#         ax.set_yticks(np.arange(0, 1.05, .10))
#         ax.set_xlabel('DWD Stations')
#         ax.legend(loc='lower left')
#         ax.set_ylabel('Spearman Correlation')
#  
#         plt.savefig((main_dir / r'ppt_cross_valid_RH_60min' / (
#              r'ppt_temporal_spr_corr_%s_events_%s.png'
#                                   % (temp_freq, used_data_acc ))),
#                     papertype='a4', bbox_inches='tight', pad_inches=.2)
#         plt.close()
# #===============================================================================
#===============================================================================
