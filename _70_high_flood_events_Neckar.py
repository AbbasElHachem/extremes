'''
Created on 23 Jan 2023

@author: hachem
'''

import pandas as pd
import glob
import os

from pathlib import Path

if __name__ == '__main__':

    path_events_bw = ()

    main_path = Path(
        r'X:\exchange\ElHachem\High_flood_events')
    os.chdir(main_path)

    path_to_events_single = glob.glob(
        '*.csv')

    dfs = []
    for _df in path_to_events_single:
        in_df = pd.read_csv(_df, index_col=0,
                            sep=';', parse_dates=True,
                            infer_datetime_format=True)
        dfs.append(in_df)
        # break
    df_comb = pd.concat(dfs, axis=0)
    unique_idf = df_comb.index.drop_duplicates()
    df_final = pd.DataFrame(index=unique_idf).sort_index()
    df_final['Idx'] = range(1, len(unique_idf) + 1)
    df_final.to_csv('comb_events.csv', sep=';')

    # path_events_bw = Path(
    # r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\kriging_ppt_netatmo")

    events_hourly = path_events_bw / r'dwd_edf_transf_vg_60min.csv'
    events_3hours = path_events_bw / r'dwd_edf_transf_vg_180min.csv'
    events_6hours = path_events_bw / r'dwd_edf_transf_vg_360min.csv'
    events_12hours = path_events_bw / r'dwd_edf_transf_vg_720min.csv'
    events_24hours = path_events_bw / r'dwd_edf_transf_vg_1440min.csv'

    dfs_bw = [events_hourly, events_3hours, events_6hours,
              events_12hours, events_24hours]
    dfs_bw_comb = []
    for _df in dfs_bw:
        in_df = pd.read_csv(_df, index_col=0,
                            sep=';', parse_dates=True,
                            infer_datetime_format=True)
        dfs_bw_comb.append(in_df)

    df_comb_bw = pd.concat(dfs_bw_comb, axis=0)
    unique_idf_bw = df_comb_bw.index.drop_duplicates()
    df_final_bw = pd.DataFrame(index=unique_idf_bw).sort_index()
    df_final_bw['Idx'] = range(1, len(unique_idf_bw) + 1)
    df_final_bw_daily = df_final_bw.resample('D').mean().dropna()
    df_final_bw_daily.to_csv('comb_events_BW.csv', sep=';')

    cmn_events = df_final_bw_daily.index.intersection(
        df_final.index)
    df_final_bw_nk = pd.DataFrame(index=cmn_events).sort_index()
    df_final_bw_nk['Idx'] = range(1, len(cmn_events) + 1)
    df_final_bw_nk.to_csv('final_events.csv', sep=';')

    cmn_events.to_csv('.csv', sep=';')
    pass
