# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
GOAL: To what areal and temporal extent are extreme precipitation
     events simultaneously occurring in NetAtmo Hourly Data

Get NetAtmo precipitation data (Hourly freq)

Select one station with an extreme event ( ppt > 24mm / hour)

Look how often does it happen that extreme values are measured
First look at exactly the same time (0min)

Second shift the second station with time intervals of +-30min
(+-60min, +-120min)

Save the result in a dataframe for every extreme event and station
"""

__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# ============================================================================
import os
import timeit
import time

import pandas as pd
import numpy as np

from pathlib import Path
from datetime import timedelta

#==============================================================================
#
#==============================================================================
main_dir = Path(os.getcwd())
os.chdir(main_dir)

path_to_netatmo_data_df = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                           r'\NetAtmo_BW'
                           r'\ppt_all_netatmo_hourly_stns_combined_.csv')

assert os.path.exists(path_to_netatmo_data_df), 'wrong data_df path'

path_to_netatmo_coords_df = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                             r'\NetAtmo_BW'
                             r'\rain_bw_1hour'
                             r'\netatmo_bw_1hour_coords.csv')
assert os.path.exists(path_to_netatmo_coords_df), 'wrong coords_df path'

out_save_dir = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\thr24NetAtmo')

if not os.path.exists(out_save_dir):
    os.mkdir(out_save_dir)


ppt_thrs_stn1 = [24]  # 20
ppt_thrs_stn2 = 0

ppt_thr_max_val = 60


time_shifts = [timedelta(minutes=int(m)) for m in np.arange(60, 121, 60)]
time_shifts_arr_floats = [i for i in np.arange(-120, 121, 60)]


#==============================================================================
#
#==============================================================================

def find_simulataneous_Netatmo_events(netatmo_ppt_df_file,
                                      ppt_thrs_lst,
                                      ppt_thrs2,
                                      ppt_thr_max,
                                      time_shifts_lst,
                                      time_shifts_arr,
                                      out_dir):

    in_netatmo_stns_df = pd.read_csv(netatmo_ppt_df_file,
                                     index_col=0, sep=';',
                                     parse_dates=True,
                                     infer_datetime_format=True,
                                     engine='c')

    for thr in ppt_thrs_lst:
        print('PPT Threshold is: ', thr, ' mm')
        for ii, stn_id in enumerate(in_netatmo_stns_df.columns):

            print('First Station ID is: ', stn_id,
                  ' Station Index is ', ii)
            try:
                # read station one
                idf1 = in_netatmo_stns_df.loc[:, stn_id]
            except Exception as msg:
                print(msg)
                continue
            # select values above threshold and drop nans
            stn1_abv_thr = idf1[(thr <= idf1) & (idf1 < ppt_thr_max)]
            stn1_abv_thr.dropna(inplace=True)

            if len(stn1_abv_thr.values) > 0:

                # start going through events of station one
                for ix, val in zip(stn1_abv_thr.index, stn1_abv_thr.values):
                    print('First time index is:', ix,
                          'Ppt Station 1 is ', val)

                    ix_time = ix.isoformat().replace(':',
                                                     '_').replace('T', '_')

                    # remove the ID of the first station from all Ids

                    stns_lst_without_stn1 = in_netatmo_stns_df.columns[
                        np.where(in_netatmo_stns_df.columns != stn_id)]

                    # create new dataframe to hold the ouput per event
                    df_result = pd.DataFrame(columns=stns_lst_without_stn1,
                                             index=time_shifts_arr)

                    # to know where the iteration is
                    count_all_stns = len(stns_lst_without_stn1)

                    for ii2, stn_id_2 in enumerate(df_result.columns):

                        print('Second Station ID is: ', stn_id_2,
                              ' index is ', ii2,
                              ' Count of Stns is :', count_all_stns)
                        # read second station
                        try:
                            idf2 = in_netatmo_stns_df.loc[:, stn_id_2]
                        except Exception as msg:
                            print(msg)
                            continue

                        # select values above threshold and drop nans
                        stn2_abv_thr = idf2[(ppt_thrs2 <= idf2) &
                                            (idf2 < ppt_thr_max)]
                        stn2_abv_thr.dropna(inplace=True)

                        if len(stn2_abv_thr.values) > 0:
                            # check if at the time event at stn2

                            if ix in stn2_abv_thr.index:
                                print('Same time in Station 2', ix)

                                val2 = np.round(stn2_abv_thr.loc[ix], 2)

                                print(' Ppt at Station 2 is', val2)
                                df_result.loc[0, stn_id_2] = val2

                            for time_shift in time_shifts_lst:
                                print('Shifting time by +- ', time_shift)

                                # get the shift as float, for index in df
                                shift_minutes = (
                                    time_shift / 60).total_seconds()

                                ix2_pos = ix + time_shift  # pos shifted
                                ix2_neg = ix - time_shift  # neg shifted

                                if ix2_pos in stn2_abv_thr.index:
                                    print('+ shifted idx present', ix2_pos)

                                    val2_pos = stn2_abv_thr.loc[ix2_pos]

                                    df_result.loc[
                                        shift_minutes,
                                        stn_id_2] = np.round(val2_pos, 2)

                                if ix2_neg in stn2_abv_thr.index:
                                    print('- shifted idx present', ix2_neg)

                                    shift_minutes_neg = - shift_minutes

                                    val2_neg = stn2_abv_thr.loc[
                                        ix2_neg]
                                    df_result.loc[
                                        shift_minutes_neg,
                                        stn_id_2] = np.round(val2_neg, 2)

                        del (idf2, stn2_abv_thr)
                        count_all_stns -= 1
                    df_result.dropna(axis=1, how='all', inplace=True)
                    # save df for every event
                    if df_result.values[0].shape[0] > 0:

                        print('Saving dataframe')
                        df_result.to_csv(
                            os.path.join(
                                out_dir,
                                'ppt_val_%0.2f_date_%s_stn_%s_thr_%s.csv'
                                % (val, ix_time, stn_id, thr)),
                            float_format='%0.2f')
                        del df_result
                    else:
                        print('df is empty will be deleted')
                        del df_result
            else:
                print('Station %s, has no data above % 0.1f mm' %
                      (stn_id, thr))
                continue

    return


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    find_simulataneous_Netatmo_events(netatmo_ppt_df_file=path_to_netatmo_data_df,
                                      ppt_thrs_lst=ppt_thrs_stn1,
                                      ppt_thrs2=ppt_thrs_stn2,
                                      ppt_thr_max=ppt_thr_max_val,
                                      time_shifts_lst=time_shifts,
                                      time_shifts_arr=time_shifts_arr_floats,
                                      out_dir=out_save_dir)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
