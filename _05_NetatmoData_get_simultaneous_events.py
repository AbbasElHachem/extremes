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


from _00_additional_functions import list_all_full_path
#==============================================================================
#
#==============================================================================
main_dir = Path(os.getcwd())
os.chdir(main_dir)

path_to_netatmo_data = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                        r'\NetAtmo_BW'
                        r'\rain_bw_1hour')
print(path_to_netatmo_data)
assert os.path.exists(path_to_netatmo_data), 'wrong data_df path'

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


def get_stns_ids_coords_file(netatmo_coords_df_file):
    ''' get all ids from coords_file'''
    in_coords_df = pd.read_csv(netatmo_coords_df_file, index_col=0, sep=',')
    ids_lst = [stn_id.replace(':', '_') for stn_id in in_coords_df.index]
    ids_lst_unique = list(set(ids_lst))
    return ids_lst_unique
# @profile

#==============================================================================
#
#==============================================================================


def split_df_file_to_get_alls_stn_ids(all_df_files_list):
    ''' get stn_id from netatmo file'''
    stns_ids = []
    for df_file in all_df_files_list:
        try:
            _, _, p3 = df_file.split('(')
            p3 = p3.replace(')).csv', '')
        except Exception as msg:
            print(msg)
        stns_ids.append(p3)
    stns_ids_unique = list(set(stns_ids))
    return stns_ids_unique

#==============================================================================
#
#==============================================================================


def split_one_df_file_to_get_stn_id(df_file):
    ''' get stn_id from netatmo file'''
    try:
        _, _, p3 = df_file.split('(')
        p3 = p3.replace(')).csv', '')
    except Exception as msg:
        print(msg)
    return p3

#==============================================================================
#
#==============================================================================


def ids_list_intersection(lst1, lst2):
    ''' intersect two lists'''
    # Use of hybrid method
    temp = set(lst2)
    lst3 = [value for value in lst1 if value in temp]
    return lst3

#==============================================================================
#
#==============================================================================


def remove_one_elt_lst(elt, lst):
    ''' remove a specific element from a list'''
    new_lst = []
    for el in lst:
        if el != elt:
            new_lst.append(el)
    return new_lst
#==============================================================================
#
#==============================================================================


def find_simulataneous_Netatmo_events(all_dfs_files,
                                      ppt_thrs_lst,
                                      ppt_thrs2,
                                      ppt_thr_max,
                                      stns_ids_lst,
                                      time_shifts_lst,
                                      time_shifts_arr,
                                      out_dir):

    for thr in ppt_thrs_lst:
        print('PPT Threshold is: ', thr, ' mm')
        for ii, df_file in enumerate(all_dfs_files):
            stn_id = split_one_df_file_to_get_stn_id(df_file)
            if stn_id in stns_ids_lst:
                print('stn in coords df', df_file)
                print('First Station ID is: ', stn_id,
                      ' Station Index is ', ii)
                try:
                    # read station one
                    idf1 = pd.read_csv(df_file, sep=';', index_col=0,
                                       parse_dates=True, engine='c')
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

                        # remove the file of the first station from all dfs
                        all_dfs_files_left = remove_one_elt_lst(
                            df_file, all_dfs_files)

                        # remove the id of the first station from all stn ids
                        stns_ids_lst_2 = remove_one_elt_lst(
                            stn_id, stns_ids_lst)

                        # create new dataframe to hold the ouput per event
                        df_result = pd.DataFrame(columns=stns_ids_lst_2,
                                                 index=time_shifts_arr)

                        # to know where the iteration is
                        count_all_stns = len(stns_ids_lst_2)

                        for ii2, df_file_2 in enumerate(all_dfs_files_left):
                            iid2 = split_one_df_file_to_get_stn_id(df_file_2)
                            print('Second Station ID is: ', iid2,
                                  ' index is ', ii2,
                                  ' Count of Stns is :', count_all_stns)
                            # read second station
                            try:
                                idf2 = pd.read_csv(df_file_2,
                                                   sep=';',
                                                   index_col=0,
                                                   parse_dates=True,
                                                   engine='c')
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

                                    val2 = np.round(
                                        stn2_abv_thr.loc[ix, :].values[0], 2)

                                    print(' Ppt at Station 2 is', val2)
                                    df_result.loc[0, iid2] = val2

                                for time_shift in time_shifts_lst:
                                    print('Shifting time by +- ', time_shift)

                                    # get the shift as float, for index in df
                                    shift_minutes = (
                                        time_shift / 60).total_seconds()

                                    ix2_pos = ix + time_shift  # pos shifted
                                    ix2_neg = ix - time_shift  # neg shifted

                                    if ix2_pos in stn2_abv_thr.index:
                                        print('+ shifted idx present', ix2_pos)

                                        val2_pos = stn2_abv_thr.loc[
                                            ix2_pos, :].values[0]

                                        df_result.loc[
                                            shift_minutes,
                                            iid2] = np.round(val2_pos, 2)

                                    if ix2_neg in stn2_abv_thr.index:
                                        print('- shifted idx present', ix2_neg)

                                        shift_minutes_neg = - shift_minutes

                                        val2_neg = stn2_abv_thr.loc[
                                            ix2_neg, :].values[0]
                                        df_result.loc[
                                            shift_minutes_neg,
                                            iid2] = np.round(val2_neg, 2)

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

    # get all Netatmo df_files and all station ids
    all_netatmo_df_files = list_all_full_path('.csv', path_to_netatmo_data)

    stn_ids_df_files = split_df_file_to_get_alls_stn_ids(all_netatmo_df_files)
    stn_ids_coords_file = get_stns_ids_coords_file(path_to_netatmo_coords_df)

    stn_ids_in_common = ids_list_intersection(stn_ids_df_files,
                                              stn_ids_coords_file)

    find_simulataneous_Netatmo_events(all_dfs_files=all_netatmo_df_files,
                                      ppt_thrs_lst=ppt_thrs_stn1,
                                      ppt_thrs2=ppt_thrs_stn2,
                                      ppt_thr_max=ppt_thr_max_val,
                                      stns_ids_lst=stn_ids_in_common,
                                      time_shifts_lst=time_shifts,
                                      time_shifts_arr=time_shifts_arr_floats,
                                      out_dir=out_save_dir)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
