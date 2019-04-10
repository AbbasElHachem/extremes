# !/usr/bin/env python.
# -*- coding: utf-8 -*-

"""
GOAL: To what areal and temporal extent are extreme precipitation
     events simultaneously occurring

Get high resolution precipitation data (5min freq)

Select 2 stations, work pairwise

Look how often does it happen that extreme values are measured
First look at exactly the same time (5min)

Second shift the second station with time intervals of +-5min
(+-5min, +-10min, +-15min,..., +-30min, ..., +-60min)

Save the result in a dataframe for every extreme event and station
"""

__author__ = "Abbas El Hachem"
__copyright__ = 'Institut fuer Wasser- und Umweltsystemmodellierung - IWS'
__email__ = "abbas.el-hachem@iws.uni-stuttgart.de"

# ============================================================================

from pathlib import Path
from datetime import timedelta

from b_get_data import HDF5

import os
import timeit
import time

import pandas as pd
import numpy as np


main_dir = Path(os.getcwd())
os.chdir(main_dir)

path_to_ppt_hdf_data = (r'X:\exchange\ElHachem'
                        r'\niederschlag_deutschland'
                        r'\1993_2016_5min_merge_nan.h5')

ppt_thrs = [4, 6, 8, 12, 14, 16, 20]  # 0.1, 1, 2, 4, 6, 8, 10 12, 14,

time_shifts = [timedelta(minutes=5), timedelta(minutes=10),
               timedelta(minutes=15), timedelta(minutes=20),
               timedelta(minutes=25), timedelta(minutes=30),
               timedelta(minutes=35), timedelta(minutes=40),
               timedelta(minutes=45), timedelta(minutes=50),
               timedelta(minutes=55), timedelta(minutes=60)]
time_shifts_arr_floats = [i for i in np.arange(-60, 61, 5)]


def find_simulataneous_events(ppt_thrs_lst, stns_ids_lst,
                              time_shifts_lst, time_shifts_arr):

    for thr in ppt_thrs_lst:
        print('PPT Threshold is: ', thr, ' mm')
        for ii, iid in enumerate(stns_ids_lst):
            print('First Station ID is: ', iid, ' Station Index is ', ii)
            try:
                # read station one
                idf1 = HDF52.get_pandas_dataframe(ids=[iid])
            except Exception as msg:
                print(msg)
                continue
            # select values above threshold and drop nans
            stn1_abv_thr = idf1[idf1 >= thr]
            stn1_abv_thr.dropna(inplace=True)

            if len(stn1_abv_thr.values) > 0:

                # start going through events of station one
                for ix, val in zip(stn1_abv_thr.index, stn1_abv_thr.values):
                    print('First time index is:', ix, 'Ppt Station 1 is ', val)

                    # remove the id of the first station from all IDS
                    ids2 = np.array(list(filter(lambda x: x != iid, ids)))

                    # create new dataframe to hold the ouput per event
                    df_result = pd.DataFrame(columns=ids2,
                                             index=time_shifts_arr)

                    # to know where the iteration is
                    count_all_stns = len(ids2)

                    for ii2, iid2 in enumerate(ids2):
                        if iid != iid2:
                            print('Second Station ID is: ', iid2,
                                  ' index is ', ii2,
                                  ' Count of Stns is :', count_all_stns)
                            # read second station
                            try:
                                idf2 = HDF52.get_pandas_dataframe(ids=[iid2])
                            except Exception as msg:
                                print(msg)
                                continue

                            # select values above threshold and drop nans
                            stn2_abv_thr = idf2[idf2 >= thr]
                            stn2_abv_thr.dropna(inplace=True)

                            if len(stn2_abv_thr.values) > 0:
                                # check if at the time event at stn2
                                if ix in stn2_abv_thr.index:
                                    print('Same time in Station 2', ix)

                                    val2 = stn2_abv_thr.loc[ix, :].values[0]

                                    print(' Ppt at Station 2 is', val2)
                                    df_result.loc[0, iid2] = val2

                                for time_shift in time_shifts_lst:
                                    print('Shifting time index by +- ',
                                          time_shift)

                                    # get the shift as float, for index in df
                                    shift_minutes = (
                                        time_shift / 60).total_seconds()

                                    ix2_pos = ix + time_shift  # pos shifted
                                    ix2_neg = ix - time_shift  # neg shifted

                                    if ix2_pos in stn2_abv_thr.index:
                                        print('+ shifted idx present', ix2_pos)

                                        val2_pos = stn2_abv_thr.loc[
                                            ix2_pos, :].values[0]

                                        df_result.loc[shift_minutes,
                                                      iid2] = val2_pos

                                    if ix2_neg in stn2_abv_thr.index:
                                        print('- shifted idx present', ix2_neg)

                                        shift_minutes_neg = - shift_minutes

                                        val2_neg = stn2_abv_thr.loc[
                                            ix2_neg, :].values[0]
                                        df_result.loc[shift_minutes_neg,
                                                      iid2] = val2_neg

                            del (idf2, stn2_abv_thr)
                        count_all_stns -= 1
                    df_result.dropna(axis=1, how='all', inplace=True)
                    df_result.to_csv(
                        'ppt_val_%0.2f_date_%s_stn_%s_thr_%s.csv'
                        % (val,
                           ix.isoformat().replace(':', '_').replace('T', '_'),
                           iid, thr))
                    # break

            else:
                print('Station %s, has no data above % 0.1f mm' % (iid, thr))
                continue
            break
        # break


if __name__ == '__main__':

    print('**** Started on %s ****\n' % time.asctime())
    START = timeit.default_timer()  # to get the runtime of the program

    HDF52 = HDF5(infile=path_to_ppt_hdf_data)
    ids = HDF52.get_all_ids()
    metadata = HDF52.get_metadata(ids=ids)

    find_simulataneous_events(ppt_thrs_lst=ppt_thrs,
                              stns_ids_lst=ids,
                              time_shifts_lst=time_shifts,
                              time_shifts_arr=time_shifts_arr_floats)

    STOP = timeit.default_timer()  # Ending time
    print(('\n****Done with everything on %s.\nTotal run time was'
           ' about %0.4f seconds ***' % (time.asctime(), STOP - START)))
