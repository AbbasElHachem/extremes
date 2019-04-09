#-------------------------------------------------------------------------
# Name:        get_data
# Purpose:     read the precipitation data from the hdf5 file
#
# Author:      mueller
#
# Created:     07.08.2013
# Copyright:   (c) mueller 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------

import os
import numpy as np
import tables
import matplotlib.pyplot as plt
import scipy
import datetime
import pandas as pd

import scipy.spatial.distance as dist

# ################### import modules as  #################################
#
# modulepath=r'X:\projects\2012Niedsim\04_Data\01_Code\00_General_Functions'
# sys.path.append(modulepath)
#
# or for shira
#basepath = os.path.dirname(os.path.dirname(os.getcwd()))
# ##############################################################################

"""
Class: Creates HDF5 object with open file and several methods
--------------------------------------------------------------------------------

Parameters:
----------series
    HDF5.plot_timeseries
    HDF5.get_dates_isoformat


    series      : string, optional, default = 'cut'
                    'all': maximum timeseries according to hdf5 is returned
                    'cut': timeseries is cut according to the very first and the
                           very last entry of all given stations. NaN values at
                           beginning/end are discarded.

    ids         : 'string' or [list of strings]
    idx         :  integer or [list of integers]

Keyword Arguments:
------------------
The following keywords can be used with:
    HDF5.get_data
    HDF5.get_pandas_time
    start       : datetime object
                    Start of the defined timeseries

    end         : datetime object
                    End of the defined timeseries

Definitions:
------------
    ids         : station ids
    idxs        : station indices in hdf5 file
    start_idx   : first index of all timeseries
    end_idx     : last index of all timeseries
"""


class HDF5(object):
    def __init__(self,
                 infile_def=None,
                 infile=None,
                 basepath=r'P:\2017_SYNOPSE_II',
                 rwmode='r',
                 datapath='all_stations'):
        """Creates an object containing the hdf5 file
        ---------
        Parameters
            infile_def: str
                default hdf5 data file
                choose time aggregation ['daily', '5' or '60']
                + version ['_op', '_re']
                e.g.: 'daily_op'

            infile: str, optional
                File and pathname of input hdf5. Default are the hdf5 files
                containing all data depending on the given agg_in

            basepath: str, optional
                Basepath of Niedsim Folder, default 'X:\\projects\\2012Niedsim'
                (ONLY for Windows systems)

        """

        # default paths and filenames of hdf5
        if 'all_stations' in datapath:
            print('all stations')
            datapath = os.path.join(basepath, '03_Data')
            if infile is None:
                agg_dic = {'daily': r'1781_2016_daily.h5',
                           '60': r'1949_2017_hourly.h5',
                           '5': os.path.join(r'04_merge_data',
                                             '1993_2016_5min_merge_nan.h5')}
                self.infile = os.path.join(datapath,
                                           agg_dic[infile_def])
            else:
                self.infile = infile
        else:
            print('without ref stations')
            datapath = os.path.join(
                basepath, '03_Data', '07_without_ref_station')
            if infile is None:
                agg_dic = {'daily': r'1781_2016_daily.h5',
                           '60': r'1949_2017_hourly.h5',
                           '5': r'1993_2016_5min_merge_nan.h5'}
                self.infile = os.path.join(datapath,
                                           agg_dic[infile_def])
            else:
                self.infile = infile

        # open file
        self.f = tables.open_file(self.infile, rwmode)

        # find aggregation of input hdf5
        # datetime object (seconds)
        iagg_dt = (pd.to_datetime(self.f.root.timestamps.isoformat[1].astype(str))
                   - pd.to_datetime(self.f.root.timestamps.isoformat[0].astype(str)))
        if iagg_dt.days == 1:
            self.agg_in = 1440
        else:
            #integer (minutes)
            if pd.__version__ == '0.15.2':
                self.agg_in = int(iagg_dt.to_pytimedelta().seconds / 60)
            else:
                self.agg_in = int(iagg_dt.seconds / 60)

    def close(self):
        """ close hdf5 file"""
        self.f.close()

    def check_idx_id(self, ids, idxs):
        """Check if input is idxs or ids and return the idxs """
        if idxs is None:
            if ids is None:
                raise ValueError('Neither ids nor idxs are specified')
            else:
                idxs = self.get_idxs_by_ids(ids)
        return idxs

    def get_all_ids(self):
        """ get all ids in the hdf5 file"""
        ids = self.f.root.id[:].astype(str)
        return ids

    def get_all_states(self):
        """ get all ids in the hdf5 file"""
        states = self.f.root.state_s[:].astype(str)
        return states

    def get_all_z(self):
        """get all z in the hdf5 file"""
        z_all = self.f.root.z[:].astype(str)
        return(z_all)

    def get_all_names(self):
        """get all names in the hdf5 file"""
        names = self.f.root.name[:].astype(str)
        return(names)

    def get_provider(self):
        """ get providers for station data"""
        provider = self.f.root.provider[:].astype(str)
        return(provider)

    def get_ids_by_coords(self, x, y, r, gk):
        """ get station ids within a given radius from the center
        -----------
        Parameters:
            x: float
                x-coordinates of the center
            y: float
                y-coordinates of the center
            r: float
                radial distances from the center in meter
            gk: str
                gk-coordinate system ('gk2','gk3','gk4')
        """

        x_coords = self.f.get_node('/coord/', 'x{:s}'.format(gk[-1]))[:]
        y_coords = self.f.get_node('/coord/', 'y{:s}'.format(gk[-1]))[:]
        coords = np.array([x_coords, y_coords]).T
        distance = dist.cdist(np.array([[x, y]]), coords, 'euclidean').ravel()

        return self.f.root.id[distance <= r].astype(str)

    def get_ids_by_idxs(self, idxs):
        """Get the ids of the stations for the given indices in the hdf5 file
        """
        # transform idxs into an array
        idxs = np.asanyarray(idxs)
        # get indices of stats in hdf5
        ids = self.f.root.id[idxs]
        return ids

    def get_idxs_by_ids(self, ids):
        """Get the indices in the hdf5 file for given ids"""

        # transform ids into an array
        ids = np.asanyarray(ids)
        # get indices of stats in hdf5
        idxs = np.where(np.in1d(self.f.root.id[:].astype(str), ids))[0]
        return idxs

    def get_idxs_by_time(self, start=None, end=None):
        """Get the indices of all stations in the hdf5 file with the given
        start and/or end indices
        Paramters: '1992-01-01'
        """
        if HDF5.agg_in == 'daily':
            freq = 'D'
        else:
            freq = '{:d}Min'.format(HDF5.agg_in)

        # get full timeseries for the given hdf5 file
        timeseries = pd.date_range(HDF5.f.root.timestamps.isoformat[0].astype(str),
                                   HDF5.f.root.timestamps.isoformat[-1].astype(
                                       str),
                                   freq=freq)

        # if start indix is given
        if start != None:
            start = pd.to_datetime(start)
            # find the correct start time index in the hdf5 file
            start_idx = np.where(timeseries == start)[0][0]
            # check which time series starts at or AFTER this time
            idx_1 = HDF5.f.root.timestamps.start_idx[...].astype(
                int) >= start_idx
            # idx is only returned if end is not given
            #(otherwise it is overwritten)
            idxs = idx_1

        # if end indix is given do same
        if end != None:
            end = pd.to_datetime(end)
            end_idx = np.where(timeseries == end)[0][0]
            # at or BEFORE this time
            idx_2 = HDF5.f.root.timestamps.end_idx[...].astype(int) <= end_idx
            idxs = idx_2

        if (start != None) & (end != None):
            idxs = idx_1 & idx_2

        if (start == None) & (end == None):
            AssertionError('Missing input parameter!')
        return idxs

    def get_start_end_idx(self,
                          ids=None, idxs=None,
                          series='cut',
                          start=None, end=None):
        """Get the start_idx and end_idx of time series for given ids"""

        if (start != None) & (end != None):
            if self.agg_in == 1440:
                freq = '1D'
            else:
                freq = '{:d}Min'.format(self.agg_in)

            dates = pd.date_range(self.f.root.timestamps.isoformat[0].astype(str),
                                  self.f.root.timestamps.isoformat[-1].astype(
                                      str),
                                  freq=freq)
            start_idx = np.where(dates == pd.to_datetime(start))[0]
            end_idx = np.where(dates == pd.to_datetime(end))[0]

        elif series == 'all':
            start_idx = 0
            end_idx = self.f.root.timestamps.isoformat.shape[0] - 1

        elif series == 'cut':
            idxs = self.check_idx_id(ids, idxs)
            start_idx = self.f.root.timestamps.start_idx[idxs].astype(int)
            end_idx = self.f.root.timestamps.end_idx[idxs].astype(int)

        else:
            raise Exception('Please check your time selection!')

        return start_idx, end_idx

    def get_min_start_max_end_idx(self, ids=None, idxs=None, **kwargs):
        """Get the minimum start_idx and max end_idx of time series
        for given ids
        """
        start_idx, end_idx = self.get_start_end_idx(ids=ids,
                                                    idxs=idxs,
                                                    **kwargs)

        # get start and end indices of stats, in order to be able to store
        # all data in one array, the length has to be the maximum of all
        # timeseries
        return np.min(start_idx), np.max(end_idx)

    def get_data(self, ids=None, idxs=None, **kwargs):
        """Get the precipitation data for given ids or idxs
        ----------
        Keyword Arguments: series, start, end

        """

        idxs = self.check_idx_id(ids, idxs)

        start_idx, end_idx = self.get_min_start_max_end_idx(
            idxs=idxs, **kwargs)

        data_org = (self.f.root.data[start_idx:end_idx + 1, idxs])

        return data_org

    def get_metadata(self, ids=None, idxs=None):
        """get metadata for given ids"""
        idxs = self.check_idx_id(ids, idxs)
        nodes_meta = ['z', 'id', 'provider',
                      'state_s', 'vers', 'name', 'comment']
        metadata = dict(zip(
            nodes_meta, [self.f.get_node('/', value)[idxs].astype(str)
                         for value in nodes_meta]))
        nodes_meta1 = ['start_idx', 'end_idx']
        metadata1 = dict(zip(
            nodes_meta1, [self.f.get_node('/timestamps', value)[idxs].astype(int)
                          for value in nodes_meta1]))
        start_time = [self.f.root.timestamps.isoformat[item].astype(str)
                      for item in metadata1['start_idx']]
        end_time = [self.f.root.timestamps.isoformat[item].astype(str)
                    for item in metadata1['end_idx']]
        metadata2 = {'start_time': start_time,
                     'end_time': end_time}
        metadata.update(metadata1)
        metadata.update(metadata2)

        return metadata

    def get_coordinates(self, ids=None, idxs=None):
        """get coordinates for given ids"""
        idxs = self.check_idx_id(ids, idxs)
        ##nodes_coord = ['x2','y2','x3','y3','x4','y4']
        nodes_coord = [a.name for a in self.f.get_node(
            '/coord/')._f_walknodes()]
        coordinates = dict(zip(nodes_coord, [self.f.get_node('/coord/', value)[idxs]
                                             for value in nodes_coord]))
        return coordinates

    def get_ids_of_closest_stations(self, x, y, gk, nstats=1, sort_dist=False):
        """Get the ids of the nstats closest stations
        ----------
        Parameters:
            x: int
                gk-coordinate in x direction
            y: int
                gk-coordinate in x direction
            gk: str
                coordinate system, e.g. 'gk3'
            nstats: int, default = 1
                number of closest stations
            sort_dist: bool, default = False
                if True, sort the values according to the closest distance
                if False, order according to the input metafile
        --------
        Return:
            ids, distances

        """
        all_ids = self.get_all_ids()
        coords = self.get_coordinates(all_ids)

        if gk == 'umt':
            raise Exception('STRANGE coordinate system!')
        if gk == 'utm':
            x_all = coords['x']
            y_all = coords['y']
        else:
            x_all = coords['x{:s}'.format(gk[-1])]
            y_all = coords['y{:s}'.format(gk[-1])]

        differences = ((x - x_all)**2. + (y - y_all)**2.)**0.5
        # create the postion number
        rank_diff = np.argsort(differences)
        # get the index that have a position number <= the number of
        # closest stations (excluding the position 0, which is the station
        # itself
        idx_close = rank_diff[:nstats]
        ids = all_ids[idx_close]
        distances = differences[idx_close]

        # Loop over gauges to get indices of the sorted gauges
        if not sort_dist:
            ids_sorted = self.get_metadata(ids)['id']
            idx_sorted = []
            for iid, id in enumerate(ids_sorted):
                idx = np.where(ids == id)[0][0]
                idx_sorted.append(idx)
            idx_sorted = np.array(idx_sorted)
            distances = distances[idx_sorted]
            ids = ids[idx_sorted]

        return ids, distances

    def get_dates_isoformat(self, ids=None, idxs=None,
                            series='cut', start=None, end=None):
        """Get dates for stations in isoformat
        ----------
        Keyword Arguments: series, start, end

        ---------
        Returns
        dates : list [dim = n] of list[dim = k]
                n: all stations
                k: for all dates
        """
        start_idx, end_idx = self.get_start_end_idx(ids=ids,
                                                    idxs=idxs,
                                                    series=series,
                                                    start=start,
                                                    end=end)

        if series == 'all':
            dates = [
                self.f.root.timestamps.isoformat[start_idx:end_idx + 1].astype(str)]
        else:
            dates = []
            for ii in range(start_idx.shape[0]):
                dates.append(
                    self.f.root.timestamps.isoformat[start_idx[ii]:end_idx[ii] + 1].astype(str))

        return dates

    def get_pandas_dataframe(self, ids=None, idxs=None, **kwargs):
        """create pandas dataframe
        ----------
        Keyword Arguments: series, start, end

        """

        data = self.get_data(ids=ids, idxs=idxs, **kwargs)
        dates = self.get_dates(ids=ids, idxs=idxs, **kwargs)

        # if station ids are not unique, colums of pandas object will
        # have same id. But as ids are unique, the id of the hdf5 have to be
        # used
        hdf5_ids = self.get_metadata(ids=ids, idxs=idxs)['id']
        if isinstance(hdf5_ids, str):
            hdf5_ids = [hdf5_ids]
        df = pd.DataFrame(data, index=dates, columns=hdf5_ids)

        return df

    def get_dates(self, ids=None, idxs=None, **kwargs):

        if self.agg_in == 1440:
            freq = '1D'
        else:
            freq = '{:d}Min'.format(self.agg_in)

        start_idx, end_idx = self.get_min_start_max_end_idx(ids=ids,
                                                            idxs=idxs,
                                                            **kwargs)

        ##start = times.iso2datetime(self.f.root.timestamps.isoformat[start_idx])
        ##end = times.iso2datetime(self.f.root.timestamps.isoformat[end_idx])
        start = pd.to_datetime(
            self.f.root.timestamps.isoformat[start_idx].astype(str))
        end = pd.to_datetime(
            self.f.root.timestamps.isoformat[end_idx].astype(str))
        dates = pd.date_range(start, end, freq=freq)

        return dates

    def plot_timeseries(self, ids=None, idxs=None, figname=None, **kwargs):
        """plot timeseries of given ids or given idxs
        ----------
        Keyword Arguments: series, start, end
        """

        # get data
        data = self.get_data(ids=ids, idxs=idxs, **kwargs)
        # get dates
        dates = self.get_dates(ids=ids, idxs=idxs, **kwargs)

        # get metadata
        metadata = self.get_metadata(ids=ids, idxs=idxs)

        # add a second axis, if data is 1D
        if np.ndim(data) == 1:
            data = data[:, np.newaxis]
        for ii in range(data.shape[1]):
            plt.plot(dates, data[:, ii], label=metadata['id'][ii])
        plt.legend()
        if figname == None:
            pass
        else:
            plt.savefig(figname)


if __name__ == '__main__':

    ##

    HDF52 = HDF5(
        infile=r'X:\exchange\ElHachem\niederschlag_deutschland\1993_2016_5min_merge_nan.h5')

    ids = HDF52.get_all_ids()
    metadata = HDF52.get_metadata(ids=ids)
    # data = HDF52.get_data(ids = ['5100', '2115'])

#     ids_bw = metadata['id'][metadata['state_s'] == 'BW']

#     df = HDF52.get_pandas_dataframe(ids=ids_bw)
    # '15Min', '120Min', '2H', '24H', 'D'
#     df_5min = df.resample('5Min', label='right', closed='right').sum()

#     for ii, iid in enumerate(ids):
    #         if ii < 10:
    #             continue
#         try:
#             idf = HDF52.get_pandas_dataframe(ids=[iid])
#             idf.to_pickle(
#                 r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\PPT_DATA_5min\%s.pkl' %
#                 str(iid))
#         except Exception as msg:
#             print(msg)
#             continue
#         idf.to_parquet(
#             r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\PPT_DATA_5min\%s.gzip' % str(
#                 iid), compression='gzip')

#         break
#         print(idf.head())
#         if ii == 15:
# break  # It terminates the current loop and resumes execution at the
# next statement, just like the traditional break statement in C
