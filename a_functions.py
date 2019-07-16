#-------------------------------------------------------------------------
# Name:        create_hdf5
# Purpose:     creates hdf5 file for the precipitation data
#
# Author:      mueller
#
# Created:     12.06.2013
# Copyright:   (c) mueller 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------
import tables
import numpy as np
import math
# import mpl_toolkits.basemap.pyproj as pyproj
import pyproj as pyproj
import datetime
import pandas as pd
import os
import shutil
from matplotlib import pyplot as plt
from scipy.spatial import distance as dist
import b_get_data as gd
plt.ioff()
# ################### import modules as  #################################
#
# modulepath=r'X:\projects\2012Niedsim\04_Data\01_Code\00_General_Functions'
# sys.path.append(modulepath)
#
# ##############################################################################


def create_hdf5(hf5_file,
                start_dt,
                end_dt,
                nstats,
                freq,
                data_title,
                in_ppt_df,
                in_stns_df,
                utm=False):
    """Creates empty hdf5 file
    -----------
    Parameters:

        hf5_file: string
            filepath and filename of hdf5 file
        start_dt: daytime object (year,month,day)
            First day of the timeseries
        end_dt: daytime object (year,month,day)
            Last date of the timeseries
        nstats: int
            Number of stations
        freq: str
            Pandas aggregation definition:
            M':Monthly 'D': Daily, 'H': Hourly, 'Min': Minutely
        utm: boolean
            if True, create entry for utm coordinates (instead of GK)

    """
    dates = pd.date_range(start_dt, end_dt, freq=freq)
    # number of maximum timesteps
    nts_max = dates.shape[0]

    # , filters=tables.Filters(complevel=6))
    hf = tables.open_file(hf5_file, 'w')

    # timestamps
    hf.create_group(where=hf.root,
                    name='timestamps',
                    title='Timestamps of respective Aggregation as Start Time')
    hf.create_carray(where=hf.root.timestamps,
                     name='isoformat',
                     atom=tables.StringAtom(19),
                     shape=(nts_max,),
                     chunkshape=(10000,),
                     title='Strings of Timestamps in Isoformat')
    hf.create_carray(where=hf.root.timestamps,
                     name='year',
                     atom=tables.IntAtom(),
                     shape=(nts_max,),
                     chunkshape=(10000,),
                     title='Yearly Timestamps')
    hf.create_carray(where=hf.root.timestamps,
                     name='month',
                     atom=tables.IntAtom(),
                     shape=(nts_max,),
                     chunkshape=(10000,),
                     title='Monthly Timestamps')
    hf.create_carray(where=hf.root.timestamps,
                     name='day',
                     atom=tables.IntAtom(),
                     shape=(nts_max,),
                     chunkshape=(10000,),
                     title='Daily Timestamps')
    hf.create_carray(where=hf.root.timestamps,
                     name='yday',
                     atom=tables.IntAtom(),
                     shape=(nts_max,),
                     chunkshape=(10000,),
                     title='Yearday Timestamps')
    hf.create_carray(where=hf.root.timestamps,
                     name='start_idx',
                     atom=tables.IntAtom(),
                     shape=(nstats,),
                     title='First Index of the Timeseries')
    hf.create_carray(where=hf.root.timestamps,
                     name='end_idx',
                     atom=tables.IntAtom(),
                     shape=(nstats,),
                     title='Last Index of the Timeseries')

    # data
    hf.create_carray(where=hf.root,
                     name='data',
                     atom=tables.FloatAtom(dflt=np.nan),
                     shape=(nts_max, nstats),
                     chunkshape=(10000, 1),
                     title=data_title)

    # data flags
    hf.create_carray(where=hf.root,
                     name='data_flags',
                     atom=tables.StringAtom(50),
                     shape=(nts_max, nstats),
                     chunkshape=(10000, 1),
                     title='Information about noticeable values')

    # coordinates
    if not utm:
        hf.create_group(where=hf.root,
                        name='coord',
                        title='Coordinates of Stations with different Gauss-Krueger'
                              + ' References')
        hf.create_carray(where=hf.root.coord,
                         name='x',
                         atom=tables.IntAtom(),
                         shape=(nstats,),
                         title='X-Coordinate (Gauss-Krueger 2)')
        hf.create_carray(where=hf.root.coord,
                         name='y',
                         atom=tables.IntAtom(),
                         shape=(nstats,),
                         title='Y-Coordinate (Gauss-Krueger 2)')
# hf.create_carray(where=hf.root.coord,
# name='x3',
# atom=tables.IntAtom(),
# shape=(nstats,),
# title='X-Coordinate (Gauss-Krueger 3)')
# hf.create_carray(where=hf.root.coord,
# name='y3',
# atom=tables.IntAtom(),
# shape=(nstats,),
# title='Y-Coordinate (Gauss-Krueger 3)')
# hf.create_carray(where=hf.root.coord,
# name='x4',
# atom=tables.IntAtom(),
# shape=(nstats,),
# title='X-Coordinate (Gauss-Krueger 4)')
# hf.create_carray(where=hf.root.coord,
# name='y4',
# atom=tables.IntAtom(),
# shape=(nstats,),
# title='Y-Coordinate (Gauss-Krueger 4)')
    else:
        hf.create_group(where=hf.root,
                        name='coord',
                        title='UTM Coordinates of Stations')
        hf.create_carray(where=hf.root.coord,
                         name='x',
                         atom=tables.IntAtom(),
                         shape=(nstats,),
                         title='X-Coordinate')
        hf.create_carray(where=hf.root.coord,
                         name='y',
                         atom=tables.IntAtom(),
                         shape=(nstats,),
                         title='Y-Coordinate')

    # metadata
    hf.create_carray(where=hf.root,
                     name='name',
                     atom=tables.StringAtom(50),
                     shape=(nstats,),
                     title='Name of Station')
    hf.create_carray(where=hf.root,
                     name='state_s',
                     atom=tables.StringAtom(200),
                     shape=(nstats,),
                     title='Abbreviation of State, where Station is situated in')
    hf.create_carray(where=hf.root,
                     name='provider',
                     atom=tables.StringAtom(50),
                     shape=(nstats,),
                     title='Provider/Operator of the Station')
    hf.create_carray(where=hf.root,
                     name='id',
                     atom=tables.StringAtom(10),
                     shape=(nstats,),
                     title='Identification Number of Station')
    hf.create_carray(where=hf.root,
                     name='z',
                     atom=tables.IntAtom(),
                     shape=(nstats,),
                     title='Elevation')
    hf.create_carray(where=hf.root,
                     name='vers',
                     atom=tables.StringAtom(50),
                     shape=(nstats,),
                     title='NiedSim Versions, where the Station is used for')
    hf.create_carray(where=hf.root,
                     name='comment',
                     atom=tables.StringAtom(200),
                     shape=(nstats,),
                     title='Individual comments for the Station')

    # convert timestamp to isoformat
    ts_iso = []
    ts_year = []
    ts_month = []
    ts_day = []
    ts_yday = []

    for ii in range(dates.shape[0]):
        ts = dates[ii]
        ts_iso.append(ts.isoformat())
        # get timestamp years
        ts_year.append(ts.year)
        # get timestamp months
        ts_month.append(ts.month)
        # get timestamp days
        ts_day.append(ts.day)
        # get timestamp year days
        ts_yday.append(ts.timetuple().tm_yday)

    print('create_hdf5 done')

    # added by Abbas
    ppt_data_all = np.array(in_ppt_df.values)

    stns_int = []
    for stn_str in in_ppt_df.columns:
        stns_int.append(int(stn_str))

    stns_int_arr = np.array(stns_int)
    stns_coords_lon = in_stns_df.loc[stns_int_arr, 'geoLaenge'].values
    stns_coords_lat = in_stns_df.loc[stns_int_arr, 'geoBreite'].values
    stns_coords_z = in_stns_df.loc[stns_int_arr, 'Stationshoehe'].values
    stns_name = in_stns_df.loc[stns_int_arr, 'Stationsname'].values
    stns_bundesland = in_stns_df.loc[stns_int_arr, 'Bundesland'].values

    start_idx_values, last_idx_values = [], []
    # TODO: MAKE ME FASTER
    for stn_id in in_ppt_df.columns:
        try:
            stn_vals = in_ppt_df.loc[:, stn_id].values.ravel()
            stn_idx = in_ppt_df.loc[:, stn_id].index.values

            idx_and_val_not_nan = next((ix, x)
                                       for ix, x in zip(stn_idx, stn_vals)
                                       if not math.isnan(x))
            first_time_ix = idx_and_val_not_nan[0]  # get time index of start
            # get arg of idx and append to list
            first_time_ix_ix = np.where(stn_idx == first_time_ix)[0][0]
            start_idx_values.append(first_time_ix_ix)

            idx_and_val_is_nan = next((ix, x)
                                      for ix, x in zip(stn_idx[::-1], stn_vals[::-1])
                                      if not math.isnan(x))
            last_time_ix = idx_and_val_is_nan[0]
            last_time_ix_ix = np.where(stn_idx == last_time_ix)[0][0]
            last_idx_values.append(last_time_ix_ix)
        except Exception as msg:
            print(msg)
            start_idx_values.append(np.nan)
            last_idx_values.append(np.nan)
    stns_name_utf8, stn_bundesland_utf8 = [], []
    for stn_name, stn_land in zip(stns_name, stns_bundesland):
        stns_name_utf8.append(str(stn_name).encode('utf-8'))
        stn_bundesland_utf8.append(str(stn_land).encode('utf-8'))

    # fill hf5 with predefined stamps
    hf.root.timestamps.isoformat[:] = ts_iso[:]
    hf.root.timestamps.year[:] = ts_year[:]
    hf.root.timestamps.month[:] = ts_month[:]
    hf.root.timestamps.day[:] = ts_day[:]
    hf.root.timestamps.yday[:] = ts_yday[:]
    hf.root.timestamps.start_idx[:] = start_idx_values[:]
    hf.root.timestamps.end_idx[:] = last_idx_values[:]

    hf.root.data[:] = ppt_data_all[:]

    hf.root.coord.x[:] = stns_coords_lon[:]
    hf.root.coord.y[:] = stns_coords_lat[:]

    hf.root.id[:] = in_ppt_df.columns.values[:]
    hf.root.name[:] = stns_name_utf8[:]
    hf.root.z[:] = stns_coords_z[:]
    hf.root.state_s[:] = stn_bundesland_utf8[:]

    hf.close()


def copy_metadata(hdf5_out, hdf5_in):
    """ Copy metadata and coordinates from hdf5_in to hdf5_out

    Parameters:
    ----------
        hdf5_out: string
            filepath and name of hdf5file in which the metadata is copied
        hdf5_in: string
            filepath and name of hdf5file from which the metadata is copied

    Remark:
    --------
        Both file must exist and have to be of the same structure!!

    """
    hf_in = tables.open_file(hdf5_in, 'r')
    hf_out = tables.open_file(hdf5_out, 'r+')

    # write metadata
    for node in ['z', 'id', 'provider', 'state_s', 'vers', 'name', 'comment']:
        metadata = hf_in.get_node('/', node)[...]
        hf_out.get_node('/', node)[...] = metadata
        print(node, ' written')
    # write coord
    for node in ['x', 'y']:
        metadata = hf_in.get_node('/coord', node)[...]
        hf_out.get_node('/coord', node)[...] = metadata

    hf_in.close()
    hf_out.close()


def copy_metadata_from_dwd(hdf5, hdf5_old,
                           meta_path_dwd=(r'\\lhg-000\abt3\2012Niedsim\03_Import' +
                                          '\dwd_stations.csv'),
                           id_tag='STAT'):
    """ Copy metadata and coordinates from official dwd list
    Parameters:
    ----------
        hdf5: string
            filepath and name of hdf5file in which the metadata is copied
        hdf5_old: string
            filepath and name of hdf5file from which the metadata is copied
        meta_path_dwd: str, optional
            Path and filename of DWD Metafile, default:
            'X:\projects\2012Niedsim\03_Import\dwd_stations.csv'
            (ONLY for Windows systems)
        id_tag: str, optional
            define wether STAT_ID oder STAT are used in hdf5 file
            default is STAT
    """

    mdwd = pd.read_csv(meta_path_dwd, sep=';\s*')

    shutil.copyfile(hdf5_old, hdf5)

    hf = tables.open_file(hdf5, 'r+')

    ids = hf.root.id[:]
    provider = hf.root.provider[:]
    vers = hf.root.provider[:]
    for ii, iid in enumerate(ids):
        # if iid != '891':continue
        # check if iid is in dwd file and if it is a rainfall station
        # biid is boolean index of the dwd file (True for the correct iid)
        ids_dwd = (mdwd[id_tag].values == int(iid))
        rain_dwd = (mdwd['KE'].values == 'RR')
        biid = ids_dwd & rain_dwd
        # if station available, get the metadata
        if biid.sum() == 1:
            print(iid, provider[ii], vers[ii], 'updated!')
            # get longitudes and latitudes in degrees and minutes
            lat_deg = mdwd['BR'][biid].values[0] / 100
            lat_min = mdwd['BR'][biid].values[0] % 100
            lon_deg = mdwd['LA'][biid].values[0] / 100
            lon_min = mdwd['LA'][biid].values[0] % 100
            # transform to decimal degrees
            lats = DMS_to_DD(lat_deg, lat_min, 0)
            lons = DMS_to_DD(lon_deg, lon_min, 0)
            # transform to gauss-krueger
            x2, y2 = lonlat_to_gk(lons, lats, 'gk2')
            x3, y3 = lonlat_to_gk(lons, lats, 'gk3')
            x4, y4 = lonlat_to_gk(lons, lats, 'gk4')

            # get other metadata
            z = mdwd['HS'].values[biid][0]
            name = mdwd['STAT_NAME'].values[biid][0]
            state = mdwd['BL'].values[biid][0]

            # change coordinates
            hf.root.coord.x2[ii] = int(x2)
            hf.root.coord.y2[ii] = int(y2)
            hf.root.coord.x3[ii] = int(x3)
            hf.root.coord.y3[ii] = int(y3)
            hf.root.coord.x4[ii] = int(x4)
            hf.root.coord.y4[ii] = int(y4)

            # change meta data
            hf.root.z[ii] = z
            hf.root.name[ii] = name
            hf.root.state_s[ii] = str.lower(state)
            hf.root.provider[ii] = 'DWD'

    hf.close()


def update_startidx_endidx(hf5_file):
    """ Updates the first and last idx in the given hdf5_file

    Parameters:
    ----------
        hf5_file: string
            filepath and name of hdf5file
    """
    # open file again
    hf = tables.open_file(hf5_file, 'r+')
    nstats = hf.root.id.shape[0]

    for istat in range(nstats):
        # get start idx and end idx of the entire timeseries
        #!= syear_idx/eyear_idx) if first or last year is not fully available
        data_all = hf.root.data[:, istat]
        # check if data is not nan (=>boolean array)
        val_data = ~np.isnan(data_all)
        # indices of start and end of timeperiod
        # indices of start and end of timeperiod
        try:
            start_idx = np.where(val_data)[0][0]
            end_idx = np.where(val_data)[0][-1]
        except:
            start_idx = np.nan
            end_idx = np.nan
            print('No Data in Station', hf.root.id[istat])
        hf.root.timestamps.start_idx[istat] = start_idx
        hf.root.timestamps.end_idx[istat] = end_idx
    hf.close()
    print('update start-end_idx finished')


def find_duplicate(HDF5, outpath, calc_data=True, plot_fig=False):
    """Find and plot all duplicate stations as well as some metainformation

    Parameters:
    ---------
    HDF5: object
        HDF5 object created in b_get_data
    outpath: string
        outpath where figures and metafiles will be saved
    calc_data: boolean, optional
        If True, differences in the dataset are calculated, default = True
    plot_fig: boolean, optional
        If True, figures of timeseries are created, default = False

    """

    figpath = os.path.join(outpath, 'figures')
    try:
        os.makedirs(figpath)
    except:
        pass

    f = open(os.path.join(outpath, 'metadata_duplicate.csv'), 'w')

    # write title of metadata textfile
    f.write('id;version;start;end;altitude;diff_alt;diff_sys;'
            + 'diff_all;diff_coord;IDX\n')

    a_id = HDF5.get_all_ids()
    # get stats without duplicates
    u_id = np.unique(a_id)

    print('number of invalid stats', (a_id.shape[0] - u_id.shape[0]))

    for ii, iid in enumerate(u_id):
        print('start', ii, iid)
        metadata = HDF5.get_metadata(iid)

        # amount of dublicates in the hdf5 with same id
        ndup = HDF5.get_metadata(iid)['id'].shape[0]

        # do nothing if station is already unique
        if ndup == 1:
            continue
        # get timeseries data of entire length (nescessary to compare different
        # timeseries)
        coords = HDF5.get_coordinates(iid)
        dates_cut = HDF5.get_dates_isoformat(iid)
        idx = HDF5.get_idxs_by_ids(iid)

        # ################### differences in the data #########################

        if calc_data == True:
            ndata = HDF5.f.root.data.shape[0]

            # in order to avoid memory errors with larger datasets,
            # these datasets are sliced
            # lists containing the differences of the slices
            n_diff_data_sys = []
            n_diff_data_all = []

            if ndata < 3000000:
                start = [0]
                end = [ndata]
            else:
                start = [0, 3000000]
                end = [3000000, ndata]
            for islice in range(len(start)):
                istart = start[islice]
                iend = end[islice]

                if ndup == 2:
                    data1 = HDF5.f.root.data[istart:iend, idx[0]]
                    data2 = HDF5.f.root.data[istart:iend, idx[1]]
                    diff_data = data1 - data2

# if ndup == 3:
##                    data1 = HDF5.f.root.data[istart:iend,idx[0]]
##                    data2 = HDF5.f.root.data[istart:iend,idx[1]]
##                    data3 = HDF5.f.root.data[istart:iend,idx[3]]
##                    diff_data1 = data1 - data2
##                    diff_data2 = data1 - data3
##                    diff_data3 = data2 - data3

                # calculate the number of differing days, which differ more than
                # +-1 -> small changes are ignored!
                # only systematic offsets are regarded
                n_diff_data_sys.append(np.sum((diff_data < -0.1) |
                                              (diff_data > .1)))
                # all differences are regarded
                n_diff_data_all.append(np.sum((diff_data < -0.01)
                                              | (diff_data > 0.01)))
            # sum over all slices
            n_diff_data_sys = str(sum(n_diff_data_sys))
            n_diff_data_all = str(sum(n_diff_data_all))

        else:
            n_diff_data_sys = '-'
            n_diff_data_all = '-'

        # ################ differences in the coordinates #####################

        diff_coord = dist.pdist(np.vstack((coords['x3'], coords['y3'])).T)

        # only for two stations (three can be observed by eye)
        diff_altitude = metadata['z'][0] - metadata['z'][1]

        # write other information to file
        for jj in range(ndup):
            f.write('{:10s};{:3s};{:19s};{:19s};'.format(
                iid,
                metadata['vers'][jj],
                dates_cut[jj][0],
                dates_cut[jj][-1]))
            f.write('{:4.0f};{:10.0f};{:10s};{:10s};'.format(
                metadata['z'][jj],
                diff_altitude,
                n_diff_data_sys,
                n_diff_data_all))
            if ndup == 2:
                f.write('{:10.0f};'.format(diff_coord[0]))
            if ndup == 3:
                f.write('{:10.0f};'.format(diff_coord[jj]))
            f.write('{:10.0f}\n'.format(idx[jj]))

        f.flush()

        # plot data
        if plot_fig == True:
            dates = HDF5.get_dates_isoformat(iid, series='all')
            fig = plt.figure(figsize=(20, 10))
            plt.plot(
                # pd.to_datetime(dates[0]),
                HDF5.f.root.data[:, idx[0]],
                label=metadata['provider'][0] + ' ' + metadata['name'][0])
            plt.plot(
                # pd.to_datetime(dates[0]),
                HDF5.f.root.data[:, idx[1]], '--',
                label=metadata['provider'][1] + ' ' + metadata['name'][1])
            plt.legend()
            plt.show()
            # fig.savefig(os.path.join(figpath,iid+'.png'))
            plt.close(fig)

    f.close()

    print('finding duplicates done')


def del_ids(HDF5, hdf5_out, id_delete, utm=False):
    """Delete all stations with the given id

    Parameters:
    ---------
    HDF5: object
        HDF5 object created in b_get_data
    hdf5_out: string
        filepath and name of new file
    id_delete: 1D-array of strings
        station id of stations to be deleted

    """
    # get ...

    #... number of new stations
    ids = HDF5.f.root.id[...]
    n_old = ids.shape[0]

    # ... the number of new ids. The ids to be deleted need not all to be
    # included in the existing hdf5.
    n_new = n_old - len(set(ids).intersection(set(id_delete)))

    print(n_old, ' reduced to ', n_new)
    #... HDF5 information
    start_dt, end_dt, freq = _read_HDF5_structure(HDF5)

    create_hdf5(hdf5_out, start_dt, end_dt, n_new, freq, utm=utm)

    hf_new = tables.open_file(hdf5_out, 'r+')

    # idx: index in the old hdf5
    # ii: index of the new hdf5
    ii = 0

    for idx in range(n_old):
        # if idx != 559: continue
        if idx % 100 == 0:
            print(idx)
        iid = HDF5.f.root.id[idx]
        # do nothing if idx is to be deleted
        if iid in id_delete:
            print('del', iid)
            continue

        # copy data
        data = HDF5.f.root.data[:, idx]
        vers = HDF5.f.root.vers[idx]
        iid = HDF5.f.root.id[idx]

        # fill data
        hf_new.root.data[:, ii] = data[:]
        data = 0

        # get and fill data flags
        # has to be done in loops due to memory errors loading whole data_flags
        # ########################################
        n_dataflags = HDF5.f.root.data_flags.shape[0]
        max_data = 1000000
        niter = n_dataflags / max_data
        nrest = n_dataflags % max_data
        istart = 0
        iend = istart + max_data
        # add the flags stepwise with an array length of max_data
        for jj in range(niter):
            data_flags = HDF5.f.root.data_flags[istart:iend, idx]
            hf_new.root.data_flags[istart:iend, ii] = data_flags[:]
            istart += max_data
            iend += max_data
        # add the rest of the flags
        data_flags = HDF5.f.root.data_flags[istart:, idx]
        hf_new.root.data_flags[istart:, ii] = data_flags[:]
        # ########################################

        # fill coordinates
        if utm:
            coords = ['x', 'y']
        else:
            coords = ['x2', 'y2', 'x3', 'y3', 'x4', 'y4']
        for node in coords:
            hf_new.get_node(
                '/coord', node)[ii] = HDF5.f.get_node('/coord', node)[idx]

        # fill meta data
        for node in ['z', 'provider', 'state_s', 'name', 'comment']:
            hf_new.get_node('/', node)[ii] = HDF5.f.get_node('/', node)[idx]
        hf_new.get_node('/', 'id')[ii] = iid
        hf_new.get_node('/', 'vers')[ii] = vers

        ii += 1

    hf_new.close()

    update_startidx_endidx(hdf5_out)

    print('cleaning data done')


def clean_stations(HDF5,
                   hdf5_out,
                   idx_delete=np.array([]),
                   idx_survive=np.array([]),
                   idx_merge=np.array([]),
                   idx_change=np.array([]),
                   nan_endtime=None,
                   nan_starttime=None,
                   utm=False):
    """Make the station ids unique using different methods
       Merging, deleting or renaming

    Parameters:
    ---------
    HDF5: object
        HDF5 object created in b_get_data
    hdf5_out: string
        filepath and name of new file
    idx_delete: 1D-array of strings
        station idx of stations to be deleted
    idx_survive: 1D-array of strings
        station idx of stations where all data and metadata is used
    idx_merge: 1D-array of strings
        station idx of stations where only data is used for which the
        corresponding values in idx_survive are nan
    idx_change: 1D-array of strings
        station idx for which id will be changed to id_1 and both stations are
        kept
    nan_endtime: string or list of strings (in datetime isoformat), optional
        All data before the given date (date including) is set to nan for
        all merged stations in order to keep the data of the surviving station.

        Either one date for all surviving stations can be given,
        or a list of dates for each surviving station individually
        (idx_survive and dates need to be of the same structure AND order)

        =>> dates               || idx_survive
            2012-12-31T23:00:00     325
            2010-12-31T23:00:00     110
                   .                 .
                   .                 .
                   .                 .
    nan_starttime: string or list of strings (in datetime isoformat), optional
        Same as nan_endtime, except all data AFTER the given date is set to nan
        in the merged product

    """
    # get ...

    #... number of new stations
    n_old = HDF5.f.root.id.shape[0]
    n_merg = idx_merge.shape[0]
    n_del = idx_delete.shape[0]

    nstats = n_old - n_del - n_merg

    #... HDF5 information
    start_dt, end_dt, freq = _read_HDF5_structure(HDF5)

    create_hdf5(hdf5_out, start_dt, end_dt, nstats, freq, utm=utm)

    hf_new = tables.open_file(hdf5_out, 'r+')

    # idx: index in the old hdf5
    # ii: index of the new hdf5
    ii = 0

    for idx in range(n_old):
        # if idx != 559: continue
        if idx % 100 == 0:
            print(idx)

        # do nothing if idx is to be deleted
        if idx in idx_delete:
            continue

        # if idx will be merged, do nothing (as the merging will be done in for
        # the idx of the surviving station)
        elif idx in idx_merge:
            continue

        elif idx in idx_change:
            # copy data and change id to id_1
            data = HDF5.f.root.data[:, idx]
            vers = HDF5.f.root.vers[idx]
            iid = HDF5.f.root.id[idx] + '_1'

        elif idx in idx_survive:
            ##raise EnvironmentError('PLEASE CHECK HOW THE data_flags can be merged ->> ask Thomas ')
            # get data of surviving idx
            data = HDF5.f.root.data[:, idx]
            vers = HDF5.f.root.vers[idx]
            iid = HDF5.f.root.id[idx]

            # find all corresponding idxs for CURRENT id
            iidx_all = HDF5.get_idxs_by_ids(HDF5.get_ids_by_idxs(idx))
            # find idx to be merged for CURRENT id
            # iidx_all and iidx_merg  != idx_merge which are ALL idx to be
            # merged
            iidx_merg = iidx_all[~np.in1d(iidx_all, idx)]

            # print a warning if idx merge is empty

            # not is used, to find if array is empty
            # any() has to be used to avoid ambiguous errors for arrays with more
            # than one entry
            a = not any(iidx_merg)
            # however, if one idx is 0, not([0]) is False, although 0 is a valid
            # id -> additionally, not(iidx_all)==0 is applied
            b = not any(iidx_merg) == 0
            if a & b:
                id_missing = HDF5.get_ids_by_idxs(idx)
                print('Warning, for {:s} the merging idx is missing '.format(
                    id_missing))
            else:
                for kidx in iidx_merg:
                    # find the corresponding data to be merged
                    data_merg = HDF5.f.root.data[:, kidx]

                    if (nan_endtime != None):
                        # set all data of the given date and earlier to NAN
                        # for the merged data
                        if type(nan_endtime) == str:
                            inan_endtime = nan_endtime
                        else:
                            # if nan_endtime is a list,
                            # get the endtime corrresponging to the idx_survive
                            inan_endtime = nan_endtime[idx_survive == idx]

                        dates = HDF5.f.root.timestamps.isoformat[...]
                        nan_eidx = np.where((dates == inan_endtime))[0][0]
                        data_merg[:nan_eidx + 1] = np.nan

                    if (nan_starttime != None):
                        # set all data of the given date and later to NAN
                        # for the merged data
                        if type(nan_starttime) == str:
                            inan_starttime = nan_starttime
                        else:
                            # if nan_endtime is a list,
                            # get the endtime corrresponging to the idx_survive
                            inan_starttime = nan_starttime[idx_survive == idx]
                        dates = HDF5.f.root.timestamps.isoformat[...]
                        nan_sidx = np.where((dates == inan_starttime))[0][0]
                        data_merg[nan_sidx:] = np.nan

                    # find indices that are nan and which are to be replace by
                    # the new data
                    isnan = np.isnan(data)
                    # fill nan values of the surviving data with merged data
                    data[isnan] = data_merg[isnan]
                    vers += '_' + HDF5.f.root.vers[kidx]

        else:
            # copy data
            data = HDF5.f.root.data[:, idx]
            vers = HDF5.f.root.vers[idx]
            iid = HDF5.f.root.id[idx]

        # fill data
        hf_new.root.data[:, ii] = data[:]
        data = 0

        # get and fill data flags
        # has to be done in loops due to memory errors loading whole data_flags
        # ########################################
        n_dataflags = HDF5.f.root.data_flags.shape[0]
        max_data = 1000000
        niter = n_dataflags / max_data
        nrest = n_dataflags % max_data
        istart = 0
        iend = istart + max_data
        # add the flags stepwise with an array length of max_data
        for jj in range(niter):
            data_flags = HDF5.f.root.data_flags[istart:iend, idx]
            hf_new.root.data_flags[istart:iend, ii] = data_flags[:]
            istart += max_data
            iend += max_data
        # add the rest of the flags
        data_flags = HDF5.f.root.data_flags[istart:, idx]
        hf_new.root.data_flags[istart:, ii] = data_flags[:]
        # ########################################

        # fill coordinates
        if utm:
            coords = ['x', 'y']
        else:
            coords = ['x2', 'y2', 'x3', 'y3', 'x4', 'y4']
        for node in coords:
            hf_new.getNode(
                '/coord', node)[ii] = HDF5.f.getNode('/coord', node)[idx]

        # fill meta data
        for node in ['z', 'provider', 'state_s', 'name', 'comment']:
            hf_new.getNode('/', node)[ii] = HDF5.f.getNode('/', node)[idx]
        hf_new.getNode('/', 'id')[ii] = iid
        hf_new.getNode('/', 'vers')[ii] = vers

        ii += 1

    hf_new.close()

    update_startidx_endidx(hdf5_out)

    print('cleaning data done')


def merge_hdf5(hdf5_out, hdf5_files):
    """
    Merges several hdf5 files to one hdf5 file -> similiar to np.append()

    Parameters:
    ---------
    hdf5_out: string
        filepath and name of new file
    hdf5_files: list of objects created in b_get_data
        hdf5 objects that are to be merged [HDF5_1, HDF5_2,...]

    Remarks:
    --------
    All HDF5 objects need to be of the same structure!!


    """
    nfiles = len(hdf5_files)

    # slices is a list with the numbers of stations of the different hdf5 files
    # this is used to identify  the first_stat and last_stat of the
    # different hdf5 files (station index of the hdf5_out changes compared to
    # the input hdf5)
    slices = []
    ndata = []

    for ifile in range(nfiles):
        slices.append(hdf5_files[ifile].f.root.id.shape[0])
        # number of timesteps of each hdf5file
        ndata.append(hdf5_files[ifile].f.root.data.shape[0])

    nstats = sum(slices)

    #... HDF5 information (use the hdf5 with the longest timeseries
    # for structure information of the new hdf5)
    #stat_max_data = np.where(ndata == np.max(ndata))[0][0]
    start_dts = []
    end_dts = []
    for ifile in hdf5_files:
        start_dt, end_dt, freq = _read_HDF5_structure(ifile)
        start_dts.append(start_dt)
        end_dts.append(end_dt)

    create_hdf5(hdf5_out, np.min(start_dts), np.max(end_dts), nstats, freq)

    hf = tables.open_file(hdf5_out, 'r+')

    first_stat = 0

    for ifile in range(nfiles):
        last_stat = first_stat + slices[ifile]
        # calculate the start and end index in the new hdf5 file
        # nescessary if the input hdf5 have different length in the timeseries
        # the output hdf5 will get the maximum length
        start_iso = hdf5_files[ifile].f.root.timestamps.isoformat[0]
        start_idx = np.where(
            hf.root.timestamps.isoformat[:] == start_iso)[0][0]
        end_idx = start_idx + ndata[ifile]

        # write data
        for ii in range(slices[ifile]):
            data = hdf5_files[ifile].f.root.data[:, ii]
            hf.getNode('/', 'data')[start_idx:end_idx, first_stat + ii] = data
            if ii % 10 == 0:
                print(ii, 'written')
        # write metadata
        for node in ['z', 'id', 'provider', 'state_s', 'vers', 'name', 'comment']:
            metadata = hdf5_files[ifile].f.getNode('/', node)[...]
            hf.getNode('/', node)[first_stat:last_stat] = metadata
            print(node, ' written')
        # write coord
        for node in ['x2', 'y2', 'x3', 'y3', 'x4', 'y4']:
            metadata = hdf5_files[ifile].f.getNode('/coord', node)[...]
            hf.getNode('/coord', node)[first_stat:last_stat] = metadata
        # change
        first_stat = last_stat

    hf.close()

    update_startidx_endidx(hdf5_out)

    print('merging files done')


def create_metadata(hf5_file, outfile, utm=False):
    """ Creates a metadata file

    Parameters:
    ----------
        hf5_file: string
            filepath and name of hdf5file
        outfile: string
            filepath and name of metafile

    """

    hf = tables.open_file(hf5_file)
    f_out = open(outfile, 'w')

    nstats = hf.root.id[:].shape[0]

    f_out.write('Rechtswert;Hochwert;Hoehe;Stationsnr.;Name;')
    f_out.write('Bundesland;Betreiber;Version;Kommentar;')
    f_out.write('Anfang;Ende;max._Werte;Fehlwerte;Prozent\n')

    for ii in range(nstats):

        data = hf.root.data[:, ii]
        # check if data is not nan (=>boolean array)
        val_data = ~np.isnan(data)
        # if no data is available
        if val_data.sum() == 0:
            n_val = 0
            n_inval_tp = '-999'
            frac = -999
            start_iso = '-999'
            end_iso = '-999'
        else:
            # indices of start and end of timeperiod
            start_idx = hf.root.timestamps.start_idx[ii]
            end_idx = hf.root.timestamps.end_idx[ii]
            start_iso = hf.root.timestamps.isoformat[start_idx]
            end_iso = hf.root.timestamps.isoformat[end_idx]
            # number of values in time period
            n_val = end_idx - start_idx + 1
            # data of timeperiod
            data_tp = data[start_idx:end_idx + 1]
            # invalid data within timeperiod
            n_inval_tp = np.sum(np.isnan(data_tp))
            # fraction of invalid data
            frac = n_inval_tp / np.float(n_val) * 100

        name = hf.root.name[ii]
        if name == '':
            name = '-'

        if not utm:
            f_out.write('%s;' % hf.root.coord.x3[ii])
            f_out.write('%s;' % hf.root.coord.y3[ii])
        else:
            f_out.write('%s;' % hf.root.coord.x[ii])
            f_out.write('%s;' % hf.root.coord.y[ii])
        f_out.write('%s;' % hf.root.z[ii])
        f_out.write('%s;' % hf.root.id[ii])
        f_out.write('%s;' % name)
        f_out.write('%s;' % hf.root.state_s[ii])
        f_out.write('%s;' % hf.root.provider[ii])
        f_out.write('%s;' % hf.root.vers[ii])
        f_out.write('%s;' % hf.root.comment[ii])
        f_out.write('%s;' % start_iso)
        f_out.write('%s;' % end_iso)
        f_out.write('%s;' % n_val)
        f_out.write('%s;' % n_inval_tp)
        f_out.write('%6.3f\n' % frac)

        if ii % 100 == 0:
            print(ii)

    f_out.close()
    hf.close()

    print('create metadata finished')


def _read_HDF5_structure(HDF5):
    # get...

    #...start time
    start_dt = pd.to_datetime(HDF5.f.root.timestamps.isoformat[0])
    #...end time
    end_dt = pd.to_datetime(HDF5.f.root.timestamps.isoformat[-1])

    #...freq
    # datetime object (seconds)
    iagg_dt = (pd.to_datetime(HDF5.f.root.timestamps.isoformat[1])
               - pd.to_datetime(HDF5.f.root.timestamps.isoformat[0]))
    if iagg_dt.days == 1:
        freq = 'D'
    elif iagg_dt.days > 1:
        raise ValueError('NOT implemented yet!!')
    else:
        #integer (minutes)
        # iagg_dt.minutes, falls verfuegbar
        freq = '{:d}Min'.format(iagg_dt.seconds / 60)

    return start_dt, end_dt, freq


def gk_to_lonlat(gk_x, gk_y, gk_choose):

    # Convert Gauss-Krueger coordinates to longitudes and latitudes
    # Output geographic coordinates: WGS84
    wgs84 = pyproj.Proj(init='epsg:4326')

    # Input kartesian coordinates: Gauss-Krueger
    # (looked up: http://www.spatialreference.org/ref/epsg/31466/)
    if gk_choose == 'gk2':
        # wgs84 bounds: 5.87 ,49.10, 7.50, 53.75
        gk = pyproj.Proj(init="epsg:31466")
    elif gk_choose == 'gk3':
        # wgs84 bounds: 7.50, 47.27, 10.50, 55.06
        gk = pyproj.Proj(init="epsg:31467")
    elif gk_choose == 'gk4':
        # wgs84 bounds: 10.50, 47.27, 13.50, 55.06
        gk = pyproj.Proj(init="epsg:31468")

    # Coordinate system transformation for the different Gauss-Krueger "fields"
    lons, lats = pyproj.transform(gk, wgs84, gk_x, gk_y)

    return lons, lats


def lonlat_to_gk(lons, lats, gk_choose):

    # Convert longitudes and latitudes to Gauss-Krueger coordinates
    # Input geographic coordinates: WGS84
    wgs84 = pyproj.Proj(init='epsg:4326')

    # Output kartesian coordinates: Gauss-Krueger
    # (looked up: http://www.spatialreference.org/ref/epsg/31466/)
    if gk_choose == 'gk2':
        # wgs84 bounds: 5.87 ,49.10, 7.50, 53.75
        gk = pyproj.Proj(init="epsg:31466")
    elif gk_choose == 'gk3':
        # wgs84 bounds: 7.50, 47.27, 10.50, 55.06
        gk = pyproj.Proj(init="epsg:31467")
    elif gk_choose == 'gk4':
        # wgs84 bounds: 10.50, 47.27, 13.50, 55.06
        gk = pyproj.Proj(init="epsg:31468")

    # Coordinate system transformation for the different Gauss-Krueger "fields"
    gk_x, gk_y = pyproj.transform(wgs84, gk, lons, lats)

    return gk_x, gk_y


def utm32n_to_gk(utm_x, utm_y, gk_choose):

    # Convert utm coordiantes to Gauss-Krueger coordinates
    # Input utm coordinates
    utm32n = pyproj.Proj(init="epsg:25832")

    # Output kartesian coordinates: Gauss-Krueger
    # (looked up: http://www.spatialreference.org/ref/epsg/23032/)
    if gk_choose == 'gk2':
        # wgs84 bounds: 5.87 ,49.10, 7.50, 53.75
        gk = pyproj.Proj(init="epsg:31466")
    elif gk_choose == 'gk3':
        # wgs84 bounds: 7.50, 47.27, 10.50, 55.06
        gk = pyproj.Proj(init="epsg:31467")
    elif gk_choose == 'gk4':
        # wgs84 bounds: 10.50, 47.27, 13.50, 55.06
        gk = pyproj.Proj(init="epsg:31468")

    # Coordinate system transformation for the different Gauss-Krueger "fields"
    gk_x, gk_y = pyproj.transform(utm32n, gk, utm_x, utm_y)

    return gk_x, gk_y


def utm32n_to_lon_lat(utm_x, utm_y):

    # Convert utm coordiantes to geographic coordinates
    # Input utm coordinates
    # (looked up: http://www.spatialreference.org/ref/epsg/25832/)
    utm32n = pyproj.Proj(init="epsg:25832")

    # Output kartesian coordinates: geographic coordinates
    wgs84 = pyproj.Proj(init='epsg:4326')

    # Coordinate system transformation for the different Gauss-Krueger "fields"
    lon, lat = pyproj.transform(utm32n, wgs84, utm_x, utm_y)

    return lon, lat


def lon_lat_to_utm32n(lon, lat):

    # Convert geographic coordinate into utm coordiantes
    # Input
    wgs84 = pyproj.Proj(init='epsg:4326')

    # Output utm32
    # (looked up: http://www.spatialreference.org/ref/epsg/25832/)
    utm32n = pyproj.Proj(init="epsg:25832")

    utm_x, utm_y = pyproj.transform(wgs84, utm32n, lon, lat)

    return utm_x, utm_y


def DMS_to_DD(degree, minutes, seconds):
    # Convert Degrees, Minutes, and Seconds to Decimal Degrees
    decimal_degrees = degree + (minutes / 60.) + (seconds / 3600.)

    return decimal_degrees


def get_datetimeseries(syear, eyear, td):
    """ Creates a timeseries for start year, end year in the given aggregation

    Parameters:
    ----------
        syear: int
            start year
        eyear: int
            end year
        td: int
            timedelta in minutes

    Returns:
    --------
        start_dt: datetime object
            first entry of the timeseries
        end_dt:
            last entry of the timeseries
        dates: DatetimeIndex
            timeseries
    """
    start_dt = datetime.datetime(syear, 1, 1)
    end_dt = datetime.datetime(eyear + 1, 1, 1) - \
        datetime.timedelta(0, int(td) * 60)
    freq = '{:d}Min'.format(td)
    dates = pd.date_range(start_dt, end_dt, freq=freq)

    return start_dt, end_dt, dates


if __name__ == '__main__':
    ##
    b = (2658500, 5539500)
    print(b)
    a = gk_to_lonlat(b[0], b[1], 'gk2')
    print(lonlat_to_gk(a[0], a[1], 'gk2'))
    print(lon_lat_to_utm32n(a[0], a[1]))
