#==============================================================================
# TODO: ADD DOCUMENATION ORGANIZE
#==============================================================================


import datetime
import os
import time
import timeit
import fnmatch
import tables
import pandas as pd
import requests
import glob
import zipfile
import numpy as np
import math

from _00_additional_functions import (list_all_full_path,
                                      resampleDf,
                                      select_df_within_period)

#==============================================================================
# 10minutenwerte_rr_00003_19930428_19991231_hist.zip
# stundenwerte_RR_00003_19950901_20110401_hist.zip
# tageswerte_RR_00001_19120101_19860630_hist.zip
#==============================================================================
download_data = True
delete_zip_files = True


build_one_df = True
delete_df_files = True

make_hdf5_dataset_new = False

temporal_freq = '1_minute'  # '1_minute' '10_minutes' 'hourly' 'daily'
time_period = 'historical'  # 'recent' 'historical'
temp_accr = 'stundenwerte'  # '1minuten' '10minuten' 'stundenwerte'  'tages'

ppt_act = 'nieder'  # nieder 'rr' 'RR' 'RR'

stn_name_len = 5  # length of dwd stns Id, ex: 00023

year_period = '2017'  # fpr 1_minute data

start_date = '%s-01-01 00:00:00' % year_period
end_date = '%s-12-30 23:00:00' % year_period

temp_freq = 'Min'  # '60Min'  # 'Min' '10Min' 'H' 'D'

# '%Y%m%d%H'  # when reading download dwd station df '%Y%m%d%H%M'
time_fmt = '%Y%m%d%H%M'  # for 10Min data use: '%Y%m%d%H%M%S'

# name of station id and rainfall column name in df
stn_id_col_name = 'STATIONS_ID'
ppt_col_name = 'RS_01'  # 'TT_TU'  # RS_01 'RWS_10' '  R1' 'RS'

freqs_list = ['60Min']  # '60Min']  # '1D' '5Min',


date_range = pd.date_range(start=start_date,
                           end=end_date, freq=temp_freq)

# years_ = [year for year in date_range.years]
#==============================================================================
#
#==============================================================================
# main_dir = os.path.join(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\download_DWD_data_recent')

main_dir = r'/run/media/abbas/EL Hachem 2019/download_DWD_data_recent/'
# r'F:\download_DWD_data_recent'
os.chdir(main_dir)

dwd_in_coords_df = pd.read_csv(
    #r"X:\staff\elhachem\Data\DWD_DE_Data\station_coordinates_names_hourly.csv",
    r"file:///run/media/abbas/EL Hachem 2019/home_office/NetAtmo_BW/station_coordinates_names_hourly_only_in_BW.csv",
    sep=',', index_col=0, encoding='latin-1')

# added by Abbas, for DWD stations
stndwd_ix = ['0' * (5 - len(str(stn_id))) + str(stn_id)
              if len(str(stn_id)) < 5 else str(stn_id)
              for stn_id in dwd_in_coords_df.index]
 
dwd_in_coords_df.index = stndwd_ix
dwd_in_coords_df.index = list(map(str, dwd_in_coords_df.index))
#==============================================================================
#
#==============================================================================


def firstNonNan(listfloats):
    '''return first non NaN value in python list'''
    idx_not_nan = np.argwhere(math.isnan(listfloats) == False)[0]
    return idx_not_nan


#==============================================================================
# # TODO:
#==============================================================================
# download station coordinates name (maybe do it seperat later)
# stn_names_files = os.path.join(main_dir, r'station_coordinates_names.csv')
# stn_names_files = os.path.join(
#     main_dir, r'station_coordinates_names_hourly_only_in_BW_utm32.csv')
#
# assert os.path.exists(stn_names_files)
# #     in_stn_coords_df_loc = os.path.join(
# #         r"F:\download_DWD_data_recent\station_coordinates_names_hourly_only_in_BW_utm32.csv")
#
# stn_names_df = pd.read_csv(
#     stn_names_files, index_col=0, sep=';', encoding='latin-1')
# stations_id = stn_names_df.index
#
# # get all station ids, make them a string for generating file_names
# stations_id_str_lst = []
# for stn_id in stations_id:
#     stn_id = str(stn_id)
#     if len(stn_id) < 5:
#         stn_id = '0' * (5 - len(stn_id)) + stn_id
#     stations_id_str_lst.append(stn_id)
# print(stations_id_str)

#==============================================================================
# DOWNLOAD ALL ZIP FILES
#==============================================================================

out_dir = os.path.join(main_dir, 'DWD_data_%s_%s_%s_%s'
                       % (temporal_freq, time_period, temp_accr, year_period))
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

os.chdir(out_dir)


# generate url
if download_data:  # download all zip files

    # precipitation
    base_url = (r'https://opendata.dwd.de/climate_environment/CDC'
                r'/observations_germany/climate/%s/precipitation/%s/%s/'
                % (temporal_freq, time_period, year_period))
    # r'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/1_minute/precipitation/historical/1993/'
    r = requests.get(base_url, stream=True)

    all_urls = r.content.split(b'\n')

#     if time_period == 'historical':
#         _zip_file = -3
#     if time_period == 'recent':
#         _zip_file = -2

    all_urls_zip = [url_zip.split(b'>')[-2].replace(b'</a', b'').decode('utf-8')
                    for url_zip in all_urls if b'.zip' in url_zip]

    if temporal_freq == '1_minute':
        all_urls_zip = [url_zip.split(b'>')[0].replace(
            b'<a', b'').replace(b' href=', b'').strip(b'"').decode('utf-8')
            for url_zip in all_urls if b'.zip' in url_zip]
    # station namea in BW
    stns_list = dwd_in_coords_df.index.to_list()
    
    # files in BW
    all_urls_zip_bw = [_file for stn in stns_list 
                       for _file in all_urls_zip
                       if str(stn) in _file]
    # len(all_urls_zip_bw)
    for file_name in all_urls_zip_bw:

        #         if time_period == 'historical':
        #             file_name = file_name.split('=')[1].replace('"', '')
        print('getting data for', file_name)

        file_url = base_url + file_name
        
        try:
            zip_url = requests.get(file_url)
        except Exception as msg:
            print(msg)
            continue
        local_filename = file_url.split('/')[-1]

        with open(local_filename, 'wb') as f:
            for chunk in zip_url.iter_content(chunk_size=1024):
                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)
        zip_url.close()

    r.close()
    print('\n####\n Done Getting all Data \n#####\n')
    #==========================================================================
    #
    #==========================================================================
    
    # get all zip files as list
    all_zip_files = list_all_full_path('.zip', out_dir)

    out_extract_df_dir = ('DWD_extracted_zip_%s_%s_%s_%s' % 
                          (temporal_freq, time_period, temp_accr, year_period))

    # extract all zip files, dfs in zip files
    for zip_file in all_zip_files:
        print('exctracting data from ', zip_file)
        with zipfile.ZipFile(zip_file, "r") as zip_ref:
            # Get a list of all archived file names from the zip
            listOfFileNames = zip_ref.namelist()

            for fileName in listOfFileNames:
                # Check filename endswith csv
                if 'produkt' in fileName and fileName.endswith('.txt'):
                    # Extract a single file from zip
                    zip_ref.extract(fileName, out_extract_df_dir)

    if delete_zip_files:
        print('deleting downloaded zip files')
        for zip_file in all_zip_files:
            os.remove(zip_file)
    #==========================================================================
    #
    #==========================================================================
if build_one_df:
    #     all_df_files = list_all_full_path('.txt', out_extract_df_dir)

    all_df_files = list_all_full_path('.txt', out_dir)

#     all_df_files_bw = [_df_f for stn_name in stations_id_str_lst
#                        for _df_f in all_df_files if stn_name in _df_f]
    # get all downloaded station ids, used as columns for df_final
    available_stns = []

    for df_txt_file in all_df_files:

        stn_name = df_txt_file.split('_')[-1].split('.')[0]
        print('getting station name from file ', stn_name)
        assert len(stn_name) == 5
        if stn_name not in available_stns:
            available_stns.append(stn_name)

    # for stations only in BW do this
#     available_stns = stations_id_str_lst

    # create daterange and df full of nans

    date_range = pd.date_range(start=start_date, end=end_date, freq=temp_freq)

    data = np.zeros(shape=(date_range.shape[0], len(available_stns)))
    data[data == 0] = np.nan

    final_df_combined = pd.DataFrame(data=data,
                                     index=date_range,
                                     columns=available_stns)
#     final_df_combined = pd.DataFrame(columns=available_stns)
    # read df stn file and fill final df for all stns
    all_files_len = len(all_df_files)

    for df_txt_file in all_df_files:
        print('\n++\n Total number of files is ++', all_files_len, '\n')
        
        stn_name = df_txt_file.split('_')[-1].split('.')[0]
#         if stn_name in stations_id_str_lst:
#             print('\n##\n Stations in BW \n##\n', stns_bw_len)
        print('\n##\n Station ID is \n##\n', stn_name)
        in_df = pd.read_csv(df_txt_file, sep=';', index_col=1, engine='c')

        assert int(stn_name) == in_df.loc[:, stn_id_col_name].values[0]
        assert stn_name in final_df_combined.columns, 'assertion error'

        if temporal_freq == '10_minutes':
            try:
                in_df.index = pd.to_datetime(in_df.index, format=time_fmt)
            except AttributeError as msg:
                print(msg)
                in_df.index = [ix.split(":")[0] for ix in in_df.index]
                in_df.index = pd.to_datetime(in_df.index, format=time_fmt)
                continue

        else:
            in_df.index = pd.to_datetime(in_df.index, format=time_fmt)

        in_df = select_df_within_period(in_df, start_date, end_date)
        in_df = in_df[in_df[ppt_col_name] >= -100]
        ppt_data = in_df[ppt_col_name].values.ravel()

        if ppt_data.shape[0] > 0:
            print('\n++\n  Data shape is \n++\n', ppt_data.shape)
#             try:
            cmn_idx = final_df_combined.index.intersection(
                in_df.index)
            try:
                final_df_combined.loc[cmn_idx, stn_name] = ppt_data
            except Exception as msg:
                print(msg)

            sum_vals = final_df_combined.loc[in_df.index, stn_name].sum()

            print('**Sum of station data \n**', sum_vals)


#             except Exception as msg:
#                 print(msg)
#                 continue
        # break
        all_files_len -= 1

        # break
        if delete_df_files:
            os.remove(df_txt_file)
    #final_df_combined = final_df_combined[final_df_combined >= 0]
#     print('Resampling Dataframe')
#     final_df_combined_resampled = resampleDf(final_df_combined, '5min')

    print('Saving Dataframe')
    final_df_combined.dropna(how='all', inplace=True)
    #final_df_combined.set_index('Time', inplace=True, drop=True)
    final_df_combined.to_csv(
        os.path.join(out_dir,
                     'BW_dwd_PPT_data_%s.csv' % year_period),
        sep=';', float_format='%0.2f')
#     final_df_combined.reset_index(inplace=True)
#     final_df_combined.rename({'index': 'Time'}, inplace=True)
# #
#     final_df_combined.to_feather(
#         os.path.join(out_dir,
#                      'all_dwd_10min_ppt_data_combined_2018_2019.fk'))

    print('done creating df')

# # get all df files as list
# all_df_files = list_all_full_path('.txt', out_dir)
#
# # get all downloaded station ids, used as columns for df_final
# available_stns = []
#
# for df_txt_file in all_df_files:
#
#     stn_name = df_txt_file.split('_')[-1].split('.')[0]
#     print('getting station name from file ', stn_name)
#     assert len(stn_name) == 5
#     if stn_name not in available_stns:
#         available_stns.append(stn_name)
#
#
# nstats = len(available_stns)
#
# hdf5_path = os.path.join(out_dir, 'DWD_DE_data_%s_.h5' % temporal_freq)
# blank_df = pd.DataFrame(index=date_range)  # , data=np.zeros(dates.shape))
#
# if make_hdf5_dataset_new:
#     # number of maximum timesteps
#     nts_max = date_range.shape[0]
#
#     hf = tables.open_file(hdf5_path, 'w', filters=tables.Filters(complevel=6))
#
#     # timestamps
#     hf.create_group(where=hf.root,
#                     name='timestamps',
#                     title='Timestamps of respective Aggregation as Start Time')
#     hf.create_carray(where=hf.root.timestamps,
#                      name='isoformat',
#                      atom=tables.StringAtom(19),
#                      shape=(nts_max,),
#                      chunkshape=(10000,),
#                      title='Strings of Timestamps in Isoformat')
#     hf.create_carray(where=hf.root.timestamps,
#                      name='year',
#                      atom=tables.IntAtom(),
#                      shape=(nts_max,),
#                      chunkshape=(10000,),
#                      title='Yearly Timestamps')
#     hf.create_carray(where=hf.root.timestamps,
#                      name='month',
#                      atom=tables.IntAtom(),
#                      shape=(nts_max,),
#                      chunkshape=(10000,),
#                      title='Monthly Timestamps')
#     hf.create_carray(where=hf.root.timestamps,
#                      name='day',
#                      atom=tables.IntAtom(),
#                      shape=(nts_max,),
#                      chunkshape=(10000,),
#                      title='Daily Timestamps')
#     hf.create_carray(where=hf.root.timestamps,
#                      name='yday',
#                      atom=tables.IntAtom(),
#                      shape=(nts_max,),
#                      chunkshape=(10000,),
#                      title='Yearday Timestamps')
#
#     # data
#     hf.create_carray(where=hf.root,
#                      name='data',
#                      atom=tables.FloatAtom(dflt=np.nan),
#                      shape=(nts_max, nstats),
#                      chunkshape=(10000, 1),
#                      title='DWD %s' % temporal_freq)
#
#     # coordinates
#     hf.create_group(where=hf.root,
#                     name='coord',
#                     title='WGS84 of Stations')
#     hf.create_carray(where=hf.root.coord,
#                      name='lon',
#                      atom=tables.FloatAtom(dflt=np.nan),
#                      shape=(nstats,),
#                      title='Longitude')
#     hf.create_carray(where=hf.root.coord,
#                      name='lat',
#                      atom=tables.FloatAtom(dflt=np.nan),
#                      shape=(nstats,),
#                      title='Latitude')
#     hf.create_carray(where=hf.root.coord,
#                      name='z',
#                      atom=tables.FloatAtom(dflt=np.nan),
#                      shape=(nstats,),
#                      title='Z-Coordinate')
#
#     # metadata
#     hf.create_carray(where=hf.root,
#                      name='name',
#                      atom=tables.StringAtom(50),
#                      shape=(nstats,),
#                      title='Name of Station')
#
#     # convert timestamp to isoformat
#     ts_iso = []
#     ts_year = []
#     ts_month = []
#     ts_day = []
#     ts_yday = []
#
#     for ii in range(date_range.shape[0]):
#         ts = date_range[ii]
#         ts_iso.append(ts.isoformat())
#         # get timestamp years
#         ts_year.append(ts.year)
#         # get timestamp months
#         ts_month.append(ts.month)
#         # get timestamp days
#         ts_day.append(ts.day)
#         # get timestamp year days
#         ts_yday.append(ts.timetuple().tm_yday)
#
#     # fill hf5 with predefined stamps
#     hf.root.timestamps.isoformat[:] = ts_iso[:]
#     hf.root.timestamps.year[:] = ts_year[:]
#     hf.root.timestamps.month[:] = ts_month[:]
#     hf.root.timestamps.day[:] = ts_day[:]
#     hf.root.timestamps.yday[:] = ts_yday[:]
#     hf.close()
#
# #==============================================================================
# #
# #==============================================================================
# # write station data
# hf = tables.open_file(hdf5_path, 'r+')
#
# # get all df files as list
# all_df_files = list_all_full_path('.txt', out_dir)
#
#
# all_files_len = len(all_df_files)
# #     stns_bw_len = len(all_df_files_bw)
#
# for df_txt_file in all_df_files:
#     print('\n++\n Total number of files is \n++\n', all_files_len)
#
#     stn_name = df_txt_file.split('_')[-1].split('.')[0]
#
#     i_stn = np.where(stn_name in available_stns)[0][0]
#     print('\n##\n Station ID is \n##\n', stn_name)
#     in_df = pd.read_csv(df_txt_file, sep=';', index_col=1, engine='c')
#
#     assert int(stn_name) == in_df.loc[:, stn_id_col_name].values[0]
#
#     in_df = in_df[in_df[ppt_col_name] >= 0]
#     in_df = in_df[~in_df.index.duplicated()]
#     assert not in_df.index.duplicated().any(), 'still duplicates in DF'
#     try:
#         in_df = in_df[ppt_col_name]
#     except Exception as msg:
#         print(msg)
#         pass
#     print('\n++\n  Data shape is \n++\n', in_df.shape)
#
#     # TODO: check resampling
#     # resample station data
# #     temp_station = resampleDf(df=temp_station, agg=freq,
# #                               closed='right', label='right',
# #                               shift=False, leave_nan=True,
# #                               label_shift=None,
# #                               temp_shift=0,
#                               max_nan=0)
    # assert (temp_station_res.values.sum() == temp_station.values.sum(),
    #        'error in resampling')
    # if temp_station.index.freq
#     hf.root.data[:, i_stn] = blank_df.join(in_df).values.flatten()
#     hf.root.coord.lon[i_stn] = dwd_in_coords_df['geoLaenge'].loc[stn_name]
#     hf.root.coord.lat[i_stn] = dwd_in_coords_df['geoBreite'].loc[stn_name]
#     hf.root.coord.z[i_stn] = dwd_in_coords_df['Stationshoehe'].loc[stn_name]
#     hf.root.name[i_stn] = np.string_(stn_name)

# hf.close()
"""

#     print('\n+++\n reading df \n+++\n')
#     final_df_combined = pd.read_csv(
#         os.path.join(
#             main_dir, r'all_dwd_hourly_ppt_data_combined_1995_2019.csv'),
#         #                             r'all_ppt_data_combined.csv'),
#         sep=';',
#         index_col=0)
#     final_df_combined.index = pd.to_datetime(final_df_combined.index,
#                                              format='%Y-%m-%d %H:%M:%S')
#     print('\n+++\n resampling df \n+++\n')
#
#     for temp_freq_resample in freqs_list:
#         if temp_freq_resample == '60Min':
#             df_resampled = final_df_combined
#         else:
#             print(temp_freq_resample)
#             df_resampled = resampleDf(
#                 final_df_combined,
#                 temp_freq_resample)
#     #     ppt_data_all = np.array(df_60min.values)
#         df_resampled = df_resampled[df_resampled >= 0]
#         print('\n+++\n creating HDF5 file \n+++\n')
#
#         str_date = str(df_resampled.index[0]).replace(
#             '-', '').replace(':', '').replace(' ', '')
#         ed_date = str(df_resampled.index[-1]).replace(
#             '-', '').replace(':', '').replace(' ', '')
#
#         output_hdf5_file = r"DWD_%s_ppt_stns_%s_%s_new.h5" % (
#             temp_freq_resample, str_date, ed_date)
#
#         create_hdf5(hf5_file=output_hdf5_file,
#                     start_dt=df_resampled.index[0],
#                     end_dt=df_resampled.index[-1],
#                     nstats=len(df_resampled.columns),
#                     freq=temp_freq_resample,
#                     data_title=(r'Precipitation DWD Data Historical and Recent'
#                                 r'Aggregation %s') % temp_freq_resample,
#                     in_ppt_df=df_resampled,
#                     in_stns_df=stn_names_df,
#                     utm=False)


"""
