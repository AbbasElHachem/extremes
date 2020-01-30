#==============================================================================
# TODO: ADD DOCUMENATION ORGANIZE
#==============================================================================


import datetime
import os
import time
import timeit
import fnmatch
import pandas as pd
import requests
import zipfile
import numpy as np
import math

from _00_additional_functions import (list_all_full_path, resampleDf,
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

temporal_freq = 'hourly'  # '1_minute' '10_minutes' 'hourly' 'daily'
time_period = 'historical'  # 'recent' 'historical'
temp_accr = 'stundenwerte'  # '1minuten' '10minuten' 'stundenwerte'  'tages'

ppt_act = 'RR'  # nieder 'rr' 'RR' 'RR'

stn_name_len = 5  # length of dwd stns Id, ex: 00023

start_date = '1995-01-01 00:00:00'
end_date = '2019-12-30 23:59:00'

temp_freq = 'H'  # '60Min'  # 'Min' '10Min' 'H' 'D'

# '%Y%m%d%H'  # when reading download dwd station df '%Y%m%d%H%M'
time_fmt = '%Y%m%d%H'  # :%S' # for 10Min data use: '%Y%m%d%H%M:%S'

# name of station id and rainfall column name in df
stn_id_col_name = 'STATIONS_ID'
ppt_col_name = '  R1'  # 'TT_TU'  # RS_01 'RWS_10' '  R1' 'RS'

freqs_list = ['60Min']  # '60Min']  # '1D' '5Min',
#==============================================================================
#
#==============================================================================
# main_dir = os.path.join(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\download_DWD_data_recent')

main_dir = r'F:\download_DWD_data_recent'
os.chdir(main_dir)


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
stn_names_files = os.path.join(
    main_dir, r'station_coordinates_names_hourly_only_in_BW_utm32.csv')

assert os.path.exists(stn_names_files)
#     in_stn_coords_df_loc = os.path.join(
#         r"F:\download_DWD_data_recent\station_coordinates_names_hourly_only_in_BW_utm32.csv")

stn_names_df = pd.read_csv(
    stn_names_files, index_col=0, sep=';', encoding='latin-1')
stations_id = stn_names_df.index

# get all station ids, make them a string for generating file_names
stations_id_str_lst = []
for stn_id in stations_id:
    stn_id = str(stn_id)
    if len(stn_id) < 5:
        stn_id = '0' * (5 - len(stn_id)) + stn_id
    stations_id_str_lst.append(stn_id)
# print(stations_id_str)

#==============================================================================
# DOWNLOAD ALL ZIP FILES
#==============================================================================
out_dir = os.path.join(main_dir, 'DWD_data_%s_%s_%s'
                       % (temporal_freq, time_period, temp_accr))
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

os.chdir(out_dir)

# generate url
if download_data:  # download all zip files
    # precipitation
    base_url = (r'https://opendata.dwd.de/climate_environment/CDC'
                r'/observations_germany/climate/%s/precipitation/%s/'
                % (temporal_freq, time_period))

    r = requests.get(base_url, stream=True)

    all_urls = r.content.split(b'\n')

#     if time_period == 'historical':
#         _zip_file = -3
#     if time_period == 'recent':
#         _zip_file = -2

    all_urls_zip = [url_zip.split(b'>')[-2].replace(b'</a', b'').decode('utf-8')
                    for url_zip in all_urls if b'.zip' in url_zip]

    for file_name in all_urls_zip:

        #         if time_period == 'historical':
        #             file_name = file_name.split('=')[1].replace('"', '')
        print('getting data for', file_name)

        file_url = base_url + file_name
        zip_url = requests.get(file_url)

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

    out_extract_df_dir = ('DWD_extracted_zip_%s_%s_%s' %
                          (temporal_freq, time_period, temp_accr))

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

# get all df files as list
if build_one_df:
    #     all_df_files = list_all_full_path('.txt', out_extract_df_dir)

    all_df_files = list_all_full_path(
        '.txt', out_dir)

#     all_df_files_bw = [_df_f for stn_name in stations_id_str_lst
#                        for _df_f in all_df_files if stn_name in _df_f]

    all_df_files = [_df_f for _df_f in all_df_files]

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
#     stns_bw_len = len(all_df_files_bw)

    for df_txt_file in all_df_files:
        print('\n++\n Total number of files is \n++\n', all_files_len)

        stn_name = df_txt_file.split('_')[-1].split('.')[0]
        if stn_name in stations_id_str_lst:
            #             print('\n##\n Stations in BW \n##\n', stns_bw_len)
            print('\n##\n Station ID is \n##\n', stn_name)
            in_df = pd.read_csv(df_txt_file, sep=';', index_col=1, engine='c')

            assert int(stn_name) == in_df.loc[:, stn_id_col_name].values[0]
            assert stn_name in final_df_combined.columns

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
#             in_df = in_df[in_df[ppt_col_name] >= 0]
            try:
                ppt_data = in_df[ppt_col_name].values.ravel()
            except Exception as msg:
                print(msg)
                pass
            print('\n++\n  Data shape is \n++\n', ppt_data.shape)

            try:
                final_df_combined.loc[in_df.index, stn_name] = ppt_data
            except Exception:
                print('Error')
                pass
#             sum_vals = final_df_combined.loc[in_df.index, stn_name].sum()
#             if sum_vals < 0:
#                 raise Exception
#             print('**Sum of station data \n**', sum_vals)
#             stns_bw_len -= 1

#             except Exception as msg:
#                 print(msg)
#                 continue
        # break
            #all_files_len -= 1

        # break
#             if delete_df_files:
#                 os.remove(df_txt_file)
    final_df_combined = final_df_combined[final_df_combined >= 0]
#     print('Resampling Dataframe')
#     final_df_combined_resampled = resampleDf(final_df_combined, '5min')
    final_df_combined.plot(legend=False)
    print('Saving Dataframe')
    final_df_combined.dropna(how='all', inplace=True)
    #final_df_combined.set_index('Time', inplace=True, drop=True)
    final_df_combined.to_csv(
        os.path.join(out_dir,
                     'all_dwd_60min_ppt_data_combined.csv'),
        sep=';')
    # all_dwd_60min_temp_data_combined_2010_2018.csv
#     final_df_combined.reset_index(inplace=True)
#     final_df_combined.rename({'index': 'Time'}, inplace=True)
# #
#     final_df_combined.to_feather(
#         os.path.join(out_dir,
#                      'all_dwd_10min_ppt_data_combined_2018_2019.fk'))

    print('done creating df')
#==============================================================================
# # In[129]:
#==============================================================================


if make_hdf5_dataset_new:
    from a_functions import create_hdf5

    print('\n+++\n reading df \n+++\n')
    final_df_combined = pd.read_csv(
        os.path.join(
            main_dir, r'all_dwd_hourly_ppt_data_combined_1995_2019.csv'),
        #                             r'all_ppt_data_combined.csv'),
        sep=';',
        index_col=0)
    final_df_combined.index = pd.to_datetime(final_df_combined.index,
                                             format='%Y-%m-%d %H:%M:%S')
    print('\n+++\n resampling df \n+++\n')

    for temp_freq_resample in freqs_list:
        if temp_freq_resample == '60Min':
            df_resampled = final_df_combined
        else:
            print(temp_freq_resample)
            df_resampled = resampleDf(
                final_df_combined,
                temp_freq_resample)
    #     ppt_data_all = np.array(df_60min.values)
        df_resampled = df_resampled[df_resampled >= 0]
        print('\n+++\n creating HDF5 file \n+++\n')

        str_date = str(df_resampled.index[0]).replace(
            '-', '').replace(':', '').replace(' ', '')
        ed_date = str(df_resampled.index[-1]).replace(
            '-', '').replace(':', '').replace(' ', '')

        output_hdf5_file = r"DWD_%s_ppt_stns_%s_%s_new.h5" % (
            temp_freq_resample, str_date, ed_date)

        create_hdf5(hf5_file=output_hdf5_file,
                    start_dt=df_resampled.index[0],
                    end_dt=df_resampled.index[-1],
                    nstats=len(df_resampled.columns),
                    freq=temp_freq_resample,
                    data_title=(r'Precipitation DWD Data Historical and Recent'
                                r'Aggregation %s') % temp_freq_resample,
                    in_ppt_df=df_resampled,
                    in_stns_df=stn_names_df,
                    utm=False)
