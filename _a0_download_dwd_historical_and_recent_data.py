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

from _00_additional_functions import resampleDf
#==============================================================================
# 10minutenwerte_rr_00003_19930428_19991231_hist.zip
# stundenwerte_RR_00003_19950901_20110401_hist.zip
# tageswerte_RR_00001_19120101_19860630_hist.zip
#==============================================================================
download_data = True
delete_zip_files = True


build_one_df = True
delete_df_files = False

make_hdf5_dataset_new = False

temporal_freq = '1_minute'  # '1_minute' '10_minutes' 'hourly' 'daily'
time_period = 'recent'  # 'recent' ''
temp_accr = '1minuten'  # '1minuten' '10minuten' 'stundenwerte'  'tages'

ppt_act = 'RR'  # nieder 'rr' 'RR' 'RR'

stn_name_len = 5  # length of dwd stns Id, ex: 00023

start_date = '2014-01-01 00:00:00'
end_date = '2019-08-20 00:00:00'

temp_freq = 'Min'  # '60Min'  # 'Min' '10Min' 'H' 'D'

# '%Y%m%d%H'  # when reading download dwd station df '%Y%m%d%H%M'
time_fmt = '%Y%m%d%H%M'

# name of station id and rainfall column name in df
stn_id_col_name = 'STATIONS_ID'
ppt_col_name = 'RS_01'  # '  R1'  # RS_01 'RWS_10' '  R1' 'RS'

freqs_list = ['5Min']  # '60Min']  # '1D' '5Min',
#==============================================================================
#
#==============================================================================
# main_dir = os.path.join(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\download_DWD_data_recent')

main_dir = r'F:\download_DWD_data_recent'
os.chdir(main_dir)


#==============================================================================
#
#==============================================================================


def list_all_full_path(ext, file_dir):
    """
    Purpose: To return full path of files in all dirs of a 
            given folder with a
            given extension in ascending order.
    Description of the arguments:
        ext (string) = Extension of the files to list
            e.g. '.txt', '.tif'.
        file_dir (string) = Full path of the folder in which the files
            reside.
    """
    new_list = []
    patt = '*' + ext
    for root, _, files in os.walk(file_dir):
        for elm in files:
            if fnmatch.fnmatch(elm, patt):
                full_path = os.path.join(root, elm)
                new_list.append(full_path)
    return(sorted(new_list))

#==============================================================================
#
#==============================================================================


# def resampleDf(data_frame,  # dataframe to resample (or series)
#                temp_freq,  # temp frequency to resample
#                df_sep_=None,  # sep used if saving dataframe
#                out_save_dir=None,  # out direcory for saving df
#                fillnan=False,  # if True, fill nans with fill_value
#                fillnan_value=0,  # if True, replace nans with 0
#                df_save_name=None  # name of outpot df
#                ):
#     ''' sample DF based on freq and time shift and label shift '''
#
#     # data_frame = data_frame[data_frame >= 0]
#     df_res = data_frame.resample(rule=temp_freq,
#                                  axis=0,
#                                  label='left',
#                                  closed='right',
#                                  convention='end').apply(lambda x: np.nansum(x.values))
#
#     if fillnan:
#         df_res.fillna(value=fillnan_value, inplace=True)
#     if df_save_name is not None and out_save_dir is not None:
#         df_res.to_csv(os.path.join(out_save_dir, df_save_name),
#                       sep=df_sep_)
#     if not isinstance(df_res, pd.DataFrame):
#         df_res = pd.DataFrame(index=df_res.index, data=df_res.values)
#     return df_res


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
    main_dir, r'station_coordinates_names_hourly.csv')
assert os.path.exists(stn_names_files)

stn_names_df = pd.read_csv(
    stn_names_files, index_col=0, sep=',', encoding='latin-1')
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

# generate url
if download_data:  # download all zip files
    base_url = (r'https://opendata.dwd.de/climate_environment/CDC'
                r'/observations_germany/climate/%s/precipitation/%s/'
                % (temporal_freq, time_period))

    out_dir = os.path.join(main_dir, 'DWD_data_%s_%s_%s'
                           % (temporal_freq, time_period, temp_accr))
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    os.chdir(out_dir)

    r = requests.get(base_url, stream=True)

    all_urls = r.content.split(b'\n')

    all_urls_zip = [url_zip.split(b'>')[-2].replace(b'</a', b'').decode('utf-8')
                    for url_zip in all_urls if b'.zip' in url_zip]

    for file_name in all_urls_zip:

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
        '.txt',
        r'F:\download_DWD_data_recent\DWD_data_1_minute_recent_1minuten')

    # get all downloaded station ids, used as columns for df_final
    available_stns = []

    for df_txt_file in all_df_files:

        stn_name = df_txt_file.split('_')[-1].split('.')[0]
        print('getting station name from file ', stn_name)
        assert len(stn_name) == 5
        if stn_name not in available_stns:
            available_stns.append(stn_name)

    # create daterange and df full of nans

    date_range = pd.date_range(start=start_date, end=end_date, freq=temp_freq)

    data = np.zeros(shape=(date_range.shape[0], len(available_stns)))
    data[data == 0] = np.nan

    final_df_combined = pd.DataFrame(data=data,
                                     index=date_range,
                                     columns=available_stns)

    # read df stn file and fill final df for all stns
    all_files_len = len(all_df_files)
    for df_txt_file in all_df_files:
        print('\n++\n Total number of files is \n++\n', all_files_len)

        stn_name = df_txt_file.split('_')[-1].split('.')[0]

        print('\n##\n Station ID is \n##\n', stn_name)
        in_df = pd.read_csv(df_txt_file, sep=';', index_col=1, engine='c')
        assert int(stn_name) == in_df.loc[:, stn_id_col_name].values[0]
        assert stn_name in final_df_combined.columns

        in_df.index = pd.to_datetime(in_df.index, format=time_fmt)
        ppt_data = in_df[ppt_col_name].values.ravel()
        print('\n++\n  Data shape is \n++\n', ppt_data.shape)

        final_df_combined.loc[in_df.index, stn_name] = ppt_data

        all_files_len -= 1

        if delete_df_files:
            os.remove(df_txt_file)
    final_df_combined = final_df_combined[final_df_combined >= 0]
    final_df_combined_resampled = resampleDf(final_df_combined, '5min')
    final_df_combined_resampled.to_csv(
        'all_dwd_5min_ppt_data_combined_2015_2019.csv',
        sep=';')
    final_df_combined.reset_index()
    final_df_combined.rename({'index': 'Time'})
    final_df_combined.to_feather(
        'all_dwd_5min_ppt_data_combined_2015_2019.csv',)

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
