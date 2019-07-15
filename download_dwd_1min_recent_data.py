
# coding: utf-8

# In[2]:


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
import h5py


#==============================================================================
# 10minutenwerte_rr_00003_19930428_19991231_hist.zip
# stundenwerte_RR_00003_19950901_20110401_hist.zip
# tageswerte_RR_00001_19120101_19860630_hist.zip
#==============================================================================
download_data = False

delete_zip_files = True
delete_df_files = True

make_hdf5_dataset = True

temporal_freq = '1_minute'  # '10_minutes' 'hourly' 'daily'
time_period = 'recent'  # 'historical'
temp_accr = '1minuten'  # '10minuten' 'stundenwerte'  'tages'

ppt_act = 'nieder'  # 'rr' 'RR' 'RR'

stn_name_len = 5  # length of dwd stns Id, ex: 00023

start_date = '2018-01-01 00:00:00'
end_date = '2019-07-11 00:00:00'

temp_freq = 'Min'  # '10Min' 'H' 'D'

time_fmt = '%Y%m%d%H%M'  # when reading download dwd station df

# name of station id and rainfall column name in df
stn_id_col_name = 'STATIONS_ID'
ppt_col_name = 'RS_01'  # 'RWS_10' 'R1' 'RS'

temp_freq_resample = '1D'
#==============================================================================
#
#==============================================================================
# main_dir = os.path.join(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\download_DWD_data_recent')

main_dir = r'E:\download_DWD_data_recent'
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


def resampleDf(data_frame,  # dataframe to resample (or series)
               temp_freq,  # temp frequency to resample
               df_sep_=None,  # sep used if saving dataframe
               out_save_dir=None,  # out direcory for saving df
               fillnan=False,  # if True, fill nans with fill_value
               fillnan_value=0,  # if True, replace nans with 0
               df_save_name=None  # name of outpot df
               ):
    ''' sample DF based on freq and time shift and label shift '''

    # data_frame = data_frame[data_frame >= 0]
    df_res = data_frame.resample(rule=temp_freq,
                                 axis=0,
                                 label='left',
                                 closed='right',
                                 convention='end').apply(lambda x: x.values.sum())

    if fillnan:
        df_res.fillna(value=fillnan_value, inplace=True)
    if df_save_name is not None and out_save_dir is not None:
        df_res.to_csv(os.path.join(out_save_dir, df_save_name),
                      sep=df_sep_)
    if not isinstance(df_res, pd.DataFrame):
        df_res = pd.DataFrame(index=df_res.index, data=df_res.values)
    return df_res


#==============================================================================
#
#==============================================================================

def firstNonNan(listfloats):
    '''return first non NaN value in python list'''
    idx_not_nan = np.argwhere(math.isnan(listfloats) == False)[0]
    return idx_not_nan
#     for item in listfloats:
#         if math.isnan(item) == False:
#             return item


#==============================================================================
# # TODO:
#==============================================================================
# download station coordinates name (maybe do it seperat later)
stn_names_files = os.path.join(main_dir, r'station_coordinates_names.csv')
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
base_url = (r'https://opendata.dwd.de/climate_environment/CDC'
            r'/observations_germany/climate/%s/precipitation/%s/'
            % (temporal_freq, time_period))


# generate url
if download_data:  # download all zip files
    urls_list = []

    # generate file_names, which files to download
    file_names = []
    for stn_id_stn in stations_id_str_lst:
        filename_template = '%swerte_%s_%s_akt.zip' % (temp_accr,
                                                       ppt_act, stn_id_stn)
    file_names.append(filename_template)
    for file_name in file_names:
        file_url = base_url + file_name
        r = requests.get(file_url, stream=True)
        local_filename = file_url.split('/')[-1]

        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024):
                if chunk:  # filter out keep-alive new chunks
                    f.write(chunk)
            r.close()
    #==========================================================================
    #
    #==========================================================================

    # get all zip files as list
    all_zip_files = list_all_full_path('.zip', main_dir)

    # extract all zip files, dfs in zip files
    for zip_file in all_zip_files:
        with zipfile.ZipFile(zip_file, "r") as zip_ref:
            zip_ref.extractall('targetdir')
        if delete_zip_files:
            os.remove(zip_file)

    #==========================================================================
    #
    #==========================================================================
    # get all df files as list
    all_df_files = list_all_full_path('.txt', main_dir)

    # get all downloaded station ids, used as columns for df_final
    available_stns = []

    for df_txt_file in all_df_files:
        stn_name = df_txt_file.split('_')[-1].split('.')[0]
        assert len(stn_name) == 5
        available_stns.append(stn_name)

    # create daterange and df full of nans

    date_range = pd.date_range(start=start_date, end=end_date, freq=temp_freq)

    data = np.zeros(shape=(date_range.shape[0], len(all_df_files)))
    data[data == 0] = np.nan

    final_df_combined = pd.DataFrame(data=data,
                                     index=date_range,
                                     columns=available_stns)

    # read df stn file and fill final df for all stns

    for df_txt_file in all_df_files:
        stn_name = df_txt_file.split('_')[-1].split('.')[0]

        print(stn_name)
        in_df = pd.read_csv(df_txt_file, sep=';', index_col=1, engine='c')
        assert int(stn_name) == in_df.loc[:, stn_id_col_name].values[0]
        assert stn_name in final_df_combined.columns

        in_df.index = pd.to_datetime(in_df.index, format=time_fmt)
        ppt_data = in_df[ppt_col_name].values.ravel()
        final_df_combined.loc[in_df.index, stn_name] = ppt_data

        if delete_df_files:
            os.remove(df_txt_file)

    final_df_combined.to_csv('all_ppt_data_combined.csv', sep=';')


#==============================================================================
# # In[129]:
#==============================================================================

if make_hdf5_dataset:
    print('\n+++\n reading df \n+++\n')
    final_df_combined = pd.read_csv(
        os.path.join(main_dir, r'1D_all_ppt_data_combined.csv'),
        #                      r'all_ppt_data_combined.csv'),
        sep=';',
        index_col=0)
    final_df_combined.index = pd.to_datetime(final_df_combined.index,
                                             format='%Y-%m-%d %H:%M:%S')
    print('\n+++\n resampling df \n+++\n')

    df_60min = resampleDf(
        final_df_combined,
        temp_freq=temp_freq_resample,
        df_sep_=';', out_save_dir=main_dir,
        df_save_name='%s_all_ppt_data_combined.csv' % temp_freq_resample)
    ppt_data_all = np.array(df_60min.values)

    print('\n+++\n creating HDF5 file \n+++\n')

    str_date = start_date.replace('-', '').replace(':', '').replace(' ', '')
    ed_date = end_date.replace('-', '').replace(':', '').replace(' ', '')
    hf = h5py.File("DWD_%s_ppt_stns_%s_%s.h5"
                   % (temp_freq_resample, str_date, ed_date), "w")

    #hf.create_dataset('data', data=ppt_data_all, dtype='f8')
    stns_int = []
    for stn_str in final_df_combined.columns:
        stns_int.append(int(stn_str))

    stns_coords_lon = stn_names_df.loc[stns_int, 'geoBreite'].values
    stns_coords_lat = stn_names_df.loc[stns_int, 'geoLaenge'].values
    stns_coords_z = stn_names_df.loc[stns_int, 'Stationshoehe'].values
    # stns_name = stn_names_df.loc[stns_int, 'Stationsname'].values
    stns_bundesland = stn_names_df.loc[stns_int, 'Bundesland'].values

    start_idx_values, last_idx_values = [], []
    for stn_id in final_df_combined.columns:
        stn_vals = final_df_combined.loc[:, stn_id].values.ravel()
        stn_idx = final_df_combined.loc[:, stn_id].index.values

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

    start_end_idx_df = pd.DataFrame(index=final_df_combined.columns,
                                    data=start_idx_values,
                                    columns=['start_idx'])
    start_end_idx_df['end_idx'] = last_idx_values

    # assign ppt data as dataset
    hf.create_dataset('data', data=ppt_data_all, dtype='f8')

    dt = h5py.special_dtype(vlen=str)
    dset = hf.create_dataset("Bundesland", data=stns_bundesland, dtype=dt)

    hf.create_dataset("id", data=final_df_combined.columns, dtype=dt)

    g2 = hf.create_group('coords')
    g2.create_dataset('lon', data=stns_coords_lon, dtype='f8')
    g2.create_dataset('lat', data=stns_coords_lat, dtype='f8')
    g2.create_dataset('z', data=stns_coords_z, dtype='f8')

    # get time index as string
    time_ix_str_lst = []
    for ix_time in df_60min.index:
        #         time_ix_str_lst.append(ix_time.strftime('%Y-%m-%d %H:%M:%S'))
        time_ix_str_lst.append(ix_time.isoformat())

#     ix0 = df_60min.index[0]
#     d = ix0.strftime('%Y-%m-%d %H:%M:%S')
    g3 = hf.create_group('timestamps')
    g3.create_dataset('day',
                      data=df_60min.index.day,
                      dtype='i4')
    g3.create_dataset('isoformat',
                      data=np.array(time_ix_str_lst).astype(np.object_),
                      dtype=dt)
    g3.create_dataset('start_idx',
                      data=start_end_idx_df.loc[:, 'start_idx'].values,
                      dtype='i4')
    g3.create_dataset('end_idx',
                      data=start_end_idx_df.loc[:, 'end_idx'].values,
                      dtype='i4')

    g3.create_dataset('month',
                      data=df_60min.index.month,
                      dtype='i4')
    g3.create_dataset('yday',
                      data=df_60min.index.dayofyear,
                      dtype='i4')
    g3.create_dataset('year',
                      data=df_60min.index.year,
                      dtype='i4')
    g3.create_dataset('Start_Date',
                      data=np.string_(df_60min.index[0]))
    g3.create_dataset('End_Date',
                      data=np.string_(df_60min.index[-1]))
    g3.create_dataset('Index_%s_Freq' % temp_freq_resample,
                      data=np.arange(0, len(df_60min.index), 1),
                      dtype='i4')
    hf.close()

    print('\n+++\n saved HDF5 file \n+++\n')
