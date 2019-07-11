
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

import h5py
#==============================================================================
#
#==============================================================================
download_data = False

delete_zip_files = True
delete_df_files = True

make_hdf5_dataset = True

temporal_freq = '1_minute'  # '10_minutes' 'hourly' 'daily'
time_period = 'recent'  # 'historical'
temp_accr = '1minuten'  # '10minuten' 'stundenwerte'  'tages'

stn_name_len = 5  # length of dwd stns Id, ex: 00023

start_date = '2018-01-01 00:00:00'
end_date = '2019-07-11 00:00:00'

temp_freq = 'Min'  # '10Min' 'H' 'D'

time_fmt = '%Y%m%d%H%M'  # when reading download dwd station df

# name of station id and rainfall column name in df
stn_id_col_name = 'STATIONS_ID'
ppt_col_name = 'RS_01'  # 'RWS_10' 'R1' 'RS'
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


# generate file_names, which files to download
file_names = []
for stn_id_stn in stations_id_str_lst:
    filename_template = '%swerte_nieder_%s_akt.zip' % (temp_accr, stn_id_stn)
    file_names.append(filename_template)

# generate url
if download_data:  # download all zip files
    urls_list = []

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

    final_df_combined = pd.read_csv(
        os.path.join(main_dir,
                     r'all_ppt_data_combined.csv'),
        sep=';',
        index_col=0)

    ppt_data_all = np.array(final_df_combined.values)

    str_date = start_date.replace('-', '').replace(':', '').replace(' ', '')
    ed_date = start_date.replace('-', '').replace(':', '').replace(' ', '')
    hf = h5py.File("DWD_ppt_stns_%s_%s.h5" % (str_date, ed_date), "w")

    #hf.create_dataset('data', data=ppt_data_all, dtype='f8')
    stns_int = []
    for stn_str in available_stns:
        stns_int.append(int(stn_str))

    stns_coords_lon = stn_names_df.loc[stns_int, 'geoBreite'].values
    stns_coords_lat = stn_names_df.loc[stns_int, 'geoLaenge'].values
    stns_coords_z = stn_names_df.loc[stns_int, 'Stationshoehe'].values
    # stns_name = stn_names_df.loc[stns_int, 'Stationsname'].values
    stns_bundesland = stn_names_df.loc[stns_int, 'Bundesland'].values

    # assign ppt data as dataset
    hf.create_dataset('data', data=ppt_data_all, dtype='f8')

    dt = h5py.special_dtype(vlen=str)
    dset = hf.create_dataset("Bundesland", data=stns_bundesland, dtype=dt)

    g2 = hf.create_group('coords')
    g2.create_dataset('lon', data=stns_coords_lon, dtype='f8')
    g2.create_dataset('lat', data=stns_coords_lat, dtype='f8')
    g2.create_dataset('z', data=stns_coords_z, dtype='f8')

    g3 = hf.create_group('Timestamps')
    g3.create_dataset('Start_Date',
                      data=np.string_(final_df_combined.index[0]))
    g3.create_dataset('End_Date',
                      data=np.string_(final_df_combined.index[-1]))
    g3.create_dataset('Index_1Min_Freq',
                      data=np.arange(0, len(final_df_combined.index), 1),
                      dtype='i4')
    hf.close()
