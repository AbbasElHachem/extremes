
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


# In[3]:


# main_dir = os.path.join(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\download_DWD_data_recent')

main_dir = r'E:\download_DWD_data_recent'
os.chdir(main_dir)


# In[4]:


base_url = r'https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/1_minute/precipitation/recent/'

filename_template = '1minutenwerte_nieder_%s_akt.zip'

stn_names_files =  os.path.join(main_dir, r'station_coordinates_names.csv')
assert os.path.exists(stn_names_files)


# In[5]:


stn_names_df = pd.read_csv(stn_names_files, index_col=0, sep=',', encoding='latin-1')
stations_id = stn_names_df.index

# get all station ids, make them a string for generating file_names
stations_id_str_lst = []
for stn_id in stations_id:
    stn_id = str(stn_id)
    if len(stn_id) < 5:
        stn_id = '0'*(5-len(stn_id)) + stn_id
    stations_id_str_lst.append(stn_id)
#print(stations_id_str)


# In[80]:


# generate file_names, which files to download
file_names = []
for stn_id_stn in stations_id_str_lst:
    filename_template = r'1minutenwerte_nieder_%s_akt.zip' % stn_id_stn
    file_names.append(filename_template)


# In[68]:


# generate url

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
        


# In[81]:


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

all_zip_files = list_all_full_path('.zip', main_dir)


# In[84]:


# extract all zip files, dfs in zip files
for zip_file in all_zip_files:
    with zipfile.ZipFile(zip_file,"r") as zip_ref:
        zip_ref.extractall('targetdir')
    


# In[127]:


all_df_files = list_all_full_path('.txt', main_dir)

available_stns = []

for df_txt_file in all_df_files:
    stn_name = df_txt_file.split('_')[-1].split('.')[0]
    assert len(stn_name) == 5
    available_stns.append(stn_name)
    
    
start_date = '2018-01-01 00:00:00'
end_date = '2019-07-11 00:00:00'

date_range = pd.date_range(start=start_date, end=end_date, freq='Min')

data = np.zeros(shape=(date_range.shape[0], len(all_df_files)))
data[data == 0] = np.nan

final_df_combined = pd.DataFrame(data=data, index=date_range,
                                columns=available_stns)


    


# In[128]:


time_fmt = '%Y%m%d%H%M'

for df_txt_file in all_df_files:
    stn_name = df_txt_file.split('_')[-1].split('.')[0]
    print(stn_name)
    in_df = pd.read_csv(df_txt_file, sep=';', index_col=1, engine='c')
    assert int(stn_name) == in_df.loc[:, 'STATIONS_ID'].values[0]
    assert stn_name in final_df_combined.columns
    
    in_df.index = pd.to_datetime(in_df.index, format=time_fmt)
    ppt_data = in_df['RS_01'].values.ravel()
    final_df_combined.loc[in_df.index, stn_name] = ppt_data
    


# In[130]:


final_df_combined.to_csv('all_ppt_data_combined.csv', sep=';')


# In[129]:


ppt_data_all = np.array(final_df_combined.values)


# In[1]:


hf = h5py.File("DWD_ppt_stns_201801010000_201907100000.h5", "w")

hf_test = h5py.File("DWD_ppt_stns_201801010000_201907100000_Test.h5", "w")


# In[152]:


#hf.create_dataset('data', data=ppt_data_all, dtype='f8')
stns_int = []
for stn_str in available_stns:
    stns_int.append(int(stn_str))

stns_coords_lon =  stn_names_df.loc[stns_int, 'geoBreite'].values
stns_coords_lat =  stn_names_df.loc[stns_int, 'geoLaenge'].values
stns_coords_z =  stn_names_df.loc[stns_int, 'Stationshoehe'].values
stns_name = stn_names_df.loc[stns_int, 'Stationsname'].values
stns_bundesland = stn_names_df.loc[stns_int, 'Bundesland'].values


# In[153]:


stns_bundesland


# In[154]:


# assign ppt data as dataset
hf.create_dataset('data', data=ppt_data_all, dtype='f8')


# In[165]:


dt = h5py.special_dtype(vlen=str)

dset = hf_test.create_dataset("Bundesland", data=stns_bundesland, dtype=dt)
#hf_test.create_dataset('Bundesland', data=stns_bundesland)


# In[166]:


hf_test.close()


# In[155]:


g2 = hf.create_group('coords')
g2.create_dataset('lon', data=stns_coords_lon, dtype='f8')
g2.create_dataset('lat', data=stns_coords_lat, dtype='f8')
g2.create_dataset('z', data=stns_coords_z, dtype='f8')


# In[148]:


np.string_(stns_bundesland)


# In[169]:


g3 = hf.create_group('Timestamps')
g3.create_dataset('Start_Date', data=np.string_(final_df_combined.index[0]))
g3.create_dataset('End_Date', data=np.string_(final_df_combined.index[-1]))
g3.create_dataset('Index_1Min_Freq',
                  data=np.arange(0, len(final_df_combined.index), 1),
                  dtype='i4')
hf.close()

