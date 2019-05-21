import os

import fnmatch
import pyproj

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

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


def convert_coords_fr_wgs84_to_utm32_(epgs_initial_str, epsg_final_str,
                                      first_coord, second_coord):
    '''fct to convert points from wgs 84 to utm32'''
    initial_epsg = pyproj.Proj(epgs_initial_str)
    final_epsg = pyproj.Proj(epsg_final_str)
    x, y = pyproj.transform(initial_epsg, final_epsg,
                            first_coord, second_coord)
    return x, y

#==============================================================================
# Get colors by name
#==============================================================================


colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

# Sort colors by hue, saturation, value and name.
by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
                for name, color in colors.items())
sorted_names_all = [name for hsv, name in by_hsv]


def pltcolor():
    cols_time_dict = {i: [] for i in np.arange(-60, 61, 5)}

    for i, l in enumerate(cols_time_dict.keys()):
        if l == -60:
            cols_time_dict[l].append(sorted_names_all[i * 12])
        elif l == -55:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -50:
            cols_time_dict[l].append(sorted_names_all[i * 8])
        elif l == -45:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -40:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -35:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -30:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -25:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -20:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -15:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -10:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == -5:
            cols_time_dict[l].append(sorted_names_all[i * 6])

        elif l == 0:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 5:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 10:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 15:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 20:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 25:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 30:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 35:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 40:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 45:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 50:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 55:
            cols_time_dict[l].append(sorted_names_all[i * 6])
        elif l == 60:
            cols_time_dict[l].append(sorted_names_all[i * 6])
    return cols_time_dict

#==============================================================================
#
#==============================================================================


def get_stns_ids_coords_file(netatmo_coords_df_file):
    ''' get all ids from coords_file'''
    in_coords_df = pd.read_csv(netatmo_coords_df_file, index_col=0, sep=';')
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


def resampleDf(data_frame,
               temp_freq,
               temp_shift=0,
               label_shift=None,
               df_sep_=None,
               out_save_dir=None,
               fillnan=False,
               df_save_name=None):
    ''' sample DF based on freq and time shift and label shift '''

    df_ = data_frame.copy()
    df_res = df_.resample(temp_freq,
                          label='right',
                          closed='right',
                          loffset=label_shift,
                          base=temp_shift).sum()

    if fillnan:
        df_res.fillna(value=0, inplace=True)
    if df_save_name is not None and out_save_dir is not None:
        df_res.to_csv(os.path.join(out_save_dir, df_save_name),
                      sep=df_sep_)
    return df_res

#==============================================================================
#
#==============================================================================


def select_df_within_period(df, start, end):
    ''' a function to select df between two dates'''
    mask = (df.index > start) & (df.index <= end)
    df_period = df.loc[mask]
    return df_period

#==============================================================================
#
#==============================================================================


def calculate_probab_ppt_below_thr(ppt_data, ppt_thr):
    ''' calculate probability of values being below threshold '''
    origin_count = ppt_data.shape[0]
    count_below_thr = ppt_data[ppt_data <= ppt_thr].shape[0]
    p0 = np.divide(count_below_thr, origin_count)
    return p0


#==============================================================================
#
#==============================================================================


def build_edf_fr_vals(ppt_data):
    # Construct EDF
    ''' construct empirical distribution function given data values '''
    data_sorted = np.sort(ppt_data, axis=0)[::-1]
    x0 = np.squeeze(data_sorted)[::-1]
    y0 = (np.arange(data_sorted.size) / len(data_sorted))
    return x0, y0


#==============================================================================
#
#==============================================================================


def get_cdf_part_abv_thr(ppt_data, ppt_thr):
    ''' select part of the CDF that is abv ppt thr '''

    p0 = calculate_probab_ppt_below_thr(ppt_data, ppt_thr)

    x0, y0 = build_edf_fr_vals(ppt_data)
    x_abv_thr = x0[x0 > ppt_thr]
    y_abv_thr = y0[np.where(x0 > ppt_thr)]
    assert y_abv_thr[0] == p0, 'something is wrong with probability cal'

    return x_abv_thr, y_abv_thr
#==============================================================================
#
#==============================================================================
