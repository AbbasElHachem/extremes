import os
import shutil
import urllib.request as urequest
import zipfile
from collections import OrderedDict
from datetime import datetime
from ftplib import FTP
from pathlib import Path
from functools import reduce
from warnings import warn

import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt
from matplotlib.transforms import offset_copy
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import cartopy.feature as cfeature


data_dir = Path("~/data/opendata_dwd").expanduser()
# ftp_url = "ftp-cdc.dwd.de"
ftp_url = "opendata.dwd.de"
ftp_root = "climate_environment/CDC/observations_germany/climate/hourly/{variable}/historical"
# ftp_root = "pub/CDC/observations_germany/climate/hourly/{variable}/historical"
# ftp_root = "pub/CDC/observations_germany/climate/daily/{variable}/historical"
meta_url = (f"ftp://{ftp_url}/{ftp_root}"
            "/{variable_short}_Stundenwerte_Beschreibung_Stationen.txt")
data_url = (f"ftp://{ftp_url}/{ftp_root}")
zip_template = "stundenwerte_{variable_short}_{Stations_id:05d}"
variable_shorts = {"cloudiness": "N",
                   # solar is f**d
                   # "solar": "ST",
                   "sun": "SD",
                   "precipitation": "RR",
                   "pressure": "P0",
                   "soil_temperature": "EB",
                   "air_temperature": "TU",
                   "relative_humidity": "TU",
                   "wind": "FF"}
variable_cols = {"wind_speed": ["F"],
                 "wind_direction": ["D"],
                 "wind": ["F", "D"],
                 "air_temperature": ["TT_TU"],
                 "relative_humidity": ["RF_TU"],
                 "sun": ["SD_SO"],
                 "precipitation": ["R1"],
                 }
cols_variable = {val[0]: key
                 for key, val in variable_cols.items()}
header_dtypes = OrderedDict((("Stations_id", int),
                             ("von_datum", int),
                             ("bis_datum", int),
                             ("Stationshoehe", float),
                             ("geoBreite", float),
                             ("geoLaenge", float),
                             ("Stationsname", str),
                             ("Bundesland", str),))
meta_header = list(header_dtypes.keys())
data_dir.mkdir(parents=True, exist_ok=True)


def load_metadata(variable):
    meta_filename = f"metadata_{variable}.txt"
    meta_filepath = data_dir / meta_filename
    if not meta_filepath.exists():
        url = meta_url.format(variable=variable,
                              variable_short=variable_shorts[variable])
        with meta_filepath.open("w") as fi:
            with urequest.urlopen(url) as req:
                content = req.read()
                fi.write(content.decode("latin-1"))

    def date_parser(stamp):
        return datetime.strptime(stamp, "%Y%m%d")

    # ok, this is brittle as the DWD might not be consistent about
    # column widths, but colspecs=infer only checks the first 100
    # lines and ends up with widths that are too small for some
    # station names.
    colspecs = [(0, 5),
                (6, 14),
                (15, 23),
                (23, 38),
                (38, 50),
                (50, 60),
                (61, 102),
                (102, 200)]
    df = pd.read_fwf(meta_filepath, skiprows=[0, 1],
                     colspecs=colspecs,
                     dtype=header_dtypes,
                     parse_dates=[1, 2],
                     date_parser=date_parser)
    df.columns = meta_header
    return df


def _load_station_one_var(station_name, variable):
    variable_col = variable_cols[variable]
    variable_short = variable_shorts[variable]
    if variable == "relative_humidity":
        rh = True
        variable = "air_temperature"
    else:
        rh = False
    metadata_df = load_metadata(variable)
    metadata_df = metadata_df[metadata_df["Stationsname"] == station_name]
    if metadata_df.size == 0:
        print(f"No {variable} for {station_name}")
        return pd.DataFrame([])
    metadata = {key: val[0] for key, val in
                metadata_df.to_dict("list").items()}
    zip_filename = zip_template.format(variable_short=variable_short,
                                       **metadata)
    (data_dir / variable).mkdir(parents=True, exist_ok=True)
    zip_filepath = data_dir / variable / zip_filename
    zip_filepath_zip = zip_filepath.with_suffix(".zip")
    if not zip_filepath_zip.exists():
        print(f"Downloading {zip_filepath_zip.name}")
        # we do not know what the filename is exactly because at least
        # the "bis_datum" is different in metadata and filename
        ftp = FTP(ftp_url)
        ftp_path = ftp_root.format(variable=variable)
        ftp.login()
        zip_filename = [name[0] for name in ftp.mlsd(ftp_path)
                        if name[0].startswith(zip_filename)]
        try:
            zip_filename = zip_filename[0]
        except IndexError:
            print(f"{zip_filename} not available online. Alert DWD.")
            return pd.DataFrame([])
        url = os.path.join(data_url.format(variable=variable),
                           zip_filename)
        with zip_filepath_zip.open("wb") as fi:
            with urequest.urlopen(url) as req:
                content = req.read()
                fi.write(content)
        ftp.close()
    with zipfile.ZipFile(str(zip_filepath_zip)) as zif:
        data_name = [name for name in zif.namelist()
                         if name.startswith("produkt_")
                     ][0]
        filepath_txt = zip_filepath_zip.with_suffix(".txt")
        if not filepath_txt.exists():
            zif.extract(data_name, path=zip_filepath.parent)
            produkt_path = zip_filepath.parent / data_name
            shutil.move(produkt_path, filepath_txt)

    def date_parser(stamp):
        # faster version than this:
        # return datetime.strptime(stamp, "%Y%m%d%H")
        return datetime(int(stamp[:4]), int(stamp[4:6]),
                   int(stamp[6:8]), int(stamp[8:10]))

    cols = ["MESS_DATUM"] + variable_col
    # colspecs = [(0, 11),
    #             (12, 22),
    #             (23, 28),
    #             (29, 35),
    #             (36, 40),
    #             (41, 44),
    #             ]
    # usecols = [1, 3]
    # if len(variable_col) == 2:
    #     usecols += [4]
    # names = ["time"] + [cols_variable[col] for col in cols[1:]]
    # df = pd.read_fwf(filepath_txt,
    #                  skiprows=1,
    #                  index_col=0,
    #                  names=names,
    #                  colspecs=colspecs,
    #                  usecols=usecols,
    #                  parse_dates=True,
    #                  date_parser=date_parser,
    #                  na_values="-999",
    #                  )
    # despite having the regex-separator, this is faster:
    df = pd.read_csv(filepath_txt, sep=r";\s*", index_col=0,
                     parse_dates=True, date_parser=date_parser,
                     engine="python", usecols=cols, na_values="-999")
    df.index.name = "time"
    df.columns = [cols_variable[col] for col in cols[1:]]
    df.name = station_name
    return df


def load_station(station_name, variables):
    if isinstance(variables, str):
        variables = variables,
    var_list = [_load_station_one_var(station_name, variable)
                for variable in variables]
    var_list = [vals for vals in var_list if vals.size > 0]
    ds = xr.merge(var_list)
    return ds.to_array(dim="variable", name=station_name)


def load_data(station_names, variables):
    if isinstance(station_names, str):
        station_names = station_names,
    if isinstance(variables, str):
        variables = variables,
    data_dict = {station_name: load_station(station_name, variables)
                 for station_name in station_names}
    return xr.Dataset(data_dict).to_array(dim="station")


def filter_metadata(metadata, lon_min=None, lat_min=None, lon_max=None,
                  lat_max=None, start=None, end=None):
    lons = metadata["geoLaenge"]
    lats = metadata["geoBreite"]
    starts = metadata["von_datum"]
    ends = metadata["bis_datum"]
    if lon_min is None:
        lon_min = lons.min()
    if lat_min is None:
        lat_min = lats.min()
    if lon_max is None:
        lon_max = lons.max()
    if lat_max is None:
        lat_max = lats.max()
    if start is None:
        start = starts.min()
    if end is None:
        end = ends.max()
    mask = ((lons >= lon_min) & (lons < lon_max) &
            (lats >= lat_min) & (lats < lat_max) &
            (starts < start) & (ends > end))
    if not mask.sum():
        warn("No stations match requirements.")
    return metadata[mask]


def get_metadata(variables):
    variable_metadata = {variable: load_metadata(variable)
                         for variable in variables}
    # find the stations that offer all variables
    station_ids = [var["Stations_id"] for var in variable_metadata.values()]
    station_ids = reduce(np.intersect1d, station_ids)
    metadata = (variable_metadata
                [variables[0]]
                .set_index("Stations_id")
                .loc[station_ids]
                )
    return metadata


def map_stations(variables, lon_min=None, lat_min=None, lon_max=None,
               lat_max=None, start=None, end=None, **skwds):
    metadata = get_metadata(variables)
    metadata = filter_metadata(metadata, lon_min, lat_min, lon_max,
                               lat_max, start, end)
    
    stamen_terrain = cimgt.StamenTerrain()
    crs = ccrs.PlateCarree()
    land_10m = cfeature.NaturalEarthFeature('cultural',
                                            'admin_0_countries',
                                            '10m',
                                            edgecolor=(0, 0, 0, 0),
                                            facecolor=(0, 0, 0, 0), )
    fig = plt.figure(**skwds)
    ax = fig.add_subplot(111, projection=stamen_terrain.crs)
    geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax)
    text_transform = offset_copy(geodetic_transform, units='dots', x=-25)
    lons = metadata["geoLaenge"]
    lats = metadata["geoBreite"]
    station_names = metadata["Stationsname"]
    starts = metadata["von_datum"]  # don't trust this too much!
    ends = metadata["bis_datum"]
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=crs)
    for i in metadata.index:
        text = f"{station_names[i]}\n{starts[i].year}-{ends[i].year}"
        ax.text(lons[i], lats[i], text, fontsize=8, transform=text_transform)
    ax.scatter(lons, lats, transform=crs)
    ax.add_feature(land_10m, alpha=.1)
    ax.add_image(stamen_terrain, 10)
    return fig, ax

def map_stations_var(variables, lon_min=None, lat_min=None,
                   lon_max=None, lat_max=None, **skwds):
    if isinstance(variables, str):
        variables = [variables]
    n_variables = len(variables)
    variable_metadata = {variable: load_metadata(variable)
                         for variable in variables}
    all_station_names = [variable_metadata
                         [variable]
                         ["Stationsname"]
                         .values
                         for variable in variables]
    station_names_inter = reduce(np.intersect1d, all_station_names)

    stamen_terrain = cimgt.StamenTerrain()
    crs = ccrs.PlateCarree()
    land_10m = cfeature.NaturalEarthFeature('cultural',
                                            'admin_0_countries',
                                            '10m',
                                            edgecolor=(0, 0, 0, 0),
                                            facecolor=(0, 0, 0, 0), )
    
    if n_variables == 1:
        fig, axs = plt.subplots()
        axs = [axs]
    else:
        ncols = 2
        nrows = n_variables // ncols + n_variables % ncols
        fig = plt.figure(**skwds)
        # fig, axs = plt.subplots(nrows=nrows, ncols=ncols, **skwds)
        # axs = np.ravel(axs)
    for ax_i, variable in enumerate(variables):
        ax = fig.add_subplot(nrows, ncols, ax_i + 1,
                             projection=stamen_terrain.crs)
        geodetic_transform = ccrs.Geodetic()._as_mpl_transform(ax)
        text_transform = offset_copy(geodetic_transform, units='dots', x=-25)
        # station_metadata = load_metadata(variable)
        station_metadata = variable_metadata[variable]
        lons = station_metadata["geoLaenge"]
        lats = station_metadata["geoBreite"]
        station_names = station_metadata["Stationsname"]
        starts = station_metadata["von_datum"]
        ends = station_metadata["bis_datum"]
        if lon_min is None:
            lon_min = lons.min()
        if lat_min is None:
            lat_min = lats.min()
        if lon_max is None:
            lon_max = lons.max()
        if lat_max is None:
            lat_max = lats.max()
        mask = ((lons >= lon_min) & (lons < lon_max) &
                (lats >= lat_min) & (lats < lat_max))
        lons = lons[mask]
        lats = lats[mask]
        station_names = station_names[mask]
        starts = starts[mask]
        ends = ends[mask]
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=crs)
        for i in station_metadata[mask].index:
            text = f"{station_names[i]}\n{starts[i].year}-{ends[i].year}"
            ax.text(lons[i], lats[i], text, fontsize=6,
                    transform=text_transform)
        ax.scatter(lons, lats, transform=crs)
        red_station_names = np.intersect1d(station_names.values,
                                           station_names_inter)
        ii = np.where(station_names.values == red_station_names[:, None])[1]
        ax.scatter(lons.iloc[ii], lats.iloc[ii], facecolor="red")
        ax.set_title(variable)
        ax.add_feature(land_10m, alpha=.1)
        ax.add_image(stamen_terrain, 10)
    # for ax in axs[n_variables:]:
    #     ax.set_axis_off()
    fig.tight_layout()
    return fig


# load_station("Friedrichshafen-Unterraderach", "relative_humidity")
# varnames = ("air_temperature", "sun", "precipitation")
# fig, axs = map_stations(varnames,
#                         llcrnrlon=8.6, llcrnrlat=47.4,
#                         urcrnrlon=10., urcrnrlat=48.0,
#                         figsize=(6, 4))
# plt.show()

if __name__ == '__main__':
    _load_station_one_var("Konstanz", "wind")
    start = datetime(1980, 1, 1)
    end = datetime(2016, 12, 31)
    # df = map_stations(["wind", "air_temperature", "sun",
    #                    "precipitation"],
    #                   lon_min=7, lat_min=47.4,
    #                   lon_max=12., lat_max=49.0,
    #                   start=start, end=end)
    df = map_stations_var(["wind", "air_temperature", "sun",
                           "precipitation"],
                          lon_min=7, lat_min=47.4,
                          lon_max=12., lat_max=49.0)

    plt.show()
