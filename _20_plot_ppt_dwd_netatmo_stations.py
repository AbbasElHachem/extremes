
# coding: utf-8

# In[1]:


import shapefile
import matplotlib.pyplot as plt
import numpy as np


# In[16]:


import pandas as pd
path_to_netatmo_coords_df = (r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes'
                             r'\NetAtmo_BW'
                             r'\rain_bw_1hour'
                             r'\netatmo_bw_1hour_coords.csv')
df_c = pd.read_csv(path_to_netatmo_coords_df, sep=';', index_col=0)


plt.ioff()
fig = plt.figure(figsize=(15, 15), dpi=200)

ax = fig.add_subplot(111)

path_to_shpfile = (
    r"X:\exchange\ElHachem\Netatmo\Landesgrenze_ETRS89\Landesgrenze_10000_ETRS89_lon_lat.shp")
shp_de = shapefile.Reader(path_to_shpfile)
# read and plot shapefile (BW or Germany) should be lon lat
for shape_ in shp_de.shapeRecords():
    lon = [i[0] for i in shape_.shape.points[:][::-1]]
    lat = [i[1] for i in shape_.shape.points[:][::-1]]

    ax.scatter(lon, lat, marker='.', c='lightgrey',
               alpha=0.25, s=1)



df_coords = pd.read_csv(
    r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\DWD_coords_BW.csv", sep=',', index_col=0)


ax.scatter(df_c[' lon'], df_c[' lat'], c='b', alpha=0.85,
           marker='o', s=25, label='Netatmo stations')
ax.scatter(df_coords['lon'], df_coords['lat'], c='r', alpha=0.85,
           s=25,  marker='d', label='DWD stations')
plt.axis('equal')
plt.grid(alpha=.05)
plt.legend(loc=0, fontsize=12)
plt.xlabel('Longitude', fontsize=10)
plt.ylabel('Latitude', fontsize=10)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.tight_layout()
plt.savefig(r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\stations.png",
            frameon=True, papertype='a4',
            bbox_inches='tight', pad_inches=.2)


# In[7]:


plt.savefig(r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\stations.png",
            frameon=True, papertype='a4',
            bbox_inches='tight', pad_inches=.2)

