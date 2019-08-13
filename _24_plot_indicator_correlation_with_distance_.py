
# coding: utf-8

# In[6]:


import pandas as pd
import os

import matplotlib.pyplot as plt


percent = '95'
time_freq = '60min'
data_source = 'netatmo'
# In[2]:


df0 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_%s_%s_netatmo_upper_%s_percent_data_considered_neighbor_0_.csv" % (
    time_freq, data_source, percent)
df1 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_%s_%s_netatmo_upper_%s_percent_data_considered_neighbor_1_.csv" % (
    time_freq, data_source, percent)
df2 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_%s_%s_netatmo_upper_%s_percent_data_considered_neighbor_2_.csv" % (
    time_freq, data_source, percent)
df3 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_%s_%s_netatmo_upper_%s_percent_data_considered_neighbor_3_.csv" % (
    time_freq, data_source, percent)
# df4 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_max_sep_dist_500000_freq_%s_dwd_netatmo_upper_%s_percent_data_considered_neighbor_4_.csv" % (
#     time_freq, percent)
#df0 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_indicator_correlations_DWD_DWD_ppt_\year_allyears_df_dwd_correlationsfreq_60min_dwd_netatmo_upper_95_percent_data_considered_neighbor_1_.csv"
#df1 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_indicator_correlations_DWD_DWD_ppt_\year_allyears_df_dwd_correlationsfreq_60min_dwd_netatmo_upper_95_percent_data_considered_neighbor_2_.csv"


# In[3]:


in_df0 = pd.read_csv(df0, index_col=0, sep=';').dropna(how='all')
in_df1 = pd.read_csv(df1, index_col=0, sep=';').dropna(how='all')
in_df2 = pd.read_csv(df2, index_col=0, sep=';').dropna(how='all')
in_df3 = pd.read_csv(df3, index_col=0, sep=';').dropna(how='all')
# in_df4 = pd.read_csv(df4, index_col=0, sep=';').dropna(how='all')

cmn_idx = (in_df0.index.intersection(in_df1.index).intersection(
    in_df2.index).intersection(in_df3.index))

# In[4]:


dfs_all = pd.DataFrame(index=cmn_idx)  # in_df0.index)


# In[5]:


# in_df0.head(5)


# In[6]:

x0 = in_df0.loc[:, 'Distance to neighbor'].values.ravel()
x1 = in_df1.loc[:, 'Distance to neighbor'].values.ravel()
x2 = in_df2.loc[:, 'Distance to neighbor'].values.ravel()
x3 = in_df3.loc[:, 'Distance to neighbor'].values.ravel()
# dfs_all.loc[in_df4.index, 'distance4'] = in_df4.loc[:, 'Distance to neighbor']

y0 = in_df0.loc[:,
                'Bool_Spearman_Correlation'].values.ravel()
y1 = in_df1.loc[:,
                'Bool_Spearman_Correlation'].values.ravel()
y2 = in_df2.loc[:,
                'Bool_Spearman_Correlation'].values.ravel()
y3 = in_df3.loc[:,
                'Bool_Spearman_Correlation'].values.ravel()

# dfs_all.loc[in_df0.index, 'distance0'] = in_df0.loc[:, 'Distance to neighbor']
# dfs_all.loc[in_df1.index, 'distance1'] = in_df1.loc[:, 'Distance to neighbor']
# dfs_all.loc[in_df2.index, 'distance2'] = in_df2.loc[:, 'Distance to neighbor']
# dfs_all.loc[in_df3.index, 'distance3'] = in_df3.loc[:, 'Distance to neighbor']
# # dfs_all.loc[in_df4.index, 'distance4'] = in_df4.loc[:, 'Distance to neighbor']
#
# dfs_all.loc[in_df0.index, 'Indic_corr0'] = in_df0.loc[:,
#                                                       'Bool_Spearman_Correlation']
# dfs_all.loc[in_df1.index, 'Indic_corr1'] = in_df1.loc[:,
#                                                       'Bool_Spearman_Correlation']
# dfs_all.loc[in_df2.index, 'Indic_corr2'] = in_df2.loc[:,
#                                                       'Bool_Spearman_Correlation']
# dfs_all.loc[in_df3.index, 'Indic_corr3'] = in_df3.loc[:,
#                                                       'Bool_Spearman_Correlation']
# dfs_all.loc[in_df4.index, 'Indic_corr4'] = in_df4.loc[:,
#                                                       'Bool_Spearman_Correlation']


# In[11]:


plt.ioff()
plt.figure(figsize=(12, 8), dpi=300)
plt.scatter(x0, y0, c='r',
            alpha=0.5, marker='x', label='First Neighbor', s=28)
plt.scatter(x1, y1, c='b',
            alpha=0.5, marker='.', label='Second Neighbor', s=28)
plt.scatter(x2, y2, c='g',
            alpha=0.5, marker='d', label='Third Neighbor', s=28)
plt.scatter(x3, y3, c='orange',
            alpha=0.5, marker='1', label='Fourth Neighbor', s=28)
# plt.scatter(dfs_all.distance4, dfs_all.Indic_corr4, c='m',
#             alpha=0.5, marker='2', label='Fifth Neighbor', s=28)

#plt.scatter([],[], label='Number of stations %d' % dfs_all.shape[0], marker='.')
plt.xlim([0, max(x3.max(), x2.max()) + 1000])
plt.ylim([-0.1, 1.1])
plt.xlabel('Distance (m)')
plt.ylabel('Indicator Spearman Correlation')
plt.legend(loc=0)
plt.grid(alpha=.25)
plt.tight_layout()
plt.title('Netatmo %s %d stations %s: Indicator correlation with distance for upper %s percent of data values'
          % (data_source, y0.shape[0], time_freq, percent))
plt.savefig(os.path.join(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature',
                         'Netatmo_%s_dwd_percent_indic_corr_freq_%s_filtred.png'
                         % (percent, time_freq)),
            frameon=True, papertype='a4',
            bbox_inches='tight', pad_inches=.2)
# plt.show()
