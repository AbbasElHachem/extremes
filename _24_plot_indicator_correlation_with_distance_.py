
# coding: utf-8

# In[6]:


import pandas as pd
import os

import matplotlib.pyplot as plt


# In[2]:


df0 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_p5_freq_60min_dwd_netatmo_upper_98_percent_data_considered_neighbor_0_.csv"
df1 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_p5_freq_60min_dwd_netatmo_upper_98_percent_data_considered_neighbor_1_.csv"
df2 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_p5_freq_60min_dwd_netatmo_upper_98_percent_data_considered_neighbor_2_.csv"
df3 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_p5_freq_60min_dwd_netatmo_upper_98_percent_data_considered_neighbor_3_.csv"
df4 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature\year_allyears_df_comparing_correlations_p5_freq_60min_dwd_netatmo_upper_98_percent_data_considered_neighbor_4_.csv"
#df0 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_indicator_correlations_DWD_DWD_ppt_\year_allyears_df_dwd_correlationsfreq_60min_dwd_netatmo_upper_95_percent_data_considered_neighbor_1_.csv"
#df1 = r"X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_indicator_correlations_DWD_DWD_ppt_\year_allyears_df_dwd_correlationsfreq_60min_dwd_netatmo_upper_95_percent_data_considered_neighbor_2_.csv"


# In[3]:


in_df0 = pd.read_csv(df0, index_col=0, sep=';')
in_df1 = pd.read_csv(df1, index_col=0, sep=';')
in_df2 = pd.read_csv(df2, index_col=0, sep=';')
in_df3 = pd.read_csv(df3, index_col=0, sep=';')
in_df4 = pd.read_csv(df4, index_col=0, sep=';')


# In[4]:


dfs_all = pd.DataFrame(index=in_df1.index)


# In[5]:


# in_df0.head(5)


# In[6]:


dfs_all.loc[in_df0.index, 'distance0'] = in_df0.loc[:, 'Distance to neighbor']
dfs_all.loc[in_df1.index, 'distance1'] = in_df1.loc[:, 'Distance to neighbor']
dfs_all.loc[in_df2.index, 'distance2'] = in_df2.loc[:, 'Distance to neighbor']
dfs_all.loc[in_df3.index, 'distance3'] = in_df3.loc[:, 'Distance to neighbor']
dfs_all.loc[in_df4.index, 'distance4'] = in_df4.loc[:, 'Distance to neighbor']

dfs_all.loc[in_df0.index, 'Indic_corr0'] = in_df0.loc[:, 'Bool_Spearman_Correlation']
dfs_all.loc[in_df1.index, 'Indic_corr1'] = in_df1.loc[:, 'Bool_Spearman_Correlation']
dfs_all.loc[in_df2.index, 'Indic_corr2'] = in_df2.loc[:, 'Bool_Spearman_Correlation']
dfs_all.loc[in_df3.index, 'Indic_corr3'] = in_df3.loc[:, 'Bool_Spearman_Correlation']
dfs_all.loc[in_df4.index, 'Indic_corr4'] = in_df4.loc[:, 'Bool_Spearman_Correlation']


# In[11]:


plt.ioff()
plt.figure(figsize=(12,8), dpi=300)
plt.scatter(dfs_all.distance0, dfs_all.Indic_corr0, c='r', alpha=0.5, marker='x', label='First Neighbor', s=28)
plt.scatter(dfs_all.distance1, dfs_all.Indic_corr1, c='orange', alpha=0.5, marker='.', label='Second Neighbor', s=28)
plt.scatter(dfs_all.distance2, dfs_all.Indic_corr2, c='g', alpha=0.5, marker='d', label='Third Neighbor', s=28)
plt.scatter(dfs_all.distance3, dfs_all.Indic_corr3, c='b', alpha=0.5, marker='1', label='Fourth Neighbor', s=28)
plt.scatter(dfs_all.distance4, dfs_all.Indic_corr4, c='m', alpha=0.5, marker='2', label='Fifth Neighbor', s=28)

#plt.scatter([],[], label='Number of stations %d' % dfs_all.shape[0], marker='.')
plt.xlim([0, max(dfs_all.distance3.max(),dfs_all.distance4.max()) + 1000])
plt.ylim([-0.1,1.1])
plt.xlabel('Distance (m)')
plt.ylabel('Indicator Spearman Correlation')
plt.legend(loc=0)
plt.grid(alpha=.25)
plt.tight_layout()
plt.title('Netatmo %d stations: Indicator correlation with distance for upper 98percent of data values' % dfs_all.shape[0])
plt.savefig(os.path.join(r'X:\hiwi\ElHachem\Prof_Bardossy\Extremes\plots_NetAtmo_ppt_NetAtmo_temperature', 'Netatmo_98percent_indic_corr.png'),
                frameon=True, papertype='a4',
                bbox_inches='tight', pad_inches=.2)
plt.show()

