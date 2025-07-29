# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 16:54:40 2025

@author: felix
"""

import pandas as pd
from functions import mean_data
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("whitegrid")

#importing splicing_secretion dataframes
splicing_secretion = pd.read_excel('Results/splicing_secretion.xlsx', index_col = [1,2,3,4,5,6,7,8,9,10])
splicing_secretion = splicing_secretion.drop('Unnamed: 0', axis = 1)




# Create a dataframe with average PSI of all secretory related splice events ofer the tissues
splicing_secretion_lolliplot = mean_data(splicing_secretion)
splicing_secretion_lolliplot = splicing_secretion_lolliplot.drop('PSI_lymphnode_1', axis = 1)

#creating and delta column for sorting by the maximum difference
splicing_secretion_lolliplot['max'] = splicing_secretion_lolliplot.max(axis=1)
splicing_secretion_lolliplot['min'] = splicing_secretion_lolliplot.min(axis=1)
splicing_secretion_lolliplot['delta'] = splicing_secretion_lolliplot['max'] - splicing_secretion_lolliplot['min']
splicing_secretion_lolliplot = splicing_secretion_lolliplot.sort_values(by = 'delta', ascending = False)
#filtering by delta PSI of >0.5, <-0.5
splicing_secretion_lolliplot = splicing_secretion_lolliplot.query('delta > .5 | delta < -.5')

splicing_secretion_lolliplot = splicing_secretion_lolliplot.rename(columns={"PSI_kidney_": "kidney", "PSI_lung_": "lung"
                , 'PSI_skin_': 'skin', 'PSI_testis_' : 'testis', 'PSI_urinarybladder_' : 'urinarybladder'
                , 'PSI_spleen_' : 'spleen', 'PSI_Liver_' : 'liver', 'PSI_smallintestine_' : 'smallintestine'
                , 'PSI_fat_' : 'fat', 'PSI_thyroid_' : 'thyroid', 'PSI_esophagus_' : 'esophagus'
                , 'PSI_brain_' : 'brain', 'PSI_heart_' : 'heart', 'PSI_lymphnode_' : 'lymphnode'
                , 'PSI_placenta_' : 'placenta', 'PSI_stomach_' : 'stomach', 'PSI_Bonemarrow_' : 'bonemarrow'
                , 'PSI_colon_' : 'colon'})
splicing_secretion_lolliplot = splicing_secretion_lolliplot.fillna(0)
splicing_secretion_lolliplot = splicing_secretion_lolliplot.reset_index()
splicing_secretion_lolliplot = splicing_secretion_lolliplot.drop_duplicates(subset='coord_1')
splicing_secretion_lolliplot['coords'] =  splicing_secretion_lolliplot['geneSymbol'] + '  ' + splicing_secretion_lolliplot['coord_chr'] + ': ' + splicing_secretion_lolliplot['coord_1'].astype(str) + '-'  + splicing_secretion_lolliplot['coord_2'].astype(str)
splicing_secretion_lolliplot = splicing_secretion_lolliplot.drop(['geneSymbol', 'strand', 'coord_chr', 'coord_1', 'coord_2',
                                                                  'coord_3', 'coord_4', 'coord_5', 'coord_6', 
                                                                  'Type', 'liver', 'lymphnode', 'brain', 'bonemarrow', 
                                                                  'testis', 'kidney', 'skin', 'stomach', 'colon',
                                                                 'placenta', 'spleen', 'esophagus', 'thyroid', 'lung',
                                                                  'urinarybladder', 'fat', 'smallintestine', 'heart'
                                                                 ],
                                                                 axis = 1)
 
# Reorder it following the values of the first value:
splicing_secretion_lolliplot_ordered = splicing_secretion_lolliplot.sort_values(by='delta')
my_range=range(1,len(splicing_secretion_lolliplot.index)+1)
 
palette_generator = sns.color_palette("PiYG_r", 10)
palette_generator2 =  sns.color_palette("cividis", 5)
# The horizontal plot is made using the hline function
plt.figure(figsize=(8,6))
plt.hlines(y=my_range, xmin=splicing_secretion_lolliplot_ordered['min'],
           xmax=splicing_secretion_lolliplot_ordered['max']
           , color='grey'
           ,alpha=0.8, linewidths = 3)

plt.scatter(splicing_secretion_lolliplot_ordered['min'], my_range, color=palette_generator2[4], alpha=1, label='Delta PSI'
           ,s = 100)
plt.scatter(splicing_secretion_lolliplot_ordered['min'], my_range, color=palette_generator[0], alpha=1, label='Min PSI'
           ,s = 100)
plt.scatter(splicing_secretion_lolliplot_ordered['max'], my_range, color=palette_generator[9], alpha=1 , label='Max PSI'
           ,s = 100)
plt.legend()
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
# Add title and axis names
plt.yticks(my_range, splicing_secretion_lolliplot_ordered['coords'])
plt.xlabel('PSI')

plt.yticks(fontsize=14)
plt.xticks(fontsize=13)


# Show the graph
plt.savefig("Figures/Tissue_secretion_splicing_MINMAX.png", dpi = 600, bbox_inches = "tight")




# Reorder it following the values of the first value:
splicing_secretion_lolliplot_ordered = splicing_secretion_lolliplot.sort_values(by='delta')
my_range=range(1,len(splicing_secretion_lolliplot.index)+1)
 
palette_generator = sns.color_palette("PiYG_r", 10)
palette_generator2 =  sns.color_palette("cividis", 5)
# The horizontal plot is made using the hline function
plt.figure(figsize=(3,6))
#plt.legend()
   
plt.scatter(splicing_secretion_lolliplot_ordered['delta'], my_range, color=palette_generator2[4], alpha=1, label='Delta PSI'
           ,s = 100)
# Add title and axis names
plt.yticks(my_range, splicing_secretion_lolliplot_ordered['coords'])
plt.xticks((0.5, 0.75, 1.0))
#plt.title("Comparison of the lowest and the highest PSI per exon", loc='left')
plt.xlabel('Delta PSI')
plt.xlim(0.4,1.05)
#plt.ylabel('coords')

plt.yticks(fontsize=14)
plt.xticks(fontsize=13)

# Show the graph
plt.savefig("Figures/Tissue_secretion_splicing_DELTA.png", dpi = 600, bbox_inches = "tight")