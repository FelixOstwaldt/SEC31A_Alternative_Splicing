# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 16:57:52 2025

@author: felix
"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("whitegrid")

#importing splicing_secretion dataframes
splicing_secretion = pd.read_excel('Results/splicing_secretion.xlsx', index_col = [1,2,3,4,5,6,7,8,9,10])
splicing_secretion = splicing_secretion.drop('Unnamed: 0', axis = 1)


#extract SEC31A exon 24c of the global splicing list
splicing_E24c  = splicing_secretion.loc['SEC31A', '-', 'chr4', 82830936, 82830975, 82828999, 82829058]
splicing_E24c['exon'] = ['E24c_1', 'E24c_3', 'PSI SEC31A Exon 24c']
splicing_E24c = splicing_E24c.reset_index()
splicing_E24c = splicing_E24c.set_index('exon')
splicing_E24c = splicing_E24c.drop(['coord_5', 'coord_6', 'Type'], axis = 1)
splicing_E24c = splicing_E24c.transpose()
splicing_E24c = splicing_E24c.reset_index()
splicing_E24c.columns = ['tissue','E24c_1', 'E24c_3', 'PSI SEC31A Exon 24c']

# remove names from sinlge samples, and replace it with an index of the tissues
string_list = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'PSI_', '_']
for x in string_list:
    splicing_E24c['tissue'] = splicing_E24c['tissue'].str.replace(x, '')

# homogenize the index    
splicing_E24c = splicing_E24c.replace('Liver', 'liver')
splicing_E24c = splicing_E24c.replace('Bonemarrow', 'bone marrow')
splicing_E24c = splicing_E24c.replace('smallintestine', 'small intestine')
splicing_E24c = splicing_E24c.replace('urinarybladder', 'urinary bladder')
splicing_E24c = splicing_E24c.replace('lymphnode', 'lymph node')

splicing_E24c = splicing_E24c.set_index('tissue')

#create column to sort the samples from mean maximum to mean minimum
list_PSI24c_mean = []
for i in splicing_E24c.transpose():
    list_PSI24c_mean.append(np.mean(splicing_E24c.loc[i]['PSI SEC31A Exon 24c']))
splicing_E24c['PSI24c_mean'] = list_PSI24c_mean
splicing_E24c = splicing_E24c.sort_values(by = 'PSI24c_mean', ascending = False)

splicing_E24c.to_excel('Results/splicing_E24c.xlsx')

# plot the bar plot
size = 4
plt.figure(figsize=(size * 1.38, size))

palette_generation = sns.color_palette("hls", 22)
plot = sns.barplot(data=splicing_E24c.reset_index()
            , x="tissue", y="PSI SEC31A Exon 24c", ci = "sd", palette=palette_generation
           , linewidth=2, edgecolor="0"
           , errorbar=("pi", 50), capsize=.3, errcolor="0")
plt.tick_params(axis='x', rotation=90)

plot.set(xlabel=None)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)

plt.ylabel('PSI SEC31A Exon 24c', fontsize=15)

plt.savefig("Figures/SEC31A_E24c_splicing.png", dpi = 600, bbox_inches = "tight")