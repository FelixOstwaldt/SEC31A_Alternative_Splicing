# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 16:41:15 2025

@author: felix
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("whitegrid")

#importing splicing and splicing_secretion dataframes
splicing = pd.read_excel('Results/splicing.xlsx', index_col = [1,2,3,4,5,6,7,8,9,10])
splicing = splicing.drop('Unnamed: 0', axis = 1)
splicing_secretion = pd.read_excel('Results/splicing_secretion.xlsx', index_col = [1,2,3,4,5,6,7,8,9,10])
splicing_secretion = splicing_secretion.drop('Unnamed: 0', axis = 1)
splicing = splicing.fillna(0)



#preparing and transposing splicing data for PCA 
plot = splicing.loc[(splicing).any(axis=1)]
plot_transposed = plot.transpose()
plot_transposed.reset_index(inplace=True)

# Separating out the tarplott
y = plot_transposed.loc[:,['index']].values

plot_transposed = plot_transposed.drop(['index'], axis=1)
features = list(plot_transposed.columns)
x = plot_transposed.loc[:, features].values

#Standardizing the features
x = StandardScaler().fit_transform(x)

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2'])

plot_transposed = plot.transpose()
plot_transposed.reset_index(inplace=True)

finalDf = pd.concat([principalDf, plot_transposed[['index']]], axis = 1)

finalDf.columns = ['principal component 1', 'principal component 2', 'index']
finalDf


# Finalising data and correct index
finalDf_split = finalDf['index'].str.split('_', n = 1, expand = True)
finalDf_split.columns = [ 'b' , 'index']
finalDf_split = finalDf_split['index'].str.split('_', n = 1, expand = True)
finalDf_split.columns = ['index', 'b']

#Add a column PSI24c_mean, as we wanted to preserve the color code for the tissues in all plots
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
finalDf['PSI24c_mean'] = list_PSI24c_mean
finalDf = finalDf.sort_values(by = 'PSI24c_mean', ascending = False)

#correct the index
finalDf['index'] = finalDf_split['index']
finalDf = finalDf.replace('Liver', 'liver')
finalDf = finalDf.replace('Bonemarrow', 'bone marrow')
finalDf = finalDf.replace('smallintestine', 'small intestine')
finalDf = finalDf.replace('urinarybladder', 'urinary bladder')
finalDf = finalDf.replace('lymphnode', 'lymph node')
finalDf


#draw PCA plot
sns.scatterplot(data = finalDf, x = 'principal component 1', y = 'principal component 2', hue = 'index' 
                , palette = sns.color_palette("hls", 22)
               , style = 'index')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.xlim([-25, 20])
plt.ylim([-40, 150])
plt.savefig("Figures/TissueData_splicing_PCA.png", dpi = 600, bbox_inches = "tight")


#draw zoomed in PCA plot
sns.scatterplot(data = finalDf, x = 'principal component 1', y = 'principal component 2', hue = 'index' 
                , palette = sns.color_palette("hls", 22)
               , style = 'index')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.xlim([20, 5])
plt.ylim([-20, 15])
plt.savefig("Figures/TissueData_splicing_PCA_zoom.png", dpi = 600, bbox_inches = "tight")