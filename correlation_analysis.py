# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 17:11:10 2025

@author: felix
"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("whitegrid")
from scipy.stats import linregress
from functions import statistic_f
from functions import cohens_d
from functions import statistic_asterisk


#importing splicing_secretion dataframes
splicing_secretion = pd.read_excel('Results/splicing_secretion.xlsx', index_col = [1,2,3,4,5,6,7,8,9,10])
splicing_secretion = splicing_secretion.drop('Unnamed: 0', axis = 1)


#preparing exon 24c alternative splicing data to mach the gene expression DataFrame
splicing_E24c  = splicing_secretion.loc['SEC31A', '-', 'chr4', 82830936, 82830975, 82828999, 82829058, 82842184, 82842481]
splicing_E24c = splicing_E24c.drop(['PSI_Bonemarrow_3', 'PSI_Bonemarrow_7', 'PSI_colon_4','PSI_heart_5', 'PSI_lymphnode_7'
                                   ,'PSI_lymphnode_14', 'PSI_testis_3'
                                    , 'PSI_thyroid_6', 'PSI_lung_7']
                                  , axis = 1)
splicing_E24c = splicing_E24c.transpose()#.reset_index().sort_values(by = 'index').set_index('index')
splicing_E24c.columns = ['SEC31A Exon 24c']
splicing_E24c = splicing_E24c.transpose()


#read in DEseq2 DataFrame with gene counts of all tissues 
genecounts_all = pd.read_csv('Data/dataframeall.csv', sep=',', index_col=0)
genecounts_all = genecounts_all.reset_index()
genecounts_all = genecounts_all.drop(['liver4', 'liver5', 'bonemarrow8', 'colon6', 'lung7','heart5', 'lymphnode'
                                      , 'ovary1', 'ovary2', 'ovary3'
                                      , 'ovary4',  'bonemarrow9', 'testis3'
                                     ], axis = 1)
genecounts_all = genecounts_all.set_index(['index'])
genecounts_all = genecounts_all.transpose().reset_index().set_index('index').transpose()


#merging SEC31A splicing PSI with gene counts
splicing_E24c.columns = genecounts_all.columns
genecounts_all_Exon24 = pd.concat([splicing_E24c, genecounts_all])


#calculating the pearson correlation between the PSI and Exon 24a and all other genes in all tissues
GE_all = genecounts_all_Exon24
corr_list = []
for i in GE_all.reset_index()['index']:
    Corr_df = pd.DataFrame([genecounts_all_Exon24.loc['SEC31A Exon 24c'], genecounts_all_Exon24.loc[f'{i}']])
    #print(Corr_df.transpose().corr()[f'A3SS'][f'{i}'])
    corr_list.append(Corr_df.transpose().corr(method='pearson')['SEC31A Exon 24c'][f'{i}'])
    
corr_list[0] = 1.0
Correlation = pd.DataFrame(corr_list)
Correlation['Gene_ID'] = GE_all.reset_index()['index']
Correlation.columns = ['pearson_24c', 'GENE_ID']
Correlation = Correlation.set_index('GENE_ID')


#calculating the spearman correlation between the PSI and Exon 24a and all other genes in all tissues
GE_all = genecounts_all_Exon24
corr_list = []
for i in GE_all.reset_index()['index']:
    Corr_df = pd.DataFrame([genecounts_all_Exon24.loc['SEC31A Exon 24c'], genecounts_all_Exon24.loc[f'{i}']])
    #print(Corr_df.transpose().corr()[f'A3SS'][f'{i}'])
    corr_list.append(Corr_df.transpose().corr(method = 'spearman')['SEC31A Exon 24c'][f'{i}'])
 
    
corr_list[0] = 1.0
Correlation['spearman_24c'] = corr_list


Correlation = Correlation.drop('SEC31A Exon 24c')
Correlation.columns = ['Pearson_correlation', 'Spearman_correlation']
Correlation = Correlation.reset_index()
Correlation.to_excel('Results/Correlation.xlsx')


#plotting the correlation with highlighted GO terms

# loading the correltaion list
#Correlation = pd.read_excel('Results/Correlation.xlsx', index_col= 0)
#Correlation = Correlation.sort_index(ascending=True)

#loading various GO temr lists 
Collagen_Chylomicrons = pd.read_excel('Helper_lists/Collagen_Chylomicrons.xlsx')

GOsecretion_human2 = pd.read_excel('Helper_lists/GOsecretion_human2.xlsx')
GOsecretion_human2.columns = ['GENE_ID', 'GO_term']

GO_term_RNA_splicing = pd.read_excel('Helper_lists/GO_term_RNA_splicing.xlsx')
GO_term_RNA_splicing.columns = ['GENE_ID', 'GO_term']

GO_terms = pd.concat([GOsecretion_human2, GO_term_RNA_splicing, Collagen_Chylomicrons])
GO_terms = GO_terms.drop_duplicates(subset = 'GENE_ID')

#merge Correlation list and GO term lists
Correlation = Correlation.merge(GO_terms, on = 'GENE_ID', how = 'outer')
Correlation = Correlation.dropna(subset=['Pearson_correlation', 'Spearman_correlation'])
Correlation = Correlation.fillna('other')

Correlation = Correlation.reset_index().drop('index', axis = 1)

#create various color palettes
palette_generation_2 = sns.color_palette("PiYG_r", 20)
palette3 = sns.color_palette("Dark2", 4)

#plot each GO term in a scatter plot
plt.figure(figsize=(5,5))
plt.scatter(data = Correlation.query('GO_term == "other"'), x = 'Pearson_correlation', y = 'Spearman_correlation'
           , color = 'grey', alpha = .4, label = 'other')
plt.scatter(data = Correlation.query('GO_term == "secretion"'), x = 'Pearson_correlation', y = 'Spearman_correlation'
           , color = palette3[0], alpha = .9, label = 'secretion')
plt.scatter(data = Correlation.query('GO_term == "mRNA splicing"'), x = 'Pearson_correlation', y = 'Spearman_correlation'
           , color = palette3[1], alpha = .9, label = 'mRNA splicing')

plt.scatter(data = Correlation.query('GO_term == "Collagens"'), x = 'Pearson_correlation', y = 'Spearman_correlation'
           , color = palette3[2], alpha = .9, label = 'Collagens')
plt.scatter(data = Correlation.query('GO_term == "Chylomicrons"'), x = 'Pearson_correlation', y = 'Spearman_correlation'
           , color = palette3[3], alpha = .9, label = 'Chylomicrons')

#plot legend
plt.legend(loc='best', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
plt.legend(title='GO terms')

#mark highlighted genes with error and gene name
text = Correlation['GENE_ID']
for i, label in enumerate(text):
    if text[i] == 'RBM47':
        plt.annotate(text[i],  xy = (Correlation['Pearson_correlation'][i]+0.00, Correlation['Spearman_correlation'][i])
                     , xytext = (Correlation['Pearson_correlation'][i]+.1, Correlation['Spearman_correlation'][i]+0.1)
                     , size=15, color = palette3[1],
                    arrowprops=dict(facecolor=palette3[1]
                                    , shrink=0.05
                                    , alpha = .8
                                   ))
    if text[i] == 'MIA2':
        plt.annotate(text[i],  xy = (Correlation['Pearson_correlation'][i]+0.00, Correlation['Spearman_correlation'][i])
                     , xytext = (Correlation['Pearson_correlation'][i]+.1, Correlation['Spearman_correlation'][i]-.15)
                     , size=15, color = palette3[0],
                    arrowprops=dict(facecolor=palette3[0]
                                    , shrink=0.05
                                    , alpha = .8
                                   ))
        
    if text[i] == 'LSR':
        plt.annotate(text[i],  xy = (Correlation['Pearson_correlation'][i]+0.005, Correlation['Spearman_correlation'][i])
                     , xytext = (Correlation['Pearson_correlation'][i]+.1, Correlation['Spearman_correlation'][i]-.15)
                     , size=15, color = palette3[3],
                    arrowprops=dict(facecolor=palette3[3]
                                    , shrink=0.05
                                    , alpha = .8
                                   ))
    if text[i] == 'APOB':
        plt.annotate(text[i],  xy = (Correlation['Pearson_correlation'][i]+0.005, Correlation['Spearman_correlation'][i])
                     , xytext = (Correlation['Pearson_correlation'][i]+.12, Correlation['Spearman_correlation'][i]-.1)
                     , size=15, color = palette3[3],
                    arrowprops=dict(facecolor=palette3[3]
                                    , shrink=0.05
                                    , alpha = .8
                                   ))
    if text[i] == 'APOA1':
        plt.annotate(text[i],  xy = (Correlation['Pearson_correlation'][i]+0.005, Correlation['Spearman_correlation'][i])
                     , xytext = (Correlation['Pearson_correlation'][i]+.12, Correlation['Spearman_correlation'][i]-.1)
                     , size=15, color = palette3[3],
                    arrowprops=dict(facecolor=palette3[3]
                                    , shrink=0.05
                                    , alpha = .8
                                   ))
    if text[i] == 'APOH':
        plt.annotate(text[i],  xy = (Correlation['Pearson_correlation'][i]+0.005, Correlation['Spearman_correlation'][i])
                     , xytext = (Correlation['Pearson_correlation'][i]-.2, Correlation['Spearman_correlation'][i]+.15)
                     , size=15, color = palette3[3],
                    arrowprops=dict(facecolor=palette3[3]
                                    , shrink=0.05
                                    , alpha = .8
                                   ))
        
    if text[i] == 'COL1A1':
        plt.annotate(text[i],  xy = (Correlation['Pearson_correlation'][i]+0.000, Correlation['Spearman_correlation'][i])
                     , xytext = (Correlation['Pearson_correlation'][i]-.5, Correlation['Spearman_correlation'][i]+.2)
                     , size=15, color = palette3[2],
                    arrowprops=dict(facecolor=palette3[2]
                                    , shrink=0.05
                                    , alpha = .8
                                   ))
        
    if text[i] == 'COL1A2':
        plt.annotate(text[i],  xy = (Correlation['Pearson_correlation'][i]+0.005, Correlation['Spearman_correlation'][i])
                     , xytext = (Correlation['Pearson_correlation'][i]+.2, Correlation['Spearman_correlation'][i]-.2)
                     , size=15, color = palette3[2],
                    arrowprops=dict(facecolor=palette3[2]
                                    , shrink=0.05
                                    , alpha = .8
                                   ))
        
    if text[i] == 'COL2A1':
        plt.annotate(text[i],  xy = (Correlation['Pearson_correlation'][i]+0.005, Correlation['Spearman_correlation'][i])
                     , xytext = (Correlation['Pearson_correlation'][i]-.5, Correlation['Spearman_correlation'][i]+.2)
                     , size=15, color = palette3[2],
                    arrowprops=dict(facecolor=palette3[2]
                                    , shrink=0.05
                                    , alpha = .8
                                   ))

plt.xlabel("Pearson Correlation", fontsize = 16)
plt.ylabel("Spearman Correlation", fontsize = 16)

plt.yticks(fontsize=15)
plt.xticks(fontsize=15)

plt.savefig("Figures/SEC31A_E24c_GE_correlation_2.png", dpi = 600, bbox_inches = "tight")


#preparing Splicing data and gene epxression data to be ploted in linear regression
genecounts_all_Exon24_transpose = genecounts_all_Exon24.transpose().reset_index()
genecounts_all_Exon24_transpose['index'] = genecounts_all_Exon24_transpose['index'].str.replace('1', '')
genecounts_all_Exon24_transpose['index'] = genecounts_all_Exon24_transpose['index'].str.replace('2', '')
genecounts_all_Exon24_transpose['index'] = genecounts_all_Exon24_transpose['index'].str.replace('3', '')
genecounts_all_Exon24_transpose['index'] = genecounts_all_Exon24_transpose['index'].str.replace('4', '')
genecounts_all_Exon24_transpose['index'] = genecounts_all_Exon24_transpose['index'].str.replace('5', '')
genecounts_all_Exon24_transpose['index'] = genecounts_all_Exon24_transpose['index'].str.replace('6', '')
genecounts_all_Exon24_transpose['index'] = genecounts_all_Exon24_transpose['index'].str.replace('7', '')
genecounts_all_Exon24_transpose['index'] = genecounts_all_Exon24_transpose['index'].str.replace('8', '')
genecounts_all_Exon24_transpose['index'] = genecounts_all_Exon24_transpose['index'].str.replace('9', '')
genecounts_all_Exon24_transpose['index'] = genecounts_all_Exon24_transpose['index'].str.replace('0', '')

genecounts_all_Exon24_transpose['index'] = genecounts_all_Exon24_transpose['index'].replace('Bonemarrow', 'bone marrow')
genecounts_all_Exon24_transpose['index'] = genecounts_all_Exon24_transpose['index'].replace('smallintestine', 'small intestine')
genecounts_all_Exon24_transpose['index'] = genecounts_all_Exon24_transpose['index'].replace('urinarybladder', 'urinary bladder')
genecounts_all_Exon24_transpose['index'] = genecounts_all_Exon24_transpose['index'].replace('lymphnode', 'lymph node')

list_PSI24c_mean = []
for i in genecounts_all_Exon24_transpose.set_index('index').transpose():
    list_PSI24c_mean.append(np.mean(genecounts_all_Exon24_transpose.set_index('index').loc[i]['SEC31A Exon 24c']))
genecounts_all_Exon24_transpose['PSI24c_mean'] = list_PSI24c_mean
genecounts_all_Exon24_transpose = genecounts_all_Exon24_transpose.sort_values(by = 'PSI24c_mean', ascending = False)


#plot linear regression
size = 5
plt.figure(figsize=(size * 1.38, size))

x = 'SEC31A Exon 24c'
y = 'MIA2'
sns.regplot(data = genecounts_all_Exon24_transpose, x = f'{x}', y = f'{y}', marker='', ci = 0, color = 'gray')
sns.scatterplot(data = genecounts_all_Exon24_transpose, x = f'{x}', y = f'{y}', hue = 'index',  style = 'index' ,
                palette = sns.color_palette("hls", 22)).set(xlabel = 'PSI SEC31A Exon 24c', 
                                                            ylabel = 'MIA2 (normalized gene counts)')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

reg_stats = linregress(genecounts_all_Exon24_transpose.dropna()[f'{x}']
                 , genecounts_all_Exon24_transpose.dropna()[f'{y}'])

plt.xlabel("PSI SEC31A Exon 24c", fontsize = 15)
plt.ylabel("MIA2 \n(normalized gene counts)", fontsize = 15)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)

print(reg_stats)                                 
plt.savefig("Figures/24c_MIA2_corr.png", dpi = 600, bbox_inches = "tight")


#plot linear regression
x = 'SEC31A Exon 24c'
y = 'RBM47'

size = 3.5
plt.figure(figsize=(size * 1.38, size))

sns.regplot(data = genecounts_all_Exon24_transpose, x = f'{x}', y = f'{y}', marker="", ci = 0, color = 'gray')
sns.scatterplot(data = genecounts_all_Exon24_transpose, x = f'{x}', y = f'{y}', hue = 'index',  style = 'index' ,
                palette = sns.color_palette("hls", 22)).set(xlabel = 'PSI SEC31A Exon 24c', 
                                                            ylabel = 'RBM47 (normalized gene counts)')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

reg_stats = linregress(genecounts_all_Exon24_transpose.dropna()[f'{x}']
                 , genecounts_all_Exon24_transpose.dropna()[f'{y}'])

plt.xlabel("PSI SEC31A Exon 24c", fontsize = 15)
plt.ylabel("RBM47 \n(normalized gene counts)", fontsize = 15)

plt.yticks(fontsize=14)
plt.xticks(fontsize=14)

print(reg_stats)                                 
plt.savefig("Figures/24c_RBM47_corr.png", dpi = 600, bbox_inches = "tight")


y_limit = .50
DATA = Correlation.query('GO_term !="other"').query('GO_term != "secretion"').query('GO_term != "mRNA splicing"')

x_data = 'GO_term'
y_data = 'Pearson_correlation'
group1 = DATA.query('GO_term == "Chylomicrons"')[y_data]
group2 = DATA.query('GO_term == "Collagens"')[y_data]


size = 5
plt.figure(figsize=(size * 0.45, size))
plot = sns.swarmplot(data = DATA, x = x_data, y = y_data,  alpha=1, color='black')
plot = sns.boxplot(data= DATA, x = x_data, y = y_data, palette = ['#7fbc41', '#de77ae']
           # ,capsize = 0.2,             
           # saturation = 8,             
           # errcolor = 'black', errwidth = 2,
           # linewidth=2, edgecolor="0"
                  )

#statistics
p_test = statistic_f(group1, group2)

#calculating Cohen´s d
cohensd = cohens_d(group1, group2)

#statistics depiction with asterisks
asterisk_text = statistic_asterisk(p_test)
    
x1, x2 = 0, 1   
y, h, col = y_limit-(y_limit*0.06), y_limit*0.015, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=2, c=col)
plt.text((x1+x2)*.5, y, asterisk_text, ha='center', va='bottom', color=col, size = 20)
plot.set_ylabel("Pearson correlation index \nto exon 24c ", fontsize=15)

plt.tick_params(axis='x', rotation=20)

plot.set(xlabel=None)
plt.yticks(fontsize=14)
plt.xticks(fontsize=12)

print('Pvalue: ' + str(p_test))
print('Cohen´s d: ' + str(cohensd))
plt.savefig("Figures/E24c_corr_APO_COL.png", dpi=600, bbox_inches='tight')
