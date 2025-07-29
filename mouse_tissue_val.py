# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 16:59:46 2025

@author: felix
"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("whitegrid")
from scipy.stats import linregress

#create DataFrame of the data from the quantification of the radioactive gel
E24csplicing_mousetissues = pd.DataFrame({'PSI 24c' : [ 0.703915241, 0.718047615, 0.710920947
                                                        , 0.432203698, 0.386756401, 0.420612286
                                                        , 0.09155929, 0.097797553, 0.151566682
                                                        , 0.040899208, 0.028932168, 0.000193505
                                                        , 0.145476394, 0.183124488, 0.22928863
                                                        , 0.020357807, 0.046647126, 0.016379024
                                                        ],
                                           'tissue' : ['liver', 'liver','liver',
                                                     'kidney', 'kidney', 'kidney',
                                                    'spleen', 'spleen', 'spleen',
                                                    'lung', 'lung', 'lung',
                                                    'heart', 'heart', 'heart',
                                                    'brain', 'brain', 'brain',
                                                    ]
                                         })
E24csplicing_mousetissues = E24csplicing_mousetissues.set_index('tissue')


#create a DataFrame assigning a color for each tissue, to make the color code consistent
palette_generation = sns.color_palette("hls", 22)
palette_dataframe = pd.DataFrame({'tissue' :  ['liver', 'kidney', 'stomach', 'thyroid', 'colon', 'smallintestine',
                                                    'bonemarrow', 'lung', 'placenta', 
                                                   'urinarybladder', 'esophagus','skin',
                                                   'testis', 'lymphnode', 'spleen','fat', 
                                                   'heart', 'brain'],
                                'color' : palette_generation[0:18]})
palette_dataframe = palette_dataframe.set_index('tissue')


# draw barplot
size = 4
plt.figure(figsize=(size * 1.38, size))

palette = [palette_dataframe.loc['liver'][0], palette_dataframe.loc['kidney'][0], palette_dataframe.loc['spleen'][0]
          ,palette_dataframe.loc['lung'][0], palette_dataframe.loc['heart'][0], palette_dataframe.loc['brain'][0]]
plot = sns.barplot(data=E24csplicing_mousetissues.reset_index(), x="tissue", y="PSI 24c", ci = "sd", palette=palette
           , linewidth=2, edgecolor="0"
           , errorbar=("pi", 50), capsize=.3, errcolor="0")
plt.tick_params(axis='x', rotation=90)

plot.set(xlabel=None)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)

plt.ylabel('PSI SEC31A Exon 24c', fontsize=15)
plt.savefig("Figures/SEC31A_E24c_splicing_mouse_validation.png", dpi = 600, bbox_inches = "tight")



#preparing data for linear regression between mouse and human data
tissue_list = ['liver', 'kidney', 'spleen', 'lung', 'heart', 'brain']
mouse_mean = []
human_mean = []
mouse_std = []
human_std = []
#creating lists for the mean and standart deviasion of exon 24c PSI in mouse and human data
splicing_E24c = pd.read_excel('Results/splicing_E24c.xlsx', index_col = 0)
for i in tissue_list:
    mouse_mean.append(np.mean(E24csplicing_mousetissues.loc[i]['PSI 24c']))
    mouse_std.append(np.std(E24csplicing_mousetissues.loc[i]['PSI 24c']))
    human_mean.append(np.mean(splicing_E24c.loc[i]['PSI SEC31A Exon 24c']))
    human_std.append(np.std(splicing_E24c.loc[i]['PSI SEC31A Exon 24c']))

size = 3
plt.figure(figsize=(size * 1.38, size))

# creating DataFrame of the mouse and human data
E24csplicing_validationplot_mouse = pd.DataFrame({'PSI Exon 24c mouse RT-PCR' : mouse_mean,
                                                 'mouse_std' : mouse_std,
                                                'tissue' : tissue_list
                                                 })    
E24csplicing_validationplot_human = pd.DataFrame({'PSI Exon 24c human RNA-Seq data' : human_mean,
                                                 'human_std' : human_std 
                                                 })    
E24csplicing_validationplot = pd.concat([E24csplicing_validationplot_human, E24csplicing_validationplot_mouse]
                                       ,axis = 1)
#plot the y and x error bars in propper color 
for n in range(0,6):
    plt.errorbar(human_mean[n], mouse_mean[n], xerr = human_std[n], yerr=mouse_std[n]
             , fmt="none"
            ,color = palette[n]
             ,capsize = 3
             ,alpha= .5
            )
    
#plot the mean values
sns.scatterplot(data = E24csplicing_validationplot, x = 'PSI Exon 24c human RNA-Seq data'
                , y = 'PSI Exon 24c mouse RT-PCR'
                , hue = 'tissue'
                , palette = palette)

#plot the regression line
sns.regplot(data = E24csplicing_validationplot, x = 'PSI Exon 24c human RNA-Seq data'
                , y = 'PSI Exon 24c mouse RT-PCR'
           , marker = ''
           , ci = 0
           , color = 'gray'
           )

plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)


plt.yticks(fontsize=14)
plt.xticks(fontsize=14)

reg_stats = linregress(E24csplicing_validationplot['PSI Exon 24c human RNA-Seq data']
                 , E24csplicing_validationplot['PSI Exon 24c mouse RT-PCR'])
plt.xlabel('PSI Exon 24c human RNA-Seq data', fontsize=15)
plt.ylabel('PSI Exon 24c \nmouse RT-PCR', fontsize=15)

plt.savefig("Figures/SEC31A_E24c_splicing_mouse_validation_reg.png", dpi = 600, bbox_inches = "tight")