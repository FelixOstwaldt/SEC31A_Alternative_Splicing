# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 22:00:58 2025

@author: felix
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("whitegrid")
from functions import statistic_f
from functions import cohens_d
from functions import statistic_asterisk


"""   Caco2 differentiation and MO treatment   """

splicing_Caco2_MO_diff = pd.DataFrame({'PSI_E24c' : [0.151119692, 0.150565958, 0.13446239,
                                                     0.542274688, 0.605173626, 0.603249424,
                                                     0.066010369, 0.043055281, 0.073993873],
            
                                       'treatment' : ['MOCtrl', 'MOCtrl', 'MOCtrl',
                                                     'diff. MOCtrl', 'diff. MOCtrl', 'diff. MOCtrl',
                                                     'diff. MO24c', 'diff. MO24c', 'diff. MO24c']
                                                      })


palette_generation = sns.color_palette("RdGy_r", 9)
DATA = splicing_Caco2_MO_diff
y_limit = .70
x_data = 'treatment'
y_data = 'PSI_E24c'
group1 = DATA.query('treatment == "MOCtrl"')[y_data]
group2 = DATA.query('treatment == "diff. MOCtrl"')[y_data]
group3 = DATA.query('treatment == "diff. MO24c"')[y_data]

size = 6
plt.figure(figsize=(size * 0.666666, size))

plot = sns.swarmplot(data = DATA, x = x_data, y = y_data,  alpha=1, color='black')
plot = sns.barplot(data= DATA, x = x_data, y = y_data
            , palette = [palette_generation[1],
                        palette_generation[2],
                        palette_generation[7]]
            ,capsize = 0.2,             
            saturation = 8,             
            errcolor = 'black', errwidth = 2,
            linewidth=2, edgecolor="0")

#statistics1
p_ttest1 = statistic_f(group1, group2)
cohensd1 = cohens_d(group1, group2)

p_ttest2 = statistic_f(group2, group3)
cohensd2 = cohens_d(group2, group3)


#statistics depiction with asterisks
asterisk_text_1 = statistic_asterisk(p_ttest1)
x1, x2 = [0, 0.9]   
y, h, col = y_limit-(y_limit*0.1), y_limit*0.015, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=2, c=col)
plt.text((x1+x2)*.5, y, asterisk_text_1, ha='center', va='bottom', color=col, size = 20)

#statistics depiction with asterisks
asterisk_text_2 = statistic_asterisk(p_ttest2)
x1, x2 = [1.1, 2]   
y, h, col = y_limit-(y_limit*0.1), y_limit*0.015, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=2, c=col)
plt.text((x1+x2)*.5, y, asterisk_text_2, ha='center', va='bottom', color=col, size = 20)


#print(stat_text)
print('Pvalue1: ' + str(p_ttest1))
print('Cohen´s d 1: ' + str(cohensd1))

print('Pvalue2: ' + str(p_ttest2))
print('Cohen´s d 2: ' + str(cohensd2))

plot.set(xlabel=None)


plot.set(xlabel=None)
plt.yticks(fontsize=14)
plt.xticks(fontsize=13)

plot.set_ylabel("PSI exon 24c", fontsize=15)
plt.savefig("Figures/E24c_splicing_Caco2diff_MO.png", dpi=600, bbox_inches='tight')

Exon24c_splicing_Caco2_MO_diff = pd.DataFrame({'Y_data' : [0.004693891,
                                                        0.003400588,
                                                        0.003594483,
                                                        0.013461659,
                                                        0.022097087,
                                                        0.015625
                                                     ],
                                       'treatment' : [
                                                      'undiff.', 'undiff.', 'undiff.',
                                                     
                                                     'diff.', 'diff.', 'diff.'
                                                     ]
                                                      })


palette_generation = sns.color_palette("RdGy_r", 9)
DATA = Exon24c_splicing_Caco2_MO_diff
y_limit = 0.025
x_data = 'treatment'
y_data = 'Y_data'
group1 = DATA.query('treatment == "undiff."')[y_data]
group2 = DATA.query('treatment == "diff."')[y_data]


size = 5.5
plt.figure(figsize=(size * 0.45, size))

plot = sns.swarmplot(data = DATA, x = x_data, y = y_data,  alpha=1, color='black')
plot = sns.barplot(data= DATA, x = x_data, y = y_data
            , palette = [palette_generation[1], palette_generation[2]]
            ,capsize = 0.2,             
            saturation = 8,             
            errcolor = 'black', errwidth = 2,
            linewidth=2, edgecolor="0")

#statistics1
p_ttest1 = statistic_f(group1, group2)
cohensd1 = cohens_d(group1, group2)


#statistics depiction with asterisks
asterisk_text_1 = statistic_asterisk(p_ttest1)
x1, x2 = [0, 1]   
y, h, col = y_limit-(y_limit*0.1), y_limit*0.015, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=2, c=col)
plt.text((x1+x2)*.5, y, asterisk_text_1, ha='center', va='bottom', color=col, size = 20)

#print(stat_text)
print('Pvalue1: ' + str(p_ttest1))
print('Cohen´s d 1: ' + str(cohensd1))

plot.set(xlabel=None)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)

plot.set_ylabel("RBM47 gene expression (rel. to GAPDH)", fontsize=16)
plt.savefig("Figures/RBM47_GE_Caco2diff.png", dpi=600, bbox_inches='tight')


"""  Caco-2 Lipid transport assay   """

NBD_basal_ratio_Caco2_MO_diff = pd.DataFrame({'y_data' : [0.387943593,
                                                        0.418247793,
                                                        0.473844901,
                                                        0.242668238,
                                                        0.173260996,
                                                        0.253714371
                                                     ],
                                       'treatment' : [#'MOCtrl', 'MOCtrl', 'MOCtrl',
                                                      'MOCtrl', 'MOCtrl', 'MOCtrl',  
                                                     #'MO24c', 'MO24c', 'MO24c', 
                                                      'MO24c', 'MO24c',  'MO24c'
                                                     ]
                                                      })


palette_generation = sns.color_palette("RdGy_r", 9)
DATA = NBD_basal_ratio_Caco2_MO_diff
y_limit = 0.55
x_data = 'treatment'
y_data = 'y_data'
group1 = DATA.query('treatment == "MOCtrl"')[y_data]
group2 = DATA.query('treatment == "MO24c"')[y_data]


size = 5.5
plt.figure(figsize=(size * 0.45, size))

plot = sns.swarmplot(data = DATA, x = x_data, y = y_data,  alpha=1, color='black')
plot = sns.barplot(data= DATA, x = x_data, y = y_data
            , palette = [palette_generation[2],
                        palette_generation[7]]
            ,capsize = 0.2,             
            saturation = 8,             
            errcolor = 'black', errwidth = 2,
            linewidth=2, edgecolor="0")

#statistics1
p_ttest1 = statistic_f(group1, group2)
cohensd1 = cohens_d(group1, group2)


#statistics depiction with asterisks
asterisk_text_1 = statistic_asterisk(p_ttest1)
x1, x2 = [0, 1]   
y, h, col = y_limit-(y_limit*0.1), y_limit*0.015, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=2, c=col)
plt.text((x1+x2)*.5, y, asterisk_text_1, ha='center', va='bottom', color=col, size = 20)

print('Pvalue1: ' + str(p_ttest1))
print('Cohen´s d 1: ' + str(cohensd1))


plot.set(xlabel=None)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)

plot.set_ylabel("Fluorescence signal ratio (Basal / Apical)", fontsize=15)
plt.savefig("Figures/Lipid_assay_rel_basal_signal.png", dpi=600, bbox_inches='tight')
