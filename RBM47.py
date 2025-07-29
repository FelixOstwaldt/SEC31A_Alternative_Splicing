# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 21:33:23 2025

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



""" RBM47 KO RNAseq"""

#import RNA sequencing data set of an RBM47 knock out
RBM47_KO_SE = pd.read_excel('Data/CTRLvsRBM47KD/SE_CTRLvsRBM47KD.xlsx',
                           index_col = [2,3,4,5,6,7,8,9,10])

#modulate dataframe
RBM47_KO_E24c = RBM47_KO_SE.loc['SEC31A', 'chr4', '-', 82830936, 82830975, 82828999, 82829058, 82842139]
RBM47_KO_E24c = RBM47_KO_E24c.drop(['Unnamed: 0', 'GeneID', 'PValueCTRLvsRBM47KD', 'DPSI_CTRLvsRBM47KD'], axis = 1)
RBM47_KO_E24c_plot = RBM47_KO_E24c.transpose().reset_index()
RBM47_KO_E24c_plot = RBM47_KO_E24c_plot.replace('PSI_CTRL_1', 'Ctrl')
RBM47_KO_E24c_plot = RBM47_KO_E24c_plot.replace('PSI_CTRL_2', 'Ctrl')
RBM47_KO_E24c_plot = RBM47_KO_E24c_plot.replace('PSI_CTRL_3', 'Ctrl')
RBM47_KO_E24c_plot = RBM47_KO_E24c_plot.replace('PSI_RBM47KD_1', 'RBM47 KD')
RBM47_KO_E24c_plot = RBM47_KO_E24c_plot.replace('PSI_RBM47KD_2', 'RBM47 KD')
RBM47_KO_E24c_plot = RBM47_KO_E24c_plot.replace('PSI_RBM47KD_3', 'RBM47 KD')
RBM47_KO_E24c_plot.columns = ['treatment', 'PSI_E24c']



#plot PSI values of SEC31A exon 24c
y_limit = .40
DATA = RBM47_KO_E24c_plot


x_data = 'treatment'
y_data = 'PSI_E24c'
group1 = RBM47_KO_E24c_plot.query('treatment == "Ctrl"')[y_data]
group2 = RBM47_KO_E24c_plot.query('treatment == "RBM47 KD"')[y_data]

size = 3.5
plt.figure(figsize=(size * 0.45, size))
plot = sns.swarmplot(data = DATA, x = x_data, y = y_data,  alpha=1, color='black')
plot = sns.barplot(data= DATA, x = x_data, y = y_data, palette = ['#7fbc41', '#de77ae']
            ,capsize = 0.2,             
            saturation = 8,             
            errcolor = 'black', errwidth = 2,
            linewidth=2, edgecolor="0")

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
plot.set_ylabel("PSI SEC31A exon 24c\n(predicted)", fontsize=15)

plot.set(xlabel=None)
plt.yticks(fontsize=14)
plt.xticks(fontsize=12)

print('Pvalue: ' + str(p_test))
print('Cohen´s d: ' + str(cohensd))
plt.savefig("Figures/E24c_splicing_RBM47KD_predicted.png", dpi=600, bbox_inches='tight')

#import gene differential expression  
genecounts_all_RBM47KD = pd.read_csv('Data/dataframe_CTRLvsRBM47KD.csv', sep=',', index_col=0)
genecounts_RBM47KD = genecounts_all_RBM47KD.loc['RBM47']
genecounts_RBM47KD = genecounts_RBM47KD.transpose().reset_index()
genecounts_RBM47KD = genecounts_RBM47KD.replace('CTRL_1', 'Ctrl')
genecounts_RBM47KD = genecounts_RBM47KD.replace('CTRL_2', 'Ctrl')
genecounts_RBM47KD = genecounts_RBM47KD.replace('CTRL_3', 'Ctrl')
genecounts_RBM47KD = genecounts_RBM47KD.replace('RBM47KD_1', 'RBM47 KD')
genecounts_RBM47KD = genecounts_RBM47KD.replace('RBM47KD_2', 'RBM47 KD')
genecounts_RBM47KD = genecounts_RBM47KD.replace('RBM47KD_3', 'RBM47 KD')
genecounts_RBM47KD.columns = ['treatment', 'GE_RBM47']

#plot RBM47 gene expression
DATA = genecounts_RBM47KD
y_limit = 11000
x_data = 'treatment'
y_data = 'GE_RBM47'
group1 = DATA.query('treatment == "Ctrl"')[y_data]
group2 = DATA.query('treatment == "RBM47 KD"')[y_data]

size = 3.5
plt.figure(figsize=(size * 0.45, size))
plot = sns.swarmplot(data = DATA, x = x_data, y = y_data,  alpha=1, color='black')
plot = sns.barplot(data= DATA, x = x_data, y = y_data, palette = ['#b8e186', '#f1b6da']
            ,capsize = 0.2,             
            saturation = 8,             
            errcolor = 'black', errwidth = 2,
            linewidth=2, edgecolor="0")
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
plot.set_ylabel("RBM47 gene expression", fontsize=15)

plot.set(xlabel=None)
plt.yticks(fontsize=14)
plt.xticks(fontsize=12)

print('Pvalue: ' + str(p_test))
print('Cohen´s d: ' + str(cohensd))
plt.savefig("Figures/RBM47_GE_RBM47KD_predicted.png", dpi=600, bbox_inches='tight')


""" RBM47 gene expression validation """

#type in RBM47 gene expression data gathered by rt-qPCR from mice tissues
RBM47GE_mousetissues = pd.DataFrame({'GE_RBM47' : [0.367292, 0.158769, 0.156583,
                                                    0.17254, 0.264255,
                                                    0.090873, 0.088083, 0.029462,
                                                    0.092783, 0.092783, 0.096723,
                                                    0.033609, 0.005263, 0.00879,
                                                    0.002971, 0.001724, 0.002036
                                                    ],
                                           'tissue' : ['liver', 'liver', 'liver',
                                                'kidney', 'kidney',
                                                'spleen', 'spleen', 'spleen',
                                                'lung', 'lung', 'lung',
                                                'heart', 'heart','heart',
                                                'brain', 'brain', 'brain'
                                                    ]
                                         })
RBM47GE_mousetissues = RBM47GE_mousetissues.set_index('tissue')

#import human RBM47 gene expression level from RNA sequencing data
genecounts_all = pd.read_csv('Data/dataframeall.csv', sep=',', index_col=0)
genecounts_all_RBM47 = pd.DataFrame(genecounts_all.loc['RBM47'])
genecounts_all_RBM47 = genecounts_all_RBM47.reset_index()
genecounts_all_RBM47['index'] = genecounts_all_RBM47['index'].str.replace('1', '')
genecounts_all_RBM47['index'] = genecounts_all_RBM47['index'].str.replace('2', '')
genecounts_all_RBM47['index'] = genecounts_all_RBM47['index'].str.replace('3', '')
genecounts_all_RBM47['index'] = genecounts_all_RBM47['index'].str.replace('4', '')
genecounts_all_RBM47['index'] = genecounts_all_RBM47['index'].str.replace('5', '')
genecounts_all_RBM47['index'] = genecounts_all_RBM47['index'].str.replace('6', '')
genecounts_all_RBM47['index'] = genecounts_all_RBM47['index'].str.replace('7', '')
genecounts_all_RBM47['index'] = genecounts_all_RBM47['index'].str.replace('8', '')
genecounts_all_RBM47['index'] = genecounts_all_RBM47['index'].str.replace('9', '')
genecounts_all_RBM47['index'] = genecounts_all_RBM47['index'].str.replace('0', '')
genecounts_all_RBM47.columns = ['tissue', 'RBM47']
genecounts_all_RBM47 = genecounts_all_RBM47.set_index('tissue')


#E24csplicing_mousetissues = E24csplicing_mousetissues.set_index('tissue') 
size = 3.5
plt.figure(figsize=(size * 1.38, size))

tissue_list = ['liver', 'kidney', 'spleen', 'lung', 'heart', 'brain']
mouse_RBM47_mean = []
human_RBM47_mean = []
mouse_RBM47_std = []
human_RBM47_std = []
for i in tissue_list:
    mouse_RBM47_mean.append(np.mean(RBM47GE_mousetissues.loc[i]['GE_RBM47']))
    mouse_RBM47_std.append(np.std(RBM47GE_mousetissues.loc[i]['GE_RBM47']))
    human_RBM47_mean.append(np.mean(genecounts_all_RBM47.loc[i]['RBM47']))
    human_RBM47_std.append(np.std(genecounts_all_RBM47.loc[i]['RBM47']))
    
RBM47GE_validationplot_mouse = pd.DataFrame({'RBM47_GE_mouse_validation' : mouse_RBM47_mean,
                                                 'mouse_std' : mouse_RBM47_std,
                                                'tissue' : tissue_list,
                                                 #'organism' : 'mouse' 
                                                 })    
RBM47GE_validationplot_human = pd.DataFrame({'RBM47_GE_human_prediction' : human_RBM47_mean,
                                                 'human_std' : human_RBM47_std,
                                               # 'tissue' : tissue_list,
                                                 #'organism' : 'human' 
                                                 })    
RBM47GE_validationplot = pd.concat([RBM47GE_validationplot_human, RBM47GE_validationplot_mouse]
                                       ,axis = 1)

palette = [(0.86, 0.3712, 0.33999999999999997),
 (0.86, 0.5130181818181818, 0.33999999999999997),
 (0.33999999999999997, 0.4033454545454544, 0.86),
 (0.35607272727272715, 0.86, 0.33999999999999997),
 (0.560290909090909, 0.33999999999999997, 0.86),
 (0.7021090909090906, 0.33999999999999997, 0.86)]

for n in range(0,6):
    plt.errorbar(human_RBM47_mean[n], mouse_RBM47_mean[n], xerr = human_RBM47_std[n], yerr=mouse_RBM47_std[n]
             , fmt="none"
            ,color = palette[n]
             ,capsize = 3
             ,alpha= .5
            )



sns.regplot(data = RBM47GE_validationplot, x = 'RBM47_GE_human_prediction'
                , y = 'RBM47_GE_mouse_validation'
           , marker = ''
           , ci = 0
           , color = 'gray'
           ).set(xlabel = 'human RBM47 from RNA-Seq data'
                                         , ylabel = 'mouse RBM47 gene expression (rel. to Hprt)')

sns.scatterplot(data = RBM47GE_validationplot, x = 'RBM47_GE_human_prediction'
                , y = 'RBM47_GE_mouse_validation'
                , hue = 'tissue'
                , palette = palette).set(xlabel = 'human RBM47 from RNA-Seq data'
                                         , ylabel = 'mouse RBM47 gene expression (rel. to Hprt)')

plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

reg_stats = linregress(RBM47GE_validationplot['RBM47_GE_human_prediction']
                 , RBM47GE_validationplot['RBM47_GE_mouse_validation'])

plt.xlabel("human RBM47 from RNA-Seq data", fontsize = 15)
plt.ylabel("mouse RBM47 gene expression \n (rel. to Hprt)", fontsize = 15)

plt.yticks(fontsize=14)
plt.xticks(fontsize=14)

print(reg_stats)
plt.savefig("Figures/RBM47_GE_mouse_validation_reg.png", dpi = 600, bbox_inches = "tight")


"""  Mouse tisse RBM47 GE E24c splicing """

#import experimental data from radioactive splice PCR on mice tissues

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

RBM47GE_mousetissues = pd.DataFrame({'GE_RBM47' : [0.367292, 0.158769, 0.156583,
                                                    0.17254, 0.264255,
                                                    0.090873, 0.088083, 0.029462,
                                                    0.092783, 0.092783, 0.096723,
                                                    0.033609, 0.005263, 0.00879,
                                                    0.002971, 0.001724, 0.002036
                                                    ],
                                           'tissue' : ['liver', 'liver', 'liver',
                                                'kidney', 'kidney',
                                                'spleen', 'spleen', 'spleen',
                                                'lung', 'lung', 'lung',
                                                'heart', 'heart','heart',
                                                'brain', 'brain', 'brain'
                                                    ]
                                         })
RBM47GE_mousetissues = RBM47GE_mousetissues.set_index('tissue')



#plot SEC31A exon 24c and RBM47 gene expression from mice tissues
#E24csplicing_mousetissues = E24csplicing_mousetissues.set_index('tissue') 
tissue_list = ['liver', 'kidney', 'spleen', 'lung', 'heart', 'brain']
splicing_mean = []
GE_mean = []
splicing_std = []
GE_std = []
for i in tissue_list:
    splicing_mean.append(np.mean(E24csplicing_mousetissues.loc[i]['PSI 24c']))
    splicing_std.append(np.std(E24csplicing_mousetissues.loc[i]['PSI 24c']))
    GE_mean.append(np.mean(RBM47GE_mousetissues.loc[i]['GE_RBM47']))
    GE_std.append(np.std(RBM47GE_mousetissues.loc[i]['GE_RBM47']))
    
E24csplicing_validationplot_mouse = pd.DataFrame({'Mouse Exon 24c PSI' : splicing_mean,
                                                 'splicing_std' : splicing_std,
                                                'tissue' : tissue_list,
                                                 #'organism' : 'mouse' 
                                                 })    
E24csplicing_validationplot_human = pd.DataFrame({'Mouse RBM47 gene expression (rel. to Hprt)' : GE_mean,
                                                 'GE_std' : GE_std,
                                               # 'tissue' : tissue_list,
                                                 #'organism' : 'human' 
                                                 })    
E24csplicing_validationplot = pd.concat([E24csplicing_validationplot_human, E24csplicing_validationplot_mouse]
                                       ,axis = 1)
for n in range(0,6):
    plt.errorbar(GE_mean[n], splicing_mean[n], xerr = GE_std[n], yerr=splicing_std[n]
             , fmt="none"
            ,color = palette[n]
             ,capsize = 3
             ,alpha= .5
            )

sns.scatterplot(data = E24csplicing_validationplot, x = 'Mouse RBM47 gene expression (rel. to Hprt)'
                , y = 'Mouse Exon 24c PSI'
                , hue = 'tissue'
                , palette = palette) 
sns.regplot(data = E24csplicing_validationplot, x = 'Mouse RBM47 gene expression (rel. to Hprt)'
                , y = 'Mouse Exon 24c PSI'
           , marker = ''
           , ci = 0
           , color = 'gray'
           )

plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

reg_stats = linregress(E24csplicing_validationplot['Mouse RBM47 gene expression (rel. to Hprt)']
                 , E24csplicing_validationplot['Mouse Exon 24c PSI'])
print(reg_stats)
plt.xlabel("Mouse RBM47 gene expression (rel. to Hprt)", fontsize = 20)
plt.ylabel("Mouse Exon 24c PSI", fontsize = 20)
plt.savefig("Figures/E24c_splicing_RBM47_GE_mouse.png", dpi = 600, bbox_inches = "tight")


"""  RBM47 over expression  """

E24c_splicing_RBM47_OE = pd.DataFrame({'PSI_E24c' : [0.0141585, 0.012621038,
                                                    0.164037232, 0.086022517, 0.138930506
                                                    ],
                                           'treatment' : ['Ctrl', 'Ctrl',
                                                'RBM47 OE', 'RBM47 OE', 'RBM47 OE' 
                                                
                                                
                                                    ]
                                         })


palette_generation_3 = sns.color_palette("RdGy_r", 5)
DATA = E24c_splicing_RBM47_OE
y_limit = 0.18
x_data = 'treatment'
y_data = 'PSI_E24c'
group1 = E24c_splicing_RBM47_OE.query('treatment == "Ctrl"')[y_data]
group2 = E24c_splicing_RBM47_OE.query('treatment == "RBM47 OE"')[y_data]

size = 3.5
plt.figure(figsize=(size * 0.45, size))
plot = sns.swarmplot(data = DATA, x = x_data, y = y_data,  alpha=1, color='black')
plot = sns.barplot(data= DATA, x = x_data, y = y_data
            , palette = [palette_generation_3[0],
                        palette_generation_3[4],
                        ]
            ,capsize = 0.2,             
            saturation = 8,             
            errcolor = 'black', errwidth = 2,
            linewidth=2, edgecolor="0")

#statistics1
p_ttest_1 = statistic_f(group1, group2)
cohensd_1 = cohens_d(group1, group2)

#statistics depiction with asterisks
asterisk_text_1 = statistic_asterisk(p_ttest_1)
x1, x2 = [0, 1]   
y, h, col = y_limit-(y_limit*0.06), y_limit*0.015, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=2, c=col)
plt.text((x1+x2)*.5, y, asterisk_text_1, ha='center', va='bottom', color=col, size = 20)

#print(stat_text)
print('Pvalue1: ' + str(p_ttest_1))
print('Cohen´s d 1: ' + str(cohensd_1))

plot.set(xlabel=None)
plt.yticks(fontsize=14)
plt.xticks(fontsize=12)

plot.set_ylabel("PSI exon 24c", fontsize=15)
plt.savefig("Figures/E24c_splicing_RBM47_OE.png", dpi=600, bbox_inches='tight')


"""  Quantification RBM47 overexpression  """

GE_RBM47_OE_quantification = pd.DataFrame({'PSI_E24c' : [0.015356572,
                                                                    0.019572882,
                                                                    0.020545981,
                                                                    124.9320654,
                                                                    121.9376637,
                                                                    96.33579183
                                                    ],
                                           'treatment' : ['Hek Ctrl', 'Hek Ctrl', 'Hek Ctrl',
                                                'Hek RBM47 OE', 'Hek RBM47 OE', 'Hek RBM47 OE']
                                         })

palette_generation_3 = sns.color_palette("RdGy_r", 5)
DATA = GE_RBM47_OE_quantification
y_limit = 140
x_data = 'treatment'
y_data = 'PSI_E24c'
group1 = DATA.query('treatment == "Hek Ctrl"')[y_data]
group2 = DATA.query('treatment == "Hek RBM47 OE"')[y_data]
plt.figure(figsize=(3,6))
plot = sns.swarmplot(data = DATA, x = x_data, y = y_data,  alpha=1, color='black')
plot = sns.barplot(data= DATA, x = x_data, y = y_data
            , palette = [palette_generation_3[0],
                        palette_generation_3[4]]
            ,capsize = 0.2,             
            saturation = 8,             
            errcolor = 'black', errwidth = 2,
            linewidth=2, edgecolor="0")



#statistics1
p_ttest_1 = statistic_f(group1, group2)
cohensd_1 = cohens_d(group1, group2)


#statistics depiction with asterisks
asterisk_text_1 = statistic_asterisk(p_ttest_1)
x1, x2 = [0, 1]   
y, h, col = y_limit-(y_limit*0.06), y_limit*0.015, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=2, c=col)
plt.text((x1+x2)*.5, y, asterisk_text_1, ha='center', va='bottom', color=col, size = 20)

#print(stat_text)
print('Pvalue1: ' + str(p_ttest_1))
print('Cohen´s d 1: ' + str(cohensd_1))

plot.set_ylabel("RBM47 expression (rel. to GAPDH)", fontsize=14)
plt.savefig("Figures/RBM47_OE_quantification_reltoGAPDH.png", dpi=600, bbox_inches='tight')


"""  RBM47 OE Westernblott   """

Westernblott_RBM47_OE_quantification = pd.DataFrame({'PSI_E24c' : [
0.000567614,
0.000886921,
0.001879443,
0.567954112,
0.680190195,
1.510938541
                                                    ],
                                           'treatment' : ['Hek Ctrl', 'Hek Ctrl', 'Hek Ctrl',
                                                'Hek RBM47 OE', 'Hek RBM47 OE', 'Hek RBM47 OE']
                                         })

palette_generation_3 = sns.color_palette("RdGy_r", 8)
DATA = Westernblott_RBM47_OE_quantification
y_limit = 1.8
x_data = 'treatment'
y_data = 'PSI_E24c'
group1 = DATA.query('treatment == "Hek Ctrl"')[y_data]
group2 = DATA.query('treatment == "Hek RBM47 OE"')[y_data]

plt.figure(figsize=(3,6))
plot = sns.swarmplot(data = DATA, x = x_data, y = y_data,  alpha=1, color='black')
plot = sns.barplot(data= DATA, x = x_data, y = y_data
            , palette = [palette_generation_3[0],
                        palette_generation_3[4]]
            ,capsize = 0.2,             
            saturation = 8,             
            errcolor = 'black', errwidth = 2,
            linewidth=2, edgecolor="0")



#statistics1
p_ttest_1 = statistic_f(group1, group2)
cohensd_1 = cohens_d(group1, group2)


#statistics depiction with asterisks
asterisk_text_1 = statistic_asterisk(p_ttest_1)
x1, x2 = [0, 1]   
y, h, col = y_limit-(y_limit*0.06), y_limit*0.015, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=2, c=col)
plt.text((x1+x2)*.5, y, asterisk_text_1, ha='center', va='bottom', color=col, size = 20)


#print(stat_text)
print('Pvalue1: ' + str(p_ttest_1))
print('Cohen´s d 1: ' + str(cohensd_1))

plot.set_ylabel("RBM47-FLAG expression (rel. to GAPDH)", fontsize=14)
plt.savefig("Figures/RBM47_OE_quantification_reltoGAPDH_westernblott.png", dpi=600, bbox_inches='tight')


"""  RBM47 knockdown  """


E24c_splicing_RBM47_KO = pd.DataFrame({'PSI_E24c' : 
                                       [0.153331631, 0.161688104, 0.101034164, 
                                       0.065283243, 0.03235441, 0.044363869]
                                      ,
                                      'treatment' : 
                                       ['siCtrl', 'siCtrl','siCtrl',
                                        'siRBM47', 'siRBM47', 'siRBM47',]
                                      })

    
palette_generation_3 = sns.color_palette("RdGy_r", 5)
DATA = E24c_splicing_RBM47_KO
y_limit = 0.18
x_data = 'treatment'
y_data = 'PSI_E24c'
group1 = DATA.query('treatment == "siCtrl"')[y_data]
group2 = DATA.query('treatment == "siRBM47"')[y_data]



size = 3.5
plt.figure(figsize=(size * 0.45, size))
plot = sns.swarmplot(data = DATA, x = x_data, y = y_data,  alpha=1, color='black')
plot = sns.barplot(data= DATA, x = x_data, y = y_data
            , palette = [palette_generation_3[0],
                        palette_generation_3[4],
                        palette_generation_3[1],
                        palette_generation_3[3]]
            ,capsize = 0.2,             
            saturation = 8,             
            errcolor = 'black', errwidth = 2,
            linewidth=2, edgecolor="0")

#statistics1
p_ttest = statistic_f(group1, group2)
cohensd = cohens_d(group1, group2)

#statistics depiction with asterisks
asterisk_text_1 = statistic_asterisk(p_ttest)
x1, x2 = [0, 1]   
y, h, col = y_limit-(y_limit*0.06), y_limit*0.015, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=2, c=col)
plt.text((x1+x2)*.5, y, asterisk_text_1, ha='center', va='bottom', color=col, size = 20)


#print(stat_text)
print('Pvalue1: ' + str(p_ttest))
print('Cohen´s d 1: ' + str(cohensd))

plot.set(xlabel=None)
plt.yticks(fontsize=14)
plt.xticks(fontsize=13)

plot.set_ylabel("PSI exon 24c", fontsize=14)
plt.savefig("Figures/E24c_splicing_RBM47_KO.png", dpi=600, bbox_inches='tight')


RBM47_GE_RBM47_KO = pd.DataFrame({'RBM47_GE' : 
                                       [0.968312644, 0.878762273, 1.171652576,
                                        0.263074777, 0.166494317, 0.355653641]
                                      ,
                                      'treatment' : 
                                       ['siCtrl', 'siCtrl','siCtrl',
                                        'siRBM47', 'siRBM47', 'siRBM47',]
                                      })

    
palette_generation_3 = sns.color_palette("RdGy_r", 5)
DATA = RBM47_GE_RBM47_KO
y_limit = 1.3
x_data = 'treatment'
y_data = 'RBM47_GE'
group1 = DATA.query('treatment == "siCtrl"')[y_data]
group2 = DATA.query('treatment == "siRBM47"')[y_data]

size = 3.5
plt.figure(figsize=(size * 0.45, size))
plot = sns.swarmplot(data = DATA, x = x_data, y = y_data,  alpha=1, color='black')
plot = sns.barplot(data= DATA, x = x_data, y = y_data
            , palette = [palette_generation_3[1],
                        palette_generation_3[3],
                        palette_generation_3[1],
                        palette_generation_3[3]]
            ,capsize = 0.2,             
            saturation = 8,             
            errcolor = 'black', errwidth = 2,
            linewidth=2, edgecolor="0")

#statistics1
p_ttest = statistic_f(group1, group2)
cohensd = cohens_d(group1, group2)

#statistics depiction with asterisks
asterisk_text_1 = statistic_asterisk(p_ttest)
x1, x2 = [0, 1]   
y, h, col = y_limit-(y_limit*0.06), y_limit*0.015, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=2, c=col)
plt.text((x1+x2)*.5, y, asterisk_text_1, ha='center', va='bottom', color=col, size = 20)


#print(stat_text)
print('Pvalue1: ' + str(p_ttest))
print('Cohen´s d 1: ' + str(cohensd))

plot.set(xlabel=None)
plt.yticks(fontsize=14)
plt.xticks(fontsize=13)


plot.set_ylabel("RBM47 gene expresison (FC)", fontsize=14)
plt.savefig("Figures/RBM47_GE_RBM47_KO.png", dpi=600, bbox_inches='tight')


"""  Minigene  """

Minigene_hek = pd.read_excel('Data/SEC31AE24c_minigene_Hek_results.xlsx')

Minigene_hek['sample'] = Minigene_hek['FLAG'] + ' ' + Minigene_hek['Minigene']

palette_generation = sns.color_palette("RdGy_r", 9)
DATA = Minigene_hek.query('sample != "CTRL SD"').query('sample != "RBM47 SD"')
x_data = 'sample'
y_data = 'E24 inclusion'
group1 = DATA.query('sample == "CTRL S"')[y_data]
group2 = DATA.query('sample == "RBM47 S"')[y_data]
group5 = DATA.query('sample == "CTRL L"')[y_data]
group6 = DATA.query('sample == "RBM47 L"')[y_data]
group7 = DATA.query('sample == "CTRL ES"')[y_data]
group8 = DATA.query('sample == "RBM47 ES"')[y_data]

size = 3.5
plt.figure(figsize=(size * 1.5, size))
plot = sns.swarmplot(data = DATA, x = x_data, y = y_data,  alpha=1, color='black')
plot = sns.barplot(data= DATA, x = x_data, y = y_data
            , palette = [palette_generation[1],
                        palette_generation[7]]
            ,capsize = 0.2,             
            saturation = 8,             
            errcolor = 'black', errwidth = 2,
            linewidth=2, edgecolor="0")

#statistics1
p_ttest1 = statistic_f(group1, group2)
cohensd1 = cohens_d(group1, group2)

#statistics depiction with asterisks
y_limit = max(group2)+0.02
asterisk_text_1 = statistic_asterisk(p_ttest1)
x1, x2 = [4, 5]   
y, h, col = y_limit-(y_limit*0.1), y_limit*0.015, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=2, c=col)
plt.text((x1+x2)*.5, y, asterisk_text_1, ha='center', va='bottom', color=col, size = 20)

#print(stat_text)
print('Pvalue1: ' + str(p_ttest1))
print('Cohen´s d 1: ' + str(cohensd1))


p_ttest3 = statistic_f(group5, group6)
cohensd3 = cohens_d(group5, group6)

#statistics depiction with asterisks
y_limit = max(group6)+0.03
asterisk_text_2 = statistic_asterisk(p_ttest3)
x1, x2 = [6, 7]   
y, h, col = y_limit-(y_limit*0.1), y_limit*0.015, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=2, c=col)
plt.text((x1+x2)*.5, y-(y_limit*0.02), asterisk_text_2, ha='center', va='bottom', color=col, size = 20)

print('Pvalue3: ' + str(p_ttest3))
print('Cohen´s d 3: ' + str(cohensd3))

p_ttest6 = statistic_f(group2, group6)
cohensd6 = cohens_d(group2, group6)

#statistics depiction with asterisks
y_limit=0.23
asterisk_text_2 = statistic_asterisk(p_ttest6)
x1, x2 = [5, 7]   
y, h, col = y_limit-(y_limit*0.1), y_limit*0.015, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=2, c=col)
plt.text((x1+x2)*.5, y+.002, asterisk_text_2, ha='center', va='bottom', color=col, size = 20)

print('Pvalue6: ' + str(p_ttest6))
print('Cohen´s d 6: ' + str(cohensd6))

p_ttest7 = statistic_f(group7, group8)
cohensd7 = cohens_d(group7, group8)

#statistics depiction with asterisks
y_limit = max(group7)+0.01
asterisk_text_2 = statistic_asterisk(p_ttest6)
x1, x2 = [2, 3]   
y, h, col = y_limit-(y_limit*0.1), y_limit*0.015, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=2, c=col)
plt.text((x1+x2)*.5, y-(y_limit*0.02), asterisk_text_2, ha='center', va='bottom', color=col, size = 20)

print('Pvalue7: ' + str(p_ttest7))
print('Cohen´s d 7: ' + str(cohensd7))

plot.set(xlabel=None)
plt.yticks(fontsize=14)
plt.xticks(fontsize=0)

plot.set_ylabel("PSI exon 24c", fontsize=15)
plt.savefig("Figures/E24c_minigene_Hek.png", dpi=600, bbox_inches='tight')