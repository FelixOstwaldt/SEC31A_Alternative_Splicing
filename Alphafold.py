# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 22:03:25 2025

@author: felix
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("whitegrid")
import glob


""" pLDDT   """

path = 'Data/Alphafold/multimer/active_fragment/SAR1A_SEC23_SEC31A_Iso1_activeFragment_multimer/'
file_path = glob.glob(path + '*ranking_confidence.npy') 
model_list = []
rankingconf_list = []
for file in file_path:
    file_name = file.replace(path[:-1] + '\\' , '')
    file_name = file_name.replace('_ranking_confidence.npy' , '')
    file_name = file_name.replace('_multimer_v3_pred' , '')
    Iso1_rankingconf = np.load(file)
    model_list.append(file_name)
    rankingconf_list.append(Iso1_rankingconf)
Iso1_ranks = pd.DataFrame({'Models' : model_list, 'ranking confidence' : rankingconf_list})
Iso1_ranks.sort_values(by = 'ranking confidence', ascending = False)




path = 'Data/Alphafold/multimer/active_fragment/SAR1A_SEC23_SEC31A_Iso8_activeFragment_multimer/'
file_path = glob.glob(path + '*ranking_confidence.npy') 
model_list = []
rankingconf_list = []
for file in file_path:
    file_name = file.replace(path[:-1] + '\\' , '')
    file_name = file_name.replace('_ranking_confidence.npy' , '')
    file_name = file_name.replace('_multimer_v3_pred' , '')
    Iso8_rankingconf = np.load(file)
    model_list.append(file_name)
    rankingconf_list.append(Iso8_rankingconf)
Iso8_ranks = pd.DataFrame({'Models' : model_list, 'ranking confidence' : rankingconf_list})
Iso8_ranks.sort_values(by = 'ranking confidence', ascending = False)




path = 'Data/Alphafold/multimer/active_fragment/SAR1A_SEC23_SEC31A_Iso1_activeFragment_multimer/'
file_path = glob.glob(path + '*plddt.npy') 
Iso1_plddt = {}
for file in file_path:
    file_name = file.replace(path[:-1] + '\\' , '')
    file_name = file_name.replace('_plddt.npy' , '')
    file_name = file_name.replace('result_model_' , '')
    file_name = file_name.replace('_multimer_v3_pred' , '')
    plddt_i = np.load(file)
    Iso1_plddt[file_name] = plddt_i
    
    
    
path = 'Data/Alphafold/multimer/active_fragment/SAR1A_SEC23_SEC31A_Iso8_activeFragment_multimer/'
file_path = glob.glob(path + '*plddt.npy') 
Iso8_plddt = {}
for file in file_path:
    file_name = file.replace(path[:-1] + '\\' , '')
    file_name = file_name.replace('_plddt.npy' , '')
    file_name = file_name.replace('result_model_' , '')
    file_name = file_name.replace('_multimer_v3_pred' , '')
    plddt_i = np.load(file)
    Iso8_plddt[file_name] = plddt_i
    
    

plt.figure(figsize=(12,5))
plt.plot(Iso1_plddt['5_0'], color = '#dd1c77', alpha = .6, label = '24c-')
plt.plot(Iso8_plddt['4_2'], color = '#2ca25f', alpha = .6, label = '24c+')
plt.legend(loc='lower right')


plt.plot([0, 0], [110, 0], linewidth=2, color = 'grey')
plt.plot([198, 198], [110, 0], linewidth=2, color = 'grey')
plt.plot([963, 963], [110, 0], linewidth=2, color = 'grey')
#plt.plot([1050, 1050], [110, 0], linewidth=2, color = 'grey')

plt.text(60, 102, 'SAR1A', fontsize=15)
plt.text(500, 102, 'SEC23A', fontsize=15)
plt.text(970, 105, 'SEC31A', fontsize=10)
plt.text(970, 101, 'active fragment', fontsize=10)

plt.xlabel("Residue number", fontsize = 20)
plt.ylabel("pLDDT", fontsize = 20)

plt.savefig("Figures/Alphafold_SAR1_SEC23_SEC31A_activeFragment_plddt_Iso1_8.png", dpi = 600, bbox_inches = "tight")


"""  PAE  """

path = 'Data/Alphafold/multimer/active_fragment/SAR1A_SEC23_SEC31A_Iso1_activeFragment_multimer/'
file_path = glob.glob(path + '*pae.npy') 
Iso1_pae = {}
for file in file_path:
    file_name = file.replace(path[:-1] + '\\' , '')
    file_name = file_name.replace('_pae.npy' , '')
    file_name = file_name.replace('result_model_' , '')
    file_name = file_name.replace('_multimer_v3_pred' , '')
    pae_i = np.load(file)
    Iso1_pae[file_name] = pd.DataFrame(pae_i)
    

plt.figure(figsize=(10,7))
plot = sns.heatmap(Iso1_pae['5_0'], cmap = 'PiYG_r', xticklabels=100, yticklabels=100, square = True
              )
#plt.xticks(ticks=np.arange(0,len(Iso8_pae['5_2']),100))
#plt.yticks(ticks=np.arange(0,len(Iso8_pae['5_2']),100))
#yticks = np.linspace(10,100,10)
#ylabels = np.linspace(100,1000,10)

plt.xlabel("Scored residue", fontsize = 20)
plt.ylabel("Aligned residue", fontsize = 20)
plt.savefig("Figures/Alphafold_SAR1_SEC23_SEC31A_activeFragment_pae_Iso1.png", dpi = 600, bbox_inches = "tight")

plt.show()
    


path = 'Data/Alphafold/multimer/active_fragment/SAR1A_SEC23_SEC31A_Iso8_activeFragment_multimer/'
file_path = glob.glob(path + '*pae.npy') 
Iso8_pae = {}
for file in file_path:
    file_name = file.replace(path[:-1] + '\\' , '')
    file_name = file_name.replace('_pae.npy' , '')
    file_name = file_name.replace('result_model_' , '')
    file_name = file_name.replace('_multimer_v3_pred' , '')
    pae_i = np.load(file)
    Iso8_pae[file_name] = pd.DataFrame(pae_i)
    
    
    
plt.figure(figsize=(10,7))
plot = sns.heatmap(Iso8_pae['5_2'], cmap = 'PiYG_r', xticklabels=100, yticklabels=100, square = True
              )
#plt.xticks(ticks=np.arange(0,len(Iso8_pae['5_2']),100))
#plt.yticks(ticks=np.arange(0,len(Iso8_pae['5_2']),100))
#yticks = np.linspace(10,100,10)
#ylabels = np.linspace(100,1000,10)

plt.xlabel("Scored residue", fontsize = 20)
plt.ylabel("Aligned residue", fontsize = 20)
plt.savefig("Figures/Alphafold_SAR1_SEC23_SEC31A_activeFragment_pae_Iso8.png", dpi = 600, bbox_inches = "tight")

plt.show()


PAE_Iso1_plot_df = pd.DataFrame({'pae' : Iso1_pae['5_1'][963:][800]})#.reset_index()
#PAE_Iso1_plot_df = PAE_Iso1_plot_df.drop('index', axis = 1)
PAE_Iso1_plot_df['index'] = range(976, 1037)
PAE_Iso1_plot_df = PAE_Iso1_plot_df.set_index('index')



size = 3
plt.figure(figsize=(size * 2.4, size))


plt.plot(PAE_Iso1_plot_df, color = '#c51b7d', alpha = .8, label = '24c-')
plt.plot(Iso8_pae['5_2'][963:975][800], color = '#4d9221', alpha = .8, label = '24c+')
plt.plot(Iso8_pae['5_2'][974:989][800], color = sns.color_palette("cividis", 10)[9], alpha = 1, label = 'Exon 24c')
plt.plot(Iso8_pae['5_2'][988:][800], color = '#4d9221', alpha = .8)
plt.legend(loc='lower right', fontsize="12")


plt.xlabel("Residue number", fontsize = 15)
plt.ylabel("Predicted Aligned Error", fontsize = 15)

plot.set(xlabel=None)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)

plt.savefig("Figures/Alphafold_pae_Iso1_8.png", dpi = 600, bbox_inches = "tight")


"""   Polypeptide alignment  """


peptide_align = pd.DataFrame({'H. sapiens' : [10,10,10,10,10,10,10,8,8,8,4,10,10,8,10,6,4,10,8,10,10,10,10,10,10,10,10,10,10,10,10,10,10,7,6,10,10,10,10]
                             ,'P. abelii' : [10,10,10,10,10,10,10,8,8,8,4,10,10,8,10,6,4,10,8,10,10,10,10,10,10,10,10,10,10,10,10,10,10,7,6,10,10,10,10]
                             ,'V. vulpes' : [10,10,10,10,10,10,10,8,8,8,4,10,10,8,10,6,4,10,8,10,10,10,10,10,10,10,10,10,10,10,10,10,10,7,6,10,10,10,10]
                             ,'M. musculus' : [10,10,10,10,10,10,10,8,8,8,4,10,10,8,10,6,4,10,8,10,10,10,10,10,10,10,10,10,10,10,10,10,10,7,6,10,10,10,10]
                             ,'R. rattus' : [10,10,10,10,10,10,10,8,8,8,4,10,10,8,10,6,4,10,8,10,10,10,10,10,10,10,10,10,10,10,10,10,10,7,6,10,10,10,10]
                             ,'T. s. elegans' : [10,10,10,10,10,10,10,8,8,8,4,10,10,8,10,6,4,10,8,10,10,10,10,10,10,10,10,10,10,10,10,10,10,7,6,10,10,10,10]
                             ,'G. gallus' : [10,10,10,10,10,10,10,8,8,8,4,10,10,8,10,6,4,10,8,10,10,10,10,10,10,10,10,10,10,10,10,10,10,7,6,10,10,10,10]
                             ,'Conservation' : [10,10,10,10,10,10,10,8,8,8,4,10,10,8,10,6,4,10,8,10,10,10,10,10,10,10,10,10,10,10,10,10,10,7,6,10,10,10,10]})

plt.figure(figsize=(10,3))
plot = sns.heatmap(peptide_align.transpose(), cmap = 'PiYG', square = False,  vmin=0, vmax=10,linewidth=.0
                   , cbar_kws={'label': 'Conservation', 'location' : "right"}
                   )
for label in plot.get_yticklabels():
    if label.get_text() != 'Conservation':
        label.set_fontstyle("italic")
        label.set_fontsize("18")
    
sequences = [['P', 'A', 'S', 'Q', 'R', 'T', 'E', 'N', 'Q', 'S', 'I', 'Q', 'D', 'Q', 'A', 'P', 'M', 'L', 'E', 'G', 'P', 'Q', 'N', 'G', 'W', 'N', 'D', 'P', 'P', 'A', 'L', 'N', 'R', 'V', 'P', 'K', 'K', 'K', 'K'],
['P', 'A', 'S', 'Q', 'R', 'T', 'E', 'N', 'Q', 'S', 'M', 'Q', 'D', 'Q', 'A', 'P', 'M', 'L', 'E', 'G', 'P', 'Q', 'N', 'G', 'W', 'N', 'D', 'P', 'P', 'A', 'L', 'N', 'R', 'V', 'P', 'K', 'K', 'K', 'K'],
['P', 'A', 'S', 'Q', 'R', 'T', 'E', 'N', 'Q', 'S', 'I', 'Q', 'D', 'Q', 'A', 'P', 'V', 'L', 'E', 'G', 'P', 'Q', 'N', 'G', 'W', 'N', 'D', 'P', 'P', 'A', 'L', 'N', 'R', 'V', 'P', 'K', 'K', 'K', 'K'],
['P', 'A', 'S', 'Q', 'R', 'T', 'E', 'N', 'Q', 'S', 'F', 'Q', 'D', 'Q', 'A', 'S', 'I', 'L', 'E', 'G', 'P', 'Q', 'N', 'G', 'W', 'N', 'D', 'P', 'P', 'A', 'L', 'N', 'R', 'V', 'P', 'K', 'K', 'K', 'K'],
['P', 'A', 'S', 'Q', 'R', 'T', 'E', 'N', 'Q', 'S', 'F', 'Q', 'D', 'Q', 'A', 'S', 'V', 'L', 'E', 'G', 'P', 'Q', 'N', 'G', 'W', 'N', 'D', 'P', 'P', 'A', 'L', 'N', 'R', 'V', 'P', 'K', 'K', 'K', 'K'],
['P', 'A', 'S', 'Q', 'R', 'T', 'E', 'N', 'Q', 'S', 'S', 'Q', 'D', 'K', 'A', 'S', 'T', 'L', 'E', 'G', 'P', 'Q', 'N', 'G', 'W', 'N', 'D', 'P', 'P', 'A', 'L', 'N', 'R', 'M', 'T', 'K', 'K', 'K', 'K'],
['P', 'A', 'S', 'Q', 'R', 'T', 'E', 'K', 'P', 'P', 'A', 'Q', 'D', 'Q', 'A', 'S', 'P', 'L', 'A', 'G', 'P', 'Q', 'N', 'G', 'W', 'N', 'D', 'P', 'P', 'A', 'L', 'N', 'R', 'A', 'A', 'K', 'K', 'K', 'K']
            , ['*','*','*','*','*','*','*','8','8','8','4','*','*','8','*','6','4','*','8','*','*','*','*','*','*','*','*','*','*','*','*','*','*','7','6','*','*','*','*']]
y = -1
for species in sequences:
    x = 0
    y = y + 1
    for aminoacid in species:
        
        plot.text(x + .05, y + .8, aminoacid , size = 14.5)
        x = x + 1

plt.savefig("Figures/Exon24c_conservation_long.png", dpi = 600, bbox_inches = "tight")

plt.show()