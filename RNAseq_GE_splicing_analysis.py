# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 15:50:20 2025

@author: felix
"""

import pandas as pd
from functions import read_splice_events


#import helper list GO term list regarding secretory genes
GO_secretion = pd.read_excel('Helper_lists/GOsecretion_human2.xlsx')
GO_secretion = GO_secretion.set_index('geneSymbol')


#read RMATS files and create dictionary with PSI of all splicing events of all conditions(tissues)
splicing_dict = {} #create empty dictionary
splicing_secretion_dict = {} #create specialized dictionary (GO term secretion) 
#list of conditions(tissues)
list_condition = ['liver', 'bonemarrow', 'brain',  'colon', 'esophagus', 'fat', 'heart', 'kidney'
            ,'lung', 'lymphnode', 'placenta', 'skin'
            , 'smallintestine', 'spleen', 'stomach', 'testis', 'thyroid'
            ,'urinarybladder']

source_path = ('Data/Splicing/')
#loop over condition in 
for x in list_condition:
    #read files
    SE = read_splice_events('SE', x, source_path)
    RI = read_splice_events('RI', x, source_path)
    A3SS = read_splice_events('A3SS', x, source_path)
    A5SS = read_splice_events('A5SS', x, source_path)
    
    #add to the dictionary a DataFrame with all splice events in the key of the condition 
    splicing_dict[f'{x}'] = pd.DataFrame(pd.concat([SE, RI,  A5SS, A3SS], axis = 0))
    #add to the specialized dicitonary with the merge with the secretion GO term
    splicing_secretion_dict[f'{x}'] = splicing_dict[f'{x}'].reset_index().merge(GO_secretion
                                        , on = 'geneSymbol')
    
    
#order the data and transform it into a ordered dataframe
for i in splicing_dict:
    splicing_dict[i] = splicing_dict[i].set_index(['geneSymbol', 'strand', 'coord_chr', 'coord_1', 'coord_2', 'coord_3'
                                                            , 'coord_4', 'coord_5', 'coord_6', 'Type'])
    
for i in splicing_secretion_dict:
    splicing_secretion_dict[i] = splicing_secretion_dict[i].drop(['index', 'GO_term'], axis = 1)
    splicing_secretion_dict[i] = splicing_secretion_dict[i].set_index(['geneSymbol', 'strand', 'coord_chr', 'coord_1' 
                                                            ,'coord_2', 'coord_3'
                                                            , 'coord_4', 'coord_5', 'coord_6', 'Type'])
    
    
splicing = splicing_dict['liver']
for i in splicing_dict:
    if i == 'liver':
        pass
    else: 
        splicing = splicing.merge(splicing_dict[i], left_index=True, right_index=True, how = 'inner')
        
    
splicing_secretion = splicing_secretion_dict['liver']
for i in splicing_secretion_dict:
    if i == 'liver':
        pass
    else: 
        splicing_secretion = splicing_secretion.merge(splicing_secretion_dict[i], 
                                                      left_index=True, right_index=True, how = 'inner')
    
#save splicing dataframes
splicing.reset_index().sort_values(by=['geneSymbol']).reset_index().drop('index' 
                                , axis = 1).to_excel('Results/splicing.xlsx')
splicing_secretion.reset_index().sort_values(by=['geneSymbol']).reset_index().drop('index'
                                , axis = 1).to_excel('Results/splicing_secretion.xlsx')

