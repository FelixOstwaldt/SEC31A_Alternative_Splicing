# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 15:45:39 2025

@author: felix
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import levene

def mean_data(df):
    items = list(set([c[:-1] for c in df.columns.to_list()]))

    
    numbers = list(range(0,len(items)))
    data = []
    mean = pd.DataFrame(data)
    for i in numbers:
        cols = [c for c in df if c.startswith(items[i])]

        mean[items[i]] = df.loc[:,cols].mean(1)
    return mean


def statistic_f(DATASET1, DATASET2):
    #testing normal distribution
    if len(DATASET1) > 5000 or len(DATASET2) > 5000:
        p_normaldist1 = stats.kstest(DATASET1, stats.norm.cdf)[1]
        p_normaldist2 = stats.kstest(DATASET2, stats.norm.cdf)[1]
    if len(DATASET1) < 3 or len(DATASET2) < 3:
        p_normaldist1 = 1
        p_normaldist2 = 1
    else:
        p_normaldist1 = stats.shapiro(DATASET1)[1]
        p_normaldist2 = stats.shapiro(DATASET2)[1]
    
    #testing variance
    stat, p_levene = levene(DATASET1, DATASET2)
    
    if p_normaldist1 < 0.05 or p_normaldist2 < 0.05:
        p_ttest = stats.mannwhitneyu(DATASET1, DATASET2)[1]  #Mann whitney U test
        print('The data is not normaly distributed, so the  Mann-Whitney U test was performed')
    else:
        if p_levene < 0.05:
            equal_True_False = False
            print('The data is normaly distributed but has no equal variance so Welch´s T-test was perfomred')
        else:
            equal_True_False = True
            print('The data is normaly distributed and has equal variance so Student´s T-test was perfomred')
        p_ttest = stats.ttest_ind(DATASET1, DATASET2 ,equal_var=equal_True_False)[1] # t-test
    return(p_ttest)


def cohens_d(DATASET1, DATASET2):
    #calculating Cohen´s d
    mean_difference = np.mean(DATASET1) - np.mean(DATASET2)
    pooled_standard_deviation = np.sqrt((np.var(DATASET1) + np.var(DATASET2)) / 2)
    cohen_d = mean_difference / pooled_standard_deviation
    return(cohen_d)

def statistic_asterisk(p_value):
    #statistics depiction with asterisks
    if p_value < 0.0001:
        p_value_plot = '****'
    elif p_value < 0.001:
        p_value_plot = '***'
    elif p_value < 0.01:
        p_value_plot = '**'
    elif p_value < 0.05:
        p_value_plot = '*'
    else:
        p_value_plot = 'ns'
    return(p_value_plot)
    


def read_splice_events(TYPE, tissue, path):
    if tissue == 'liver':
        splicing = pd.read_excel(f'{path}Livervsbrain/{TYPE}_Livervsbrain.xlsx'
                                ,usecols = [2,3,4,5,6,7,8,9,10,13,14,17,18,19])
        #combine coordinates in one "coord" column
        splicing['coord_chr'] = splicing.iloc[:, 1]
        splicing['coord_1'] = splicing.iloc[:, 3]
        splicing['coord_2'] = splicing.iloc[:, 4]
        splicing['coord_3'] = splicing.iloc[:, 5]
        splicing['coord_4'] = splicing.iloc[:, 6]
        splicing['coord_5'] = splicing.iloc[:, 7]
        splicing['coord_6'] = splicing.iloc[:, 8]
        
        coverage_list = []
        for i in range(len(splicing)):
            coverage_list.append(np.mean(list(map(int, splicing['IJC_SAMPLE_1'][i].split(",")))) + np.mean(list(map(int, splicing['SJC_SAMPLE_1'][i].split(",")))))
        splicing['mean_coverage'] = coverage_list
        splicing = splicing.query('mean_coverage > 10')
        
        #drop the unused columns
        splicing = splicing.drop(splicing.columns[[1,3,4,5,6,7,8,9,10,21]],axis = 1)
        #include a splice type column
        splicing['Type'] = TYPE 
        return splicing
    
    else:
        
        splicing = pd.read_excel(f'{path}Livervs{tissue}/{TYPE}_Livervs{tissue}.xlsx')
        #combine coordinates in one "coord" column
        splicing['coord_chr'] = splicing.iloc[:, 3]
        splicing['coord_1'] = splicing.iloc[:, 5]
        splicing['coord_2'] = splicing.iloc[:, 6]
        splicing['coord_3'] = splicing.iloc[:, 7]
        splicing['coord_4'] = splicing.iloc[:, 8]
        splicing['coord_5'] = splicing.iloc[:, 9]
        splicing['coord_6'] = splicing.iloc[:, 10]
        
        coverage_list = []
        for i in range(len(splicing)):
            coverage_list.append(np.mean(list(map(int, splicing['IJC_SAMPLE_1'][i].split(",")))) + np.mean(list(map(int, splicing['SJC_SAMPLE_1'][i].split(",")))))
        splicing['mean_coverage'] = coverage_list
        splicing = splicing.query('mean_coverage > 10')
        
        #drop the unused columns
        splicing = splicing.drop(splicing.columns[[0,1, 3, 5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]],axis = 1)
        splicing = splicing.drop('mean_coverage',axis = 1)
        #include a splice type column
        splicing['Type'] = TYPE 
        return splicing
