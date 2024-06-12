# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 15:26:09 2023

@author: kayve
"""
import json
import os
import re
import seaborn as sns
import pandas as pd
#import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

#will filter the bins which has 60% or more completion and get the mOTUs and their taxa
def filter_good_bins(file_50, complete_binset):
    #complete_binset variable -- in varable space
    
    #motustats import
    
    #file_50 is the df with module completeness and the count of each KO in the bins
    #file_50 = file_50.reset_index()
    file_50 = file_50.rename(columns={'members': 'bin_names'})
    file_name = file_50.columns[-1]
    
    #rounding off the decimals 
    file_50.iloc[:, -1] = file_50.iloc[:, -1].round(decimals=2)
    
    
#    file_50 = file_50[file_50[file_name] >= 0.60]
    '''
    This part does the following
    1. obtains a set of bins(non redundant) that has the completion percentage of more than 60%
    2. gets the set of mOTUs for the respective bins from complete_binset from the mOTU column
    3. removes the ones without mOTUs
    4. for the select mOTUs in mOTUs_of_interest, get the bin names from complete_binset
    5. check for unbinned bins and remove them 
    
    '''
    
    hit_bins = set(file_50.loc[file_50[file_name] >= 0.60, "bin_names"])
    mOTUs_of_interest = set(complete_binset.mOTU[complete_binset.bin_names.apply( lambda x : x in hit_bins)])
    mOTUs_of_interest.remove("NoMOTU")
    bins_of_interest = set(complete_binset.bin_names[complete_binset.mOTU.apply( lambda x : x in mOTUs_of_interest)])
    bins_of_interest = {b for b in bins_of_interest if "unbinned" not in b}
    
    '''
    Steps followed
    1. Make the bin name column as index.
    2. remove the bin_name column from the df(its already the index)
    3. mark everything greater than 0 as 1 
    4. keep the rpws where the bin names in file_50 intersect with the bins of interest that was consolidated above
    '''
    
    file_50.index = file_50.bin_names
    file_50 = file_50[file_50.columns[1:-1]]
    file_50 = file_50.applymap(lambda x : int(x > 0))
    file_50 = file_50.loc[bins_of_interest.intersection(file_50.index)]

    #hget other information from complete_binset
    file_50 = file_50.merge(complete_binset[['bin_names', 'mOTU', 'gtdbtk_classif_r207', 'percent_completion']], on = 'bin_names', how = 'left')
    
    mOTU_bins = complete_binset['mOTU'].value_counts()
    #print(f"The count for every module that present in the dataset for the module {file_name} is in the file")
    #file_50_select = file_50[file_50.iloc[:, -4].mean() > 1.05]
    
    
    print(f"the file {file_name} is saved for heatmap")
    file_50.to_csv(f"{file_name}_contents_fin.csv")
    print(f"{file_50}")
    
    file_50 = file_50.drop(columns=['bin_names', 'gtdbtk_classif_r207', 'percent_completion'])
    
    
    '''
    now groupby the column for each mOTU
    take the frequency of occurence of the KOs for the mOTUs based on number of bins they occur in the table
    that is, take the KO for each mOTU, divide them respectively by the number of bins and transpose
    
    add the number of bins from the mOTU_bins df for indices of each freq_50
    mean bin completion and highest bin completion are recorded
    
    maybe also add the lowest level of completeness for the bins that at least has 60 % of the pathways ccomplete(aim of the filtering)
    '''
    file_50 = file_50.groupby(['mOTU'], as_index=True).sum()
    print(f"the number of KOs is {len(file_50.columns)}")
    freq_50 = (file_50.transpose()/mOTU_bins[file_50.index]).transpose()
    freq_50['nb_bins'] = mOTU_bins[freq_50.index]
    freq_50['mean_bin_completion'] = complete_binset.groupby('mOTU').percent_completion.mean()[freq_50.index]
    freq_50['highest_bin_completion'] = complete_binset.groupby('mOTU').percent_completion.max()[freq_50.index]
    #maybe add a column for gtdbtk classification to put it on the plot
    print("the file with mOTU completeness is stored for the frequency of occurence of KOs in the bins of a mOTU")
    freq_50.to_csv(f"{file_name}_motu_complete_perct.csv")
    
    return

#to make heatmaps with percent completion for each mOTU - this is for selecting taxa
def replace_with_column_name(row):
    return [col if val >= 1 else val for col, val in zip(df.columns, row)]
m173_c = pd.read_csv('M00375_pathway_presence.csv', index_col='mOTU')
m173_c = m173_c.applymap(lambda x : int(x > 0))
m173_c_c = m173_c.apply(replace_with_column_name, axis=1)
df_with_KO = pd.Series.to_frame(m173_c_c, name='KO').applymap(lambda lst: list(filter(lambda x: x != 0, lst)) if isinstance(lst, list) else lst)
mod_complete_g = {k : { i : perct_comp_pathway(row.KO, v) for i, row in df_with_KO.iterrows() } for k, v in tqdm(ko_list_defs.items()) if k in c_pathways}
mod_complete_dfk = pd.DataFrame.from_dict(mod_complete_g)
mod_complete_dfk.to_csv('M00375_module_completeness_select_motus.csv')
m173_c = m173_c.join(mod_complete_dfk['M00375'])
m173_c.to_csv('M00375_for_heatmap_select.csv')


global motustats
motustats = pd.read_csv('motustats.csv')
global complete_binset
complete_binset = pd.read_csv('COMPLETE-SET-AND-BULK_basics.csv')    