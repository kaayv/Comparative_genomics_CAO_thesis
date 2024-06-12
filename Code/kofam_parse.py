# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 12:47:16 2023

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


def parse_kofam(txt_file, ko_list_complete):
    kofam_text = pd.read_csv(txt_file, comment="#", header= None, delimiter=r"\s{2,}", engine ="c", usecols = [0, 1, 2, 3, 4]) #just take a few columns and engine is the way the csv reader reads file, c is robust and faster and python is more comprehensive and can work with more parameters
    kofam_text.columns = ['cluster_id', 'KO', 'treshhold', 'eval', 'def']
    #grep "^\*" kofamscan.txt to get the ones that have high threshhold
    kf_dict = {'cluster_id': [row.cluster_id.replace('* ','') for i, row in kofam_text.iterrows() if row.KO in ko_list_complete],'KO': [row.KO for i, row in kofam_text.iterrows() if row.KO in ko_list_complete]}
    kf_df = pd.DataFrame(kf_dict)
    return kf_df

def perct_comp_pathway(keggs : list, def_line : str, andline = True): #keggs found in the sample, module names
    def_line = def_line.replace("+", " ")
    def_line = def_line.split() if andline else def_line.split(",")
    i = 0
    blocks = []
    fwd = None
    rev = None
    for v in def_line:
        blocks += [i]
        if v.startswith("("):
            fwd = v.count("(")
            rev = v.count(")")
        elif fwd:
            fwd += v.count("(")
            rev += v.count(")")
        if fwd == rev:
            i += 1
            fwd = None
            rev = None
    operation = " " if andline else ","
    def big_operation(x, b) :
        value = sum(x)/len(x) if andline else max(x)
        #print(f" {b} : {value}")
        return value
    
    and_block = [operation.join([b for i,b in enumerate(def_line) if blocks[i] == block]) for block in set(blocks)]
    
    
    and_block = [b[1:-1] if b.startswith("(") else b for b in and_block]    
    return big_operation([perct_comp_pathway(keggs, b, not andline) if ("," in b or " " in b) else int((b in keggs) if not b.startswith("-") else True) for b in and_block], and_block)


def heatmap_simple(df_from_json, key_n): #will this be computationally efficient or dataframe from previous step
        plt.figure(figsize=(40, 25))
        plt.clf()
        ax = plt.axes()
        sns.clustermap(df_from_json, vmin=0, vmax=1, cbar=False, cmap="Blues", ax = ax)
        ax.set_title(key_n)
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=4)
        file_name = f"{key_n}_kofam.png"
        fig = plt.savefig(file_name, dpi=400, bbox_inches='tight', pad_inches=0.5)
        return fig


def main(): #cwd for json file, text file with ko, parses them, clusters them and gives a heatmap 
    path_to_files = os.getcwd()

    for file in os.listdir(path_to_files): #gets all the files in the current working directory
        if file.endswith('.json'): #maybe use for loop   
            json_annot = file
    '''
    modules = {}    
    c_pathways = ['M00173','M00376','M00375','M00374','M00377','M00579', 'M00620']
    for f in os.listdir("modules"):
        with open(f"modules/{f}") as handle:
            modules[f] = [" ".join(l.split()[1:]) for l in handle.readlines() if l.startswith("DEFINITION") and l in c_pathways(logic)][0]
   ''' 
   
    c_pathways = ['M00173','M00376','M00375','M00374','M00377','M00579', 'M00620']
    ko_list_defs = {'M00173': '(K00169+K00170+K00171+K00172,K03737) ((K01007,K01006) K01595,K01959+K01960,K01958) K00024 (K01676,K01679,K01677+K01678) (K00239+K00240-K00241-K00242,K00244+K00245-K00246-K00247,K18556+K18557+K18558+K18559+K18560) (K01902+K01903) (K00174+K00175-K00177-K00176) K00031 (K01681,K01682) (K15230+K15231,K15232+K15233 K15234)', 'M00376': '(K02160+K01961+K01962+K01963) K14468 K14469 K15052 K05606 (K01847,K01848+K01849) (K14471+K14472) (K00239+K00240+K00241) K01679', 'M00375': '(K01964+K15037+K15036) K15017 K15039 K15018 K15019 K15020 K05606 (K01848+K01849) (K15038,K15017) K14465 (K14466,K18861) K14534 K15016 K00626', 'M00374': '(K00169+K00170+K00171+K00172) K01007 K01595 K00024 (K01677+K01678) (K00239+K00240-K00241-K18860) (K01902+K01903) (K15038,K15017) K14465 (K14467,K18861) K14534 K15016 K00626', 'M00377': 'K00198 K05299-K15022 K01938 K01491 K00297 K15023 K14138+K00197+K00194', 'M00579': '(K00625,K13788,K15024) K00925', 'M00620': '(K00169+K00170+K00171+K00172) (K01959+K01960) K00024 (K01677+K01678) (K18209+K18210) (K01902+K01903) (K00174+K00175+K00176+K00177)'}    

           
    #list of all KOs for c metabolic pathways - as list of lists
    #separates by space, all the symbols and stuff - doesnt have empty string as list elemnet and such bc it was a serious problem
    ko_list = {i: [x for x in re.split("(K\d+|[ \ \ \,\+\-\(\)])", ko_list_defs[i]) if x != "" and x != '+' and x != '(' and x != ')' and  x != '-' and x != ',' and x != ' ' and x != '\n'] for i in ko_list_defs}
    ko_complete_list = [s for v in ko_list.values() for s in v]
    
    #parse texte file to get dataframe with clusterID and KO
    kofam_parsed = parse_kofam('sub_fk.txt', ko_complete_list)
     
    #add module name to the dataframe based on the KO terms - for grouping and plotting
    res = []
    for key, value_list in ko_list.items():
        key_matches = kofam_parsed[kofam_parsed['KO'].isin(value_list)]
        if not key_matches.empty:
            key_results = key_matches.assign(key=key)
            res.append(key_results)

    kf_modules = pd.concat(res)
   
    
   #BRAND NEW!!! FOR TEXT FILE to get the sample names based on the cluster ID
    with open(json_annot) as json_file:
         json_dict = []
         current_dd = ""
         for line in json_file:
                 current_dd += line
                 if "}" in line:
                     current_dd = current_dd.strip().strip(",").strip("[")
                     dd = json.loads(current_dd)
                     current_dd = ""
                     for i, row in kofam_parsed.iterrows():
                         if dd['cluster_id'][4:] in row.cluster_id:
                             json_dict.append(dd)
     
    key_w = ['cluster_id', 'members']
    
    ref_list = [[d.get(key) for key in key_w] for d in json_dict ]
    
    refined_df = pd.DataFrame.from_records({key_w[1]: [i for k, v in ref_list for i in v], key_w[0]: [k for k, v in ref_list for i in v]})
    
    refined_df['members'] = refined_df.members.apply(lambda y: "_".join(y.split("_")[:-1]))
    
    #here i make two dataframes
    refined_df = refined_df.merge(kofam_parsed, on = 'cluster_id', how = 'inner')
    refr_df = refined_df.merge(kf_modules, on = 'cluster_id', how = 'inner')
    
    refined_df = tqdm(refined_df.explode('members'))
    refr_df = tqdm(refr_df.explode('members'))
    
    refined_df.to_csv('all_selected_members.csv')
    refr_df.to_csv('with_modules_selected_members.csv')
    
    #grouped data and modules for adding perctent completion to each table
    grouped_df = refined_df.groupby(['members'], as_index = False).agg({'KO': lambda x :  {xx for xx in x}} )
    grouped_df.to_csv('Grouped.csv')
    
    #all definitions of modules
    #ko_list_defs = {i: "".join(x for x in modules[i]) for i in c_pathways if i in modules.keys()}
    
    #look for the KOs in grouped and assign them to respective modules
    mod_complete = {k : { row.members : perct_comp_pathway(row.KO, v) for i, row in grouped_df.iterrows() } for k, v in tqdm(ko_list_defs.items()) if k in c_pathways}
    mod_complete_df = pd.DataFrame.from_dict(mod_complete)
    
       
    #for every unique key, plot heatmap and also display how many elements - this will have a function that executes R script in place of heatmap
    for key_value in refined_df['key'].unique():
        subset_df = refined_df[refined_df['key'] == key_value]
        pivot_df = subset_df.pivot_table(index='KO', columns='members', values='KO', aggfunc="size", fill_value=0).T
        pivot_df = pd.merge(pivot_df, mod_complete_df[[key_value]], left_on='members', right_index=True)
        pivot_df.to_csv(f"{key_value}_table_members_with_perct.csv")
        heatmap_simple(pivot_df, key_value)
    
    
if __name__ == '__main__':
    main()
    