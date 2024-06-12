# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 13:39:46 2023

@author: kayve
"""

import json
import os
import re
import seaborn as sns
import pandas as pd
#import numpy as np
import matplotlib.pyplot as plt
import tqdm
#import subprocess
#python setup.py install

def parse_json(json_filename, ko_list):
    key_w = ['KO', 'members']
    with open(json_filename) as file:
    
    #parse the json file and find the object that has the KO from the list of KOs 
        json_array =[]
        #ko2bin = dict()
        current_dd = ""
        for line in file: 
            current_dd += line
            if "}" in line:
                current_dd = current_dd.strip().strip(",").strip("[")
                dd = json.loads(current_dd)
                current_dd = ""
                if dd['KO'][3:] in ko_list:
                    json_array.append(dd)
                    
    #                 if dd['KO'] not in ko2members:
    #                     ko2members[dd['KO']] = set() #will give unique members from one annotation
    #                     #ko2bin[dd['KO']] = set()
    #                 ko2members[dd['KO']].update(dd['members'])
    #                 #ko2bin[dd['KO']].update(["_".join(d.split("_")[:-1]) for d in dd['members']])
                
    # with open(json_filename) as file:
    #     json_array = json.load(file)    
    
    #this will give a list of key values that are specified according to the key_list from the list of dictionaries
    #ref_list = [[(key, d[key]) for key in key_w if key in d] for d in json_array ]
    
    ref_list = [[d.get(key) for key in key_w] for d in json_array ]
    
    #ref_list = [[(key, d[key]) for key in key_w if key in d] for d in json_array ]
    
    #new & final update
    refined_df = pd.DataFrame.from_records({key_w[1]: [i for k, v in ref_list for i in v], key_w[0]: [k for k, v in ref_list for i in v]})
    
    #refined_df['members'] = refined_df.groupby("KO").sum().members  #will give the unique members
    
    
    #unique members grouped by KO without gene cluster name
    refined_df['members'] = refined_df.members.apply(lambda y: "_".join(y.split("_")[:-1])) 
    
    #remove KOs
    refined_df['KO'] = refined_df['KO'].str[3:]
    
    #group the KOs by members to get what KOs each member(bin) has
    refined_grouped_df = refined_df.groupby(['members'], as_index = False).agg({'KO': lambda x :  {xx for xx in x} })
    #refined_df = refined_df.merge(kofam, on = 'cluster_id', how = 'left')
    
    #refined_dict = refined_df.set_index('KO')['members'].to_dict() #apparently easier to navigate within dictionary
    #tt['bins'] = tt.members.apply(lambda y: set(["_".join(d.split("_")[:-1]) for d in y]))
    refined_df.to_csv('set_of_all_candidates.csv')
    refined_grouped_df.to_csv('set_of_all_candidates_grouped.csv')
    #mem = sum(refined_df.memory_usage(index=True, deep=True))
    #print(f'Memory used, hopefully this doesnt crash {mem:} B or {mem/1000000:.2f}' MB)
    
    return refined_df  #make them as list for easier access or make them dictionary


def table_parse(file, ko_complete_list):
     kofam_text = pd.read_csv(file, comment="#", header= None, delimiter=r"\s{2,}", engine ="python")
     kofam_text.columns = ['gene', 'KO', 'treshhold', 'eval', 'def']
    #grep "^\*" kofamscan.txt to get the ones that have high threshhold
     dd = {'members': [row.gene.replace('* ','') for i, row in kk.iterrows() if row.KO in ko_complete_list],'KO': [row.KO for i, row in kk.iterrows() if row.KO in ko_complete_list]}



#this function will take the kegg module definitions -- that has one or several steps of reactions and s
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
        sns.heatmap(df_from_json, vmin=0, vmax=1, cbar=False, cmap="Blues", ax = ax)
        ax.set_title(key_n)
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=4)
        file_name = f"{key_n}_kofam.png"
        fig = plt.savefig(file_name, dpi=400, bbox_inches='tight', pad_inches=0.5)
        return fig


'''
#for heatmap, pheatmap package in R 
def R_script_runner(modules_df, mod_name):
    

'''

def main(): #cwd for json file, text file with ko, parses them, clusters them and gives a heatmap 
    path_to_files = os.getcwd()
    #list_of_ko_files = []
    for file in os.listdir(path_to_files): #gets all the files in the current working directory
        #if file.endswith('.txt'):     #gets all the text files
            #list_of_ko_files.append(file)
        if file.endswith('.json'): #maybe use for loop   
            json_annot = file
    
    modules = {}    
    c_pathways = ['M00173','M00376','M00375','M00374','M00377','M00579', 'M00620']
    for f in os.listdir("modules"):
        with open(f"modules/{f}") as handle:
            modules[f] = [" ".join(l.split()[1:]) for l in handle.readlines() if l.startswith("DEFINITION")][0]
    
    '''
    result_ko = []
    for ko_file in list_of_ko_files: #produces the list of KOs for each file
        pathway, list_ko = get_KOs(ko_file)
        result_ko.append( (pathway, list_ko) )
    #The above will be instead used on all modules, look for the ones that are specific for carbon fixation, get the kegg & definition and parse through clustering_KO function

    defs_ko = {}
    for k in c_pathways:
        if k in modules.keys():
            for x, y in modules.items():
                defs_ko.update({x : y})
    '''            
    
    #list of all KOs for c metabolic pathways - as list of lists
    #separates by space, all the symbols and stuff - doesnt have empty string as list elemnet and such bc it was a serious problem
    ko_list = {i: [x for x in re.split("(K\d+|[ \ \ \,\+\-\(\)])", modules[i]) if x != "" and x != '+' and x != '(' and x != ')' and  x != '-' and x != ','] for i in c_pathways if i in modules.keys()}
   
    ko_complete_list = [s for v in ko_list.values() for s in v]
    
    #uses the KO list from modules and plots the heatmap for 
    
    json_comp_df = parse_json(json_annot, ko_complete_list)
    
    
    
    
    #BRAND NEW!!! FOR TEXT FILE to get the sample names based on the cluster ID
    with open('arctic_gene_atlas.json') as json_file:
         json_dict = []
         current_dd = ""
         for line in json_file:
                 current_dd += line
                 if "}" in line:
                     current_dd = current_dd.strip().strip(",").strip("[")
                     dd = json.loads(current_dd)
                     current_dd = ""
                     for i, row in final_df.iterrows():
                         if dd['cluster_id'][4:] in row.members:
                             json_dict.append(dd)
    
    
    
    
    
    
    
    
    #grouped data and modules for adding perctent completion to each table
    grouped_df = json_comp_df.groupby(['members'], as_index = False).agg({'KO': lambda x :  {xx for xx in x}} )
    ko_list_defs = {i: "".join(x for x in modules[i]) for i in c_pathways if i in modules.keys()}
    
    #look for the KOs in grouped and assign them to respective modules
    mod_complete = {k : { row.members : perct_comp_pathway(row.KO, v) for i, row in grouped_df.iterrows() } for k, v in tqdm(ko_list_defs.items()) if k in c_pathways}
    mod_complete_df = pd.DataFrame.from_dict(mod_complete)
    
    #add module name to the dataframe based on the KO terms - for grouping and plotting
    res = []
    for key, value_list in ko_list.items():
        key_matches = json_comp_df[json_comp_df['KO'].isin(value_list)]
        if not key_matches.empty:
            key_results = key_matches.assign(key=key)
            res.append(key_results)

    final_df = pd.concat(res)
    
    
    #for every unique key, plot heatmap and also display how many elements - this will have a function that executes R script in place of heatmap
    for key_value in final_df['key'].unique():
        subset_df = final_df[final_df['key'] == key_value]
        pivot_df = subset_df.pivot_table(index='KO', columns='members', values='KO', aggfunc="size", fill_value=0).T
        pivot_df = pd.merge(pivot_df, mod_complete_df[[key_value]], left_on='members', right_index=True)
        pivot_df.to_csv(f"{key_value}_table_members_with_perct.csv")
        heatmap_simple(pivot_df, key_value)
    
    
    '''
    #assign module names to the members based on the 
    grouped_df['Module']=0
    for i, row in grouped_df.iterrows():
        for k, v in ko_list.items():
            ll = [x for x in list(row['KO']) if x in v]
            for ss in ll:
                grouped_df.loc[i, 'Module'] = k
    grouped_df['Perct_Comp'] = 0
    for jj in grouped_df['Module'].unique():
        g_sb = grouped_df[grouped_df['Module'] == jj]
        for mod, defo in ko_list_defs.items():
            if mod in g_sb['Module'].unique():
                percts = perct_comp_pathway(g_sb['members'], defo)
                grouped_df.loc[grouped_df['Module'] == jj, 'Perct_Comp'] = percts
    '''
    #mod_merged = pd.merge(grouped_df, mod_complete_df, left_on='members', right_index=True)
    
    
    '''
    all_grouped_df = all_grouped_df.rename(columns={all_grouped_df.columns[0]: all_grouped_df.iloc[0,0]})
    all_grouped_df = all_grouped_df.drop(all_grouped_df.index[0])
    all_grouped_df.tail(-1)
    '''

    
    #to compute the percentage completion of each bin for each module and store it in a dataframe
    #table_pect_comp = pd.DataFrame.from_dict({k : {tt[1] : perct_comp_pathway(tt[2] , m) for tt in json_result_group.itertuples()} for k,m in modules.items()})

    #here i shall use the output in json_parsed_for_pathways, that has tuples that has two dataframes, use the second
    
    #next will be heatmap with all the info

if __name__ == '__main__':
    main()