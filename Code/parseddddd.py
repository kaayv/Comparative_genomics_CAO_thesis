# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 10:48:36 2023

@author: kayve
"""

import pandas as pd 
import json
import re
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

def json_parse(json_f, kofam_df):
        with open(json_f) as json_file:
            json_dict = []
            current_dd = ""
            for line in json_file:
                current_dd += line
                if "}" in line:
                    current_dd = current_dd.strip().strip(",").strip("[")
                    dd = json.loads(current_dd)
                    current_dd = ""
                    for i, row in kofam_df.iterrows():
                        if dd['cluster_id'][4:] in row.cluster_id:
                            json_dict.append(dd)
        
        #get only clusterID and members
        key_w = ['cluster_id', 'members']
        
        #in a list
        ref_list = [[d.get(key) for key in key_w] for d in json_dict ]
        
        #Make it a df
        member_df = pd.DataFrame.from_records({key_w[1]: [i for k, v in ref_list for i in v], key_w[0]: [k for k, v in ref_list for i in v]})
        
        #remove the gene cluster part in members
        member_df['members'] = member_df.members.apply(lambda y: "_".join(y.split("_")[:-1]))
        
        #add this to the kofam dataframe for all the clusters(only on clusters)
        member_df = member_df.merge(kofam_df, on = 'cluster_id', how = 'inner')
        
        #drop clusterID column ( as we only want members)
        member_wo_cluster = member_df.drop(columns = ['cluster_id'])
        
        member_df.to_csv('all_members_clusters.csv')

        #explode members for all
        #member_df = tqdm(members_df.explode('members')) 
        #member_wo_cluster = tqdm(member_wo_cluster.explode('members'))

        return member_df, member_wo_cluster

        
                            
                            
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


def main():
    c_pathways = ['M00173','M00376','M00375','M00374','M00377','M00579', 'M00620']
    ko_list_defs = {'M00173': '(K00169+K00170+K00171+K00172,K03737) ((K01007,K01006) K01595,K01959+K01960,K01958) K00024 (K01676,K01679,K01677+K01678) (K00239+K00240-K00241-K00242,K00244+K00245-K00246-K00247,K18556+K18557+K18558+K18559+K18560) (K01902+K01903)(K00174+K00175-K00177-K00176) K00031 (K01681,K01682) (K15230+K15231,K15232+K15233 K15234)', 'M00376': '(K02160+K01961+K01962+K01963) K14468 K14469 K15052 K05606 (K01847,K01848+K01849) (K14471+K14472) (K00239+K00240+K00241) K01679 K08691 K14449 K14470 K09709', 'M00375': '(K01964+K15037+K15036) K15017 K15039 K15018 K15019 K15020 K05606 (K01848+K01849) (K15038,K15017) K14465 (K14466,K18861) K14534 K15016 K00626', 'M00374': '(K00169+K00170+K00171+K00172) K01007 K01595 K00024 (K01677+K01678) (K00239+K00240-K00241-K18860) (K01902+K01903) (K15038,K15017) K14465 (K14467,K18861) K14534 K15016 K00626', 'M00377': 'K00198 (K05299-K15022,K22015+K25123+K25124) K01938 K01491-K01500 K00297-K25007-K25008 K15023 K14138+K00197+K00194', 'M00579': '(K00625,K13788,K15024) K00925', 'M00620': '(K00169+K00170+K00171+K00172) (K01959+K01960) K00024 (K01677+K01678) (K18209+K18210) (K01902+K01903) (K00174+K00175+K00176+K00177)'}    
    missing_KO = []
    ko_list = {i: [x for x in re.split("(K\d+|[ \ \ \,\+\-\(\)])", ko_list_defs[i]) if x != "" and x != '+' and x != '(' and x != ')' and  x != '-' and x != ',' and x != ' ' and x != '\n'] for i in ko_list_defs}
    ko_complete_list = {s for v in ko_list.values() for s in v}



#the kofamscan result will be parsed and put in a dataframe
    print("Loading kofamscan output")

    k_text = pd.read_csv('kofamscan.txt', comment="#", header= None, delimiter=r"\s+", engine ="c", usecols=[0,1,2,3,4])
    k_text.columns = ['cluster_id', 'KO', 'Threshold', 'Score', 'Description']
    print("Filtering kofamscan output")
    
    '''
    #this was for the missed out KOs 
    set_new = list((set(ko_complete_list_f) | set(ko_complete_list)) - (set(ko_complete_list_f) & set(ko_complete_list)))
    rem = ['K18559', 'K18558', 'K15233', 'K18860', 'K15017', 'K18861', 'K14469', 'K00194', 'K00198', 'K14138', 'K00197']
    rem_KO = set_new+rem
    '''
    
    
    filter_df = k_text.iloc[ [k.KO in ko_complete_list and float(k.Score)/float(k.Threshold) > 0.5 for i, k in tqdm(k_text.iterrows())]]
#    kf_dict = {'cluster_id': [row.cluster_id.replace('* ','') for i, row in k_text.iterrows() if row.KO in ko_complete_list],'KO': [row.KO for i, row in k_text.iterrows() if row.KO in ko_complete_list]}
#    kf_df = pd.DataFrame(kf_dict)

    #filter the hits whose score to threshold ratio is bigger than 0.5
#    filter_df = pd.DataFrame([(row.cluster_id, row.KO) for i, row in if row.score/row.threshold > 0.5], columns=['cluster_id', 'KO'])
    filter_df.to_csv('all_filt_clusters.csv')
    #drop the threshold and score columns
    filter_df = filter_df.drop(columns = ['Threshold','Score'])
    
    print("Parsing the JSON")
    
    df_members, df_members_wo_cluster = json_parse('arctic_gene_atlas.json', filter_df)
    print("Cpmpiling stuff")

    #add module name to the dataframe based on the KO terms - for grouping and plotting
    res = []
    count_KO = 0
    for key, value_list in tqdm(ko_list.items()):
        key_matches = filter_50[filter_50['KO'].isin(value_list)]
        count_KO += key_matches.count()
        if not key_matches.empty:
            key_results = key_matches.assign(key=key)
            res.append(key_results)

    kf_with_modules = pd.concat(res)
    kf_with_modules.to_csv('kf_with_modules.csv')
    
    
    #grouped data and modules for adding perctent completion to each table
    grouped_df = df_members_wo_cluster.groupby(['members'], as_index = False).agg({'KO': lambda x :  {xx for xx in x}} )
    grouped_df.to_csv('Grouped_kofam.csv')
    
    #look for the KOs in grouped and assign them to respective modules
    mod_complete = {k : { row.members : perct_comp_pathway(row.KO, v) for i, row in grouped_df.iterrows() } for k, v in tqdm(ko_list_defs.items()) if k in c_pathways}
    mod_complete_df = pd.DataFrame.from_dict(mod_complete)               

    #for every unique key, plot heatmap and also display how many elements - this will have a function that executes R script in place of heatmap
    for key_value in kf_with_modules['key'].unique():
        subset_df = kf_with_modules[kf_with_modules['key'] == key_value]
        pivot_df = subset_df.pivot_table(index='KO', columns='members', values='KO', aggfunc="size", fill_value=0).T
        pivot_df = pd.merge(pivot_df, mod_complete_df[[key_value]], left_on='members', right_index=True)
        pivot_df.to_csv(f"{key_value}_table_members_with_perct.csv")


if __name__ == '__main__':
    main()
