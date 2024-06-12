# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 10:16:58 2023

@author: kayve
"""
###execute it in ipython


c_pathways = ['M00173','M00376','M00375','M00374','M00377','M00579', 'M00620']
ko_list_defs = {'M00173': '(K00169+K00170+K00171+K00172,K03737) ((K01007,K01006) K01595,K01959+K01960,K01958) K00024 (K01676,K01679,K01677+K01678) (K00239+K00240-K00241-K00242,K00244+K00245-K00246-K00247,K18556+K18557+K18558+K18559+K18560) (K01902+K01903)(K00174+K00175-K00177-K00176) K00031 (K01681,K01682) (K15230+K15231,K15232+K15233 K15234)', 'M00376': '(K02160+K01961+K01962+K01963) K14468 K14469 K15052 K05606 (K01847,K01848+K01849) (K14471+K14472) (K00239+K00240+K00241) K01679 K08691 K14449 K14470 K09709', 'M00375': '(K01964+K15037+K15036) K15017 K15039 K15018 K15019 K15020 K05606 (K01848+K01849) (K15038,K15017) K14465 (K14466,K18861) K14534 K15016 K00626', 'M00374': '(K00169+K00170+K00171+K00172) K01007 K01595 K00024 (K01677+K01678) (K00239+K00240-K00241-K18860) (K01902+K01903) (K15038,K15017) K14465 (K14467,K18861) K14534 K15016 K00626', 'M00377': 'K00198 (K05299-K15022,K22015+K25123+K25124) K01938 K01491-K01500 K00297-K25007-K25008 K15023 K14138+K00197+K00194', 'M00579': '(K00625,K13788,K15024) K00925', 'M00620': '(K00169+K00170+K00171+K00172) (K01959+K01960) K00024 (K01677+K01678) (K18209+K18210) (K01902+K01903) (K00174+K00175+K00176+K00177)'}
ko_list = {i: [x for x in re.split("(K\d+|[ \ \ \,\+\-\(\)])", ko_list_defs[i]) if x != "" and x != '+' and x != '(' and x != ')' and  x != '-' and x != ',' and x != ' ' and x != '\n'] for i in ko_list_defs}
ko_complete_list = {s for v in ko_list.values() for s in v}


#the df name in the place of ktext fin
filter_df = k_text.iloc[ [k.KO in ko_complete_list and float(k.Score)/float(k.Threshold) > 0.5 for i, k in tqdm(k_text.iterrows())]]

#get mem for clusterID (running)

res = []

for key, value_list in tqdm(ko_list.items()):
        key_matches = df_members_wo_cluster[df_members_wo_cluster['KO'].isin(value_list)]
        #count_KO += key_matches.count()
        if not key_matches.empty:
            key_results = key_matches.assign(key=key)
            res.append(key_results)
kf_with_modules = pd.concat(res)

#add to the big dataframe
#extracting from the somethi

grouped_df = df_members_wo_cluster.groupby(['members'], as_index = False).agg({'KO': lambda x :  {xx for xx in x}} )
grouped_df.to_csv('Grouped_kofam.csv')


mod_complete = {k : { row.members : perct_comp_pathway(row.KO, v) for i, row in grouped_df.iterrows() } for k, v in tqdm(ko_list_defs.items()) if k in c_pathways}
mod_complete_df = pd.DataFrame.from_dict(mod_complete) 


for key_value in kf_with_modules['key'].unique():
        subset_df = kf_with_modules[kf_with_modules['key'] == key_value]
        pivot_df = subset_df.pivot_table(index='KO', columns='members', values='KO', aggfunc="size", fill_value=0).T
        pivot_df = pd.merge(pivot_df, mod_complete_df[[key_value]], left_on='members', right_index=True)
        pivot_df.to_csv(f"{key_value}_table_members_with_perct.csv")
        heatmap_simple(pivot_df, key_value)