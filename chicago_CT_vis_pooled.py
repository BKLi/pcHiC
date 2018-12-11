# -*- coding: utf-8 -*-

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

'''test script, EDA, individual replicates against pooled replicates'''

cell_types = ["MS001","MS002","MS003","JC005","JC006","MS127","MS128"]
pooled_cell_types = ["cortical_1", "astrocyte_1+2"]

for cell_type1 in cell_types:
    
    hicpro_df_ct1 = pd.read_table("C:\\Users\libin\\UCSF\Chicago\\hicpro\{}.ibed".format(cell_type1))
    hicpro_df_ct1.sort_values("score", inplace=True, ascending=False)
   
    plt.figure(figsize=(17,11))

    for cell_type2 in pooled_cell_types:
  
        hicpro_df_ct2 = pd.read_table("C:\\Users\libin\\UCSF\Chicago\\hicpro\pooled\{}.ibed".format(cell_type2))      
        hicpro_df_ct2.sort_values("score", inplace=True, ascending=False)
        
        hicpro_joined = pd.merge(hicpro_df_ct1, hicpro_df_ct2, how="outer", on=["bait_chr", "bait_start", "bait_end", "otherEnd_chr", "otherEnd_start", "otherEnd_end"])
        hicpro_joined.fillna(0, inplace=True)
        hicpro_ct1_filtered_above = hicpro_joined[(hicpro_joined["score_x"] >= 5)]

        sns.distplot(hicpro_ct1_filtered_above["score_y"], hist=False, label="{}-{}".format(cell_type1,cell_type2)).set(xlim=(0,20))
        
        print(cell_type1, hicpro_ct1_filtered_above['score_x'].mean())
        print(cell_type2, hicpro_ct1_filtered_above['score_y'].mean())

    plt.show()


