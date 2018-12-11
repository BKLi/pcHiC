# -*- coding: utf-8 -*-

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

''' testing 'number of significant interaction' approach'''

cell_types = ["MS001", "MS002", "MS003", "JC005", "JC006", "MS127", "MS128"]

for cell_type1 in cell_types:
    
    hicpro_df_ct1 = pd.read_table("C:\\Users\libin\\UCSF\Chicago\\hicpro\{}.ibed".format(cell_type1))
    hicpro_df_ct1.sort_values("score", inplace=True, ascending=False)

    ct1_index = cell_types.index(cell_type1)
    
    cell_type2_list = cell_types[:ct1_index]+cell_types[ct1_index+1:]
    for cell_type2 in cell_type2_list:
  
        hicpro_df_ct2 = pd.read_table("C:\\Users\libin\\UCSF\Chicago\\hicpro\{}.ibed".format(cell_type2))      
        hicpro_df_ct2.sort_values("score", inplace=True, ascending=False)
        
        hicpro_joined = pd.merge(hicpro_df_ct1, hicpro_df_ct2, how="outer", on=["bait_chr", "bait_start", "bait_end", "otherEnd_chr", "otherEnd_start", "otherEnd_end"])
        hicpro_joined.fillna(0, inplace=True)

        hicpro_ct1_filtered_above = hicpro_joined[(hicpro_joined["score_x"] >= 3)]  
        hicpro_ct2_filtered_above = hicpro_ct1_filtered_above[(hicpro_ct1_filtered_above["score_y"] >= 3)]
        
        print(cell_type1, hicpro_ct1_filtered_above.shape[0])
        print(cell_type2, hicpro_ct2_filtered_above.shape[0])
