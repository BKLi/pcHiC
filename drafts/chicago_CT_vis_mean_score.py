# -*- coding: utf-8 -*-
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

''' 'mean score' approach to assess reproducibility'''


cell_types = ["JC005","JC006","MS001","MS002","MS003","MS006","MS032","MS126","MS127",
              "MS128","MS133","MS134","MS136","MS137","MS140","MS141","MS142"]
cutoff = 5

# store mean scores for all cell types
mean_scores_matrix = {}

for cell_type1 in cell_types:
    # create nested dictionary for each cell type
    mean_scores_matrix[cell_type1] = {}
    
    hicpro_df_ct1 = pd.read_table("C:\\Users\libin\\UCSF\Chicago\\hicpro\{}.ibed".format(cell_type1))
    hicpro_df_ct1.sort_values("score", inplace=True, ascending=False)
    print(cell_type1, "full data", hicpro_df_ct1["score"].mean())
    
    # plt.figure(figsize=(17,11))
    ct1_index = cell_types.index(cell_type1)
    
    cell_type2_list = cell_types[:ct1_index]+cell_types[ct1_index+1:]
    for cell_type2 in cell_type2_list:

        hicpro_df_ct2 = pd.read_table("C:\\Users\libin\\UCSF\Chicago\\hicpro\{}.ibed".format(cell_type2))      
        hicpro_df_ct2.sort_values("score", inplace=True, ascending=False)
        
        hicpro_joined = pd.merge(hicpro_df_ct1, hicpro_df_ct2, how="outer", on=["bait_chr", "bait_start", "bait_end", "otherEnd_chr", "otherEnd_start", "otherEnd_end"])
        hicpro_joined.fillna(0, inplace=True)
        hicpro_ct1_filtered_above = hicpro_joined[(hicpro_joined["score_x"] >= cutoff)]
        
        ct1_mean_score = hicpro_ct1_filtered_above['score_x'].mean()
        ct2_mean_score = hicpro_ct1_filtered_above['score_y'].mean()
        mean_scores_matrix[cell_type1][cell_type1] = ct1_mean_score/ct1_mean_score
        mean_scores_matrix[cell_type1][cell_type2] = ct2_mean_score/ct1_mean_score

        print(cell_type1, cell_type2, ct1_mean_score, ct2_mean_score)
        score_df = pd.DataFrame(mean_scores_matrix)
        score_df_trans = score_df.transpose()
        score_df_sum = score_df.add(score_df_trans, fill_value=0)
        score_df_final = score_df_sum/2

        # plot heatmap
        plt.figure(figsize=(20,20))
        score_mean = (score_df_final.values.sum() - score_df_final.shape[0])/(score_df_final.shape[0]*(score_df_final.shape[0]-1))
        sns.heatmap(score_df_final, score_df_final.values.min(), vmax=score_df_final.values.max(), center=score_mean, cmap="RdBu_r", annot=True)
        plt.title("mean_score correlation matrix", y=1.01, x=1.1, fontsize=13)
        for item in score_df_final.get_yticklabels():
            item.set_rotation(45)


