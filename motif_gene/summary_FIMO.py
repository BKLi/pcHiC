# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 15:45:22 2019

@author: bingkun
"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


cell_type = "cortical"
interactions_final_df = pd.read_table(
        "C:\\Users\libin\\UCSF\motif_gene\\{}_pro_peak_final_df".format(cell_type),sep="\t")

fimo_out = pd.read_table(
        "C:\\Users\libin\\UCSF\motif_gene\\{}_pro_peak_2fimo.tsv".format(cell_type), sep="\t",comment='#')
fimo_out = fimo_out[
        ['motif_id',
         'motif_alt_id',
         'sequence_name',
         'score',
         'p-value']]

# aggregate rows with same sequence names(collapse on position) 
# for merging with interactions(promoter and peak info)
condition={"motif_id":";".join, "motif_alt_id":";".join, "score":"mean", "p-value":"mean"}
fimo_out_agg = (fimo_out
                .groupby('sequence_name', as_index=False)
                .agg(condition))

interaction_motif_merged = pd.merge(interactions_final_df,fimo_out_agg,on=['sequence_name'], how='inner')
interaction_motif_merged = (interaction_motif_merged[
        ['sequence_name',
         'atac_seq_peak_ID',
         'promoter_ids',
         'motif_alt_id',
         #'motif_id',
         'score',
         'p-value',
         'chrom',
         'chromStart',
         'chromEnd',
         'interaction_ID']]
.rename(index=str, columns={"sequence_name": 'atac_seq_peak_pos'}))

# split on motif id
# result would have one motif per line
interaction_motif_merged_split = \
(interaction_motif_merged.set_index(
        interaction_motif_merged.columns.drop('motif_alt_id',1).tolist())
        ['motif_alt_id'].str.split(";", expand=True)
        .stack()
        .reset_index()
        .rename(columns={0: "motif_alt_id"})
        .loc[:,interaction_motif_merged.columns])

# group by motif id
# each motif are associated with multiple promoters
interaction_motif_merged_split_cut = interaction_motif_merged_split[
        ['promoter_ids',
         'motif_alt_id',
         'atac_seq_peak_ID']]

# group by distal peaks
interaction_motif_merged_split_cut_group_by_peak = interaction_motif_merged_split_cut.groupby("atac_seq_peak_ID")
peak_dict = {}
for name, group in interaction_motif_merged_split_cut_group_by_peak:
    group_split = (group.set_index("atac_seq_peak_ID")
                ['motif_alt_id'].str.split(",", expand=True)
                .stack()
                .reset_index()
                .rename(columns={0:"motif_alt_id"})
                .loc[:,group.columns])
    peak_dict[name] = len(group_split["motif_alt_id"].value_counts())
peak_dict_df = pd.DataFrame.from_dict(peak_dict, orient="index").reset_index()\
               .rename(columns={"index": "atac_seq_peak_ID", 0: "motif number"})
plt.figure()
sns.distplot(peak_dict_df["motif number"],color="red", kde=False)    

# not sure what format like "TFAP2B(var.2)" means -- ignoring for now; treating as different motifs
# dictionary to store how many types of promoters connect with each motif
interaction_motif_merged_split_cut_grouped = interaction_motif_merged_split_cut.groupby('motif_alt_id')
motif_dict = {} 
for name, group in interaction_motif_merged_split_cut_grouped:    
    # expand on promoter id
    group_split = (group.set_index('motif_alt_id')
                ['promoter_ids'].str.split(",", expand=True)
                .stack()
                .reset_index()
                .rename(columns={0: 'promoter_ids'})
                .loc[:,group.columns])
    motif_dict[name] = len(group_split["promoter_ids"].value_counts())
    # print(len(group_split["promoter_ids"].value_counts()))
    # print(group_split["promoter_ids"].shape[0])
motif_dict_df = (pd.DataFrame.from_dict(motif_dict, orient="index").reset_index()
                .rename(columns={"index":"motif_id", 0:"promoter_occurence"}))

plt.figure()
sns.distplot(motif_dict_df["promoter_occurence"],color="green", kde=False)
plt.title('{}: # of associated promoters per motif'.format(cell_type))
plt.savefig('C:\\Users\libin\\UCSF\motif_gene\\{}_promoter_per_motif.pdf'.format(cell_type), 
            transparent=True, bbox_inches = 'tight',
            pad_inches = 0.1)

# ------- experimenting downstream analysis --------
# count the presence of each motif
# somehow the output df is already sorted
fimo_out_motif_count = fimo_out["motif_alt_id"].value_counts().to_dict()
fimo_out_motif_count_df = pd.DataFrame.from_dict(fimo_out_motif_count, orient="index") \
                                                .reset_index() \
                                                .rename(columns={"index":"motif_alt_id", 0:"occurence"})
# Todo: collapse on motif ID
motif_interactions_collapse_on_motif = \
interaction_motif_merged_split_cut \
.groupby("motif_alt_id", as_index=False) \
.agg(",".join)
                                              
# select the first 100 most presented motifs
most_occured_motif = fimo_out_motif_count_df.loc[0:99, :]   
most_occured_motif = pd.merge(most_occured_motif, motif_interactions_collapse_on_motif)                                        
