# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 15:45:22 2019

@author: bingkun
"""
import pandas as pd


cell_type = "cortical"
interactions_final_df = pd.read_table("C:\\Users\libin\\UCSF\motif_gene\\{}_pro_peak_final_df".format(cell_type),sep="\t")

fimo_out = pd.read_table("C:\\Users\libin\\UCSF\motif_gene\\{}_pro_peak_2fimo.tsv".format(cell_type), sep="\t",comment='#')
fimo_out = fimo_out[
        ['motif_id',
         'motif_alt_id',
         'sequence_name',
         'score',
         'p-value']]

# aggregate rows with same sequence names(collapse on position)
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
interaction_motif_merged_split = \
(interaction_motif_merged.set_index(
        interaction_motif_merged.columns.drop('motif_alt_id',1).tolist())
        ['motif_alt_id'].str.split(";", expand=True)
        .stack()
        .reset_index()
        .rename(columns={0: "motif_alt_id"})
        .loc[:,interaction_motif_merged.columns])

# group by motif id
interaction_motif_merged_split_cut = interaction_motif_merged_split[
        ['promoter_ids',
         'motif_alt_id',]]
interaction_motif_merged_split_cut_grouped = interaction_motif_merged_split_cut.groupby('motif_alt_id')
# not sure what format like TFAP2B(var.2) means -- ignoring for now; treating as different motifs
# dictionary to store how many types of promoters connect with each motif
motif_dict = {} 
for name, group in interaction_motif_merged_split_cut_grouped:    
    # split on promoter id
    group_split = (group.set_index('motif_alt_id')
                ['promoter_ids'].str.split(",", expand=True)
                .stack()
                .reset_index()
                .rename(columns={0: 'promoter_ids'})
                .loc[:,group.columns])
    motif_dict[name] = len(group_split["promoter_ids"].value_counts())
    # print(len(group_split["promoter_ids"].value_counts()))
    # print(group_split["promoter_ids"].shape[0])
