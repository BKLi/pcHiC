# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 14:30:38 2019

@author: bingkun
"""

import pandas as pd

cortical_peaks = pd.read_table("C:\\Users\libin\\UCSF\motif_gene\\atac.seq.peaks.cortical")
cortical_peaks = cortical_peaks[
 ['chrom',
 'chromStart',
 'chromEnd',
 'atac_seq_peak_ID']]

cortical_interactions = pd.read_table("C:\\Users\libin\\UCSF\motif_gene\\interactions.sig.res.cortical")
cortical_interactions_ann = pd.read_table("C:\\Users\libin\\UCSF\motif_gene\\interactions.sig.ann.cortical")
cortical_interactions_ann = \
cortical_interactions_ann[['promoter_lhs', 'promoter_lhs_ids', 'promoter_other_lhs', 'promoter_other_lhs_ids',
                           'promoter_rhs', 'promoter_rhs_ids', 'promoter_other_rhs', 'promoter_other_rhs_ids',
                           'distal_ATAC.seq_lhs', 'distal_ATAC.seq_lhs_ids', 
                           'distal_ATAC.seq_rhs', 'distal_ATAC.seq_rhs_ids']]
cortical_interactions_ann_lpromotor_only = \
cortical_interactions_ann.loc[((cortical_interactions_ann["promoter_lhs"] > 0) 
                          | (cortical_interactions_ann["promoter_other_lhs"] > 0))
                          & (cortical_interactions_ann['promoter_rhs'] == 0) 
                          & (cortical_interactions_ann['promoter_other_rhs'] == 0)]
cortical_interactions_ann_rpromotor_only = \
cortical_interactions_ann.loc[((cortical_interactions_ann['promoter_rhs'] > 0)  
                          | (cortical_interactions_ann['promoter_other_rhs'] > 0)) 
                          & (cortical_interactions_ann['promoter_lhs'] == 0) 
                          & (cortical_interactions_ann['promoter_other_lhs'] == 0)]

cortical_interactions_ann_lpromotor_rdpeaks = \
cortical_interactions_ann_lpromotor_only.loc[cortical_interactions_ann_lpromotor_only["distal_ATAC.seq_rhs"]>0]
cortical_interactions_ann_lpromotor_rdpeaks_full = \
pd.merge(cortical_interactions, cortical_interactions_ann_lpromotor_rdpeaks, left_index=True, right_index=True)

# split rows with multiple distal peaks into multiple lines
# https://riptutorial.com/pandas/example/25462/split--reshape--csv-strings-in-columns-into-multiple-rows--having-one-element-per-row
cortical_interactions_ann_lpromotor_rdpeaks_full_split = \
(cortical_interactions_ann_lpromotor_rdpeaks_full.set_index(
cortical_interactions_ann_lpromotor_rdpeaks_full.columns.drop("distal_ATAC.seq_rhs_ids",1).tolist())
["distal_ATAC.seq_rhs_ids"].str.split(",", expand=True)
.stack()
.reset_index()
.rename(columns={0:"distal_ATAC.seq_rhs_ids"}) 
.loc[:,cortical_interactions_ann_lpromotor_rdpeaks_full.columns]                                                       
)

# discard unimportant columns
cortical_interactions_ann_lpromotor_rdpeaks_full_split = \
(cortical_interactions_ann_lpromotor_rdpeaks_full_split[[
'interaction_ID',
 'promoter_lhs_ids',
 'promoter_other_lhs_ids',
 'distal_ATAC.seq_rhs_ids']])

# remove duplicates. I don't think I need to take into account how many 
# interactions support certain motif-gene correlations, at least for now.
cortical_interactions_ann_lpromotor_rdpeaks_full_split_dedup = \
cortical_interactions_ann_lpromotor_rdpeaks_full_split.drop_duplicates(
subset=['promoter_lhs_ids',
        'promoter_other_lhs_ids',
        'distal_ATAC.seq_rhs_ids'],
keep='first')

# merge rows with same peak_ids
cortical_interactions_ann_lpromotor_rdpeaks_full_split_dedup_merge = \
(cortical_interactions_ann_lpromotor_rdpeaks_full_split_dedup
.fillna("0")
.groupby('distal_ATAC.seq_rhs_ids', as_index=False)
.agg(",".join))

# to simplify the model, combining lhs/lhs_otherend promoters
cortical_interactions_ann_lpromotor_rdpeaks_full_split_dedup_merge["promoter_all_lhs_ids"] = \
(cortical_interactions_ann_lpromotor_rdpeaks_full_split_dedup_merge["promoter_lhs_ids"]
+ "," + cortical_interactions_ann_lpromotor_rdpeaks_full_split_dedup_merge["promoter_other_lhs_ids"])
cortical_interactions_ann_lpromotor_rdpeaks_full_split_dedup_merge = \
(cortical_interactions_ann_lpromotor_rdpeaks_full_split_dedup_merge[
['distal_ATAC.seq_rhs_ids',
 'interaction_ID',
 'promoter_all_lhs_ids']])


#cortical_interactions_ann_promoter = cortical_interactions_ann.loc[(cortical_interactions_ann["promoter_lhs"] >0) | (cortical_interactions_ann["promoter_other_lhs"] >0)\
#                                                                 | (cortical_interactions_ann['promoter_rhs'] >0) | (cortical_interactions_ann['promoter_other_rhs']>0)]
#cortical_interactions_ann_promoter_bothhand = cortical_interactions_ann_promoter.loc[((cortical_interactions_ann_promoter["promoter_lhs"] >0) | (cortical_interactions_ann_promoter["promoter_other_lhs"] >0))\
#                                                                                   & ((cortical_interactions_ann_promoter['promoter_rhs'] >0) | (cortical_interactions_ann_promoter['promoter_other_rhs'] >0))]

# interactions with promoters on one end and distal peaks on the other
#cortical_interactions_ann_promoter_dpeak = cortical_interactions_ann_promoter.loc[(((cortical_interactions_ann_promoter["promoter_lhs"] >0) | (cortical_interactions_ann_promoter["promoter_other_lhs"] >0))\
#                                                                                 & (cortical_interactions_ann_promoter["distal_ATAC.seq_rhs"]>0))\
#                                                                                 |(((cortical_interactions_ann_promoter['promoter_rhs'] >0) | (cortical_interactions_ann_promoter['promoter_other_rhs'] >0))\
#                                                                                 & (cortical_interactions_ann_promoter["distal_ATAC.seq_lhs"]>0))]

