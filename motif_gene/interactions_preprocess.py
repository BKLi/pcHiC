# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 14:30:38 2019

@author: bingkun
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

cell_type="cortical"

ctype_peaks = pd.read_table("C:\\Users\libin\\UCSF\motif_gene\\atac.seq.peaks.{}".format(cell_type))
ctype_peaks = ctype_peaks[
                 ['chrom',
                 'chromStart',
                 'chromEnd',
                 'atac_seq_peak_ID']]

ctype_interactions = pd.read_table("C:\\Users\libin\\UCSF\motif_gene\\interactions.sig.res.{}".format(cell_type))
ctype_interactions_ann = pd.read_table("C:\\Users\libin\\UCSF\motif_gene\\interactions.sig.ann.{}".format(cell_type))
# promoter_other_lhs/rhs means other types of promoters which we are not interested in at the moment
ctype_interactions_ann = \
ctype_interactions_ann[['promoter_lhs', 'promoter_lhs_ids',
                        'promoter_rhs', 'promoter_rhs_ids',
                        'distal_ATAC.seq_lhs', 'distal_ATAC.seq_lhs_ids', 
                        'distal_ATAC.seq_rhs', 'distal_ATAC.seq_rhs_ids']]
ctype_interactions_ann_lpromotor_only = \
ctype_interactions_ann.loc[(ctype_interactions_ann["promoter_lhs"] > 0)
                           & (ctype_interactions_ann['promoter_rhs'] == 0)]
ctype_interactions_ann_rpromotor_only = \
ctype_interactions_ann.loc[(ctype_interactions_ann['promoter_rhs'] > 0)
                           & (ctype_interactions_ann['promoter_lhs'] == 0)]

ctype_interactions_ann_lpromotor_rdpeaks = \
ctype_interactions_ann_lpromotor_only.loc[ctype_interactions_ann_lpromotor_only["distal_ATAC.seq_rhs"] > 0]
ctype_interactions_ann_lpromotor_rdpeaks_full = \
pd.merge(ctype_interactions, ctype_interactions_ann_lpromotor_rdpeaks, left_index=True, right_index=True)

ctype_interactions_ann_rpromotor_ldpeaks = \
ctype_interactions_ann_rpromotor_only.loc[ctype_interactions_ann_rpromotor_only["distal_ATAC.seq_lhs"] > 0]
ctype_interactions_ann_rpromotor_ldpeaks_full = \
pd.merge(ctype_interactions, ctype_interactions_ann_rpromotor_ldpeaks, left_index=True, right_index=True)

# split rows with multiple distal peaks into multiple lines
# https://riptutorial.com/pandas/example/25462/split--reshape--csv-strings-in-columns-into-multiple-rows--having-one-element-per-row
ctype_interactions_ann_lpromotor_rdpeaks_full_split = \
(ctype_interactions_ann_lpromotor_rdpeaks_full.set_index(
ctype_interactions_ann_lpromotor_rdpeaks_full.columns.drop("distal_ATAC.seq_rhs_ids", 1).tolist())
["distal_ATAC.seq_rhs_ids"].str.split(",", expand=True)
.stack()
.reset_index()
.rename(columns={0:"distal_ATAC.seq_rhs_ids"}) 
.loc[:, ctype_interactions_ann_lpromotor_rdpeaks_full.columns]
)

ctype_interactions_ann_rpromotor_ldpeaks_full_split = \
(ctype_interactions_ann_rpromotor_ldpeaks_full.set_index(
ctype_interactions_ann_rpromotor_ldpeaks_full.columns.drop("distal_ATAC.seq_lhs_ids", 1).tolist())
["distal_ATAC.seq_lhs_ids"].str.split(",", expand=True)
.stack()
.reset_index()
.rename(columns={0:"distal_ATAC.seq_lhs_ids"}) 
.loc[:, ctype_interactions_ann_rpromotor_ldpeaks_full.columns]
)

# discard unimportant columns
ctype_interactions_ann_lpromotor_rdpeaks_full_split = \
(ctype_interactions_ann_lpromotor_rdpeaks_full_split[[
 'interaction_ID',
 'promoter_lhs_ids',
 'distal_ATAC.seq_rhs_ids']])

ctype_interactions_ann_rpromotor_ldpeaks_full_split = \
(ctype_interactions_ann_rpromotor_ldpeaks_full_split[[
 'interaction_ID',
 'promoter_rhs_ids',
 'distal_ATAC.seq_lhs_ids']])

# remove duplicates. I don't think I need to take into account how many 
# interactions support certain motif-gene correlations, at least for now.
# also rename to concat later
ctype_interactions_ann_lpromotor_rdpeaks_full_split_dedup = \
(ctype_interactions_ann_lpromotor_rdpeaks_full_split.drop_duplicates(
subset=['promoter_lhs_ids',
        'distal_ATAC.seq_rhs_ids'],
keep='first')
.rename(index=str, columns={"distal_ATAC.seq_rhs_ids": 'distal_ATAC.seq_ids',
                            "promoter_lhs_ids": "promoter_ids"}))

ctype_interactions_ann_rpromotor_ldpeaks_full_split_dedup = \
(ctype_interactions_ann_rpromotor_ldpeaks_full_split.drop_duplicates(
subset=['promoter_rhs_ids',
        'distal_ATAC.seq_lhs_ids'],
keep='first')
.rename(index=str, columns={"distal_ATAC.seq_lhs_ids": 'distal_ATAC.seq_ids',
                            "promoter_rhs_ids": "promoter_ids"}))

# merge promoters on lhs and rhs
ctype_interactions_ann_promotor_dpeaks_full_split_dedup = \
pd.concat([ctype_interactions_ann_lpromotor_rdpeaks_full_split_dedup,
          ctype_interactions_ann_rpromotor_ldpeaks_full_split_dedup],
          ignore_index=True)

# merge rows with same peak_ids
ctype_interactions_ann_promotor_dpeaks_full_split_dedup_collapse = \
(ctype_interactions_ann_promotor_dpeaks_full_split_dedup
 .fillna("0")
 .groupby('distal_ATAC.seq_ids', as_index=False)
 .agg(",".join))

# to simplify the model, combining lhs/lhs_otherend promoters
# cortical_interactions_ann_lpromotor_rdpeaks_full_split_dedup_merge["promoter_all_lhs_ids"] = \
# (cortical_interactions_ann_lpromotor_rdpeaks_full_split_dedup_merge["promoter_lhs_ids"]
# + "," + cortical_interactions_ann_lpromotor_rdpeaks_full_split_dedup_merge["promoter_other_lhs_ids"])

# discard unimportant coulmns & rename column for merging peak file
ctype_interactions_ann_promotor_dpeaks_full_split_dedup_collapse = \
(ctype_interactions_ann_promotor_dpeaks_full_split_dedup_collapse[
['distal_ATAC.seq_ids',
 'interaction_ID',
 'promoter_ids']]
.rename(index=str, columns={"distal_ATAC.seq_ids": 'atac_seq_peak_ID'}))

ctype_pro_peak_final_df = pd.merge(
        ctype_interactions_ann_promotor_dpeaks_full_split_dedup_collapse,
        ctype_peaks, on=["atac_seq_peak_ID"], how="inner")

# add seqname column to merge with FIMO putput later
ctype_pro_peak_final_df['sequence_name'] = \
(ctype_pro_peak_final_df["chrom"].map(str) +
 ":" +
 ctype_pro_peak_final_df["chromStart"].map(str) +
 "-" +
 ctype_pro_peak_final_df["chromEnd"].map(str))
ctype_pro_peak_final_df.to_csv("C:\\Users\libin\\UCSF\motif_gene\\{}_pro_peak_final_df".format(cell_type),
                                 sep="\t", index=False, header=True)

# output chr,start,end as FIMO input
ctype_pro_peak_2fimo = (ctype_pro_peak_final_df[
        ['chrom',
         'chromStart',
         'chromEnd']]
.to_csv("C:\\Users\libin\\UCSF\motif_gene\\{}_pro_peak_2fimo.bed".format(cell_type), sep="\t", index=False, header=False))

# next step: bedtools getfasta

# split & group by promoter to get basic statistics
ctype_interactions_ann_promotor_dpeaks_full_dedup_doub_split = \
(ctype_interactions_ann_promotor_dpeaks_full_split_dedup.set_index(
ctype_interactions_ann_promotor_dpeaks_full_split_dedup.columns.drop("promoter_ids", 1).tolist())
["promoter_ids"].str.split(",", expand=True)
.stack()
.reset_index()
.rename(columns={0:"promoter_ids"}) 
.loc[:, ctype_interactions_ann_promotor_dpeaks_full_split_dedup.columns]
)

ctype_interactions_ann_promotor_dpeaks_full_dedup_doub_split_grouped = \
(ctype_interactions_ann_promotor_dpeaks_full_dedup_doub_split
.groupby('promoter_ids', as_index=False))
promoter_dict = {}
for name, group in ctype_interactions_ann_promotor_dpeaks_full_dedup_doub_split_grouped:
    promoter_dict[name] = len(group["distal_ATAC.seq_ids"].value_counts())
promoter_dict_df = pd.DataFrame.from_dict(promoter_dict, orient='index').reset_index()
promoter_dict_df = promoter_dict_df.rename(columns={"index":"promoter_id", 0:"atac_seq_peak_occurence"})

data_range = np.arange(promoter_dict_df["atac_seq_peak_occurence"].min(), 
                       promoter_dict_df["atac_seq_peak_occurence"].max(),1) - 0.5
sns.distplot(promoter_dict_df["atac_seq_peak_occurence"],
             bins= data_range, color="blue",
             kde=False, hist_kws={'edgecolor':'black'}).set_xticks(data_range - 0.5)
#cortical_interactions_ann_promoter = cortical_interactions_ann.loc[(cortical_interactions_ann["promoter_lhs"] >0) | (cortical_interactions_ann["promoter_other_lhs"] >0)\
#                                                                 | (cortical_interactions_ann['promoter_rhs'] >0) | (cortical_interactions_ann['promoter_other_rhs']>0)]
#cortical_interactions_ann_promoter_bothhand = cortical_interactions_ann_promoter.loc[((cortical_interactions_ann_promoter["promoter_lhs"] >0) | (cortical_interactions_ann_promoter["promoter_other_lhs"] >0))\
#                                                                                   & ((cortical_interactions_ann_promoter['promoter_rhs'] >0) | (cortical_interactions_ann_promoter['promoter_other_rhs'] >0))]

# interactions with promoters on one end and distal peaks on the other
#cortical_interactions_ann_promoter_dpeak = cortical_interactions_ann_promoter.loc[(((cortical_interactions_ann_promoter["promoter_lhs"] >0) | (cortical_interactions_ann_promoter["promoter_other_lhs"] >0))\
#                                                                                 & (cortical_interactions_ann_promoter["distal_ATAC.seq_rhs"]>0))\
#                                                                                 |(((cortical_interactions_ann_promoter['promoter_rhs'] >0) | (cortical_interactions_ann_promoter['promoter_other_rhs'] >0))\
#                                                                                 & (cortical_interactions_ann_promoter["distal_ATAC.seq_lhs"]>0))]