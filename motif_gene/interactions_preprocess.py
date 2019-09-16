# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 14:30:38 2019

@author: bingkun
Step one of motif-gene analysis
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

cell_type = "cortical"

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

# process lh and rh separately and merge later for clearance of code(I think)
ctype_interactions_ann_lpromotor_only = \
ctype_interactions_ann.loc[(ctype_interactions_ann["promoter_lhs"] > 0)
                           & (ctype_interactions_ann['promoter_rhs'] == 0)]
ctype_interactions_ann_rpromotor_only = \
ctype_interactions_ann.loc[(ctype_interactions_ann['promoter_rhs'] > 0)
                           & (ctype_interactions_ann['promoter_lhs'] == 0)]

# I don't think it matters whether there are lhs dpeaks when considering lhs promoters
ctype_interactions_ann_lpromotor_rdpeaks = \
ctype_interactions_ann_lpromotor_only.loc[ctype_interactions_ann_lpromotor_only["distal_ATAC.seq_rhs"] > 0]
# merge by index
ctype_interactions_ann_lpromotor_rdpeaks_full = \
pd.merge(ctype_interactions, ctype_interactions_ann_lpromotor_rdpeaks, left_index=True, right_index=True)

ctype_interactions_ann_rpromotor_ldpeaks = \
ctype_interactions_ann_rpromotor_only.loc[ctype_interactions_ann_rpromotor_only["distal_ATAC.seq_lhs"] > 0]
ctype_interactions_ann_rpromotor_ldpeaks_full = \
pd.merge(ctype_interactions, ctype_interactions_ann_rpromotor_ldpeaks, left_index=True, right_index=True)

# split rows with multiple distal peaks into multiple lines
# The results have one distal peak per line which may correspond to multiple promoters
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

# remove unimportant columns
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

# discard unimportant columns & rename column for merging peak file
ctype_interactions_ann_promotor_dpeaks_full_split_dedup_collapse = \
(ctype_interactions_ann_promotor_dpeaks_full_split_dedup_collapse[
['distal_ATAC.seq_ids',
 'interaction_ID',
 'promoter_ids']]
.rename(index=str, columns={"distal_ATAC.seq_ids": 'atac_seq_peak_ID'}))

ctype_pro_peak_final_df = pd.merge(
        ctype_interactions_ann_promotor_dpeaks_full_split_dedup_collapse,
        ctype_peaks, on=["atac_seq_peak_ID"], how="inner")

# write input for next step: add seqname column for merging with FIMO output later
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

## split & group by promoter to get basic statistics
# one promoter-dpeaks correspondence for each line
# one promoter per sub-group
ctype_interactions_ann_promotor_dpeaks_full_dedup_doub_split = \
(ctype_interactions_ann_promotor_dpeaks_full_split_dedup.set_index(
ctype_interactions_ann_promotor_dpeaks_full_split_dedup.columns.drop("promoter_ids", 1).tolist())
["promoter_ids"].str.split(",", expand=True)
.stack()
.reset_index()
.rename(columns={0:"promoter_ids"}) 
.loc[:, ctype_interactions_ann_promotor_dpeaks_full_split_dedup.columns]
)

ctype_interactions_ann_promotor_dpeaks_full_dedup_doub_split_grouped_by_promoter = \
(ctype_interactions_ann_promotor_dpeaks_full_dedup_doub_split
.groupby('promoter_ids', as_index=False))
# count how many distal peaks are each promoter related with
promoter_dict = {}
for name, group in ctype_interactions_ann_promotor_dpeaks_full_dedup_doub_split_grouped_by_promoter:
    promoter_dict[name] = len(group["distal_ATAC.seq_ids"].value_counts())
promoter_dict_df = pd.DataFrame.from_dict(promoter_dict, orient='index').reset_index()
promoter_dict_df = promoter_dict_df.rename(columns={"index":"promoter_id", 0:"atac_seq_peak_occurence"})

# https://stackoverflow.com/questions/27083051/matplotlib-xticks-not-lining-up-with-histogram
data_range = np.arange(promoter_dict_df["atac_seq_peak_occurence"].min(), 
                       promoter_dict_df["atac_seq_peak_occurence"].max(),1) - 0.5
                       
plt.figure()
sns.distplot(promoter_dict_df["atac_seq_peak_occurence"],
             bins= data_range, color="blue",
             kde=False, hist_kws={'edgecolor':'black'}).set_xticks(data_range - 0.5)
plt.title('{}: # of associated distal peaks per promoter'.format(cell_type))
plt.savefig('C:\\Users\libin\\UCSF\motif_gene\\{}_peak_per_promoter.pdf'.format(cell_type), transparent=True, bbox_inches = 'tight',
    pad_inches = 0.1)

print(cell_type)
print("atac_seq_peak_occurence", promoter_dict_df["atac_seq_peak_occurence"].describe())

# To do : collapse on promoter for cross_ctype comparison later

# group by distal peaks to get basic statistics
ctype_interactions_ann_promotor_dpeaks_full_dedup_doub_split_grouped_by_peak = \
(ctype_interactions_ann_promotor_dpeaks_full_dedup_doub_split
.groupby('distal_ATAC.seq_ids', as_index=False))
# count how many promoters are each distal peaks related with
distal_dict = {}
for name, group in ctype_interactions_ann_promotor_dpeaks_full_dedup_doub_split_grouped_by_peak:
    distal_dict[name] = len(group["promoter_ids"].value_counts())
distal_dict_df = pd.DataFrame.from_dict(distal_dict, orient='index').reset_index()
distal_dict_df = distal_dict_df.rename(columns={"index": "peak_id", 0:"promoter_occurence"})
data_range2 = np.arange(distal_dict_df["promoter_occurence"].min(),
                        distal_dict_df["promoter_occurence"].max(),1) - 0.5
plt.figure()
sns.distplot(distal_dict_df["promoter_occurence"],
             bins= data_range, color="orange",
             kde=False, hist_kws={'edgecolor':'white'}).set_xticks(data_range - 0.5)
plt.title('{}: # of associated promoters per distal peak'.format(cell_type))
plt.savefig('C:\\Users\libin\\UCSF\motif_gene\\{}_promoter_per_peak.pdf'.format(cell_type), transparent=True, bbox_inches = 'tight',
    pad_inches = 0.1)

# deprecated code below
# -------------------------------
#cortical_interactions_ann_promoter = \
# cortical_interactions_ann.loc[(cortical_interactions_ann["promoter_lhs"] >0) | (cortical_interactions_ann["promoter_other_lhs"] >0)\
#                               | (cortical_interactions_ann['promoter_rhs'] >0) | (cortical_interactions_ann['promoter_other_rhs']>0)]
#cortical_interactions_ann_promoter_bothhand = cortical_interactions_ann_promoter.loc[((cortical_interactions_ann_promoter["promoter_lhs"] >0) | (cortical_interactions_ann_promoter["promoter_other_lhs"] >0))\
#                                                                                   & ((cortical_interactions_ann_promoter['promoter_rhs'] >0) | (cortical_interactions_ann_promoter['promoter_other_rhs'] >0))]

# interactions with promoters on one end and distal peaks on the other
#cortical_interactions_ann_promoter_dpeak = cortical_interactions_ann_promoter.loc[(((cortical_interactions_ann_promoter["promoter_lhs"] >0) | (cortical_interactions_ann_promoter["promoter_other_lhs"] >0))\
#                                                                                 & (cortical_interactions_ann_promoter["distal_ATAC.seq_rhs"]>0))\
#                                                                                 |(((cortical_interactions_ann_promoter['promoter_rhs'] >0) | (cortical_interactions_ann_promoter['promoter_other_rhs'] >0))\
#                                                                                 & (cortical_interactions_ann_promoter["distal_ATAC.seq_lhs"]>0))]

# to simplify the model, combining lhs/lhs_otherend promoters
# cortical_interactions_ann_lpromotor_rdpeaks_full_split_dedup_merge["promoter_all_lhs_ids"] = \
# (cortical_interactions_ann_lpromotor_rdpeaks_full_split_dedup_merge["promoter_lhs_ids"]
# + "," + cortical_interactions_ann_lpromotor_rdpeaks_full_split_dedup_merge["promoter_other_lhs_ids"])