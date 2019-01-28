# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 09:00:43 2019

@author: libin
"""

import sys
import pandas as pd
import numpy as np
from statistics import mean


def eqtl_sample(sig_pair, all_pair, power, outfile):

    # count number of eQTL pairs in each distance bin
    cell_sig_pairs = pd.read_table(sig_pair, sep="\t")
    # hippocampus_sig_pairs = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\Brain_Hippocampus.v7.signif_variant_gene_pairs.txt", sep="\t")
    cell_sig_pairs["tss_distance"] = cell_sig_pairs["tss_distance"].abs()
    distance_bins = list(np.concatenate(
        (np.arange(0, 100000, 20000), np.arange(100000, 200000, 50000), np.arange(200000, 500000, 100000), np.arange(500000, 1000001, 125000)
         )))
    
    cell_sig_pairs_binned = cell_sig_pairs.groupby(pd.cut(cell_sig_pairs.tss_distance, distance_bins))
    print(cell_sig_pairs_binned.count()["tss_distance"])
    # number of samples in each group of full dataset
    sample_size = [i*int(power) for i in cell_sig_pairs_binned.count()["tss_distance"].tolist()]
    print(sample_size)
    print("total: ", sum(sample_size))
    # hippocampus_all_pairs = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\Brain_Hippocampus.v7.signif_variant_gene_pairs.txt", sep="\t")
    cell_all_pairs = pd.read_table(all_pair, sep="\t")
    cell_all_pairs_binned = cell_all_pairs.groupby(pd.cut(cell_all_pairs.tss_distance, distance_bins))
    print(cell_all_pairs_binned.count()["tss_distance"])
    # return groupby object : each sub group correspond to a df with certain distance range
    
    # turn groupby object to list for convenience
    cell_all_pairs_binned_df_list = [group.reset_index().drop(["index"], axis=1) for _, group in cell_all_pairs_binned]
    
    sampled_list = [cell_all_pairs_binned_df_list[i].sample(n=sample_size[i]) for i in range(len(sample_size))]
    
    sampled_concat = pd.concat(sampled_list).reset_index().drop(["index"], axis=1)
    sampled_binned = sampled_concat.groupby(pd.cut(sampled_concat.tss_distance, distance_bins))
    print(sampled_binned.count()["tss_distance"])
    
    sampled_concat.to_csv(outfile, sep="\t", index=False)


eqtl_sample(sig_pair=sys.argv[1], all_pair=sys.argv[2], power=sys.argv[3], outfile=sys.argv[4])
