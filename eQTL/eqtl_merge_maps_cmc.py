# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 13:02:03 2019

@author: libin
"""

import pandas as pd
import sys


def eQTL_analysis(left_hand, right_hand, merged_file):
    
    header = ["eqtl_chr","eqtl_start","eqtl_end","snp_id","gene_id","eqrl_fdr","eqtl_id",
              "inter_chr","inter_start","inter_end","inter_ID","fdr","logp"]

    lh_interactions = pd.read_table(left_hand, header=None, names=header, sep="\t")
    print("left hand read")
    rh_interactions = pd.read_table(right_hand, header=None, names=header, sep="\t")
    print("right hand read")
    
    intersect_both = pd.merge(lh_interactions, rh_interactions, how="inner", 
                              on=["eqtl_chr","gene_id","eqtl_id","eqrl_fdr", "fdr","inter_chr", "logp","inter_ID"])
    intersect_both.to_csv(merged_file, sep="\t", index=False)


eQTL_analysis(left_hand=sys.argv[1], right_hand=sys.argv[2], merged_file=sys.argv[3])