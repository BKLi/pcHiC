# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 15:33:07 2018

@author: libin
"""
import pandas as pd
import numpy as np

infile = "C:\\Users\\libin\\R_projects\\eQTL\\eqtl_test"

import sys
import pandas as pds
import numpy as np

def reforme_eQTL(infile, inter_file, outfile1, outfile2):
    # infile: raw eQTL file
    # inter_file: significant interaction
    # outfile1: bed format of eQTL file
    # outfile2: bed format of merged intersections
    
    interactions = pd.read_table(inter_file, sep="\t")

    eQTL_raw = pd.read_table(infile, delim_whitespace=True)
    eQTL_raw["chr"] = "chr" + eQTL_raw["variant_id"].str.split("_", expand=True)[0]
    eQTL_raw["end"] = eQTL_raw["variant_id"].str.split("_", expand=True)[1]
    eQTL_raw["start"] = eQTL_raw["end"].apply(int) - eQTL_raw["tss_distance"]
    eQTL_raw["start"],eQTL_raw["end"] = np.where(eQTL_raw["start"] > eQTL_raw["end"].apply(int), 
            [eQTL_raw["end"],eQTL_raw["start"]], [eQTL_raw["start"], eQTL_raw["end"]])
    eQTL_reformed = eQTL_raw[["chr","start","end","gene_id"]]

            
    eQTL_reformed = eQTL_raw[["chr", "start", "end", "gene_id", "pval_nominal"]]
    eQTL_reformed.to_csv(outfile1, sep="\t", index=False)
    
    interaction_lh = interactions[["bait_chr"," bait_start"," bait_end", "score", "interaction_ID"]]
    interaction_rh = interactions[["otherEnd_chr", "otherEnd_start", "otherEnd_end", "score", "interaction_ID"]]
    
    eQTL_bed = pybedtools.BedTool.from_dataframe(eQTL_reformed)
    lh_bed = pybedtools.BedTool.from_dataframe(interaction_lh)
    rh_bed = pybedtools.BedTool.from_dataframe(interaction_rh)
    
    eQTL_and_lh = eQTL_bed.intersect(lh_bed, wa=True, wb=True)
    eQTL_and_rh = eQTL_bed.intersect(rh_bed, wa=True, wb=True)
    eQTL_and_lh_df = pybedtools.BedTool.to_dataframe(eQTL_and_lh, names=["bait_chr"," bait_start"," bait_end", "score", "interaction_ID"])
    eQTL_and_rh_df = pybedtools.BedTool.to_dataframe(eQTL_and_rh, names=["otherEnd_chr", "otherEnd_start", "otherEnd_end", "score", "interaction_ID"])
    eQTL_both = pd.merge(eQTL_and_lh_df, eQTL_and_rh_df, on=["score", "interaction_ID"], how="inner")
    
    eQTL_both.to_csv(outfile2, sep="\t", index=False)
    
    
reforme_eQTL(infile=sys.argv[1], inter_file=sys.argv[2], outfile1=sys.argv[3], outfile2=sys.argv[4])