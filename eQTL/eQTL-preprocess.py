# -*- coding: utf-8 -*-

import sys
import pandas as pd
import numpy as np
import pybedtools


def reform_eQTL(infile, out_file):
    # infile: raw eQTL file
    # inter_file: significant interaction
    # outfile1: bed format of eQTL file
    # outfile2: bed format of merged intersections
    
    eQTL_raw = pd.read_table(infile, delim_whitespace=True)
    
    # keep only eQTL-TSS pairs with a tss_distance of at least 10kb
    eQTL_raw = eQTL_raw[eQTL_raw["tss_distance"].abs() > 10000]
    eQTL_raw["chr"] = "chr" + eQTL_raw["variant_id"].str.split("_", expand=True)[0]
    eQTL_raw["end"] = eQTL_raw["variant_id"].str.split("_", expand=True)[1]
    eQTL_raw["start"] = eQTL_raw["end"].apply(int) - eQTL_raw["tss_distance"]
    # adjust order so "start" is always smaller than "end" -- prevent bedtools error 
    eQTL_raw["start"],eQTL_raw["end"] = np.where(eQTL_raw["start"] > eQTL_raw["end"].apply(int), 
        [eQTL_raw["end"],eQTL_raw["start"]], [eQTL_raw["start"], eQTL_raw["end"]])
    eQTL_reformed = eQTL_raw[["chr","start","end","gene_id","pval_nominal"]]
    # add unique ID for each pair
    eQTL_ID = [i+1 for i in range(eQTL_reformed.shape[0])]  
    eQTL_reformed["ID"] = eQTL_ID
    eQTL_reformed.to_csv(out_file, sep="\t", index=False)
    
reform_eQTL(infile=sys.argv[1], out_file=sys.argv[2])