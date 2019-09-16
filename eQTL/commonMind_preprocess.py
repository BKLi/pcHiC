# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 22:40:23 2019

@author: libin
"""

import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None
pd.set_option('display.max_columns', 40)

cmc_eqtl = pd.read_table(r'C:\Users\libin\UCSF\eQTL\hfb\CMC_SVA\CMC_MSSM-Penn-Pitt_DLPFC_mRNA_eQTL-adjustedSVA-binned.txt',
                         delim_whitespace=True, 
                         names= ['SNP',
                                 'Gene',
                                 'Gene_symbol',
                                 "<",
                                 'FDR',
                                 'SNP_location',
                                 'Expression_Increasing_Allele',
                                 'Expression_Decreasing_Allele',
                                 'Gene_location',
                                 'Gene_strand',
                                 'eQTL_type',
                                 'mRNA_normalization' ])[1:]
cmc_eqtl = cmc_eqtl[['SNP',
                     'Gene',
                     'FDR',
                     'SNP_location',
                     'Gene_location',
                     'Gene_strand',
                     'eQTL_type']]
cmc_eqtl = cmc_eqtl[cmc_eqtl["eQTL_type"] == "cis"]
cmc_eqtl = cmc_eqtl.drop(columns=["eQTL_type"])
cmc_eqtl_pos = cmc_eqtl[cmc_eqtl["Gene_strand"] == "+"]
cmc_eqtl_neg = cmc_eqtl[cmc_eqtl["Gene_strand"] == "-"]

cmc_eqtl_pos["TSS"] = cmc_eqtl_pos.loc[:, "Gene_location"].str.extract(r'chr\d{1,2}:(\d+)..')
cmc_eqtl_neg["TSS"] = cmc_eqtl_neg.loc[:, "Gene_location"].str.extract(r'chr\d{1,2}:\d+..(\d+)')

cmc_eqtl_pos["chr"] = cmc_eqtl_pos["SNP_location"].str.extract(r'(chr\d{1,2}):\d+')
cmc_eqtl_neg["chr"] = cmc_eqtl_neg["SNP_location"].str.extract(r'(chr\d{1,2}):\d+')

cmc_eqtl_pos["SNP_pos"] = cmc_eqtl_pos["SNP_location"].str.extract(r'chr\d{1,2}:(\d+)')
cmc_eqtl_neg["SNP_pos"] = cmc_eqtl_neg["SNP_location"].str.extract(r'chr\d{1,2}:(\d+)')

cmc_eqtl_pos = \
cmc_eqtl_pos[
['chr',
 'TSS',
 'SNP_pos',
 'SNP',
 'Gene',
 'FDR',
 'Gene_strand']]\
.rename(columns={"TSS":"start", "SNP_pos":"end"})

cmc_eqtl_neg = \
cmc_eqtl_neg[
['chr',
 'SNP_pos',
 'TSS',
 'SNP',
 'Gene',
 'FDR',
 'Gene_strand']]\
.rename(columns={"SNP_pos":"start", "TSS":"end"})

cmc_eqtl_all = pd.concat([cmc_eqtl_pos, cmc_eqtl_neg], ignore_index=True)
cmc_eqtl_all = cmc_eqtl_all.dropna()
cmc_eqtl_all["start"], cmc_eqtl_all["end"] = \
np.where(cmc_eqtl_all["start"].apply(int) > cmc_eqtl_all["end"].apply(int),
         [cmc_eqtl_all["end"],cmc_eqtl_all["start"]], [cmc_eqtl_all["start"], cmc_eqtl_all["end"]])

cmc_eqtl_all["distance"] = cmc_eqtl_all["end"].apply(int) - cmc_eqtl_all["start"].apply(int)
cmc_eqtl_all = cmc_eqtl_all[cmc_eqtl_all["distance"] > 10000]
cmc_eqtl_sig = cmc_eqtl_all[(cmc_eqtl_all["FDR"] == "0.1") | (cmc_eqtl_all["FDR"] == "0.05")]
eQTL_ID = [i+1 for i in range(cmc_eqtl_sig.shape[0])]  
cmc_eqtl_sig["ID"] = eQTL_ID
cmc_eqtl_sig.to_csv(r"C:\Users\libin\UCSF\eQTL\hfb\CMC_SVA\cmc_eqtl_sig", sep="\t", index=False)
