# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 20:45:18 2019

@author: libin
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as ss
import numpy as np
from statistics import mean 

'''liver_sig_pairs = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\Liver.sig.bed", sep="\t")
lung_sig_pairs = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\Lung.sig.bed", sep="\t")
blood_sig_pairs = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\Whole_Blood.sig.bed", sep="\t")
lympho_sig_pairs = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\Cells_EBV-transformed_lymphocytes.sig.bed", sep="\t")
cor_sig_pairs = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\Brain_Cortex.sig.bed", sep="\t")
hippo_sig_pairs = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\Brain_Hippocampus.sig.bed", sep="\t")
shared_hippo_liver = pd.merge(liver_sig_pairs, hippo_sig_pairs, on=["chr","start","end","gene_id"])
shared_cor_liver = pd.merge(liver_sig_pairs, cor_sig_pairs, on=["chr","start","end","gene_id"])
shared_hippo_lympho = pd.merge(lympho_sig_pairs, hippo_sig_pairs, on=["chr","start","end","gene_id"])
shared_hippo_lung = pd.merge(lung_sig_pairs, hippo_sig_pairs, on=["chr","start","end","gene_id"])
shared_hippo_blood = pd.merge(blood_sig_pairs, hippo_sig_pairs, on=["chr","start","end","gene_id"])
shared_cor_blood = pd.merge(blood_sig_pairs, cor_sig_pairs, on=["chr","start","end","gene_id"])'''

Hippocampus_sig_QTL = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\Brain_Hippocampus.v7.signif_variant_gene_pairs.txt", sep="\t")
'''names=["variant_id", "gene_id", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se", "pval_nominal_threshold", "min_pval_nominal", "pval_beta"]'''
Hippocampus_rand_QTL = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\Brain_Hippocampus.rand.full.txt", sep="\t", names=["variant_id", "gene_id", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se", "pval_nominal_threshold", "min_pval_nominal", "pval_beta"])