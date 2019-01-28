# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 21:59:38 2018

@author: libin
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as ss
from statistics import mean 

'''merged_df = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\cortical_chr21_intersect_merged", sep="\t")
merged_df_above5 = merged_df[merged_df["score"] >= 5]
merged_df_below5 = merged_df[merged_df["score"] < 5]'''

'''cortical_self_overlapped_eQTL = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\chr21_pairs.uniq.bed", sep="\t")
cortical_all_eQTL = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\Cortex_Chr21.bed", sep="\t", names=["SNP_chr","SNP_start","SNP_end","gene_ID","p-value"])
cortical_self_merged_eQTL = pd.merge(cortical_self_overlapped_eQTL, cortical_all_eQTL, on=["SNP_chr","SNP_start","SNP_end"], how="inner")
'''

hippo_all_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippocampus_lh_interactions",sep="\t",names=["chr","start","end","score","interactions_ID"])
cortical_all_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\cortical_lh_interactions", sep="\t", names=["chr","start","end","score","interactions_ID"])
motor_all_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\motor_lh_interactions", sep="\t", names=["chr","start","end","score","interactions_ID"])

'''
cortical_self_overlapped_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\chr21_interactions.uniq.bed", sep="\t", names=["start_chr","interactions_ID"])
cortical_self_merged_interactions = pd.merge(cortical_self_overlapped_interactions, cortical_all_interactions, on=["interactions_ID"],how="inner")
'''

cortical_self_overlapped_sig_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\cor_self_sig_overlapped_interactions.ID", names=["interactions_ID"])
cortical_self_merged_sig_interactions = pd.merge(cortical_self_overlapped_sig_interactions, cortical_all_interactions, on=["interactions_ID"], how="inner")

cortical_blood_overlapped_sig_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\cor_blood_sig_overlapped_interactions.ID", names=["interactions_ID"])
cortical_blood_merged_sig_interactions = pd.merge(cortical_blood_overlapped_sig_interactions, cortical_all_interactions, on=["interactions_ID"], how="inner")

cortical_liver_overlapped_sig_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\cor_liver_sig_overlapped_interactions.ID", names=["interactions_ID"])
cortical_liver_merged_sig_interactions = pd.merge(cortical_liver_overlapped_sig_interactions, cortical_all_interactions, on=["interactions_ID"], how="inner")

cortical_rand_overlapped_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\cor_intersect_rand.ID", names=["interactions_ID"])
cortical_rand_merged_interactions = pd.merge(cortical_rand_overlapped_interactions, cortical_all_interactions, on=["interactions_ID"], how="inner")

hippo_self_overlapped_sig_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippo_self_sig_overlapped_interactions.ID", names=["interactions_ID"])
hippo_self_merged_sig_interactions = pd.merge(hippo_self_overlapped_sig_interactions, hippo_all_interactions, on=["interactions_ID"], how="inner")

hippo_blood_overlapped_sig_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippo_blood_sig_overlapped_interactions.ID", names=["interactions_ID"])
hippo_blood_merged_sig_interactions = pd.merge(hippo_blood_overlapped_sig_interactions, hippo_all_interactions, on=["interactions_ID"], how="inner")

hippo_liver_overlapped_sig_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippo_liver_sig_overlapped_interactions.ID", names=["interactions_ID"])
hippo_liver_merged_sig_interactions = pd.merge(hippo_liver_overlapped_sig_interactions, hippo_all_interactions, on=["interactions_ID"], how="inner")

hippo_rand_overlapped_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippo_rand_overlapped_interactions.ID", names=["interactions_ID"])
hippo_rand_merged_interactions = pd.merge(hippo_rand_overlapped_interactions, hippo_all_interactions, on=["interactions_ID"], how="inner")

hippo_lymph_overlapped_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippo_lymph_sig_overlapped_interactions.ID", names=["interactions_ID"])
hippo_lymph_merged_interactions = pd.merge(hippo_lymph_overlapped_interactions, hippo_all_interactions, on=["interactions_ID"], how="inner")

motor_self_overlapped_sig_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\mot_self_sig_overlapped_interactions.ID", names=["interactions_ID"])
motor_self_merged_sig_interactions = pd.merge(motor_self_overlapped_sig_interactions, motor_all_interactions, on=["interactions_ID"], how="inner")

'''
cortical_blood_overlapped_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\chr21_cortical_blood_interactions.uniq.bed", sep="\t", names=["interactions_ID"])
cortical_blood_merged_interactions = pd.merge(cortical_blood_overlapped_interactions, cortical_all_interactions, on=["interactions_ID"], how="inner")
'''

'''plt.figure(figsize=(10,6))
sns.violinplot(cortical_self_merged_interactions["score"]).set(xlim=(4,20))
plt.figure(figsize=(10,6))
sns.violinplot(cortical_blood_merged_interactions["score"]).set(xlim=(4,20))'''
'''plt.figure(figsize=(6,9))
sns.boxplot(y=cortical_self_merged_interactions["score"]).set(ylim=(4,20))'''



plt.figure(figsize=(6,9))
sns.violinplot(y=cortical_self_merged_sig_interactions["score"]).set(ylim=(4,20))
cor_self = cortical_self_merged_sig_interactions["score"].tolist()

plt.figure(figsize=(6,9))
sns.violinplot(y=cortical_blood_merged_sig_interactions["score"]).set(ylim=(4,20))
cor_blood = cortical_blood_merged_sig_interactions["score"].tolist()

plt.figure(figsize=(6,9))
sns.violinplot(y=cortical_liver_merged_sig_interactions["score"]).set(ylim=(4,20))
cor_liver = cortical_liver_merged_sig_interactions["score"].tolist()

plt.figure(figsize=(6,9))
sns.violinplot(y=hippo_self_merged_sig_interactions["score"]).set(ylim=(4,20))
hippo_self = hippo_self_merged_sig_interactions["score"].tolist()

plt.figure(figsize=(6,9))
sns.violinplot(y=hippo_blood_merged_sig_interactions["score"]).set(ylim=(4,20))
hippo_blood = hippo_blood_merged_sig_interactions["score"].tolist()

plt.figure(figsize=(6,9))
sns.violinplot(y=hippo_liver_merged_sig_interactions["score"]).set(ylim=(4,20))
hippo_liver = hippo_liver_merged_sig_interactions["score"].tolist()

'''
print ("cor_self, cor_liver", ss.anderson_ksamp([cor_self, cor_liver]), ss.ks_2samp(cor_self,cor_liver))
print ("hippo_self, hippo_liver", ss.anderson_ksamp([hippo_self, hippo_liver]), ss.ks_2samp(hippo_self, hippo_liver))
'''

'''cortical_all_subset_liver = cortical_all_interactions.sample(n=29969,random_state=10,axis=0)
cortical_all_subset_self = cortical_all_interactions.sample(n=41707,random_state=10,axis=0)
cortical_forPlot = pd.DataFrame()
cortical_forPlot["cor_self"] = cortical_self_merged_sig_interactions["score"]
cortical_forPlot["cor_all"] = cortical_all_subset_self["score"]

sns.violinplot(x="variable",y="value",data=cortical_forPlot.melt(), palette="Pastel1", inner="box").set(ylim=(2,25))

cor_self = cortical_self_merged_sig_interactions["score"].tolist()
cor_all_self = cortical_all_subset_self["score"].tolist()
ss.ks_2samp(cor_self,cor_all_self)'''


hippo_self_single_end_sig_ID = pd.read_table( "C:\\Users\\libin\\UCSF\\eQTL\\hippocampus_intersect_sig_single_end.uniq.ID", sep="\t", names=["interactions_ID"])
hippo_all_sig_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippocampus_lh_interactions",sep="\t",names=["chr","start","end","score","interactions_ID"])
hippo_self_single_end_sig = pd.merge(hippo_self_single_end_sig_ID, hippo_all_sig_interactions, on=["interactions_ID"], how="inner")
sns.violinplot(hippo_self_single_end_sig).set(xlim=(3,20))
sns.violinplot(hippo_self_single_end_sig["score"]).set(xlim=(3,20))
print(hippo_self_single_end_sig["score"].mean())
print(hippo_self_single_end_sig["score"].median())
hippo_self_single_end_sig_score = pd.read_table( "C:\\Users\\libin\\UCSF\\eQTL\\hippocampus_intersect_sig_single_end.score", sep="\t", names=["score"])
sns.violinplot(hippo_self_single_end_sig_score["score"]).set(xlim=(3,20))
print(hippo_self_single_end_sig_score["score"].median())
print(hippo_self_single_end_sig_score["score"].mean())
hippocampus_sig_pairs = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\Brain_Hippocampus.v7.signif_variant_gene_pairs.txt", sep="\t")
sns.distplot(hippocampus_sig_pairs["score"])
sns.distplot(hippocampus_sig_pairs["tss_distance"])
print(hippocampus_sig_pairs["tss_distance"].max())
print(hippocampus_sig_pairs["tss_distance"].min())
pd.head(hippocampus_sig_pairs)
hippocampus_sig_pairs.head(n=2)


sns.boxplot(x="variable",y="value",data=hippo_plot1.melt(), palette="Set2")