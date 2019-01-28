# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 21:31:05 2019

@author: libin
"""

import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


cortical_all_sig_interactions = pd.read_table("C://Users//libin//UCSF//eQTL/cortical/cortical_lh_interactions", sep="\t", names=["chr","start","end","score","interactions_ID"])

cor_sig_intersect_cor_sig_ID = pd.read_table("C://Users//libin//UCSF//eQTL/cortical/cortical_intersect_sig_cortex_merged.ID", names=["interactions_ID"])
cor_sig_intersect_cor_sig = pd.merge(cor_sig_intersect_cor_sig_ID, cortical_all_sig_interactions, on=["interactions_ID"], how="inner")

'''
cor_sig_intersect_cor_rand_ID = pd.read_table("C://Users//libin//UCSF//eQTL/cortical/cortical_intersect_rand_cortex_merged.ID", names=["interactions_ID"])
cor_sig_intersect_cor_rand = pd.merge(cortical_all_sig_interactions, cor_sig_intersect_cor_rand_ID, on=["interactions_ID"], how="inner")
cor_sig_intersect_cor_nonsig_ID = pd.read_table("C://Users//libin//UCSF//eQTL/cortical/cortical_intersect_nonsig_cortex_merged.ID", names=["interactions_ID"])
cor_sig_intersect_cor_nonsig = pd.merge(cor_sig_intersect_cor_nonsig_ID, cortical_all_sig_interactions, on=["interactions_ID"])
'''
cor_sig_intersect_cor_rand_ID = pd.read_table("C://Users//libin//UCSF//eQTL/cortical/cortical_sig.intersect.cortex_tss_dist_num_filtered.sampled_merged.ID", names=["interactions_ID"])
cor_sig_intersect_cor_rand = pd.merge(cortical_all_sig_interactions, cor_sig_intersect_cor_rand_ID, on=["interactions_ID"], how="inner")

cor_sig_intersect_cor_nonsig_ID = pd.read_table("C://Users//libin//UCSF//eQTL/cortical/cortical_sig_intersect_tss_dist_num_nonsig_sampled_cortex_merged.ID", names=["interactions_ID"])
cor_sig_intersect_cor_nonsig = pd.merge(cor_sig_intersect_cor_nonsig_ID, cortical_all_sig_interactions, on=["interactions_ID"])

'''cor_sig_intersect_cor_nonsig_ID = pd.read_table("C://Users//libin//UCSF//eQTL/cortical/cortical_intersect_nonsig_cortex_merged.ID", names=["interactions_ID"])
cor_sig_intersect_cor_nonsig = pd.merge(cor_sig_intersect_cor_nonsig_ID, cortical_all_sig_interactions, on=["interactions_ID"])
'''

cor_sig_intersect_sig_liver_ID = pd.read_table("C://Users//libin//UCSF//eQTL/cortical/cortical_intersect_sig_liver_merged.ID", names=["interactions_ID"])
cor_sig_intersect_sig_liver = pd.merge(cor_sig_intersect_sig_liver_ID, cortical_all_sig_interactions, on=["interactions_ID"])
cor_sig_intersect_sig_blood_ID = pd.read_table("C://Users//libin//UCSF//eQTL/cortical/cortical_sig_intersect_sig_blood_merged.ID", names=["interactions_ID"])
cor_sig_intersect_sig_blood = pd.merge(cor_sig_intersect_sig_blood_ID, cortical_all_sig_interactions, on=["interactions_ID"])


'''
cortical_all_interactions = pd.read_table("C://Users//libin//UCSF//eQTL/cortical/cortical_full_lh_interactions", sep="\t", names=["chr","start","end","score","interactions_ID"])
cor_full_intersect_sig_cortex_ID = pd.read_table("C://Users//libin//UCSF//eQTL/cortical/cortical_full_intersect_sig_cortex_merged.ID", names=["interactions_ID"])
cor_full_intersect_sig_cortex = pd.merge(cor_full_intersect_sig_cortex_ID, cortical_all_interactions, on=["interactions_ID"])
cor_full_intersect_sig_liver_ID = pd.read_table("C://Users//libin//UCSF//eQTL/cortical/cortical_full_intersect_sig_liver_merged.ID", names=["interactions_ID"])
cor_full_intersect_sig_liver = pd.merge(cor_full_intersect_sig_liver_ID, cortical_all_interactions, on=["interactions_ID"])
'''

cortical_final_plot = pd.DataFrame()
cortical_final_plot["random"] = cor_sig_intersect_cor_rand["score"]
cortical_final_plot["non-signif"] = cor_sig_intersect_cor_nonsig["score"]
cortical_final_plot["cortical"] = cor_sig_intersect_cor_sig["score"]
cortical_final_plot["liver"] = cor_sig_intersect_sig_liver["score"]
cortical_final_plot["blood"] = cor_sig_intersect_sig_blood["score"]
# cortical_final_plot = cortical_final_plot[["cortical","random","liver", "blood"]]
cortical_final_plot = cortical_final_plot[["cortical","random", "non-signif","liver", "blood"]]


plt.figure()
sns.boxplot(x="variable",y="value",data=cortical_final_plot.melt(), palette="Pastel2").set(ylim=(3,20))
plt.savefig('C://Users//libin//UCSF/eQTL/cortical_final-new.pdf', transparent=True)

plt.figure()
sns.boxplot(x="variable",y="value",data=cortical_final_plot.melt(), palette="Pastel2", showfliers=False).set(ylim=(3,20))
plt.savefig('C://Users//libin//UCSF/eQTL/cortical_final_Foutlier_new.pdf', transparent=True)

