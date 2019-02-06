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
cor_sig_intersect_cor_rand_ID = pd.read_table("C://Users//libin//UCSF//eQTL/cortical/Brain_Hippocampus.allpairs.tss_filt.3xgene.3x.sampled_merged.ID", names=["interactions_ID"])
cor_sig_intersect_cor_rand = pd.merge(cortical_all_sig_interactions, cor_sig_intersect_cor_rand_ID, on=["interactions_ID"], how="inner")
cor_sig_intersect_cor_nonsig_ID = pd.read_table("C://Users//libin//UCSF//eQTL/cortical/cortical_intersect_nonsig_cortex_merged.ID", names=["interactions_ID"])
cor_sig_intersect_cor_nonsig = pd.merge(cor_sig_intersect_cor_nonsig_ID, cortical_all_sig_interactions, on=["interactions_ID"])
'''
cor_sig_intersect_cor_rand_ID = pd.read_table("C://Users//libin//UCSF//eQTL/cortical/Brain_Cortex.allpairs.tss_filt.3xgene.3x.sampled_merged.ID", names=["interactions_ID"])
cor_sig_intersect_cor_rand = pd.merge(cortical_all_sig_interactions, cor_sig_intersect_cor_rand_ID, on=["interactions_ID"], how="inner")

'''cor_sig_intersect_cor_nonsig_ID = pd.read_table("C://Users//libin//UCSF//eQTL/cortical/cortical_intersect_nonsig_egenes_cortex_merged.ID", names=["interactions_ID"])
cor_sig_intersect_cor_nonsig = pd.merge(cor_sig_intersect_cor_nonsig_ID, cortical_all_sig_interactions, on=["interactions_ID"])
'''

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

cor_IQR = cor_sig_intersect_cor_sig["score"].describe()["75%"] - cor_sig_intersect_cor_sig["score"].describe()["25%"]
cor_upper = cor_sig_intersect_cor_sig["score"].describe()["75%"] + 1.5*cor_IQR
random_IQR = cor_sig_intersect_cor_rand["score"].describe()["75%"] - cor_sig_intersect_cor_rand["score"].describe()["25%"]
random_upper = cor_sig_intersect_cor_rand["score"].describe()["75%"] + 1.5*random_IQR

print(cor_upper, random_upper)

cor_sig_intersect_cor_sig_fout = cor_sig_intersect_cor_sig.loc[cor_sig_intersect_cor_sig["score"]<cor_upper] 
cor_sig_intersect_cor_rand_fout = cor_sig_intersect_cor_rand.loc[cor_sig_intersect_cor_rand["score"]<random_upper]


cortical_final_plot = pd.DataFrame()
cortical_final_plot["random"] = cor_sig_intersect_cor_rand["score"]
# cortical_final_plot["non-signif"] = cor_sig_intersect_cor_nonsig["score"]
cortical_final_plot["cortical"] = cor_sig_intersect_cor_sig["score"]
# cortical_final_plot["liver"] = cor_sig_intersect_sig_liver["score"]
# cortical_final_plot["blood"] = cor_sig_intersect_sig_blood["score"]
cortical_final_plot = cortical_final_plot[["cortical","random"]]
# cortical_final_plot = cortical_final_plot[["cortical","random", "non-signif","liver", "blood"]]


cortical_final_plot_fout = pd.DataFrame()
cortical_final_plot_fout["random"] = cor_sig_intersect_cor_rand_fout["score"]
# cortical_final_plot["non-signif"] = cor_sig_intersect_cor_nonsig["score"]
cortical_final_plot_fout["cortical"] = cor_sig_intersect_cor_sig_fout["score"]
# cortical_final_plot["liver"] = cor_sig_intersect_sig_liver["score"]
# cortical_final_plot["blood"] = cor_sig_intersect_sig_blood["score"]
cortical_final_plot_fout = cortical_final_plot_fout[["cortical","random"]]

'''
plt.figure()
sns.boxplot(x="variable",y="value",data=cortical_final_plot.melt(), palette="Pastel2").set(ylim=(3,20))
plt.savefig('C://Users//libin//UCSF/eQTL/cortical_final_1X_scontrol.pdf', transparent=True)

plt.figure()
sns.boxplot(x="variable",y="value",data=cortical_final_plot.melt(), palette="Pastel2", showfliers=False).set(ylim=(3,20))
plt.savefig('C://Users//libin//UCSF/eQTL/cortical_final_1X_scontrol_Foutlier.pdf', transparent=True)
'''

plt.figure()
sns.violinplot(x="variable",y="value",data=cortical_final_plot.melt(), palette="Pastel2")
plt.savefig('C://Users//libin//UCSF/eQTL/cortical_final_3X_scontrol_violin.pdf', transparent=True)

plt.figure()
sns.violinplot(x="variable",y="value",data=cortical_final_plot.melt(), palette="Pastel2", cut=0)
plt.savefig('C://Users//libin//UCSF/eQTL/cortical_final_3X_scontrol_cut_violin.pdf', transparent=True)

plt.figure()
sns.violinplot(x="variable",y="value",data=cortical_final_plot_fout.melt(), palette=["#5ebeeb", "#aedef5"], cut=0).set(ylim=(2,20))
sns.boxplot(x="variable",y="value",data=cortical_final_plot_fout.melt(), 
            showfliers=False, showbox=False, zorder=1,
            boxprops = {'color': 'black', 'linestyle': '-'},
            capprops = {'color': 'grey', 'linestyle': '-'},
            medianprops = {'color': 'grey', 'linestyle': '--'},
            width=0.35, whiskerprops = dict(linestyle='-', linewidth=3, color="grey"))
plt.savefig('C://Users//libin//UCSF/eQTL/cortical_final_3X_scontrol_fout_cut_violin.pdf', transparent=True)