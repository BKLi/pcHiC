# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 22:52:02 2019

@author: libin
"""

import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


cortical_all_sig_interactions = pd.read_table(r"C:\Users\libin\UCSF\eQTL\pcHiC\cortical\cortical_lh_interactions", sep="\t", names=["chr","start","end","score","interactions_ID"])

cor_sig_intersect_cor_sig_ID = pd.read_table(r"C:\Users\libin\UCSF\eQTL\pcHiC\cortical\cortical_intersect_sig_cortex_merged.ID", names=["interactions_ID"])
cor_sig_intersect_cor_sig = pd.merge(cor_sig_intersect_cor_sig_ID, cortical_all_sig_interactions, on=["interactions_ID"], how="inner")
cor_sig_intersect_cor_rand_ID = pd.read_table(r"C:\Users\libin\UCSF\eQTL\pcHiC\cortical\cortical_intersect_1x_rand_egenes_cortex_merged.ID", names=["interactions_ID"])
cor_sig_intersect_cor_rand = pd.merge(cortical_all_sig_interactions, cor_sig_intersect_cor_rand_ID, on=["interactions_ID"], how="inner")

cor_IQR = cor_sig_intersect_cor_sig["score"].describe()["75%"] - cor_sig_intersect_cor_sig["score"].describe()["25%"]
cor_upper = cor_sig_intersect_cor_sig["score"].describe()["75%"] + 1.5*cor_IQR
random_IQR = cor_sig_intersect_cor_rand["score"].describe()["75%"] - cor_sig_intersect_cor_rand["score"].describe()["25%"]
random_upper = cor_sig_intersect_cor_rand["score"].describe()["75%"] + 1.5*random_IQR

print(cor_upper, random_upper)

cor_sig_intersect_cor_sig_fout = cor_sig_intersect_cor_sig.loc[cor_sig_intersect_cor_sig["score"]<cor_upper] 
cor_sig_intersect_cor_rand_fout = cor_sig_intersect_cor_rand.loc[cor_sig_intersect_cor_rand["score"]<random_upper]

cortical_final_plot = pd.DataFrame()
cortical_final_plot["random"] = cor_sig_intersect_cor_rand["score"]
cortical_final_plot["cortical"] = cor_sig_intersect_cor_sig["score"]
cortical_final_plot = cortical_final_plot[["cortical","random"]]

cortical_final_plot_fout = pd.DataFrame()
cortical_final_plot_fout["random"] = cor_sig_intersect_cor_rand_fout["score"]
cortical_final_plot_fout["cortical"] = cor_sig_intersect_cor_sig_fout["score"]
cortical_final_plot_fout = cortical_final_plot_fout[["cortical","random"]]

hippocampus_all_sig_interactions = pd.read_table(r"C:\Users\libin\UCSF\eQTL\pcHiC\hippocampus_lh_interactions",sep="\t",names=["chr","start","end","score","interactions_ID"])
hippo_sig_intersect_hippo_sig_ID = pd.read_table(r'C:\Users\libin\UCSF\eQTL\pcHiC\hippocampus_intersect_self_sig_merged.ID', sep="\t", names=["interactions_ID"])
hippo_sig_intersect_hippo_sig = pd.merge(hippocampus_all_sig_interactions, hippo_sig_intersect_hippo_sig_ID, on=["interactions_ID"], how="inner")
hippo_sig_intersect_hippo_rand_ID = pd.read_table(r'C:\Users\libin\UCSF\eQTL\pcHiC\hippo_intersect_1x_egenes_rand_merged.ID', names=["interactions_ID"])
hippo_sig_intersect_hippo_rand = pd.merge(hippocampus_all_sig_interactions, hippo_sig_intersect_hippo_rand_ID, on=["interactions_ID"], how="inner")

hippo_IQR = hippo_sig_intersect_hippo_sig["score"].describe()["75%"] - hippo_sig_intersect_hippo_sig["score"].describe()["25%"]
hippo_upper = hippo_sig_intersect_hippo_sig["score"].describe()["75%"] + 1.5*hippo_IQR
random_IQR = hippo_sig_intersect_hippo_rand["score"].describe()["75%"] - hippo_sig_intersect_hippo_rand["score"].describe()["25%"]
random_upper = hippo_sig_intersect_hippo_rand["score"].describe()["75%"] + 1.5*random_IQR

print(hippo_upper, random_upper)

hippo_sig_intersect_hippo_sig_fout = hippo_sig_intersect_hippo_sig.loc[hippo_sig_intersect_hippo_sig["score"]<hippo_upper] 
hippo_sig_intersect_hippo_rand_fout = hippo_sig_intersect_hippo_rand.loc[hippo_sig_intersect_hippo_rand["score"]<random_upper]


hippo_final_plot = pd.DataFrame()
hippo_final_plot["random"] = hippo_sig_intersect_hippo_rand["score"]
hippo_final_plot["hippocampus"] = hippo_sig_intersect_hippo_sig["score"]
hippo_final_plot = hippo_final_plot[["hippocampus","random"]]

hippo_final_plot_fout = pd.DataFrame()
hippo_final_plot_fout["random"] = hippo_sig_intersect_hippo_rand_fout["score"]
hippo_final_plot_fout["hippocampus"] = hippo_sig_intersect_hippo_sig_fout["score"]
hippo_final_plot_fout = hippo_final_plot_fout[["hippocampus","random"]]



both_final_plot_fout = pd.DataFrame()
both_final_plot_fout["random_hippo"] = hippo_sig_intersect_hippo_rand_fout["score"]
both_final_plot_fout["hippocampus"] = hippo_sig_intersect_hippo_sig_fout["score"]
both_final_plot_fout["random_cor"] = cor_sig_intersect_cor_rand_fout["score"]
both_final_plot_fout["cortical"] = cor_sig_intersect_cor_sig_fout["score"]
both_final_plot_fout = both_final_plot_fout[["cortical", "random_cor", "hippocampus", "random_hippo"]]

both_final_plot = pd.DataFrame()
both_final_plot["random_hippo"] = hippo_sig_intersect_hippo_rand["score"]
both_final_plot["hippocampus"] = hippo_sig_intersect_hippo_sig["score"]
both_final_plot["random_cor"] = cor_sig_intersect_cor_rand["score"]
both_final_plot["cortical"] = cor_sig_intersect_cor_sig["score"]
both_final_plot = both_final_plot[["cortical", "random_cor", "hippocampus", "random_hippo"]]



'''plt.figure(figsize=(10,5))
sns.violinplot(x="variable",y="value",data=both_final_plot_fout.melt(), palette=["#5ebeeb", "#aedef5", "#eb8b5e","#f7d0be"], cut=0).set(ylim=(2,20))
sns.boxplot(x="variable",y="value",data=both_final_plot.melt(), 
            showfliers=False, showbox=True,
            boxprops = {'color': 'black', 'linestyle': '-'},
            capprops = {'color': 'grey', 'linestyle': '-'},
            medianprops = {'color': 'grey', 'linestyle': '--'},
            width=0.1, whiskerprops = dict(linestyle='-', linewidth=3, color="grey"))

plt.savefig('C://Users//libin//UCSF/eQTL/both_final_1X_egenes_scontrol_fout_cut_violin.pdf', transparent=True)
'''