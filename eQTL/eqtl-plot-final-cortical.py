# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 21:31:05 2019

@author: libin
"""

import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


cortical_all_sig_interactions = pd.read_table(r"C:\Users\libin\UCSF\eQTL\pcHiC\cortical\cortical_lh_interactions", sep="\t", names=["chr","start","end","score","interactions_ID"])

cor_sig_intersect_cor_sig_ID = pd.read_table(r"C:\Users\libin\UCSF\eQTL\pcHiC\cortical\cortical_intersect_sig_cortex_merged.ID", names=["interactions_ID"])
cor_sig_intersect_cor_sig = pd.merge(cor_sig_intersect_cor_sig_ID, cortical_all_sig_interactions, on=["interactions_ID"], how="inner")
cor_sig_intersect_cor_rand_ID = pd.read_table(r"C:\Users\libin\UCSF\eQTL\pcHiC\cortical\Brain_Cortex.allpairs.tss_filt.3xgene.3x.sampled_merged.ID", names=["interactions_ID"])
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
print (stats.ks_2samp(cor_sig_intersect_cor_rand["score"],cor_sig_intersect_cor_sig["score"]))

cortical_final_plot_fout = pd.DataFrame()
cortical_final_plot_fout["random"] = cor_sig_intersect_cor_rand_fout["score"]
cortical_final_plot_fout["cortical"] = cor_sig_intersect_cor_sig_fout["score"]
cortical_final_plot_fout = cortical_final_plot_fout[["cortical","random"]]


plt.figure()
sns.violinplot(x="variable",y="value",data=cortical_final_plot.melt(), palette="Pastel2")
# plt.savefig('C://Users//libin//UCSF/eQTL/cortical_final_1X_egene_scontrol_violin.pdf', transparent=True)

plt.figure()
sns.violinplot(x="variable",y="value",data=cortical_final_plot.melt(), palette="Pastel2", cut=0)
# plt.savefig('C://Users//libin//UCSF/eQTL/cortical_final_1X_egene_scontrol_cut_violin.pdf', transparent=True)


'''plt.figure()
sns.violinplot(x="variable",y="value",data=cortical_final_plot_fout.melt(), palette=["#5ebeeb", "#aedef5"], cut=0, zorder=0).set(ylim=(2,20))
# ax2 = ax.twinx()
sns.boxplot(x="variable",y="value",data=cortical_final_plot.melt(), 
            showfliers=False, showbox=True,
            boxprops = {'color': 'black', 'linestyle': '-'},
            capprops = {'color': 'grey', 'linestyle': '-'},
            medianprops = {'color': 'grey', 'linestyle': '--'},
            width=0.1, whiskerprops = dict(linestyle='-', linewidth=3, color="grey"))
plt.savefig('C://Users//libin//UCSF/eQTL/cortical_final_1X_egenes_scontrol_fout_cut_violin.pdf', transparent=True)
'''
# fig1, ax1 = plt.subplots()
# ax1.boxplot(data=cortical_final_plot_fout.melt(),
                             # showmeans=False, showcaps=True, showbox=True, showfliers=False, zorder=10)

# plt.show()

