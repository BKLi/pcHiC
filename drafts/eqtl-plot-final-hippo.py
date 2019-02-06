# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 20:00:08 2019

@author: libin
"""

import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sns.set_style("white")


hippocampus_all_sig_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippocampus_lh_interactions",sep="\t",names=["chr","start","end","score","interactions_ID"])
hippo_sig_intersect_hippo_sig_ID = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippocampus_intersect_self_sig_merged.ID", sep="\t", names=["interactions_ID"])
hippo_sig_intersect_hippo_sig = pd.merge(hippocampus_all_sig_interactions, hippo_sig_intersect_hippo_sig_ID, on=["interactions_ID"], how="inner")
hippo_sig_intersect_hippo_rand_ID = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\Brain_Hippocampus.allpairs.tss_filt.3xgene.3x.sampled_merged.ID", names=["interactions_ID"])
hippo_sig_intersect_hippo_rand = pd.merge(hippocampus_all_sig_interactions, hippo_sig_intersect_hippo_rand_ID, on=["interactions_ID"], how="inner")
# hippo_sig_intersect_hippo_nonsig_ID = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippo_intersect_egenes_nonsig_merged.ID",names=["interactions_ID"])
# hippo_sig_intersect_hippo_nonsig = pd.merge(hippocampus_all_sig_interactions, hippo_sig_intersect_hippo_nonsig_ID, on=["interactions_ID"], how="inner")
hippo_sig_intersect_liver_sig_ID = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippo_intersect_Liver.signif_merged.ID", names=["interactions_ID"])
hippo_sig_intersect_liver_sig = pd.merge(hippo_sig_intersect_liver_sig_ID, hippocampus_all_sig_interactions, on=["interactions_ID"], how="inner")
hippo_sig_intersect_blood_sig_ID = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippo_sig_intersect_sig_blood_merged.ID", names=["interactions_ID"])
hippo_sig_intersect_blood_sig = pd.merge(hippo_sig_intersect_blood_sig_ID, hippocampus_all_sig_interactions, on=["interactions_ID"], how="inner")


'''hippocampus_all_interactions = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippocampus_full_lh_interactions",sep="\t",names=["chr","start","end","score","interactions_ID"])
hippo_full_intersect_hippo_sig_ID = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippocampus_full_intersect_hippo_sig_merged.ID", sep="\t", names=["interactions_ID"])
hippo_full_intersect_hippo_sig = pd.merge(hippocampus_all_interactions, hippo_full_intersect_hippo_sig_ID, on=["interactions_ID"], how="inner")
hippo_full_intersect_hippo_rand_ID = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippocampus_full_intersect_tss_filtered.sampled_merged.ID", names=["interactions_ID"])
hippo_full_intersect_hippo_rand = pd.merge(hippocampus_all_interactions, hippo_full_intersect_hippo_rand_ID, on=["interactions_ID"], how="inner")
hippo_full_intersect_hippo_nonsig_ID = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\Hippocampus_full_intersect_tss_filtered_nonsig_sampled_merged.ID",names=["interactions_ID"])
hippo_full_intersect_hippo_nonsig = pd.merge(hippocampus_all_interactions, hippo_full_intersect_hippo_nonsig_ID, on=["interactions_ID"], how="inner")
hippo_full_intersect_liver_sig_ID = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippo_full_intersect_Liver.signif_merged.ID", names=["interactions_ID"])
hippo_full_intersect_liver_sig = pd.merge(hippo_full_intersect_liver_sig_ID, hippocampus_all_interactions, on=["interactions_ID"], how="inner")
hippo_full_intersect_blood_sig_ID = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippo_full_intersect_sig_blood_merged.ID", names=["interactions_ID"])
hippo_full_intersect_blood_sig = pd.merge(hippo_full_intersect_blood_sig_ID, hippocampus_all_interactions, on=["interactions_ID"], how="inner")
'''


hippo_IQR = hippo_sig_intersect_hippo_sig["score"].describe()["75%"] - hippo_sig_intersect_hippo_sig["score"].describe()["25%"]
hippo_upper = hippo_sig_intersect_hippo_sig["score"].describe()["75%"] + 1.5*hippo_IQR
random_IQR = hippo_sig_intersect_hippo_rand["score"].describe()["75%"] - hippo_sig_intersect_hippo_rand["score"].describe()["25%"]
random_upper = hippo_sig_intersect_hippo_rand["score"].describe()["75%"] + 1.5*random_IQR

print(hippo_upper, random_upper)

hippo_sig_intersect_hippo_sig_fout = hippo_sig_intersect_hippo_sig.loc[hippo_sig_intersect_hippo_sig["score"]<hippo_upper] 
hippo_sig_intersect_hippo_rand_fout = hippo_sig_intersect_hippo_rand.loc[hippo_sig_intersect_hippo_rand["score"]<random_upper]


hippo_final_plot = pd.DataFrame()
hippo_final_plot["random"] = hippo_sig_intersect_hippo_rand["score"]
# hippo_final_plot["non-signif"] = hippo_sig_intersect_hippo_nonsig["score"]
hippo_final_plot["hippocampus"] = hippo_sig_intersect_hippo_sig["score"]
# hippo_final_plot["liver"] = hippo_sig_intersect_liver_sig["score"]
# hippo_final_plot["blood"] = hippo_sig_intersect_blood_sig["score"]
hippo_final_plot = hippo_final_plot[["hippocampus","random"]]

hippo_final_plot_fout = pd.DataFrame()
hippo_final_plot_fout["random"] = hippo_sig_intersect_hippo_rand_fout["score"]
# hippo_final_plot["non-signif"] = hippo_sig_intersect_hippo_nonsig["score"]
hippo_final_plot_fout["hippocampus"] = hippo_sig_intersect_hippo_sig_fout["score"]
# hippo_final_plot["liver"] = hippo_sig_intersect_liver_sig["score"]
# hippo_final_plot["blood"] = hippo_sig_intersect_blood_sig["score"]
hippo_final_plot_fout = hippo_final_plot_fout[["hippocampus","random"]]


'''plt.figure()
sns.boxplot(x="variable",y="value",data=hippo_final_plot.melt(), palette="Pastel2").set(ylim=(3,20))
plt.savefig('C://Users//libin//UCSF/eQTL/hippo_final.pdf', transparent=True)

plt.figure()
sns.boxplot(x="variable",y="value",data=hippo_final_plot.melt(), palette="Pastel2", showfliers=False).set(ylim=(3,20))
plt.savefig('C://Users//libin//UCSF/eQTL/hippo_final_Foutlier.pdf', transparent=True)
'''

# pval_bin = [1e-40, 1e-30, 1e-20, 1e-10, 1e-5, 0.1, 1]
# print(hippocampus_significant_eQTL.groupby(pd.cut(hippocampus_significant_eQTL.pval_nominal, pval_bin)).count()["pval_nominal"])

'''hippo_final_plot_2 = pd.DataFrame()
hippo_final_plot_2["random"] = hippo_full_intersect_hippo_rand["score"]
hippo_final_plot_2["non-signif"] = hippo_full_intersect_hippo_nonsig["score"]
hippo_final_plot_2["hippocampus"] = hippo_full_intersect_hippo_sig["score"]
hippo_final_plot_2["liver"] = hippo_full_intersect_liver_sig["score"]
hippo_final_plot_2["blood"] = hippo_full_intersect_blood_sig["score"]
hippo_final_plot_2 = hippo_final_plot_2[["hippocampus","random", "non-signif","liver", "blood"]]'''

'''
plt.figure()
sns.boxplot(x="variable",y="value",data=hippo_final_plot.melt(), palette="Pastel1").set(ylim=(3,20))
plt.savefig('C://Users//libin//UCSF/eQTL/hippo_final_1x_scontrol.pdf', transparent=True)

plt.figure()
sns.boxplot(x="variable",y="value",data=hippo_final_plot.melt(), palette="Pastel1", showfliers=False).set(ylim=(3,20))
plt.savefig('C://Users//libin//UCSF/eQTL/hippo_final_1x_scontrol_Foutlier.pdf', transparent=True)
'''

'''plt.figure()
sns.violinplot(x="variable",y="value",data=hippo_final_plot.melt(), palette="Pastel1")
plt.savefig('C://Users//libin//UCSF/eQTL/hippo_final_3x_scontrol_violin.pdf', transparent=True)

plt.figure()
sns.violinplot(x="variable",y="value",data=hippo_final_plot.melt(), palette="Pastel1", cut=0)
plt.savefig('C://Users//libin//UCSF/eQTL/hippo_final_3x_cut_scontrol_violin.pdf', transparent=True)
'''


'''plt.figure()
# fig, axes = plt.subplots(1,2)
sns.violinplot(x="variable",y="value",data=hippo_final_plot.melt(), palette="Pastel2").set(ylim=(2,20))
sns.boxplot(x="variable",y="value",data=hippo_final_plot.melt(), 
            showfliers=False, showbox=False, width=0.35,
            boxprops = {'color': 'black', 'linestyle': '-'},
            capprops = {'color': 'grey', 'linestyle': '-'},
            whiskerprops = dict(linestyle='-', linewidth=3, color="grey"))

plt.figure()
sns.set_style("whitegrid")
# fig, axes = plt.subplots(1,2)
# sns.violinplot(x="variable",y="value",data=hippo_final_plot_fout.melt(), palette="Pastel2", cut=0).set(ylim=(2,20))
sns.boxplot(x="variable",y="value",data=hippo_final_plot.melt(), palette="Pastel2",
            showfliers=False, showbox=True, zorder=1,
            capprops = {'color': 'grey', 'linestyle': '-'},
            medianprops = {'color': 'grey', 'linestyle': '--'},
            whiskerprops = dict(linestyle='-', linewidth=3, color="grey")).set(ylim=(4,20))'''


plt.figure()
ax = sns.boxplot(x="variable",y="value",data=hippo_final_plot.melt(), 
            showfliers=False, showbox=True,
            boxprops = {'color': 'black', 'linestyle': '-'},
            capprops = {'color': 'grey', 'linestyle': '-'},
            medianprops = {'color': 'grey', 'linestyle': '--'},
            width=0.1, whiskerprops = dict(linestyle='-', linewidth=3, color="grey"))
plt.setp(ax.lines, zorder=100)
plt.setp(ax.collections, zorder=100, label="")
sns.violinplot(x="variable",y="value",data=hippo_final_plot_fout.melt(), palette=["#eb8b5e","#f7d0be"], cut=0).set(ylim=(2,20))


plt.savefig('C://Users//libin//UCSF/eQTL/hippo_final_3x_cut_scontrol_fout_cut_violin.pdf', transparent=True)


# plt.figure()
# plt.savefig('C://Users//libin//UCSF/eQTL/hippo_final_3x_scontrol_violin.pdf', transparent=True)
# sns.set_style("whitegrid")
# tips = sns.load_dataset("tips")
# ax = sns.boxplot(x="day", y="total_bill", data=tips)



'''
sns.boxplot(x="variable",y="value",data=hippo_plot3.melt(), palette=["#98eff9","#ffffc2"]).set(ylim=(2,20))
sns.boxplot(x="variable",y="value",data=hippo_plot1.melt(), palette=["#d6b4fc","#ffb7ce"]).set(ylim=(2,20))
sns.boxplot(x="variable",y="value",data=hippo_plot1.melt(), palette=["#d6b4fc","#d1ffbd"]).set(ylim=(2,20))
sns.boxplot(x="variable",y="value",data=hippo_plot1.melt(), palette=["#feb209","#d1ffbd"]).set(ylim=(2,20))
sns.boxplot(x="variable",y="value",data=hippo_plot1.melt(), palette=["#ffb19a","#d1ffbd"]).set(ylim=(2,20))
sns.violinplot(x="variable",y="value",data=hippo_plot1.melt(), palette=["#ffb19a","#d1ffbd"]).set(ylim=(2,20))
sns.violinplot(x="variable",y="value",data=hippo_plot1.melt(), palette=["#a7dbf4","#b0a7d7"]).set(ylim=(2,20))
sns.boxplot(x="variable",y="value",data=hippo_plot1.melt(), palette=["#a7dbf4","#b0a7d7"]).set(ylim=(2,20))
'''


