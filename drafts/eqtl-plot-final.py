# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 13:32:42 2019

@author: libin
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


motor_all_interactions = pd.read_table("C://Users//libin//UCSF//eQTL/motor/motor_lh_interactions", sep="\t", names=["chr","start","end","score","interactions_ID"])
motor_intersect_sig_spinal_ID = pd.read_table("C://Users//libin//UCSF//eQTL/motor/motor_intersect_sig_spinal_merged.ID", sep="\t", names=["interactions_ID"])
motor_intersect_sig_Liver_ID = pd.read_table("C://Users//libin//UCSF//eQTL/motor/motor_intersect_sig_Liver_merged.ID", sep="\t", names=["interactions_ID"])
motor_intersect_rand_spinal_ID = pd.read_table("C://Users//libin//UCSF//eQTL/motor/motor_tss_filt_rand_spinal_merged.ID", sep="\t", names=["interactions_ID"])
# motor_intersect_nonsig_spinal_ID = pd.read_table("C://Users//libin//UCSF//eQTL/motor/motor_intersect_nonsig_spinal_merged.ID", sep="\t", names=["interactions_ID"])
motor_intersect_sig_blood_ID = pd.read_table("C://Users//libin//UCSF//eQTL/motor/motor_sig_intersect_sig_blood_merged.ID", sep="\t", names=["interactions_ID"])

motor_intersect_sig_spinal = pd.merge(motor_all_interactions,motor_intersect_sig_spinal_ID, how="inner", on=["interactions_ID"])
motor_intersect_sig_Liver = pd.merge(motor_all_interactions,motor_intersect_sig_Liver_ID, how="inner", on=["interactions_ID"])
motor_intersect_rand_spinal = pd.merge(motor_all_interactions,motor_intersect_rand_spinal_ID, how="inner", on=["interactions_ID"])
# motor_intersect_nonsig_spinal = pd.merge(motor_all_interactions,motor_intersect_nonsig_spinal_ID, how="inner", on=["interactions_ID"])
motor_intersect_sig_blood = pd.merge(motor_all_interactions, motor_intersect_sig_blood_ID, how="inner", on=["interactions_ID"])


motor_final_plot = pd.DataFrame()
motor_final_plot["random"] = motor_intersect_rand_spinal["score"]
# motor_final_plot["non-signif"] = motor_intersect_nonsig_spinal["score"]
motor_final_plot["motor"] = motor_intersect_sig_spinal["score"]
motor_final_plot["liver"] = motor_intersect_sig_Liver["score"]
motor_final_plot["blood"] = motor_intersect_sig_blood["score"]
motor_final_plot = motor_final_plot[["motor","random","liver", "blood"]]

plt.figure()
sns.boxplot(x="variable",y="value",data=motor_final_plot.melt(), palette="Pastel1").set(ylim=(3,20))
plt.savefig('C://Users//libin//UCSF/eQTL/motor_final_1x.pdf', transparent=True)

plt.figure()
sns.boxplot(x="variable",y="value",data=motor_final_plot.melt(), palette="Pastel1", showfliers=False).set(ylim=(3,20))
plt.savefig('C://Users//libin//UCSF/eQTL/motor_final_1x_Foutlier.pdf', transparent=True)
