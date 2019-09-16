# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 00:31:34 2019

@author: libin
"""

import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import sys
from scipy import stats

# check if cis-interactions

cell_type = "interneurons"

sig_interactions = pd.read_table(
        r"C:\Users\libin\UCSF\eQTL\hfb\{}\{}.bedpe".format(cell_type,cell_type), sep="\t")
sig_interactions_trans = sig_interactions[sig_interactions["chr1"] != sig_interactions["chr2"]]
if sig_interactions_trans.shape[0] > 0:
    print("check failed, please re-process")
    sys.exit(0)
sig_interactions = pd.read_table(
        r"C:\Users\libin\UCSF\eQTL\hfb\{}\{}_rh_interactions".format(cell_type,cell_type),
        names = ["chr2", "start2", "end2", "ID", "fdr", "ClusterNegLog10P", "count", "expected"], sep="\t")
sig_interactions["normalized"] = sig_interactions["count"] / sig_interactions["expected"]
# significant pairs
cell_intersect_sig_cortex_ID = pd.read_table(
        r"C:\Users\libin\UCSF\eQTL\hfb\{}\{}_intersect_sig_cortex_merged.ID".format(cell_type,cell_type),
        names = ["ID"], sep="\t")
cell_intersect_sig_cortex = pd.merge(sig_interactions, cell_intersect_sig_cortex_ID,
                                     on = ["ID"], how = "inner")
# 1x pairs
# cell_intersect_1xcontrol_cortex_ID = pd.read_table(
# r"C:\Users\libin\UCSF\eQTL\hfb\{}\{}_intersect_1xcontrol_cortex_merged.ID".format(cell_type,cell_type),
# names = ["ID"])
# cell_intersect_1xcontrol_cortex = pd.merge(sig_interactions, cell_intersect_1xcontrol_cortex_ID,
#                                         on=["ID"], how="inner")'''
# 3x pairs
cell_intersect_3xcontrol_cortex_ID = pd.read_table(
        r"C:\Users\libin\UCSF\eQTL\hfb\{}\{}_intersect_3xcontrol_cortex_merged.ID".format(cell_type,cell_type),
        names = ["ID"])
cell_intersect_3xcontrol_cortex = pd.merge(sig_interactions, cell_intersect_3xcontrol_cortex_ID,
                                           on=["ID"], how = "inner")


sig_IQR = cell_intersect_sig_cortex["normalized"].describe()["75%"] - cell_intersect_sig_cortex["normalized"].describe()["25%"]
sig_upper =  cell_intersect_sig_cortex["normalized"].describe()["75%"] + 1.5*sig_IQR
# control1x_IQR = cell_intersect_1xcontrol_cortex["normalized"].describe()["75%"] - cell_intersect_1xcontrol_cortex["normalized"].describe()["25%"]
# control1x_upper =  cell_intersect_1xcontrol_cortex["normalized"].describe()["75%"] + 1.5*control1x_IQR
control3x_IQR = cell_intersect_3xcontrol_cortex["normalized"].describe()["75%"] - cell_intersect_3xcontrol_cortex["normalized"].describe()["25%"]
control3x_upper =  cell_intersect_3xcontrol_cortex["normalized"].describe()["75%"] + 1.5*control3x_IQR

cell_intersect_sig_cortex_IQR = cell_intersect_sig_cortex[cell_intersect_sig_cortex["normalized"] < sig_upper]
# cell_intersect_1xcontrol_cortex_IQR = cell_intersect_1xcontrol_cortex[cell_intersect_1xcontrol_cortex["normalized"] < control1x_upper]
cell_intersect_3xcontrol_cortex_IQR = cell_intersect_3xcontrol_cortex[cell_intersect_3xcontrol_cortex["normalized"] < control3x_upper]

# ---------output----------
print(cell_type, cell_intersect_sig_cortex["normalized"].mean())
# print("1X control", cell_intersect_1xcontrol_cortex["normalized"].mean())
print("3x control", cell_intersect_3xcontrol_cortex["normalized"].mean())

# print(cell_intersect_sig_cortex["normalized"].median())
# print(cell_intersect_1xcontrol_cortex["normalized"].median())
# print(cell_intersect_3xcontrol_cortex["normalized"].median())

eqtl_for_plot = pd.DataFrame()
eqtl_for_plot["3x_control"] = cell_intersect_3xcontrol_cortex["normalized"]
# eqtl_for_plot["1x_control"] = cell_intersect_1xcontrol_cortex["normalized"]
eqtl_for_plot["{}".format(cell_type)] = cell_intersect_sig_cortex["normalized"]

plt.figure(figsize=(11,8))
sns.violinplot(x="variable", y="value", data=eqtl_for_plot.melt(), palette="Set2", cut=0)

# print(stats.ks_2samp(cell_intersect_sig_cortex["normalized"], cell_intersect_1xcontrol_cortex["normalized"]))
print(stats.ks_2samp(cell_intersect_sig_cortex["normalized"], cell_intersect_3xcontrol_cortex["normalized"]))

eqtl_for_plot_IQR = pd.DataFrame()
eqtl_for_plot_IQR["3x_control"] = cell_intersect_3xcontrol_cortex_IQR["normalized"]
# eqtl_for_plot_IQR["1x_control"] = cell_intersect_1xcontrol_cortex_IQR["normalized"]
eqtl_for_plot_IQR["{}".format(cell_type)] = cell_intersect_sig_cortex_IQR["normalized"]

print (cell_intersect_sig_cortex["normalized"].describe())
# print (cell_intersect_1xcontrol_cortex["normalized"].describe())
print (cell_intersect_3xcontrol_cortex["normalized"].describe())

plt.figure(figsize=(11,8))
sns.violinplot(x="variable", y="value", data=eqtl_for_plot_IQR.melt(), palette="Pastel2", cut=0).set(ylim=(0,10))