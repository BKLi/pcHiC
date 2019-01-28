# -*- coding: utf-8 -*-
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2



'''compare 0/5 cutoff of fastp'''

replicate = "hippocampal_2"

zero_cutoff = pd.read_table("C:\\Users\\libin\\UCSF\\Chicago\\12.4\\fastp\\250_0\\{}.ibed".format(replicate))
five_cutoff = pd.read_table("C:\\Users\\libin\\UCSF\\Chicago\\12.4\\fastp\\250_5\\{}.ibed".format(replicate))
zero_cutoff_above5 = zero_cutoff[zero_cutoff["score"] >= 5]
five_cutoff_cross_zero_cutoff_above5 = pd.merge(five_cutoff,zero_cutoff_above5,how="inner", on=["bait_chr", "bait_start", "bait_end", "otherEnd_chr", "otherEnd_start", "otherEnd_end"])

# cutoff_merged = pd.merge(five_cutoff, zero_cutoff, how="outer", on=["bait_chr", "bait_start", "bait_end", "otherEnd_chr", "otherEnd_start", "otherEnd_end"])
cutoff_intersect = pd.merge(five_cutoff, zero_cutoff, how="inner", on=["bait_chr", "bait_start", "bait_end", "otherEnd_chr", "otherEnd_start", "otherEnd_end"])
cutoff_intersect["difference"] = cutoff_intersect["score_x"]-cutoff_intersect["score_y"]

corr_intersect_spearman = cutoff_intersect["score_x"].corr(cutoff_intersect["score_y"], method="spearman")
corr_intersect_pearson = cutoff_intersect["score_x"].corr(cutoff_intersect["score_y"], method="pearson")

plt.figure(figsize=(15,9))
# sns.distplot(cutoff_intersect["score_x"], hist=False, label="cutoff=5").set(xlim=(0,20))
# sns.distplot(cutoff_intersect["score_y"], hist=False, label="cutoff=0").set(xlim=(0,20))
sns.distplot(cutoff_intersect["difference"], hist=False, label="difference").set(xlim=(0,7))


plt.figure(figsize=(15,9))
venn2(subsets = (five_cutoff.shape[0]-five_cutoff_cross_zero_cutoff_above5.shape[0], zero_cutoff_above5.shape[0]-five_cutoff_cross_zero_cutoff_above5.shape[0], five_cutoff_cross_zero_cutoff_above5.shape[0]), set_labels=("five_cutoff", "zero_cutoff_above5"))
plt.show()



'''compare fastp/untrimmed '''

'''replicate = "MS001"
# cutoff = 0
fastp_co0 = pd.read_table("C:\\Users\\libin\\UCSF\\Chicago\\12.4\\fastp\\250_0\\{}.ibed".format(replicate))
untrimmed_co0 = pd.read_table("C:\\Users\\libin\\UCSF\\Chicago\\12.4\\untrimmed\\250_0\\{}.ibed".format(replicate))

trimmed_merged = pd.merge(fastp_co0, untrimmed_co0, how="outer", on=["bait_chr", "bait_start", "bait_end", "otherEnd_chr", "otherEnd_start", "otherEnd_end"])

fastp_above_5 = fastp_co0[fastp_co0["score"] >= 5]
fastp_above_5_cross_untrimmed = pd.merge(fastp_above_5, untrimmed_co0, how="inner", on=["bait_chr", "bait_start", "bait_end", "otherEnd_chr", "otherEnd_start", "otherEnd_end"])
plt.figure(figsize=(15,9))
sns.distplot(fastp_above_5_cross_untrimmed["score_x"],hist=False, label="fastp").set(xlim=(0,20))
sns.distplot(fastp_above_5_cross_untrimmed["score_y"],hist=False, label="untrimmed").set(xlim=(0,20))
plt.title("score distribution (fastp > 5)", y=1.01, x=1.1, fontsize=13)
print("fastp>5 spearman", fastp_above_5_cross_untrimmed["score_x"].corr(fastp_above_5_cross_untrimmed["score_y"], method="spearman"))
print("fastp>5 pearson", fastp_above_5_cross_untrimmed["score_x"].corr(fastp_above_5_cross_untrimmed["score_y"], method="pearson"))

untrimmed_above_5 = untrimmed_co0[untrimmed_co0["score"] >= 5]
untrimmed_above5_cross_fastp = pd.merge(untrimmed_above_5, fastp_co0, how="inner", on=["bait_chr", "bait_start", "bait_end", "otherEnd_chr", "otherEnd_start", "otherEnd_end"])
plt.figure(figsize=(15,9))
sns.distplot(untrimmed_above5_cross_fastp["score_x"],hist=False, label="untrimmed").set(xlim=(0,20))
sns.distplot(untrimmed_above5_cross_fastp["score_y"],hist=False, label="fastp").set(xlim=(0,20))
plt.title("score distribution (untrimmed > 5)", y=1.01, x=1.1, fontsize=13)
print("untrimmed>5 spearman", untrimmed_above5_cross_fastp["score_x"].corr(untrimmed_above5_cross_fastp["score_y"], method="spearman"))
print("untrimmed>5 pearson", untrimmed_above5_cross_fastp["score_x"].corr(untrimmed_above5_cross_fastp["score_y"], method="pearson"))

trimmed_intersect = pd.merge(fastp_co0, untrimmed_co0, how="inner", on=["bait_chr", "bait_start", "bait_end", "otherEnd_chr", "otherEnd_start", "otherEnd_end"])
trimmed_above5_intersect = pd.merge(fastp_above_5, untrimmed_above_5, how="inner", on=["bait_chr", "bait_start", "bait_end", "otherEnd_chr", "otherEnd_start", "otherEnd_end"])

plt.figure(figsize=(15,9))
venn2(subsets = (fastp_above_5.shape[0]-trimmed_above5_intersect.shape[0], untrimmed_above_5.shape[0]-trimmed_above5_intersect.shape[0], trimmed_above5_intersect.shape[0]), set_labels=("fastp", "untrimmed"))
plt.show()'''