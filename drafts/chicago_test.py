import pandas as pd
import os
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import re

'''
test script for calculating correlations between significant interactions called using hicpro/hicup output
'''
# print(cell_type)
cell_type = "JC005"
hicpro_df = pd.read_table("C:\\Users\libin\\UCSF\Chicago\\hicpro\{}.ibed".format(cell_type))
hicpro_df.sort_values("score", inplace=True, ascending=False)
# print(hicpro_df.head())

hicup_df = pd.read_table("C:\\Users\libin\\UCSF\Chicago\\hicup\{}.ibed".format(cell_type))
hicup_df.sort_values("score", inplace=True, ascending=False)
# print(hicup_df.head())

common_set = pd.merge(hicup_df, hicpro_df, how="inner", on=["bait_chr", "bait_start", "bait_end", "otherEnd_chr", "otherEnd_start", "otherEnd_end"])
# print(common_set.head())
# print(common_set.columns)
cols_to_keep = ['bait_chr', 'bait_start', 'bait_end', 'otherEnd_chr',
                'otherEnd_start', 'otherEnd_end', 'score_x', 'score_y']

common_set_reformed = common_set[cols_to_keep]
common_set_reformed.to_csv('{}.csv'.format(cell_type), sep="\t", index=False)
# print(common_set_reformed.head())

s_corr = common_set_reformed["score_x"].corr(common_set_reformed["score_y"], method="spearman")
p_corr = common_set_reformed["score_x"].corr(common_set_reformed["score_y"], method="")

# print("{} hicpro mean score: ".format(cell_type), hicpro_df["score"].mean())
# print("{} hicup mean score: ".format(cell_type), hicup_df["score"].mean())
print("{} spearman correlation: ".format(cell_type), s_corr)
print("{} hicpro: ".format(cell_type), hicpro_df.shape[0])
print("{} hicup; ".format(cell_type), hicup_df.shape[0])
print("{} common: ".format(cell_type), common_set_reformed.shape[0])
