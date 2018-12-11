# -*- coding: utf-8 -*-

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

score_df_final = pd.read_csv("C:\\Users\\libin\\UCSF\Chicago\hicpro\correlation matrix.csv", index_col= 0, sep=",")

plt.figure(figsize=(12,12))
score_mean = (score_df_final.values.sum() - score_df_final.shape[0])/(score_df_final.shape[0]*(score_df_final.shape[0]-1))
sns.clustermap(score_df_final, vmin=score_df_final.values.min(), vmax=score_df_final.values.max(), center=score_mean, cmap="RdBu_r", annot=True)

