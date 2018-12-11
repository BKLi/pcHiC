# -*- coding: utf-8 -*-

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

ctypes = ["MS001","MS002","MS003","MS137","MS142","JC005","JC006","MS127","MS128","MS136","MS140","MS141"]
hicrep_final = pd.read_csv("C:\\Users\\libin\\Desktop\\tmp\\HiCRep-10kb_all.csv", index_col= 0, sep=",")

for rep1 in ctypes[:-1]:
    for i in range(1,len(ctypes)-ctypes.index(rep1)):
        rep2 = ctypes[ctypes.index(rep1)+i]
        hicrep_final[rep1][rep2] = hicrep_final[rep2][rep1]
    
plt.figure(figsize=(10,10))
score_mean = (hicrep_final.values.sum() - hicrep_final.shape[0])/(hicrep_final.shape[0]*(hicrep_final.shape[0]-1))
sns.clustermap(hicrep_final, vmin=hicrep_final.values.min(), vmax=hicrep_final.values.max(), center=score_mean, cmap="RdBu_r", annot=True)
# sns.heatmap(hicrep_final)

    


