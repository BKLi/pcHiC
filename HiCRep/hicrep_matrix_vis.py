# -*- coding: utf-8 -*-
'''
visualizing SCC matrix of HiCRep
'''

import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

ctypes = ["IJ050","IJ051","IJ052","IJ053","MS051","MS052","MS053","MS081","MS082","MS083","MS084"]
hicrep_final = pd.read_csv(r"C:\Users\libin\UCSF\hfb\reg\hicrep_result\plac_seq_hicrep_march10th.csv", index_col= 0, sep=",")

for rep1 in ctypes[:-1]:
    for i in range(1,len(ctypes)-ctypes.index(rep1)):
        rep2 = ctypes[ctypes.index(rep1)+i]
        hicrep_final[rep1][rep2] = hicrep_final[rep2][rep1]
    

score_mean = (hicrep_final.values.sum() - hicrep_final.shape[0])/(hicrep_final.shape[0]*(hicrep_final.shape[0]-1))

plt.figure(figsize=(10,10))
sns.clustermap(hicrep_final, vmin=hicrep_final.values.min(), vmax=hicrep_final.values.max(), center=score_mean, cmap="RdBu_r", annot=True)
# sns.heatmap(hicrep_final)
plt.yticks(rotation=0)
# plt.savefig(r'C:\Users\libin\UCSF\hfb\reg\hicrep_result\PLAC_cluster_5kb_April8th_dsp_smt.pdf', transparent=True)
    


