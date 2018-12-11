# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 14:42:15 2018

@author: libin
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

depth = ["1", "0.8", "0.6", "0.4", "0.2", "0.1"]
cutoff = "0"

hippo_ctypes = ["MS137", "MS142"]
cortical_ctypes = ["MS001","MS002","MS003"]
astro_ctypes = ["JC005","JC006","MS127","MS128"]
motor_ctypes = ["MS136","MS140","MS141"]
corr_all = []
corr_names = []

peakMatrixList = []
MS001_df = pd.DataFrame()
MS137_df = pd.DataFrame()
JC005_df = pd.DataFrame()
MS127_df = pd.DataFrame()
MS136_df = pd.DataFrame()

for d in depth:
    peak_matrix = pd.read_table("C:\\Users\\libin\\UCSF\\Chicago\\saturation\\cutoff_{}_{}.txt".format(cutoff,d))
    peakMatrixList.append(peak_matrix)
    MS001_df[d] = peak_matrix["MS001"]
    MS127_df[d] = peak_matrix["MS127"]
    MS136_df[d] = peak_matrix["MS136"]
    MS137_df[d] = peak_matrix["MS137"]
    JC005_df[d] = peak_matrix["JC005"]

print(JC005_df.corr())
    