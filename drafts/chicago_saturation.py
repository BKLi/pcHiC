# -*- coding: utf-8 -*-

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

depth = ["0.1", "0.2", "0.4", "0.6", "0.8", "1"]
cutoff = "0"

hippo_ctypes = ["MS137", "MS142"]
cortical_ctypes = ["MS001","MS002","MS003"]
astro_ctypes = ["JC005","JC006","MS127","MS128"]
motor_ctypes = ["MS136","MS140","MS141"]
corr_all = []
corr_names = []

for d in depth:
    peakMatrix = pd.read_table("C:\\Users\\libin\\UCSF\\Chicago\\saturation\\cutoff_{}_{}.txt".format(cutoff,d))
    corrs = []
    
    for hc1 in hippo_ctypes[:-1]:
        for i in range(1,len(hippo_ctypes)-hippo_ctypes.index(hc1)):
            hc2 = hippo_ctypes[hippo_ctypes.index(hc1)+i]
            corr_names.append("{}-{}".format(hc1,hc2))
            hc_corr = peakMatrix[hc1].corr(peakMatrix[hc2], method="pearson")
            print (hc1, hc2, hc_corr)
            corrs.append(hc_corr)
            
            
    for cc1 in cortical_ctypes:
        for i in range(1,len(cortical_ctypes)-cortical_ctypes.index(cc1)):
            cc2 = cortical_ctypes[cortical_ctypes.index(cc1)+i]
            corr_names.append("{}-{}".format(cc1,cc2))
            cc_corr = peakMatrix[cc1].corr(peakMatrix[cc2], method="pearson")
            print (cc1, cc2, cc_corr)
            corrs.append(cc_corr)
            
    for ac1 in astro_ctypes:
        for i in range(1,len(astro_ctypes)-astro_ctypes.index(ac1)):
            ac2 = astro_ctypes[astro_ctypes.index(ac1)+i]
            corr_names.append("{}-{}".format(ac1,ac2))
            ac_corr = peakMatrix[ac1].corr(peakMatrix[ac2], method="pearson")
            print (ac1, ac2, ac_corr)
            corrs.append(ac_corr)
            
    for mc1 in motor_ctypes:
        for i in range(1,len(motor_ctypes)-motor_ctypes.index(mc1)):
            mc2 = motor_ctypes[motor_ctypes.index(mc1)+i]
            corr_names.append("{}-{}".format(mc1,mc2))
            mc_corr = peakMatrix[mc1].corr(peakMatrix[mc2], method="pearson")
            print (mc1, mc2, mc_corr)
            corrs.append(mc_corr)
            
    corr_all.append(corrs)
# print(corr_names)  
#print(corr_all)
corr_sorted = []
for j in range(0,13):
    scorr = [i[j] for i in corr_all]
    corr_sorted.append(scorr)
corr_pd = pd.DataFrame(corr_sorted, columns=depth)
print(corr_pd)
# print(corr_names[0:13])
corr_pd["rep"] = corr_names[0:13]
corr_pd = corr_pd.set_index("rep")
corr_pd = corr_pd.T

print(corr_pd)
plt.figure(figsize=(15,15))        
for column in corr_pd:
    plt.plot(corr_pd[column], label=column) 
    plt.legend(loc=2, ncol=2)             
    
