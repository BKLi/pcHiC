# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 17:08:19 2019

@author: libin
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


hippo_signif_distance = pd.read_table("C:\\Users\libin\\UCSF\eQTL\distance\\Brain_Hippocampus.signif.tss_filtered.txt.dist")
hippo_egenes_sampled_distance = pd.read_table("C:\\Users\libin\\UCSF\eQTL\distance\\Brain_Hippocampus.allpairs.egenes.sampled.dist")
hippo_signif_distance["tss_distance"] = hippo_signif_distance["tss_distance"].abs()

plt.figure()
sns.distplot(hippo_signif_distance["tss_distance"],color="#a7dbf4")
sns.distplot(hippo_egenes_sampled_distance["tss_distance"],color="#b0a7d7")
plt.savefig('C://Users//libin//UCSF/eQTL/distance/hippo_distance.pdf', transparent=True)

print(hippo_signif_distance["tss_distance"].min())

hippo_signif_chrom = pd.read_table("C:\\Users\libin\\UCSF\eQTL\distance\\Brain_Hippocampus.sig.chr", names=['chr'])
hippo_egenes_sampled_chrom = pd.read_table("C:\\Users\libin\\UCSF\eQTL\distance\\Brain_Hippocampus.allpairs.egenes.sampled.chr", names=['chr'])

with PdfPages("C://Users//libin//UCSF/eQTL/distance/hippo_chromosome.pdf") as pdf:
    plt.figure()                                                                                                                        
    hippo_signif_chrom["chr"].value_counts().plot(kind='bar')
    plt.title('hippo_signif_chrom')
    pdf.savefig()
    plt.figure()
    hippo_egenes_sampled_chrom["chr"].value_counts().plot(kind='bar')
    plt.title('hippo_egenes_sampled_chrom')
    pdf.savefig()