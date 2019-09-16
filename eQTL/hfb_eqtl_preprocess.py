# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 11:27:57 2019

@author: libin
"""

import pandas as pd
import numpy as np

eqtl = pd.read_table(r'C:\Users\libin\UCSF\eQTL\hfb\top_eqtls_gene.txt', sep="\t")
eqtl = eqtl[eqtl["pval_nominal_threshold"] <= 1e-5]
