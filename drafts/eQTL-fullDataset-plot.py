# -*- coding: utf-8 -*-
"""
Created on Mon Dec 31 19:28:38 2018

@author: libin
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

hippo_self_fullSet = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippocampus_intersect_sig_merged", sep="\t")
hippo_blood_fullSet = pd.read_table("C:\\Users\\libin\\UCSF\\eQTL\\hippocampus_blood_intersect_sig_merged", sep="\t")

