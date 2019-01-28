# -*- coding: utf-8 -*-
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

''' test, visualization'''

cell_type = "hippocampal_2"
common_set_file = "C:\\Users\libin\\UCSF\Chicago\{}.csv".format(cell_type)
common_set = pd.read_table(common_set_file, sep="\t")
common_set_cleaned = common_set[(common_set["score_x"] >= 5) & (common_set["score_y"] >= 5)]
print(common_set_cleaned.head())

# plt.figure(figsize=(15,9))
# sns.distplot(common_set_cleaned["score_x"])

plt.figure(figsize=(15,9))
sns.distplot(common_set_cleaned["score_y"])

plt.figure(figsize=(20,20))
sns.regplot(x="score_x", y="score_y", data=common_set, fit_reg=False)

sns.kdeplot(common_set_cleaned.score_x,common_set_cleaned.score_y)
