# -*- coding: utf-8 -*-
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# prefix_1 = "cutoff_5_wiJY"
prefix_2 = "cutoff_0"

'''peakMatrix_CO5 = pd.read_table("C:\\Users\\libin\\UCSF\\Chicago\\12.4\\{}.txt".format(prefix_1))
peakMatrix_CO5_cis = peakMatrix_CO5[peakMatrix_CO5["baitChr"] == peakMatrix_CO5["oeChr"]]
# peakMatrix_CO5_trans = peakMatrix_CO5[peakMatrix_CO5["baitChr"] != peakMatrix_CO5["oeChr"]]
peakMatrix_CO5 = peakMatrix_CO5.drop(["baitChr", "baitStart","baitEnd","baitID","baitName","oeChr","oeStart","oeEnd","oeID","oeName"	,"dist"], axis=1)
peakMatrix_CO5_cis = peakMatrix_CO5_cis.drop(["baitChr", "baitStart","baitEnd","baitID","baitName","oeChr","oeStart","oeEnd","oeID","oeName"	,"dist"], axis=1)
# peakMatrix_CO5_trans = peakMatrix_CO5_trans.drop(["baitChr", "baitStart","baitEnd","baitID","baitName","oeChr","oeStart","oeEnd","oeID","oeName"	,"dist"], axis=1)

# test normaliztion
peakMatrix_CO5_above5 = peakMatrix_CO5[(peakMatrix_CO5["MS001"] >= 5) & (peakMatrix_CO5["MS002"] >= 5) & (peakMatrix_CO5["MS003"] >= 5) & (peakMatrix_CO5["MS006"] >= 5)
& (peakMatrix_CO5["MS032"] >= 5) & (peakMatrix_CO5["MS127"] >= 5) & (peakMatrix_CO5["MS128"] >= 5) & (peakMatrix_CO5["MS136"] >= 5) & (peakMatrix_CO5["MS140"] >= 5)
& (peakMatrix_CO5["MS141"] >= 5) & (peakMatrix_CO5["JC005"] >= 5) & (peakMatrix_CO5["JC006"] >= 5) & (peakMatrix_CO5["MS137"] >= 5) & (peakMatrix_CO5["MS142"] >= 5)]

peakMatrix_CO5_above5_cis = peakMatrix_CO5_cis[(peakMatrix_CO5_cis["MS001"] >= 5) & (peakMatrix_CO5_cis["MS002"] >= 5) & (peakMatrix_CO5_cis["MS003"] >= 5) & (peakMatrix_CO5_cis["MS006"] >= 5)
& (peakMatrix_CO5_cis["MS032"] >= 5) & (peakMatrix_CO5_cis["MS127"] >= 5) & (peakMatrix_CO5_cis["MS128"] >= 5) & (peakMatrix_CO5_cis["MS136"] >= 5) & (peakMatrix_CO5_cis["MS140"] >= 5)
& (peakMatrix_CO5_cis["MS141"] >= 5) & (peakMatrix_CO5_cis["JC005"] >= 5) & (peakMatrix_CO5_cis["JC006"] >= 5) & (peakMatrix_CO5_cis["MS137"] >= 5) & (peakMatrix_CO5_cis["MS142"] >= 5)]


peakMatrix_CO5_above5_cis = peakMatrix_CO5_cis[(peakMatrix_CO5_cis["cortical"] >= 5) & (peakMatrix_CO5_cis["astrocyte"] >= 5) & (peakMatrix_CO5_cis["motor"] >= 5) & (peakMatrix_CO5_cis["hippocampal"] >= 5)]
peakMatrix_CO5_above5 = peakMatrix_CO5[(peakMatrix_CO5["cortical"] >= 5) & (peakMatrix_CO5["astrocyte"] >= 5) & (peakMatrix_CO5["motor"] >= 5) & (peakMatrix_CO5["hippocampal"] >= 5)]
'''

'''peakMatrix_CO5_above0_cis = peakMatrix_CO5_cis[(peakMatrix_CO5_cis["MS001"] > 0) & (peakMatrix_CO5_cis["MS002"] > 0) & (peakMatrix_CO5_cis["MS003"] > 0) & (peakMatrix_CO5_cis["MS006"] > 0)
& (peakMatrix_CO5_cis["MS032"] > 0) & (peakMatrix_CO5_cis["MS127"] > 0) & (peakMatrix_CO5_cis["MS128"] > 0) & (peakMatrix_CO5_cis["MS136"] > 0) & (peakMatrix_CO5_cis["MS140"] > 0)
& (peakMatrix_CO5_cis["MS141"] > 0) & (peakMatrix_CO5_cis["JC005"] > 0) & (peakMatrix_CO5_cis["JC006"] > 0) & (peakMatrix_CO5_cis["MS137"] > 0) & (peakMatrix_CO5_cis["MS142"] > 0)]

peakMatrix_CO5_above10_cis = peakMatrix_CO5_cis[(peakMatrix_CO5_cis["MS001"] >= 10) & (peakMatrix_CO5_cis["MS002"] >= 10) & (peakMatrix_CO5_cis["MS003"] >= 10) & (peakMatrix_CO5_cis["MS006"] >= 10)
& (peakMatrix_CO5_cis["MS032"] >= 10) & (peakMatrix_CO5_cis["MS127"] >= 10) & (peakMatrix_CO5_cis["MS128"] >= 10) & (peakMatrix_CO5_cis["MS136"] >= 10) & (peakMatrix_CO5_cis["MS140"] >= 10)
& (peakMatrix_CO5_cis["MS141"] >= 10) & (peakMatrix_CO5_cis["JC005"] >= 10) & (peakMatrix_CO5_cis["JC006"] >= 10) & (peakMatrix_CO5_cis["MS137"] >= 10) & (peakMatrix_CO5_cis["MS142"] >= 10)]'''


peakMatrix_CO0 = pd.read_table("C:\\Users\\libin\\UCSF\\Chicago\\12.4\\{}.txt".format(prefix_2))
peakMatrix_CO0_cis = peakMatrix_CO0[peakMatrix_CO0["baitChr"] == peakMatrix_CO0["oeChr"]]

peakMatrix_CO0 = peakMatrix_CO0.drop(["baitChr", "baitStart","baitEnd","baitID","baitName","oeChr","oeStart","oeEnd","oeID","oeName"	,"dist"], axis=1)
peakMatrix_CO0_cis = peakMatrix_CO0_cis.drop(["baitChr", "baitStart","baitEnd","baitID","baitName","oeChr","oeStart","oeEnd","oeID","oeName"	,"dist"], axis=1)

'''peakMatrix_CO0_above5 = peakMatrix_CO0[(peakMatrix_CO0["MS001"] >= 5) & (peakMatrix_CO0["MS002"] >= 5) & (peakMatrix_CO0["MS003"] >= 5) & (peakMatrix_CO0["MS006"] >= 5)
& (peakMatrix_CO0["MS032"] >= 5) & (peakMatrix_CO0["MS127"] >= 5) & (peakMatrix_CO0["MS128"] >= 5) & (peakMatrix_CO0["MS136"] >= 5) & (peakMatrix_CO0["MS140"] >= 5)
& (peakMatrix_CO0["MS141"] >= 5) & (peakMatrix_CO0["JC005"] >= 5) & (peakMatrix_CO0["JC006"] >= 5) & (peakMatrix_CO0["MS137"] >= 5) & (peakMatrix_CO0["MS142"] >= 5)]

peakMatrix_CO0_above5_cis = peakMatrix_CO0_cis[(peakMatrix_CO0_cis["MS001"] >= 5) & (peakMatrix_CO0_cis["MS002"] >= 5) & (peakMatrix_CO0_cis["MS003"] >= 5) & (peakMatrix_CO0_cis["MS006"] >= 5)
& (peakMatrix_CO0_cis["MS032"] >= 5) & (peakMatrix_CO0_cis["MS127"] >= 5) & (peakMatrix_CO0_cis["MS128"] >= 5) & (peakMatrix_CO0_cis["MS136"] >= 5) & (peakMatrix_CO0_cis["MS140"] >= 5)
& (peakMatrix_CO0_cis["MS141"] >= 5) & (peakMatrix_CO0_cis["JC005"] >= 5) & (peakMatrix_CO0_cis["JC006"] >= 5) & (peakMatrix_CO0_cis["MS137"] >= 5) & (peakMatrix_CO0_cis["MS142"] >= 5)]
'''
# without JunYao
peakMatrix_CO0_above5 = peakMatrix_CO0[(peakMatrix_CO0["MS001"] >= 5) & (peakMatrix_CO0["MS002"] >= 5) & (peakMatrix_CO0["MS003"] >= 5) & (peakMatrix_CO0["MS127"] >= 5) & (peakMatrix_CO0["MS128"] >= 5) & (peakMatrix_CO0["MS136"] >= 5) & (peakMatrix_CO0["MS140"] >= 5)
& (peakMatrix_CO0["MS141"] >= 5) & (peakMatrix_CO0["JC005"] >= 5) & (peakMatrix_CO0["JC006"] >= 5) & (peakMatrix_CO0["MS137"] >= 5) & (peakMatrix_CO0["MS142"] >= 5)]

peakMatrix_CO0_above5_cis = peakMatrix_CO0_cis[(peakMatrix_CO0_cis["MS001"] >= 5) & (peakMatrix_CO0_cis["MS002"] >= 5) & (peakMatrix_CO0_cis["MS003"] >= 5) & (peakMatrix_CO0_cis["MS127"] >= 5) & (peakMatrix_CO0_cis["MS128"] >= 5) & (peakMatrix_CO0_cis["MS136"] >= 5) & (peakMatrix_CO0_cis["MS140"] >= 5)
& (peakMatrix_CO0_cis["MS141"] >= 5) & (peakMatrix_CO0_cis["JC005"] >= 5) & (peakMatrix_CO0_cis["JC006"] >= 5) & (peakMatrix_CO0_cis["MS137"] >= 5) & (peakMatrix_CO0_cis["MS142"] >= 5)]

'''plt.figure(figsize=(15,9))
sns.clustermap(peakMatrix_CO5.corr(), vmin=0, vmax=1, annot=True, cmap="RdBu_r")
plt.title("cutoff=5", y=1.1, x=1.3, fontsize=15)

plt.figure(figsize=(15,9))
sns.clustermap(peakMatrix_CO5_above5.corr(), vmin=0, vmax=1, annot=True, cmap="RdBu_r")
plt.title("cutoff=5, all_above_five", y=1.1, x=1.3, fontsize=15)

plt.figure(figsize=(15,9))
sns.clustermap(peakMatrix_CO5_above5_cis.corr(), vmin=0, vmax=1, annot=True, cmap="RdBu_r")
plt.title("cutoff=5, all_above_five, cis", y=1.1, x=1.3, fontsize=15)

plt.figure(figsize=(15,9))
sns.clustermap(peakMatrix_CO5_above0_cis.corr(), vmin=0, vmax=0.5, annot=True, cmap="RdBu_r")
plt.title("cutoff=5, all_above_zero, cis", y=1.1, x=1.3, fontsize=15)'''

'''plt.figure(figsize=(15,9))
sns.clustermap(peakMatrix_CO5_cis.corr(), vmin=0, vmax=1, annot=True, cmap="RdBu_r")
plt.title("cutoff=5, cis", y=1.1, x=1.3, fontsize=15)'''

'''plt.figure(figsize=(15,9))
sns.clustermap(peakMatrix_CO5_above5_cis,vmin=0, vmax=40, annot=False)
plt.title("cutoff=5, cis, above5", y=1.1, x=1.3, fontsize=15)'''

'''plt.figure(figsize=(15,9))
sns.clustermap(peakMatrix_CO0.corr(), vmin=peakMatrix_CO0.corr().values.min(), vmax=peakMatrix_CO0.corr().values.mean(), annot=True, cmap="RdBu_r")
plt.title("{}".format(prefix_2), y=1.1, x=1.3, fontsize=15)'''

plt.figure(figsize=(15,9))
sns.clustermap(peakMatrix_CO0_above5_cis, vmin=0, vmax=40, annot=False)
plt.title("{}".format(prefix_2), y=1.1, x=1.3, fontsize=15)


'''plt.figure(figsize=(15,9))
sns.clustermap(peakMatrix_CO0_cis.corr(), vmin=peakMatrix_CO0_cis.corr().values.min(), vmax=peakMatrix_CO0_cis.corr().values.mean(), annot=True, cmap="RdBu_r")
plt.title("{}, cis".format(prefix_2), y=1.1, x=1.3, fontsize=15)'''

'''plt.figure(figsize=(15,9))
sns.distplot(peakMatrix_CO5["MS142"], hist=False, label="MS142").set(xlim=(0,15))
sns.distplot(peakMatrix_CO5["MS137"], hist=False, label="MS137").set(xlim=(0,15))
sns.distplot(peakMatrix_CO5["MS141"], hist=False, label="MS141").set(xlim=(0,15))
sns.distplot(peakMatrix_CO5["MS006"], hist=False, label="MS136").set(xlim=(0,15))
sns.distplot(peakMatrix_CO0_cis["MS127"], hist=False, label="MS127").set(xlim=(0,15))
sns.distplot(peakMatrix_CO0_cis["MS128"], hist=False, label="MS128").set(xlim=(0,15))
sns.distplot(peakMatrix_CO0_cis["MS001"], hist=False, label="MS001").set(xlim=(0,15))
sns.distplot(peakMatrix_CO0_cis["MS002"], hist=False, label="MS002").set(xlim=(0,15))
sns.distplot(peakMatrix_CO0_cis["MS003"], hist=False, label="MS003").set(xlim=(0,15))'''

