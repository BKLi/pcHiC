import pandas as pd
import os
import glob
import seaborn as sns
import matplotlib.pyplot as plt

''' ATAC-Seq peaks motif enrichment analysis'''

# read in all results files in folder
# def motif_enrichment(known_results_folder):

known_results_folder = "C:\\Users\\libin\\UCSF\\motif_analysis\\200_human"
size = 200

list_of_result = glob.glob(known_results_folder + "\\*.txt")
cell_types = [os.path.split(path)[-1][13:-4] for path in list_of_result]
# 1e-749 will be read in as 0 if not read in as str
results = [pd.read_table(result, sep="\t", dtype=str) for result in list_of_result]

#data cleaning; rename motif name
## for hocomoco dataset
for res in results:
    res["Motif Name"] = res["Motif Name"].str.extract(r"([A-Za-z0-9]+)_.+", expand=False)
## for combined hocomoco
results_cleaned = [result.filter(["Motif Name", "Log P-value"], axis=1).drop_duplicates(subset="Motif Name") for result in results]
# results_cleaned = [result.filter(["Motif Name", "Log P-value"], axis=1) for result in results]
# extract top motifs
sig_results_cleaned = [result.iloc[0:101,:] for result in results_cleaned]
for (c, r) in zip(cell_types, sig_results_cleaned):
    ## for homer dataset
    # r["Motif Name"] = r["Motif Name"].str.extract(r"^([^\/]+)", expand=False)
    r.rename(columns={"Log P-value": c}, inplace=True)
    r.to_csv(known_results_folder + "\\" + c + str(size) + ".csv", sep=",", index=False)

# find common sig motifs
[sig_result.set_index("Motif Name", inplace=True) for sig_result in sig_results_cleaned]
common_motif = pd.concat(sig_results_cleaned, axis=1, join="inner")
common_motif.to_csv(known_results_folder + "\\common_motif.csv")

# find cell-type specific motifs
#[re.set_index("Motif Name", inplace=True) for re in results_cleaned]
#for (c1, r1) in zip(cell_types, results_cleaned):
# common_motif = common_motif.reset_index()
# drop motif name for plotting
# common_motif_1 = common_motif.drop(["Motif Name"], axis = 1)

common_motif = common_motif[common_motif.columns].astype(float)
# norm_common_motif = ((common_motif - common_motif.min())/ (common_motif.max() - common_motif.min())) * -100
# average = common_motif.mean(axis=0).mean()
# midpoint = 0 - (common_motif.values.max() - common_motif.values.min()) / 2

## Below are ploting parts

# normalize column wide
norm_col = common_motif.sub(common_motif.mean(axis=0),axis=1)
norm_col = norm_col.div(common_motif.std(axis=0),axis=1)

norm_row = common_motif.sub(common_motif.mean(axis=1),axis=0)
norm_row = norm_row.div(common_motif.std(axis=1),axis=0)

average_row = norm_row.mean(axis=0).mean()
average_col = norm_col.mean(axis=0).mean()
average_all = common_motif.mean(axis=0).mean()

plt.figure(figsize=(11,15))
ax = sns.heatmap(norm_col, vmin=norm_col.values.min(), vmax=norm_col.values.max(), cmap="RdBu", center=average_col, annot=True)
plt.title("Column-Normalized Log P-value", y =1.01, x=1.1, fontsize =13)
for item in ax.get_xticklabels():
    item.set_rotation(45)
    
plt.figure(figsize=(11,15))
ax1 = sns.heatmap(common_motif, vmin=common_motif.values.min(), vmax=common_motif.values.max(), cmap = "RdBu", center = average_all, annot=True)
plt.title("Log P-value", y =1.01, x=1.1, fontsize =13)
for item in ax.get_xticklabels():
    item.set_rotation(45)
    
plt.figure(figsize=(11,15))
ax2 = sns.heatmap(norm_row, vmin=norm_row.values.min(), vmax=norm_row.values.max(), center = average_row, cmap="RdBu", annot=True)
plt.title("Row-Normalized Log P-value", y =1.01, x=1.1, fontsize =13)
for item in ax.get_xticklabels():
    item.set_rotation(45)   

# plt.figure(figsize=(5,7))
cluster_row = sns.clustermap(common_motif, cmap = "RdBu", standard_scale=0, figsize=(11,15)) 
plt.title("Row-Normalized cluster", y =1.2, x=1.1, fontsize =13) 

cluster_column = sns.clustermap(common_motif, standard_scale=1, figsize=(11,15))
plt.title("Col-Normalized cluster", y =1.2, x=1.1, fontsize =13)


# motif_enrichment("C:\\Users\libin\Desktop\motif_analysis")

