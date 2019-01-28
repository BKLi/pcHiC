import pandas as pd
import glob
import os

'''
calc corr between hicpro/hicup output
'''

pd.set_option('display.max_columns', 30)

list_of_cell_types = glob.glob("C:\\Users\libin\\UCSF\Chicago\\\hicpro\*.ibed")
# print(list_of_cell_types)
for cell in list_of_cell_types:
    cell_type = "".join(os.path.split(cell)[-1])[:-5]
    print(cell_type)

    hicpro_df = pd.read_table("C:\\Users\libin\\UCSF\Chicago\\hicpro\{}.ibed".format(cell_type))
    hicpro_df.sort_values("score", inplace=True, ascending=False)
    # print(hicpro_df.head())

    hicup_df = pd.read_table("C:\\Users\libin\\UCSF\Chicago\\hicup\{}.ibed".format(cell_type))
    hicup_df.sort_values("score", inplace=True, ascending=False)
    # print(hicup_df.head())

    common_set = pd.merge(hicup_df, hicpro_df, how="outer", on=["bait_chr", "bait_start", "bait_end", "otherEnd_chr", "otherEnd_start", "otherEnd_end"])
    common_set.fillna(0, inplace=True)
    corr_spearman = common_set["score_x"].corr(common_set["score_y"], method="spearman")
    corr_pearson = common_set["score_x"].corr(common_set["score_y"], method="pearson")
    print(cell_type, "spearman", corr_spearman)
    print(cell_type, "pearson", corr_pearson)
    # print(common_set.shape[0])
