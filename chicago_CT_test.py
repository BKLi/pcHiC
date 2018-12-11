import pandas as pd
import glob
import os
import sys

'''
test script, calc pairwise corr between every pair of replicates
'''


def chicago_cross_CT(hicpro_folder, hicup_folder):
    
    cell_type_file1 = glob.glob("{}*.ibed".format(hicpro_folder))
    print(cell_type_file1)
    cell_type_file2 = glob.glob("{}*.ibed".format(hicup_folder))
    cell_type1 = ["".join(os.path.split(cell1)[-1])[:-5] for cell1 in cell_type_file1]
    cell_type2 = ["".join(os.path.split(cell2)[-1])[:-5] for cell2 in cell_type_file2]

    for cell1 in cell_type1[:-1]:
        for cell2 in cell_type2[cell_type1.index(cell1)+1:len(cell_type2)+1]:
            print(cell1, cell2)
            
            hicpro_df_ctype1 = pd.read_table("{}{}.ibed".format(hicpro_folder, cell1))
            hicpro_df_ctype1.sort_values("score", inplace=True, ascending=False)

            hicpro_df_ctype2 = pd.read_table("{}{}.ibed".format(hicpro_folder, cell2))
            hicpro_df_ctype2.sort_values("score", inplace=True, ascending=False)
            
            hicup_df_ctype1 = pd.read_table("{}{}.ibed".format(hicup_folder, cell1))
            hicup_df_ctype1.sort_values("score", inplace=True, ascending=False)
            
            hicup_df_ctype2 = pd.read_table("{}{}.ibed".format(hicup_folder, cell2))
            hicup_df_ctype2.sort_values("score", inplace=True, ascending=False)

            common_set_hicpro = pd.merge(hicpro_df_ctype1, hicpro_df_ctype2, how="outer", on=["bait_chr", "bait_start", "bait_end", "otherEnd_chr", "otherEnd_start", "otherEnd_end"])
            common_set_hicpro.fillna(0, inplace=True)
            
            common_set_hicup = pd.merge(hicup_df_ctype1, hicup_df_ctype2, how="outer", on=["bait_chr", "bait_start", "bait_end", "otherEnd_chr", "otherEnd_start", "otherEnd_end"])
            common_set_hicup.fillna(0, inplace=True)

            s_corr_hicpro = common_set_hicpro["score_x"].corr(common_set_hicpro["score_y"], method="spearman")
            s_corr_hicup = common_set_hicup["score_x"].corr(common_set_hicup["score_y"], method="spearman")
            p_corr_hicpro = common_set_hicpro["score_x"].corr(common_set_hicpro["score_y"], method="pearson")
            p_corr_hicup = common_set_hicpro["score_x"].corr(common_set_hicup["score_y"], method="pearson")
            
            print("hicpro spearman correlation: ", s_corr_hicpro)
            print("hicup pearson correlation: ", p_corr_hicpro)
            print("hicup spearman correlation: ", s_corr_hicup)
            print("hicup pearson correlation: ", p_corr_hicup)
            print("---"*5)
            

chicago_cross_CT(hicpro_folder=sys.argv[1], hicup_folder=sys.argv[2])
# chicago_cross_CT(hicpro_folder=sys.argv[1], hicup_folder=sys.argv[2])
# outfile.write("{}{}\t{}\t{}\t{}\t{}\n".format(cell1,cell2,s_corr_hicpro,p_corr_hicpro,s_corr_hicup,p_corr_hicup))
# print("hicup common: ", common_set_reformed_hicup.shape[0])
