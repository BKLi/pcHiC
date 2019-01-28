import pandas as pd
import sys

def eQTL_analysis(left_hand, right_hand, merged_file):
    
    header = ["SNP_chr","SNP_start","SNP_end","gene_ID","p-value","eqtl_id","inter_chr","inter_start","inter_end","score","inter_ID"]

    lh_interactions = pd.read_table(left_hand, header=None, names=header, sep="\t")
    print("left hand read")
    rh_interactions = pd.read_table(right_hand, header=None, names=header, sep="\t")
    print("right hand read")
    
    intersect_both = pd.merge(lh_interactions, rh_interactions, how="inner", on=["SNP_chr","gene_ID","p-value","eqtl_id","score","inter_ID"])
    intersect_both.to_csv(merged_file, sep="\t", index=False)
    
eQTL_analysis(left_hand=sys.argv[1], right_hand=sys.argv[2], merged_file=sys.argv[3])