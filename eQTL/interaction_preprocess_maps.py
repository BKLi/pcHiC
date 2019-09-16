# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 21:30:17 2019

@author: libin
"""

""" pre_process of MAPS output"""


import pandas as pd
import sys


def interaction_preprocess(interaction_file, cell_type):
    interaction = pd.read_table(interaction_file, sep="\t")
    interacion_ID = [cell_type+str(i+1) for i in range(interaction.shape[0])]
    interaction["ID"] = interacion_ID
    
    lhs_interaction = interaction[["chr1", "start1", "end1", "ID", "fdr", "ClusterNegLog10P", "count", "expected"]]
    rhs_interaction = interaction[["chr2", "start2", "end2", "ID", "fdr", "ClusterNegLog10P", "count", "expected"]]

    outfile_lhs = "{}_lh_interactions".format(cell_type)
    outfile_rhs = "{}_rh_interactions".format(cell_type) 
    
    lhs_interaction.to_csv(outfile_lhs,sep="\t", index=False, header=False)
    rhs_interaction.to_csv(outfile_rhs,sep="\t", index=False, header=False)

    
interaction_preprocess(interaction_file=sys.argv[1], cell_type=sys.argv[2])