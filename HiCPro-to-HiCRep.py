#!/netapp/home/yinshen/miniconda3/bin/python

import sys
import os
import numpy as np
import time
import re

start = time.clock()

def integrate_matrix(raw_bed_file, raw_matrix_file, output_bed_file, output_matrix_file, integrated_file, chr):
    # extract .bed records
    # print ("Input: raw bed file, raw matrix file, extracted bed file, extracted matrix file, integrated file")
    line_count = 0
    chr_number = "chr" + chr + "\t"
    with open(raw_bed_file) as bed_file:
        outfile1 = open(output_bed_file, "w+")
        for line in bed_file:
            if chr_number in line and "_" not in line:
                outfile1.write(line)
        outfile1.close()

    # extract .matrix records
    with open(raw_matrix_file) as raw_matrix:
        with open(output_bed_file) as bed_file:
            with open(output_matrix_file, "w+") as outfile2:
                ID = {}
                for line in bed_file:
                    sline = line.strip().split()
                    ID[sline[3]] = ''
                for entry in raw_matrix:
                    sentry = entry.strip().split()
                    if sentry[0] in ID and sentry[1] in ID:
                        outfile2.write("\t".join(sentry))
                        outfile2.write("\n")

    # initialize zero matrix
    with open(output_bed_file) as bed_file:
        bed_file_lst = bed_file.readlines()
        initial_ID = int(bed_file_lst[0].strip().split()[3])    # Starting ID of indicated chromosome
        for i in bed_file_lst:
            line_count += 1
    with open(output_matrix_file) as matrix_file:
        initial_matrix = np.zeros((line_count, line_count), dtype=int)
        matrix_listed = initial_matrix.tolist()
    # fill matrix with contact numbers
        for entry in matrix_file:
            sentry = entry.strip().split()
            bin_1 = int(sentry[0]) - initial_ID
            bin_2 = int(sentry[1]) - initial_ID
            contact_number = int(sentry[2])
            matrix_listed[bin_1][bin_2] = contact_number
    with open(integrated_file, "w+") as outfile3:
        for lst,line in zip(matrix_listed,bed_file_lst):
            outfile3.write("\t".join(line.strip().split()[0:3]) + "\t" +"\t".join(str(i) for i in lst) + "\n")

    end = time.clock()
    print(end - start)

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print ("usage: HiCrep_prepare.py")
        sys.exit()

    integrate_matrix(raw_bed_file = sys.argv[1],raw_matrix_file= sys.argv[2],
                     output_bed_file = sys.argv[3], output_matrix_file = sys.argv[4],
                     integrated_file = sys.argv[5], chr = sys.argv[6])

