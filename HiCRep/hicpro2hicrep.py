#!/usr/bin/env python3.5

"""Take output matrices from HiC-Pro and convert them into HiCRep input format."""

import sys
import numpy as np
import time
import argparse

"""
    Parameters
    -----------

    input_bed_file: string
    Filename of binned bed file from HiC-Pro output

    input_matrix_file: string
    Filename of matrix file from HiC-Pro output

    output_bed_file: string
    Filename of binned bed file of given chromosome.
    Safe to delete after script finishes running.

    output_matrix_file: string
    Filename of matrix file of given chromosome.
    Safe to delete after script finishes running.

    integrated_file: string
    Filename of final output file.

    chrom: string (single number)
    Chromosome number.

    Notes
    ------
    This script is for generating chromosome-wide matrix only.
"""

start = time.clock()

parser = argparse.ArgumentParser()
parser.add_argument("-bed", "--inputBed", type=str, help="Filename of binned bed file from HiC-Pro output", required=True)
parser.add_argument("-matrix", "--inputMatrix", type=str, help="Filename of matrix file from HiC-Pro output", required=True)
parser.add_argument("-ob", "--outputBed", type=str, help="Filename of binned bed file of given chromosome", required=True)
parser.add_argument("-om", "--outputMatrix", type=str, help="Filename of matrix file of given chromosome", required=True)
parser.add_argument("-o", "--outputFile", type=str, help="Filename of final output file", required=True)
parser.add_argument("-chr", "--chromosome", type=str, help="Chromosome number", required=True)
args = parser.parse_args()


def integrate_matrix():

    input_bed_file = args.inputBed
    input_matrix_file = args.input_matrix_file
    output_bed_file = args.outputBed
    output_matrix_file = args.outputMatrix
    integrated_file = args.outputFile
    chrom = args.chromosome

    line_count = 0
    chr_number = "chr" + chrom + "\t"
    with open(input_bed_file) as bed_file:
        outfile1 = open(output_bed_file, "w+")
        for line in bed_file:
            if chr_number in line and "_" not in line:
                outfile1.write(line)
        outfile1.close()

    # extract .matrix records
    with open(input_matrix_file) as raw_matrix:
        with open(output_bed_file) as bed_file:
            with open(output_matrix_file, "w+") as outfile2:
                bed_id = {}
                for line in bed_file:
                    sline = line.strip().split()
                    bed_id[sline[3]] = ''
                for entry in raw_matrix:
                    sentry = entry.strip().split()
                    if sentry[0] in bed_id and sentry[1] in bed_id:
                        outfile2.write("\t".join(sentry))
                        outfile2.write("\n")

    # initialize zero matrix
    with open(output_bed_file) as bed_file:
        bed_file_lst = bed_file.readlines()
        # Starting ID of indicated chromosome
        initial_id = int(bed_file_lst[0].strip().split()[3])
        for i in bed_file_lst:
            line_count += 1
    with open(output_matrix_file) as matrix_file:
        initial_matrix = np.zeros((line_count, line_count), dtype=int)
        matrix_listed = initial_matrix.tolist()
        # fill matrix with contact numbers
        for entry in matrix_file:
            sentry = entry.strip().split()
            bin_1 = int(sentry[0]) - initial_id
            bin_2 = int(sentry[1]) - initial_id
            contact_number = int(sentry[2])
            matrix_listed[bin_1][bin_2] = contact_number
    with open(integrated_file, "w+") as outfile3:
        for lst, line in zip(matrix_listed,bed_file_lst):
            outfile3.write("\t".join(line.strip().split()[0:3]) + "\t" + "\t".join(str(i) for i in lst) + "\n")

    end = time.clock()
    print("program finished running in: ", end - start)


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("usage: hicpro2hicrep.py")
        sys.exit()

    integrate_matrix()

