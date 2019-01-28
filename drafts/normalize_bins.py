import sys


def normalize_peaks(infile, bin, outfile):
    with open(infile) as f:
        with open(outfile,"w+") as out_f:
            for line in f:
                sline = line.strip().split()
                if int(sline[2]) - int(sline[1]) < int(bin):
                    midpoint_1 = (int(sline[2]) + int(sline[1]))/2
                    sline[1] = str(round(midpoint_1 - int(bin)/2))
                    sline[2] = str(round(midpoint_1 + int(bin)/2))
                out_f.write("\t".join(sline))
                out_f.write("\n")


normalize_peaks(infile=sys.argv[1],bin=sys.argv[2],outfile=sys.argv[3])