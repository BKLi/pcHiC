import sys


def get_promotor(infile, bin, outfile):
    with open(infile) as f:
        with open(outfile,"w+") as out_f:
            for line in f:
                sline = line.strip().split()
                if sline[0].startswith("#"):
                    continue
                else:
                    newLine = []
                    TSS = int(sline[1])
                    start = str(round(TSS - int(bin)/2))
                    end = str(round(TSS + int(bin)/2))
                    newLine.extend((sline[0], start, end))
                    out_f.write("\t".join(newLine))
                    out_f.write("\n")


get_promotor(infile=sys.argv[1], bin=sys.argv[2], outfile=sys.argv[3])