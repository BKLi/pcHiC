import numpy as np
import sys


def rebinning(chr, start, end, binsize, outfile):

    """
    parameters
    -----------

    chr: name of chromosome, format as "chr15"
    start: start position
    end: end position
    binsize: size of bin
    outfile:  output file
    """

    start = int(start)
    end = int(end) + 1
    binsize = int(binsize)

    bins = np.arange(start, end, binsize)
    with open(outfile, "w+") as out:
        for i in range(len(bins)-1):
            sublist = bins[i:i+2]
            out.write("{} {}\n".format(chr, " ".join([str(i) for i in sublist])))
        # print(sublist)


if __name__ == "__main__":

    # rebinning(chr=sys.argv[1], start=sys.argv[2], end=sys.argv[3], binsize=sys.argv[4], outfile=sys.argv[5])
    rebinning("chr2", 46400000, 48810000, 1000, "C:\\Users\libin\Desktop\\MSH2_MSH6.bed")
