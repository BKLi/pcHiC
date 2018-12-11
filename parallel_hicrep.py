import sys
import itertools


def parallel_hicrep_jobs(sample_list, captureFlag, res, chr, h_train, cutoff, outdir):

    combinations = list(itertools.combinations(sample_list, 2))

    for x in combinations:
        sample = "".join(list(x))
        sample1 = sample[0:5] + captureFlag
        sample2 = sample[5:10] + captureFlag

        path = outdir+"\\"+sample+".R"
        print(path)
        outfile = open(path, "w+")
        outfile.write("library(hicrep)\n\n")
        outfile.write("args <- commandArgs(trailingOnly=TRUE)\n\n")
        trackLine = "print ('res=" + str(res) + "," + "chr=" + str(chr) + "," + "h=" + str(
            h_train) + "," + "cutoff=" + str(cutoff) + "')\n"
        outfile.write(trackLine)
        outfile.write("\n")
        readSampleLine1 = sample1 + "<-read.delim2(args[1]," + "header =FALSE,sep='\\t')\n"
        readSampleLine2 = sample2 + "<-read.delim2(args[2]," + "header =FALSE,sep='\\t')\n"
        outfile.write(readSampleLine1)
        outfile.write(readSampleLine2)
        outfile.write("\n")

        line1 = sample + " <- prep(" + sample1 + "," + sample2 + "," + str(res) + "," + str(h_train) + "," + str(
            cutoff) + ")\n"
        k = "SCC." + sample
        line2 = k + " <- get.scc(" + sample + "," + str(res) + "," + str(cutoff) + ")\n"
        line3 = "print (" + "'" + k + "'" + ")\n"
        line4 = k + "[[3]]\n"
        line5 = k + "[[4]]\n"
        outfile.write(line1)
        outfile.write(line2)
        outfile.write(line3)
        outfile.write(line4)
        outfile.write(line5)
        outfile.write("\n")


'''if __name__ == "__main__":
    if len(sys.argv) != 8:
        print ("usage: generate-hicrep-script.py")
        sys.exit()'''

'''parallel_hicrep_jobs(sample_list=sys.argv[1],captureFlag=sys.argv[2],res=sys.argv[3],
                    chr=sys.argv[4],h_train=sys.argv[5],cutoff=sys.argv[6],out_script=sys.argv[7])'''

parallel_hicrep_jobs(["MS001", "MS002", "MS003",
                      "JC005", "JC006"],
                      "_0.6", 20000, 1, 11, 5000000, "C:\\Users\libin\\UCSF\PCHiC")





