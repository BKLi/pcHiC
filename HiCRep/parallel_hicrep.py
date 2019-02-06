import os
import itertools


def parallel_hicrep(sample_list, capture_flag, res, chrom, h_train, cutoff, outdir, time, memory):

    combinations = list(itertools.combinations(sample_list, 2))

    for x in combinations:
        sample = "".join(list(x))
        sample1 = sample[0:5] + capture_flag
        sample2 = sample[5:10] + capture_flag

        path = outdir+"\\"+sample+".R"
        print(path)
        outfile = open(path, "w+")
        outfile.write("library(hicrep)\n\n")
        outfile.write("args <- commandArgs(trailingOnly=TRUE)\n\n")
        trackLine = "print ('res=" + str(res) + "," + "chr=" + str(chrom) + "," + "h=" + str(
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

    # write R scripts and bash script for running
    files = [i for i in os.listdir(outdir)]
    configure = '''#!/bin/bash
                    #$ -l arch=linux-x64
                    #$ -l h_rt={}:0:0
                    #$ -l mem_free={}G
                    #$ -S /bin/bash
                    #$ -cwd
                    #$ -j y
                    #$ -r y
                    #$ -m ae
                    #$ -M libingkun1997@gmail.com\n\n'''.format(time, memory)

    for f in files:
        fname = outdir + "\\" + f[:-2] + ".sh"
        matrix1 = f[0:5] + capture_flag + ".matrix"
        matrix2 = f[5:10] + capture_flag + ".matrix"
        print(fname)
        outfile = open(fname, "w+")
        outfile.write(configure)
        outfile.write("Rscript " + f + " " + matrix1 + " " + matrix2)


if __name__ == "__main__":

    parallel_hicrep(["MS001", "MS002", "MS003",
                     "JC005", "JC006",
                     "MS006", "MS032",
                     "MS128", "MS127",
                     "MS137", "MS142",
                     "MS136", "MS140", "MS141"],
                    "_0.8", 10000, 1, 20, 5000000, "C:\\Users\libin\\UCSF\PCHiC", 128, 40)

'''parallel_hicrep_jobs(sample_list=sys.argv[1],captureFlag=sys.argv[2],res=sys.argv[3],
                    chr=sys.argv[4],h_train=sys.argv[5],cutoff=sys.argv[6],out_script=sys.argv[7])'''



