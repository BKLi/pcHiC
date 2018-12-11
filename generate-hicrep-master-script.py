import sys
import itertools

'''master-script to generate master-script for each HiCRep job - without parallelization'''


def generate_matrix(sample_list, captureFlag, res, chr, h_train, cutoff, out_script):

    flagSampleList = [i + "_" + captureFlag for i in sample_list]
    print(flagSampleList)
    combinations = list(itertools.combinations(sample_list, 2))
    
    with open(out_script, "w+") as outfile:
        outfile.write("library(hicrep)\n\n")
        outfile.write("args <- commandArgs(trailingOnly=TRUE)\n\n")
        len_sample = len(sample_list)
        outfile.write("if (length(args)!=" + str(len_sample)+ ") "+"{\n\tstop('usage:hicrep',call.=FALSE)\n}\n")
        trackLine = "print ('res="+str(res)+","+"chr="+str(chr)+ ","+"h="+str(h_train)+","+"cutoff="+str(cutoff)+"')\n"
        outfile.write(trackLine)
        outfile.write("\n")

        for sample in flagSampleList:
            readSampleLine = sample+"<-read.delim2(args["+ str(flagSampleList.index(sample) + 1)+"],"+"header =FALSE,sep='\\t')\n"
            outfile.write(readSampleLine)
        outfile.write("\n")

        for i in combinations:
            j = "".join(list(i))
            sample1 = j[0:5] + "_" + captureFlag
            sample2 = j[5:10] + "_" + captureFlag
            line1 = j+" <- prep("+sample1+","+sample2+","+str(res)+","+str(h_train)+","+str(cutoff)+")\n"
            k = "SCC." + j
            line2 = k + " <- get.scc(" +j+","+str(res)+","+str(cutoff)+")\n"
            line3 = "print (" + "'" + k + "'" + ")\n"
            line4 = k + "[[3]]\n"
            line5 = k + "[[4]]\n"
            outfile.write("{}{}{}{}{}\n".format(line1,line2,line3,line4,line5))


if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("usage: generate-hicrep-script.py")
        sys.exit()

    generate_matrix(sample_list=sys.argv[1],captureFlag=sys.argv[2],res=sys.argv[3],
                    chr=sys.argv[4],h_train=sys.argv[5],cutoff=sys.argv[6],out_script=sys.argv[7])
# generate_matrix(["MS051", "MS082", "IJ052",
# "MS052", "MS083", "IJ050"],
# "", 5000, 1, 40, 5000000, "C:\\Users\libin\Desktop\hfb_NeuIPC.R")





