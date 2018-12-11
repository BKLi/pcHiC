import os
import sys


def generate_qsub(inDir, captureFlag):
    files = [i for i in os.listdir(inDir)]
    configure = '''#!/bin/bash
#$ -l arch=linux-x64
#$ -l h_rt=256:0:0
#$ -l mem_free=168G
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -r y
#$ -m ae
#$ -M libingkun1997@gmail.com\n\n'''

    for f in files:
        fname = inDir + "\\" + f[:-2] + ".sh"
        matrix1 = f[0:5] + captureFlag + ".matrix"
        matrix2 = f[5:10] + captureFlag + ".matrix"
        print(fname)
        outfile = open(fname, "w+")
        outfile.write(configure)
        outfile.write("Rscript " + f + " " + matrix1 + " " + matrix2)


generate_qsub(inDir="C:\\Users\libin\\UCSF\PCHiC", captureFlag="_0.6")
