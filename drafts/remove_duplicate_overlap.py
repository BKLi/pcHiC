import sys


def remove_overlap(zfile, duplicate_list, rmdup_zfile):
    with open(duplicate_list) as duplicate:
        dup = {}
        for i in duplicate:
            si = i.strip().split()[0]
            dup[si] = ""
    with open(zfile) as Zfile:
        with open(rmdup_zfile, 'w+') as outfile:
            for line in Zfile:
                sline = line.strip().split()
                if sline[0] in dup:
                    outfile.write("\t".join(sline))
                    outfile.write('\n')


remove_overlap(zfile=sys.argv[1], duplicate_list=sys.argv[2], rmdup_zfile=sys.argv[3])

