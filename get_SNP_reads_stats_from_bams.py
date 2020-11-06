#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 10.10.2020
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

    
import argparse
import os
import re

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Get all.')
    parser.add_argument('-i','--input', help='Folder with sorted BAM files.', required=True)
    parser.add_argument('-o','--output', help='Output file with results', required=True)
    parser.add_argument('-c','--coordintes', help='SNP coordinates in 5:72589809 format', required=True)
    parser.add_argument('-l','--leftflank', help='Left flank for SNP, 15-23 nucleotides', required=True)
    parser.add_argument('-r','--rightflank', help='Right flank for snp, 15-23 nucleotides', required=True)
    args = vars(parser.parse_args())
    folder = args["input"]
    output_file_table = args["output"]
    output_temp_file = output_file_table + ".temp.recallx"
    left  = args["leftflank"]
    right = args["rightflank"]

    chrom, pos = args["coordintes"].split(":")
    
    command = "samtools view -S %s %s:%s-%s | cut -f 10,11 > %s"

    with open(output_file_table, "w") as fw:
        for root, dirs, files in os.walk(folder, topdown=False):
            for name in files:
                if re.search(".bam$", name):
                    apath = os.path.join(root, name)
                    c = command % (apath, chrm, pos, pos, output_temp_file)
                    os.system(c)
                    print("Computing:", apath)

                    d = [
                        name.split(".")[0], # sample
                        0, # coverage
                        0, # T
                        0, # G
                        0, # C
                        0, # A
                        0, # uncallable
                        [], # Q support T
                        [], # Q support G
                        [], # Q support C
                        [], # Q support A
                        
                    ]
                    
                    with open(output_temp_file) as fh:
                        for line in fh:
                        
                            d[1] += 1 ### increase coverage
                        
                            read, Q = line.strip().split()
                            lpos = read.find(left)
                            if lpos >= 0:
                                pos = lpos + len(left)
                            else:
                                rpos = read.find(right)
                                if not rpos >= 0:
                                    d[6] += 1
                                    continue
                                pos = rpos - 1
                                
                            variant = read[pos]
                            Q = ord(Q[pos])-33
                            if variant == "T":
                                d[2] += 1
                                d[7].append(Q)
                            elif variant == "G":
                                d[3] += 1
                                d[8].append(Q)
                            elif variant == "С":
                                d[4] += 1
                                d[9].append(Q)
                            else:
                               d[5] += 1
                               d[10].append(Q)
                        
                    d[7] = ",".join(map(str,d[7]))
                    d[8] = ",".join(map(str,d[8]))
                    d[9] = ",".join(map(str,d[9]))
                    d[10] = ",".join(map(str,d[10]))
                    
                    fw.write("%s\n" % "\t".join(map(str, d)))
    os.unlink(output_temp_file)
    print("Done.")
