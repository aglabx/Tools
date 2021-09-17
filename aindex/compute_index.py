#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 01.01.2017
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

import argparse
from PyExp import runner
import os

'''
~/Dropbox/Ariadna/Stirka/bin/R2F.exe pe500.mtDNA.reads pe500.mtDNA fasta
jellyfish count -m 23 -C -t 60 -s 20G -o pe500.mtDNA.23.jf2 pe500.mtDNA.fa 
jellyfish dump -c -t -o pe500.mtDNA.23.dat pe500.mtDNA.23.jf2
cut -f1 pe500.mtDNA.23.dat > pe500.mtDNA.23.kmers
~/Dropbox/Stirka/emphf/compute_mphf_seq.exe pe500.mtDNA.23.kmers pe500.mtDNA.23.pf
~/Dropbox/Ariadna/Stirka/bin/LU_index.exe pe500.mtDNA.23.dat pe500.mtDNA.23.pf pe500.mtDNA.23
~/Dropbox/Ariadna/Stirka/bin/LU_aindex.exe pe500.mtDNA.reads pe500.mtDNA.23.pf pe500.mtDNA.23 pe500.mtDNA.aindex 60 23 pe500.mtDNA.23.tf.bin 
'''

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Compute index from jf2 file.')
    parser.add_argument('-j', help='JF2 file if exists (None)', required=True)
    parser.add_argument('-o', help='Output prefix', required=True)
    parser.add_argument('--lu', help='-L for jellyfish', required=False, default=0)
    parser.add_argument('-P', help='Threads (12)', required=False, default=60)
    parser.add_argument('-M', help='JF2 memory in Gb (5)', required=False, default=1)

    args = vars(parser.parse_args())

    prefix = args["o"]
    threads = args["P"]
    jf2_file = args["j"]
    memory = args["M"]
    lu = args["lu"]

    commands = [
            "jellyfish dump -t -c -L %s -o %s.23.dat %s" % (lu, prefix, jf2_file),
        ]

    runner.run(commands)


    commands = [
        "cut -f1 %s.23.dat > %s.23.kmers" % (prefix, prefix),
        "~/Dropbox/Stirka/emphf/compute_mphf_seq.exe %s.23.kmers %s.23.pf" % (prefix, prefix),
        "~/Dropbox/Ariadna/Stirka/bin/LU_index.exe %s.23.dat %s.23.pf %s.23 %s" % (prefix, prefix, prefix, threads),
        "rm %s.23.dat %s.23.kmers" % (prefix, prefix),
    ]
    runner.run(commands)

   

