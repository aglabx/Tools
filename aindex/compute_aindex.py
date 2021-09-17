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

    parser = argparse.ArgumentParser(description='Compute index.')
    parser.add_argument('-i', help='Fasta, comma separated fastqs or reads', required=True)
    parser.add_argument('-j', help='JF2 file if exists (None)', required=False, default=None)
    parser.add_argument('-t', help='Reads type (reads)', required=False, default="fastq")
    parser.add_argument('-o', help='Output prefix', required=True)
    parser.add_argument('-H', help='Build header file for fasta', required=False, default=False)
    parser.add_argument('--lu', help='-L for jellyfish', required=False, default=0)
    parser.add_argument('-s', help='Sort dat file (None)', required=False, default=None)
    parser.add_argument('--interactive', help='Interactive (False)', required=False, default=None)
    parser.add_argument('-P', help='Threads (12)', required=False, default=12)
    parser.add_argument('-M', help='JF2 memory in Gb (5)', required=False, default=5)
    parser.add_argument('--onlyindex', help='Compute only index', required=False, default=False)

    args = vars(parser.parse_args())

    reads_file = args["i"]
    prefix = args["o"]
    reads_type = args["t"]
    threads = args["P"]
    jf2_file = args["j"]
    sort_dat_file = args["s"]
    memory = args["M"]
    interactive = args["interactive"]
    only_index = args["onlyindex"]

    lu = args["lu"]


    build_header = bool(args["H"])


    if not jf2_file:
        if reads_type == "reads":
            commands = [
                "python ~/Dropbox/workspace/PyBioSnippets/hiseq/reads_to_fasta.py -i %s -o %s.fa" % (reads_file, prefix),
                "jellyfish count -m 23 -t %s -s %sG -C -L %s -o %s.23.jf2 %s.fa" % (threads, memory, lu, prefix, prefix),
            ]
            runner.run(commands)
            if interactive:
                input("Continue?")
        elif reads_type == "fasta" or reads_type == "fastq":
            commands = [
                "jellyfish count -m 23 -t %s -s %sG -C -L %s -o %s.23.jf2 %s" % (threads, memory, lu, prefix, reads_file.replace(",", " ")),
            ]
            runner.run(commands)
            if interactive:
                input("Continue?")

        jf2_file = "%s.23.jf2" % prefix

    commands = []
    if reads_type == "fasta":
        if not build_header:
            commands = [
                "python ~/Dropbox/workspace/PyBioSnippets/hiseq/fasta_to_reads.py -i %s -o %s.reads" % (reads_file, prefix),   
            ]                
        else:
            commands = [
                "python ~/Dropbox/workspace/PyBioSnippets/hiseq/fasta_to_reads.py -i %s -o %s.reads -H %s.header" % (reads_file, prefix, prefix),   
            ]                
    if reads_type == "fastq":
        commands = [
            "~/Dropbox/workspace/Ariadna/Stirka/bin/V2_converter.exe %s fastq %s.reads" % (reads_file.replace(",", " "), prefix),   
        ]

    runner.run(commands)


    if lu:    
        commands = [
            "jellyfish dump -t -c -L %s -o %s.23.dat %s" % (lu, prefix, jf2_file),
        ]
    else:
        commands = [
            "jellyfish dump -t -c -o %s.23.dat %s" % (prefix, jf2_file),
        ]

    runner.run(commands)

    if sort_dat_file:
        commands = [
            "sort -k2nr %s.23.dat > %s.23.sdat" % (prefix, prefix),
        ]
        runner.run(commands)
    
    
    commands = [
        "jellyfish histo -o %s.23.histo %s" % (prefix, jf2_file),
        "cut -f1 %s.23.dat > %s.23.kmers" % (prefix, prefix),
        "~/Dropbox/workspace/Stirka/emphf/compute_mphf_seq.exe %s.23.kmers %s.23.pf" % (prefix, prefix),
        "~/aindex/bin/compute_index.exe %s.23.dat %s.23.pf %s.23 %s 0" % (prefix, prefix, prefix, threads),
    ]
    runner.run(commands)

    if not only_index:
        commands = [
            "~/aindex/bin/compute_aindex.exe %s.reads %s.23.pf %s.23 %s.23 %s 23 %s.23.tf.bin" % (prefix, prefix, prefix, prefix, threads, prefix),   
        ]
        runner.run(commands)

    commands = [
        "rm %s.23.dat %s.23.kmers" % (prefix, prefix),
    ]

        

