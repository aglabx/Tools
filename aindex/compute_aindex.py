#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 01.01.2017
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

import argparse
from PyExp import runner
import os, sys


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Compute index.')
    parser.add_argument('-i', help='Fasta, comma separated fastqs or reads', required=True)
    parser.add_argument('-j', help='JF2 file if exists (None)', required=False, default=None)
    parser.add_argument('-t', help='Reads type (reads|fasta|fastq)', required=False, default="reads")
    parser.add_argument('--aindex', help='Build aindex after index', required=False, default=False)
    parser.add_argument('-o', help='Output prefix', required=True)
    parser.add_argument('-s', help='Sort dat file (None)', required=False, default=None)
    parser.add_argument('--interactive', help='Interactive (False)', required=False, default=None)
    parser.add_argument('-P', help='Threads (12)', required=False, default=12)
    parser.add_argument('-M', help='JF2 memory in Gb (5)', required=False, default=5)

    args = vars(parser.parse_args())

    ### Checking installed tools:
    ####### 1) jellyfish
    ####### 2) compute_mphf_seq.exe
    ####### 3) compute_index.exe
    ####### 4) compute_aindex.exe
    ####### 5) convert_reads_to_fasta.exe
    ####### 5) convert_fasta_q_to_reads.exe


    reads_file = args["i"]
    prefix = args["o"]
    reads_type = args["t"].lower()
    build_aindex = bool(args["aindex"])

    if not reads_type in ["reads","fastq","fasta"]:
        print("Reads type not reds or fastq or fasta")
        sys.exit(1)

    threads = args["P"]
    jf2_file = args["j"]
    sort_dat_file = args["s"]
    memory = args["M"]
    interactive = args["interactive"]

    if not jf2_file:
        if reads_type == "reads":
            print("Converting reads to fasta...")
            commands = [
                "convert_reads_to_fasta.exe %s %s fasta" % (reads_file, prefix),
            ]
            runner.run(commands)
            print("Computing jf2 file from fasta or fastq...")
            commands = [
                "jellyfish count -m 23 -t %s -s %sG -C -o %s.23.jf2 %s.fa" % (threads, memory, prefix, prefix),
            ]
            runner.run(commands)
            if interactive:
                raw_input("Continue?")
        elif reads_type == "fasta" or reads_type == "fastq":
            print("Computing jf2 file from fasta or fastq...")
            commands = [
                "jellyfish count -m 23 -t %s -s %sG -C -o %s.23.jf2 %s" % (threads, memory, prefix, reads_file.replace(",", " ")),
            ]
            runner.run(commands)
            if interactive:
                raw_input("Continue?")
    if build_aindex and reads_type == "fasta":
        print("Computing fasta to reads...")
        commands = [
            "convert_fasta_q_to_reads.exe %s fasta %s.reads" % (reads_file, prefix),   
        ]                
        runner.run(commands)
    if build_aindex and reads_type == "fastq":
        print("Computing fastq to reads...")
        commands = [
            "convert_fasta_q_to_reads.exe %s fastq %s.reads" % (reads_file.replace(",", " "), prefix),   
        ]
        runner.run(commands)

    commands = [
        "jellyfish dump -t -c -o %s.23.dat %s.23.jf2" % (prefix, prefix),
        "jellyfish histo -o %s.23.histo %s.23.jf2" % (prefix, prefix),
        "cut -f1 %s.23.dat > %s.23.kmers" % (prefix, prefix),
        "compute_mphf_seq.exe %s.23.kmers %s.23.pf" % (prefix, prefix),
        "compute_index.exe %s.23.dat %s.23.pf %s.23 %s" % (prefix, prefix, prefix, threads),
    ]
    steps_info = [
        "Compute dat file from jf2 database...",
        "Compute histo file from jf2 database...",
        "Compute kmers from dat file...",
        "Compute mphf from kmers...",
        "Compute build index...",
        
    ]
    for i, command in commands:
        print(steps_info[i])
        runner.run([commands[i]])


    if build_aindex:
        commands = [
            "compute_aindex.exe %s.reads %s.23.pf %s.23 %s.aindex %s 23 %s.23.tf.bin" % (prefix, prefix, prefix, prefix, threads, prefix),   
        ]
        print("Compute build aindex...")
        runner.run(commands)

    commands = [
        "rm %s.23.dat %s.23.kmers" % (prefix, prefix),
    ]
    print("Remove old files build aindex...")
    runner.run(commands)

    if sort_dat_file:
        print("Sort dat file...")
        commands = [
            "sort -k2nr %s.23.dat > %s.23.sdat" % (prefix, prefix),
        ]
        runner.run(commands)
        

