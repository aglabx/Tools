#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 01.01.2017
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
import subprocess
import os


def runner(commands):
    for command in commands:
        print(command)
        subprocess.run(command, shell=True)


"""
~/Dropbox/Ariadna/Stirka/bin/R2F.exe pe500.mtDNA.reads pe500.mtDNA fasta
jellyfish count -m 23 -C -t 60 -s 20G -o pe500.mtDNA.23.jf2 pe500.mtDNA.fa 
jellyfish dump -c -t -o pe500.mtDNA.23.dat pe500.mtDNA.23.jf2
cut -f1 pe500.mtDNA.23.dat > pe500.mtDNA.23.kmers
~/Dropbox/Stirka/emphf/compute_mphf_seq.exe pe500.mtDNA.23.kmers pe500.mtDNA.23.pf
~/Dropbox/Ariadna/Stirka/bin/LU_index.exe pe500.mtDNA.23.dat pe500.mtDNA.23.pf pe500.mtDNA.23
~/Dropbox/Ariadna/Stirka/bin/LU_aindex.exe pe500.mtDNA.reads pe500.mtDNA.23.pf pe500.mtDNA.23 pe500.mtDNA.aindex 60 23 pe500.mtDNA.23.tf.bin 
"""

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compute index from jf2 file.")
    parser.add_argument("-j", help="JF2 file if exists (None)", required=True)
    parser.add_argument("-o", help="Output prefix", required=True)
    parser.add_argument("--lu", help="-L for jellyfish [0]", required=False, default=0)
    parser.add_argument("-P", help="Threads (12)", required=False, default=60)
    parser.add_argument("-M", help="JF2 memory in Gb (5)", required=False, default=1)
    parser.add_argument(
        "--path_to_aindex",
        help="Path to aindex folder inside Tools [~/Dropbox/workspace/PyBioSnippets/Tools/aindex]",
        required=False,
        default="~/Dropbox/workspace/PyBioSnippets/Tools/aindex",
    )

    args = vars(parser.parse_args())

    prefix = args["o"]
    threads = args["P"]
    jf2_file = args["j"]
    memory = args["M"]
    lu = args["lu"]
    path_to_aindex = args["path_to_aindex"]

    ### here we expect that jf2 file is created
    if not os.path.exists(jf2_file):
        print("JF2 file is missing, please, check jellyfish command")
        exit(1)

    commands = [
        f"jellyfish dump -t -c -L {lu} -o {prefix}.23.dat {jf2_file}",
        f"jellyfish histo -o {prefix}.23.histo {jf2_file}",
    ]

    runner.run(commands)

    commands = [
        f"cut -f1 {prefix}.23.dat > {prefix}.23.kmers",
        f"{path_to_aindex}/compute_mphf_seq.exe {prefix}.23.kmers {prefix}.23.pf",
        f"{path_to_aindex}/compute_index.exe {prefix}.23.dat {prefix}.23.pf {prefix}.23 {threads} 0",
        f"rm {prefix}.23.dat {prefix}.23.jf2",
    ]
    runner.run(commands)
