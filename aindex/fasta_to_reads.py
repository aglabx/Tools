#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 10.10.2013
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

import sys
import argparse
# from trseeker.seqio.fasta_file import sc_iter_fasta


def sc_iter_fasta_brute(file_name, inmem=False):
    """ Iter over fasta file."""
    
    header = None
    seq = []
    with open(file_name) as fh:
        if inmem:
            data = fh.readlines()
        else:
            data = fh
        for line in data:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    sequence = "".join(seq)
                    yield (header, sequence)
                header = line
                seq = []
                continue
            seq.append(line)
        if seq or header:
            sequence = "".join(seq)
            yield (header, sequence)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Fasta to reads.')
    parser.add_argument('-i','--input', help='Fasta input', required=True)
    parser.add_argument('-o','--output', help='Reads output', required=True)
    parser.add_argument('-H','--header', help='Headers output ()', required=False)
    args = vars(parser.parse_args())
    fasta_file = args["input"]
    output = args["output"]
    header = args["header"]
    
    assert fasta_file != output

    assert fasta_file != header
    
    pos = 0

    if header is None:
        with open(output, "w") as fh:
            for header, seq in sc_iter_fasta_brute(fasta_file):
                fh.write("%s\n" % seq.upper())
    else:
        with open(output, "w") as fh:
            with open(header, "w") as fh_header:
                for header, seq in sc_iter_fasta_brute(fasta_file):
                    fh_header.write("%s\t%s\n" % (header[1:], pos))
                    fh.write("%s\n" % seq.upper())
                    pos += len(seq)+1
                


