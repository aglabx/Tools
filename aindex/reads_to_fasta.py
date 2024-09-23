#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 10.10.2013
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reads to fasta.")
    parser.add_argument("-i", "--input", help="Reads input", required=True)
    parser.add_argument("-o", "--output", help="Fasta output", required=True)
    args = vars(parser.parse_args())
    reads_file = args["input"]
    output = args["output"]

    assert reads_file != output

    with open(output, "w") as fw:
        with open(reads_file) as fh:
            for i, read in enumerate(fh):
                fw.write(">%s\n%s" % (i, read))
