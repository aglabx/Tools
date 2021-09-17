


import aindex
import time

if __name__ == "__main__":

    settings = {
        "index_prefix": "tests/kmers.23",
        "aindex_prefix": "tests/kmers.23",
        "reads_file": "tests/reads.reads",
        "max_tf": 1000000,
    }

    # kmer2tf = aindex.AIndex(settings["index_prefix"])

    kmer2tf = aindex.load_aindex(settings, skip_aindex=False, skip_reads=False)

    time.sleep(10)

    print("P1DONE")

    s = "TAAGTTATTATTTAGTTAATACTTTTAACAATATTATTAAGGTATTTAAAAAATACTATTATAGTATTTAACATAGTTAAATACCTTCCTTAATACTGTTA"
    print(s)
    for i in range(len(s)-23+1):
        kmer = s[i:i+23]
        print(i, kmer, kmer2tf[kmer])


