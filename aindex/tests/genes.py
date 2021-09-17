

import sys
sys.path.append("/mnt/guatemala/akomissarov/Boechera_spatifolia/aindex")

from aindex import *
from collections import Counter
from trseeker.seqio.fasta_file import sc_iter_fasta


if __name__ == '__main__':
    

    settings = {
      "index_prefix": "/mnt/guatemala/akomissarov/Boechera_spatifolia/raw.23.L3",
      "aindex_prefix": "/mnt/guatemala/akomissarov/Boechera_spatifolia/raw.23.L3",
      "reads_file": "/mnt/guatemala/akomissarov/Boechera_spatifolia/raw.reads",
      "gene_fasta": "/mnt/guatemala/akomissarov/Boechera_spatifolia/apr1.fa",
    }

    k = 23
    index = load_aindex(settings)

    
    used_reads = set()
    results = []
    for seq_obj in sc_iter_fasta(settings["gene_fasta"]):

        for i in xrange(seq_obj.length-k+1):
            kmer = seq_obj.sequence[i:i+k]
            tf = index[kmer]
            if not tf:
                continue

            print i, kmer, tf

            hits = []
            for data in get_reads_se_by_kmer(kmer, index, used_reads):
                start, next_read_start, subread, pos, spring_pos, was_reversed, poses_in_read = data
                used_reads.add((start, spring_pos))
                hits.append([pos, 0, subread, poses_in_read, was_reversed])
            if not hits:
                continue
            max_pos = max([x[0] for x in hits])
            for hid, (pos, nnn, subread, poses_in_read, was_reversed) in enumerate(hits):
                if i-pos < 0:
                    hits[hid][1] = 0
                    hits[hid][2] = subread[abs(i-pos):]
                else:
                    hits[hid][1] = i-pos

                results.append(hits[hid])
    results.sort(key=lambda x: x[1]) 

    for i, (pos, nnn, subread, poses_in_read, was_reversed) in enumerate(results):
        results[i].append("N"*nnn+subread)

    for i in xrange(seq_obj.length):
        nucleotides = [x[-1][i] for x in results if i<len(x[-1]) and x[-1][i] != 'N']
        c = Counter(nucleotides)
        variants = [x for x in c.most_common() if x[1] > 1]
        print i, set(nucleotides), variants, seq_obj.sequence[i]

        if len(variants) > 1 or (len(variants) == 1 and variants[0][0] != seq_obj.sequence[i]):
            raw_input("?")



    with open("/home/akomissarov/Dropbox/PySatDNA/temp.layout", "w") as fh:
        fh.write(seq_obj.sequence)
        fh.write("\n")
        for pos, nnn, subread, poses_in_read, was_reversed in results:
                fh.write("%s\n" % ("N"*nnn+subread))



