# aindex
Perfect hash based Index for text data

## Installation

```
cd external
git clone https://github.com/ot/emphf.git
cd emphf
cmake .
make
```

```

mkdir bin

g++ Compute_index.cpp  -std=c++11 -pthread -O3 debrujin.hpp debrujin.cpp read.hpp read.cpp kmers.hpp kmers.cpp settings.cpp settings.hpp hash.hpp hash.cpp && mv a.out bin/compute_index.exe

g++ Compute_aindex.cpp  -std=c++11 -pthread -O3 debrujin.hpp debrujin.cpp read.hpp read.cpp kmers.hpp kmers.cpp settings.cpp settings.hpp hash.hpp hash.cpp && mv a.out bin/compute_aindex.exe

g++ Compute_reads.cpp  -std=c++11 -pthread -O3 debrujin.hpp debrujin.cpp read.hpp read.cpp kmers.hpp kmers.cpp settings.cpp settings.hpp hash.hpp hash.cpp && mv a.out bin/compute_reads.exe
```

Lunux compilation commadn:

```
g++ -c -std=c++11 -fPIC python_wrapper.cpp -o python_wrapper.o && g++ -c -std=c++11 -fPIC kmers.cpp kmers.hpp debrujin.cpp debrujin.hpp hash.cpp hash.hpp read.cpp read.hpp settings.hpp settings.cpp && g++ -shared -Wl,-soname,python_wrapper.so -o python_wrapper.so python_wrapper.o kmers.o debrujin.o hash.o read.o settings.o
```

Mac compilation command:

```
g++ -c -std=c++11 -fPIC python_wrapper.cpp -o python_wrapper.o && g++ -c -std=c++11 -fPIC kmers.cpp kmers.hpp debrujin.cpp debrujin.hpp hash.cpp hash.hpp read.cpp read.hpp settings.hpp settings.cpp && g++ -shared -Wl,-install_name,python_wrapper.so -o python_wrapper.so python_wrapper.o kmers.o debrujin.o hash.o read.o settings.o
```

## Usage


Compute all binary arrays:

```bash
FASTQ1=raw_reads.101bp.IS350bp25_1.fastq
FASTQ2=raw_reads.101bp.IS350bp25_2.fastq

jellyfish count -m 23 -s 2G -t 4 -C -o kmers.23.jf2 $FASTQ1 $FASTQ2
jellyfish dump -c -t -o kmers.23.dat kmers.23.jf2
cut -f1 kmers.23.dat > kmer.23.kmers
../external/emphf/compute_mphf_seq kmer.23.kmers kmer.23.pf

../bin/compute_index.exe kmers.23.dat kmers.23.pf kmers.23 4 0

../bin/compute_reads.exe raw_reads.101bp.IS350bp25_1.fastq raw_reads.101bp.IS350bp25_2.fastq fastq reads.reads

../bin/compute_aindex.exe reads.reads kmers.23.pf kmers.23 kmers.23 40 23 kmers.23.tf.bin 

```

Usage from Python:

You can simply run **demo.py** or:

```python


from aindex import *

settings = {
  "index_prefix": "tests/kmers.23",
  "aindex_prefix": "tests/kmers.23",
  "reads_file": "tests/reads.reads",
}

index = load_aindex(settings)

k = 23
sequence = "TAAGTTATTATTTAGTTAATACTTTTAACAATATTATTAAGGTATTTAAAAAATACTATTATAGTATTTAACATAGTTAAATACCTTCCTTAATACTGTTAAATTATATTCAATCAATACATATATAATATTATTAAAATACTTGATAAGTATTATTTAGATATTAGACAAATACTAATTTTATATTGCTTTAATACTTAATAAATACTACTTATGTATTAAGTAAATATTACTGTAATACTAATAACAATATTATTACAATATGCTAGAATAATATTGCTAGTATCAATAATTACTAATATAGTATTAGGAAAATACCATAATAATATTTCTACATAATACTAAGTTAATACTATGTGTAGAATAATAAATAATCAGATTAAAAAAATTTTATTTATCTGAAACATATTTAATCAATTGAACTGATTATTTTCAGCAGTAATAATTACATATGTACATAGTACATATGTAAAATATCATTAATTTCTGTTATATATAATAGTATCTATTTTAGAGAGTATTAATTATTACTATAATTAAGCATTTATGCTTAATTATAAGCTTTTTATGAACAAAATTATAGACATTTTAGTTCTTATAATAAATAATAGATATTAAAGAAAATAAAAAAATAGAAATAAATATCATAACCCTTGATAACCCAGAAATTAATACTTAATCAAAAATGAAAATATTAATTAATAAAAGTGAATTGAATAAAATTTTGAAAAAAATGAATAACGTTATTATTTCCAATAACAAAATAAAACCACATCATTCATATTTTTTAATAGAGGCAAAAGAAAAAGAAATAAACTTTTATGCTAACAATGAATACTTTTCTGTCAAATGTAATTTAAATAAAAATATTGATATTCTTGAACAAGGCTCCTTAATTGTTAAAGGAAAAATTTTTAACGATCTTATTAATGGCATAAAAGAAGAGATTATTACTATTCAAGAAAAAGATCAAACACTTTTGGTTAAAACAAAAAAAACAAGTATTAATTTAAACACAATTAATGTGAATGAATTTCCAAGAATAAGGTTTAATGAAAAAAACGATTTAAGTGAATTTAATCAATTCAAAATAAATTATTCACTTTTAGTAAAAGGCATTAAAAAAATTTTTCACTCAGTTTCAAATAATCGTGAAATATCTTCTAAATTTAATGGAGTAAATTTCAATGGATCCAATGGAAAAGAAATATTTTTAGAAGCTTCTGACACTTATAAACTATCTGTTTTTGAGATAAAGCAAGAAACAGAACCATTTGATTTCATTTTGGAGAGTAATTTACTTAGTTTCATTAATTCTTTTAATCCTGAAGAAGATAAATCTATTGTTTTTTATTACAGAAAAGATAATAAAGATAGCTTTAGTACAGAAATGTTGATTTCAATGGATAACTTTATGATTAGTTACACATCGGTTAATGAAAAATTTCCAGAGGTAAACTACTTTTTTGAATTTGAACCTGAAACTAAAATAGTTGTTCAAAAAAATGAATTAAAAGATGCACTTCAAAGAATTCAAACTTTGGCTCAAAATGAAAGAACTTTTTTATGCGATATGCAAATTAACAGTTCTGAATTAAAAATAAGAGCTATTGTTAATAATATCGGAAATTCTCTTGAGGAAATTTCTTGTCTTAAATTTGAAGGTTATAAACTTAATATTTCTTTTAACCCAAGTTCTCTATTAGATCACATAGAGTCTTTTGAATCAAATGAAATAAATTTTGATTTCCAAGGAAATAGTAAGTATTTTTTGATAACCTCTAAAAGTGAACCTGAACTTAAGCAAATATTGGTTCCTTCAAGATAATGAATCTTTACGATCTTTTAGAACTACCAACTACAGCATCAATAAAAGAAATAAAAATTGCTTATAAAAGATTAGCAAAGCGTTATCACCCTGATGTAAATAAATTAGGTTCGCAAACTTTTGTTGAAATTAATAATGCTTATTCAATATTAAGTGATCCTAACCAAAAGGAAAAATATGATTCAATGCTGAAAGTTAATGATTTTCAAAATCGCATCAAAAATTTAGATATTAGTGTTAGATGACATGAAAATTTCATGGAAGAACTCGAACTTCGTAAGAACTGAGAATTTGATTTTTTTTCATCTGATGAAGATTTCTTTTATTCTCCATTTACAAAAA"
test_kmer = "TAAGTTATTATTTAGTTAATACT"
right_kmer = "AGTTAATACTTTTAACAATATTA"


print "Task 1. Get kmer frequency"
raw_input("\nReady?")
for i in xrange(len(sequence)-k+1):
    kmer = sequence[i:i+k]
    print "Position %s kmer %s freq = %s" % (i, kmer, index[kmer])

print "Task 2. Iter read by read, print the first 20 reads"
raw_input("\nReady?")
for i, read in enumerate(index.iter_reads()):
    if i == 20:
        break
    print i, read

print "Task 3. Iter reads by kmer, returs (start, next_read_start, read, pos_if_uniq|None, all_poses)"
raw_input("\nReady?")
for read in iter_reads_by_kmer(test_kmer, index):
    print read


print "Task 4. Get distances in reads for two kmers, returns a list of (rid, left_kmer_pos, right_kmer_pos) tuples."
raw_input("\nReady?")
print get_left_right_distances(test_kmer, right_kmer, index)


print "Task 5. Get layout for kmer, returns (max_pos, reads, lefts, rights, rids, starts), for details see source code"
raw_input("\nReady?")
max_pos, reads, lefts, rights, rids, starts = get_layout_for_kmer(right_kmer, index)
print "Central layout:"
for read in reads:
    print read
print "Left flanks:"
print lefts
print "Right flanks:"
print rights

print "Task 6. Iter reads by sequence, returns (start, next_read_start, read, pos_if_uniq|None, all_poses)"
raw_input("\nReady?")
sequence = "AATATTATTAAGGTATTTAAAAAATACTATTATAGTATTTAACATA"
for read in iter_reads_by_sequence(sequence, index):
    print read


```










