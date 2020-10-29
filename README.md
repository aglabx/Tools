# Trimming raw Illumina reads with V2_trim.exe



```shell script

./V2_trim.exe <reads_file or fastq_prefix like SRRXXXX> <output_prefix> <P> <bool_save_fasta> <format: reads|fasta|fastq|fastq_se|fastq_raw> <kmer_file>

```

# Removing optical duplicates with rmdup.exe



```shell script

./rmdup.exe <fastq_prefix like SRRXXXX> <output_prefix> <distance> <hd_cutoff 1-2 per 100 bp, for PE 5 the best> <shift=0,30>

```

If you have dataset_1.fastq and dataset_2.fastq then prefix will be <em>dataset</em>.

```shell script

./V2_trim.exe dataset dataset 32 0 fastq illumina_ext.dat

./rmdup.exe dataset.trim dataset.rmdup 3000 5 0,30

```
