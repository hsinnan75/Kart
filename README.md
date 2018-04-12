Kart: A divide-and-conquer algorithm for NGS read mapping with high error tolerance
===================

Developers: Dr. Hsin-Nan Lin and Dr. Wen-Lian Hsu Institute of Information Science, Academia Sinica, Taiwan.

# Introduction

With the success of NGS technology, more and more biologists use NGS data to analyze genomic variations in a population. This technology can produce sequence reads on the order of million/billion base-pairs in a single day. Therefore, NGS applications require very fast alignments of produced reads to a reference genome.

Kart is an ultra-efficient NGS read aligner. It can rocess long reads as fast as short reads. More importantly, Kart can tolerate much higher error rates, and is suitable for long reads produced, e.g., by the PacBio system.

Kart adopts a divide-and-conquer strategy, which separates a read into regions that are easy to align and regions that require gapped alignment, and align each region independently to compose the final alignment. In our experiments, the average size of fragments requiring gapped alignment is around 20 regardless of the original read length. The experiments on synthetic datasets show that Kart spends much less time on longer reads (150bp to 7000bp) than most aligners do, and still produces reliable alignments when the error rate is as high as 15%. The experiments on real datasets fur-ther demonstrate that Kart can handle reads with poor sequencing quality.

# Download

Please use the command 
  ```
  $ git clone https://github.com/hsinnan75/Kart.git
  ```
to download the package of Kart.

# Compiling

To compile kart and the index tool, please change to kart's folder and just type 'make' to compile kart and bwt_index. If the compilation or the program fails, please contact me (arith@iis.sinica.edu.tw), Thanks.

# Changes
version 2.4.4: Fix a bug on resuce mode

version 2.4.3: Fix a bug on reading paired-end inputs.

version 2.4.2: Modify the MAPQ measure for the PacBio read mapping. 

version 2.4.1: Fix a bug on reporting alignment coordinate with soft clipping.

version 2.4.0: Use ksw algorithm for PacBio long read alignment.

version 2.3.9: Fix a bug in SAM output.

version 2.3.8: Fix the alignment of segment pairs with poor sequence identity.

version 2.3.7: Add an update command.

version 2.3.6: Fix the alignment of segment pairs with multiple mismatches.

version 2.3.5: Fix the alignment when DNA sequences are shown in lower case.

version 2.3.4: Fix the alignment of poor-quality sequences at both ends.

version 2.3.3: Remove edlib and ksw algorithms.

version 2.3.2: Allow multiple read files as the input.

version 2.3.1: Fix the buf when read number exceeds 2^32.

version 2.3.0: Add ksw2 and edlib alignment method to replace the Needleman-Wunsch algorithm.

version 2.2.6: Fix a bug in the alignment report.

version 2.2.5: Fix a bug in the alignment report.

# Get updates

We update Kart from time to time, please check if new version is available by using the following commands.

with Kart version 2.3.7 up
  ```
  $ ./kart update 
  ```
or
  ```
  $ git fetch
  $ git merge origin/master master
  $ make
  ```

# Instructions

We provide the executable file, please type 

  ```
  $ ./kart [options]
  ```
to run the program. Or you can type 'make' to build the executable file.

# Usage

To index a reference genome, Kart requires the target genome file (in fasta format) and the prefix of the index files (including the directory path).

  ```
  $ ./bwt_index ref_file[ex.ecoli.fa] index_prefix[ex. Ecoli]
  ```
The above command is to index the genome file Ecoli.fa and store the index files begining with ecoli.

Please note that if you find bwt_index does not work in your computer system, you may use bwa (http://bio-bwa.sourceforge.net/) to build the index files.
  ```
  $ ./bwa index -p index_prefix xxxx.fa
  ```

To map short reads, Kart requires the the index files of the reference genome and at least one read file (two read files for the separated paired-end reads). Users should use -i to specify the prefix of the index files (including the directory path).

 case 1: standard sam output
  ```
 $ ./kart -i ecoli -f ReadFile1.fa -f2 ReadFile2.fa -o out.sam
  ```

 case 2: multiple input 
  ```
 $ ./kart -i ecoli -f ReadFileA_1.fq ReadFileB_1.fq ReadFileC_1.fq -f2 ReadFileA_2.fq ReadFileB_2.fq ReadFileC_2.fq -o out.sam
  ```

 case 3: gzip compressed output
  ```
 $ ./kart -i ecoli -f ReadFile1.fa -f2 ReadFile2.fa -o out.sam.gz
  ```

 case 4: bam output
  ```
 $ ./kart -i ecoli -f ReadFile1.fa -f2 ReadFile2.fa | samtools view -bo out.bam
  ```

The above commands are to run Kart to align the paired-end reads in ReadFile1.fa and ReadFile2.fa with index files of ecoli.

We also provide a script (run_test.sh) to test the software. It will index a reference genome and align 2000 paired-end reads.

# File formats

- Reference genome files

    All reference genome files should be in FASTA format.

- Read files

    All reads files should be in FASTA or FASTQ format. FASTQ files can be compressed with gzip format. We do not support FASTA files with gzip compression.
    Read sequences should be capital letters. The quality scores in FASTQ are not considered in the alignments. The alignment result will not be different in either format.

    If paired-end reads are separated into two files, use -f and -f2 to indicate the two filenames. The i-th reads in the two files are paired. If paired-end reads are in the same file, use -p. The first and second reads are paired, the third and fourth reads are paired, and so on. For the latter case, use -p to indicate the input file contains paired-end reads.

- Output file

    Output is in standard SAM format. For reads aligned with reverse strand of reference genome, they are converted into obverse strand. More detailed information about SAM format, please refer to the SAMtools documents.
    Kart can also generate compressed output with gzip algorithm. 

# Parameter setting

 ```
-t INT number of threads [16]

-i STR index prefix [BWT based (BWA), required]

-f STR read filename [required, fasta or fastq or fq.gz]

-f2 STR read filename2 [optional, fasta or fastq or fq.gz], f and f2 are files with paired reads

-p the input read file consists of interleaved paired-end sequences

-o STR alignment output

-m output multiple alignments

  ```
