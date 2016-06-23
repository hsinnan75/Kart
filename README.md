Kart: A divide-and-conquer algorithm for NGS read mapping with high error tolerance

Developers: Dr. Hsin-Nan Lin and Dr. Wen-Lian Hsu Institute of Information Science, Academia Sinica, Taiwan.

Introduction
Kart is a new aligner for NGS read alignment developed by Dr. Hsin-Nan Lin and Dr. Wen-Lian Hsu. It is designed for both short and long reads and allows gapped alignment with high sensitivity and accuracy. Kart supports single-end and paired-end reads and multi-thread alignments. We describe the installation and Running instructions of Kart below.

Instructions
1.Installation

We provide the executable file, please type './kart' to run the program. Or you can type 'make' to build the executable file.

2.Usage

For indexing a reference genome, Kart requires the target genome file (in fasta format) and the prefix of the index files (including the directory path).

Ex. #./kart index -p ecoli Ecoli.fa

The above command is to index the genome file Ecoli.fa and store the index files begining with ecoli.

For mapping short reads, Kart requires the the index files of the reference genome and at least one read file (two read files for the separated paired-end reads). Users should use -i to specify the prefix of the index files (including the directory path).

Ex. #./kart aln -i ecoli -f ReadFile1.fa -f2 ReadFile2.fa -o out.sam

The above command is to run Kart to align the paired-end reads in ReadFile1.fa and ReadFile2.fa with index files of ecoli. The output is redirected to out.sam.

We also provide a script (run_test.sh) to test the software. It will index a reference genome and align 2000 paired-end reads.

3.File formats

a.Reference genome files

All reference genome files should be in FASTA format.

b.Read files

All reads files should be in FASTA or FASTQ format. Read sequences should be capital letters. The quality scores in FASTQ are not considered in the alignments. The alignment result will not be different in either format.

If paired-end reads are separated into two files, use -f and -f2 to indicate the two filenames. The i-th reads in the two files are paired. If paired-end reads are in the same file, use -p. The first and second reads are paired, the third and fourth reads are paired, and so on. For the latter case, use -p to indicate the input file contains paired-end reads.

c.Output file

Output is in standard SAM format. For reads aligned with reverse strand of reference genome, they are converted into obverse strand. More detailed information about SAM format, please refer to the SAMtools documents.

Parameter setting

-t INT number of threads [16]

-i STR index prefix [BWT based (BWA), required]

-f STR read filename [required, fasta or fastq]

-f2 STR read filename2 [optional, fasta or fastq], f and f2 are files with paired reads

-p the input read file consists of interleaved paired-end sequences

-o STR alignment output

-m output multiple alignments
