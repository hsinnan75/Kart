#!/bin/bash
#test indexing
echo
echo "Test1 -- Generate index files with a E.Coli reference file"
echo "Command=bin/bwa_index test/ecoli.fa test/EcoliIdx"
echo
bin/bwa_index test/ecoli.fa test/EcoliIdx

FILESIZE=$(du -sb test/EcoliIdx.bwt | awk '{ print $1 }')
ACT_SIZE=4639752
if [ $FILESIZE == 4639752 ];
then
	echo
	echo "[Making the index files successfully!]"
	echo
else
	echo "[Failed to generate the index files!]"
fi

#test alignment
echo
echo "Test2 -- Align 2000 E.Coli reads with 4 threads"
echo "Command=./kart -i test/EcoliIdx -t 4 -f test/r1.fq -f2 test/r2.fq -o test/alignment.sam"
echo
bin/kart -t 4 -i test/EcoliIdx -f test/r1.fq -f2 test/r2.fq -o test/alignment.sam

echo
echo "[End of test]"
