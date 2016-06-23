#!/bin/bash
#test indexing
echo
echo "Test1 -- Generate index files with a E.Coli reference file"
echo "Command=./kart index -p test/EcoliIdx test/ecoli.fa"
echo
./kart index -p test/EcoliIdx test/ecoli.fa

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
echo "Command=./kart aln -t 4 -i test/EcoliIdx -q test/r1.fq -q2 test/r2.fq -o test/alignment.sam"
echo
./kart aln -i test/EcoliIdx -f test/r1.fq -f2 test/r2.fq -o test/alignment.sam

echo
echo "[End of test]"
