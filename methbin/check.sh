#!/bin/bash

datadir=/uru/Data/Nanopore/projects/methbin

if [ $1 == align_check ] ; then

    for i in $datadir/align/* ;
    do
	prefix=`echo $i | rev | cut -d / -f 1 | rev`
	echo $i/${prefix}_40k.sorted.bam >> aligncheck.txt
	samtools flagstat $i/${prefix}_40k.sorted.bam >> aligncheck.txt
    done
fi

if [ $1 == align_check_taiyaki ] ; then

    for i in $datadir/read_ref/* ;
    do
	prefix=`echo $i | rev | cut -d / -f 1 | rev`
	echo $i/${prefix}_100k.bam >> aligncheck_taiyaki.txt
	samtools flagstat $i/${prefix}_100k.bam >> aligncheck_taiyaki.txt
    done
fi

	
