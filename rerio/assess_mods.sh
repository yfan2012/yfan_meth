#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin
ref=$datadir/reference/allsamps.fa

if [ $1 == Amod ] ; then
    ##for both 6mA samples, inspect GANTC and GATC
    for i in neb19 neb17;
    ##for i in neb17 ;
    do
	mkdir -p $datadir/rerio/$i/assess
	##can loop because mod is in the same location for both motifs (position 1)
	for motif in GANTC GATC ;
	##for motif in GANTC ;
	do
            python ~/Code/methylation/utils/methcall_check.py \
		   -r $ref \
		   -b $datadir/align/$i/${i}_rerio.md.sorted.bam \
		   -f $datadir/rerio/$i/called/workspace \
		   -o $datadir/rerio/$i/assess/$i.$motif.csv \
		   -m $motif \
		   -p 1 \
		   -s 0 \
		   -t 54
	done
    done
fi

if [ $1 == Cmod ] ; then
    ##look at CCWGG compared to GATC
    for i in neb14 nebdcm ;
    do
	mkdir -p $datadir/rerio/$i/assess
        python ~/Code/methylation/utils/methcall_check.py \
	       -r $ref \
	       -b $datadir/align/$i/${i}_rerio.md.sorted.bam \
	       -f $datadir/rerio/$i/called/workspace \
	       -o $datadir/rerio/$i/assess/$i.GATC.csv \
	       -m GATC \
	       -p 3 \
	       -s 0 \
	       -t 54
	
        python ~/Code/methylation/utils/methcall_check.py \
	       -r $ref \
	       -b $datadir/align/$i/${i}_rerio.md.sorted.bam \
	       -f $datadir/rerio/$i/called/workspace \
	       -o $datadir/rerio/$i/assess/$i.CCWGG.csv \
	       -m CCWGG \
	       -p 1 \
	       -s 0 \
	       -t 54
    done
fi
