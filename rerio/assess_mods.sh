#!/bin/bash

datadir=/mithril/Data/Nanopore/projects/methbin
ref=$datadir/reference/allsamps.fa

if [ $1 == fix ] ; then

    mkdir -p $datadir/rerio/$i/assess
    for i in neb14 neb15 nebdcm;
    do
	python ~/Code/yfan_meth/utils/methcall_check.py \
	       -r $ref \
	       -b $datadir/align/$i/${i}_rerio.md.sorted.bam \
	       -f $datadir/rerio/$i/called/workspace \
	       -o $datadir/rerio/$i/assess/$i.GCNGC.1.csv \
	       -m GCNGC \
	       -p 1 \
	       -t 54
    done
fi


if [ $1 == Amod ] ; then
    ##for both 6mA samples, inspect GANTC and GATC
    for i in neb19 neb17;
    ##for i in neb17 ;
    do
	mkdir -p $datadir/rerio/$i/assess
	##can loop because mod is in the same location for both motifs (position 1)
	for motif in GANTC GATC CCWGG GCNGC;
	do
            python ~/Code/yfan_meth/utils/methcall_check.py \
		   -r $ref \
		   -b $datadir/align/$i/${i}_rerio.md.sorted.bam \
		   -f $datadir/rerio/$i/called/workspace \
		   -o $datadir/rerio/$i/assess/$i.$motif.1.csv \
		   -m $motif \
		   -p 1 \
		   -t 20
	done
        python ~/Code/yfan_meth/utils/methcall_check.py \
	       -r $ref \
	       -b $datadir/align/$i/${i}_rerio.md.sorted.bam \
	       -f $datadir/rerio/$i/called/workspace \
	       -o $datadir/rerio/$i/assess/$i.GATC.3.csv \
	       -m $motif \
	       -p 3 \
	       -t 20
    done
fi

if [ $1 == Cmod ] ; then
    ##look at CCWGG compared to GATC
    for i in neb14 nebdcm neb15 ;
    do
	mkdir -p $datadir/rerio/$i/assess
        python ~/Code/yfan_meth/utils/methcall_check.py \
	       -r $ref \
	       -b $datadir/align/$i/${i}_rerio.md.sorted.bam \
	       -f $datadir/rerio/$i/called/workspace \
	       -o $datadir/rerio/$i/assess/$i.GATC.3.csv \
	       -m GATC \
	       -p 3 \
	       -t 54
        python ~/Code/yfan_meth/utils/methcall_check.py \
	       -r $ref \
	       -b $datadir/align/$i/${i}_rerio.md.sorted.bam \
	       -f $datadir/rerio/$i/called/workspace \
	       -o $datadir/rerio/$i/assess/$i.GCNGC.1.csv \
	       -m GCNGC \
	       -p 3 \
	       -t 54
	python ~/Code/yfan_meth/utils/methcall_check.py \
	       -r $ref \
	       -b $datadir/align/$i/${i}_rerio.md.sorted.bam \
	       -f $datadir/rerio/$i/called/workspace \
	       -o $datadir/rerio/$i/assess/$i.CCWGG.1.csv \
	       -m CCWGG \
	       -p 1 \
	       -t 54
    done
fi

if [ $1 == control ] ; then
    i=neb11
    mkdir -p $datadir/rerio/$i/assess
    ##motif postion 1
    for motif in GANTC GATC CCWGG GCNGC;
    do
        python ~/Code/yfan_meth/utils/methcall_check.py \
	       -r $ref \
	       -b $datadir/align/$i/${i}_rerio.md.sorted.bam \
	       -f $datadir/rerio/$i/called/workspace \
	       -o $datadir/rerio/$i/assess/$i.$motif.1.csv \
	       -m $motif \
	       -p 1 \
	       -t 54
    done

    ##motif position 3
    for motif in GATC ;
    do
	python ~/Code/yfan_meth/utils/methcall_check.py \
	       -r $ref \
	       -b $datadir/align/$i/${i}_rerio.md.sorted.bam \
	       -f $datadir/rerio/$i/called/workspace \
	       -o $datadir/rerio/$i/assess/$i.$motif.3.csv \
	       -m $motif \
	       -p 3 \
	       -t 54
    done
fi
