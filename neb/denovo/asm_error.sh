#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/meth
gDNA="171020_neb11 171005_neb12 171012_neb13 170906_neb14 171020_neb15 171003_neb16 171019_neb17 171019_neb19 180628_neb_dcm"
ref=$datadir/refs/er2796.fasta

##this is forward looking to make it easy to switch between gDNA and plasmids using $2
samps=$gDNA

if [ $1 == mummer ] ; then
    ##gotta mummer first
    for i in $samps ;
    do
	##check to make sure the assembly is there
	if [ -f $datadir/$i/canu_assembly/$i.contigs.fasta ] ; then
	    echo mummer for $i
	    module load gnuplot
	    mkdir -p $datadir/$i/mummer
	    nucmer -p $datadir/$i/mummer/$i $ref $datadir/$i/canu_assembly/$i.contigs.fasta
	    mummerplot --png --fat --filter -p $datadir/$i/mummer/$i.layout $datadir/$i/mummer/$i.delta -R $ref -Q $datadir/$i/canu_assembly/$i.contigs.fasta
	    dnadiff -p $datadir/$i/mummer/$i -d $datadir/$i/mummer/$i.delta
	fi
    done
fi

if [ $1 == diffs ] ; then
    ml python/2.7
    ##take mummer output and start counting which 6mers are most common
    for i in $samps ;
    do
	if [ -f $datadir/$i/mummer/$i.snps ] ; then
	    echo finding enriched motifs for $i
	    mkdir -p $datadir/$i/asm_diffs
	    rm $datadir/$i/asm_diffs/*
	    python ~/Code/utils/meth/motif_enrich.py -s $datadir/$i/mummer/$i.snps -r $ref -m 6 -o $datadir/$i/asm_diffs/$i.ref.6mer.csv
	    python ~/Code/utils/meth/motif_enrich.py -s $datadir/$i/mummer/$i.snps -r $ref -m 5 -o $datadir/$i/asm_diffs/$i.ref.5mer.csv
	    python ~/Code/utils/meth/motif_enrich.py -s $datadir/$i/mummer/$i.snps -r $ref -m 4 -o $datadir/$i/asm_diffs/$i.ref.4mer.csv
	fi
    done
fi


if [ $1 == sortdiffs ] ; then
    for i in $samps ;
    do
	if [ -d $datadir/$i/asm_diffs ] ; then
	    for file in $datadir/$i/asm_diffs/*.csv ;
	    do
		prefix=`basename $file .csv`
		sort -gr --field-separator=',' --key=4 $file > $datadir/$i/asm_diffs/$prefix.sorted.csv
		grep -vE "(AAAA|TTTT|CCCC|GGGG)" $datadir/$i/asm_diffs/$prefix.sorted.csv > $datadir/$i/asm_diffs/$prefix.hpfilt.sorted.csv
	    done
	fi
    done
fi

	
if [ $1 == listdiffs ] ; then
    ml python/2.7
    ##print out a txt file listing seqs of length 16 which are involved in some kind of sindel business
    ##for EM analysis
    for i in $samps ;
    do
	if [ -f $datadir/$i/mummer/$i.snps ] ; then
	    mkdir -p $datadir/asm_diffs
	    python findmotif.py -s $datadir/$i/mummer/$i.snps -r $ref -l 6 -o $datadir/$i/asm_diffs/$i.12merlist.txt
	fi
    done
fi

	    
	    
