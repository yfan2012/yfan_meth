#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/meth
gDNA="171020_neb11 171005_neb12 171012_neb13 170906_neb14 171020_neb15 171003_neb16 171019_neb17 171019_neb19 180628_neb_dcm"
ref=$datadir/refs/er2796.fasta

if [ $1 == mummer ] ; then
    ##gotta mummer first
    for i in $gDNA ;
    do
	##check to make sure the assembly is there
	if [ -f $datadir/$i/canu_assembly/$i.contigs.fasta ] ; then
	    echo mummer for $i
	    mkdir -p $datadir/$i/mummer
	    nucmer -p $datadir/$i/mummer/$i $ref $datadir/$i/canu_assembly/$i.contigs.fasta
	    mummerplot --png $datadir/$i/mummer/$i.layout $datadir/$i/mummer/$i.delta -R $ref -Q $datadir/$i/canu_assembly/$i.contigs.fasta
	    dnadiff -p $datadir/$i/mummer/$i -d $datadir/$i/mummer/$i.delta
	fi
    done
fi

	
	    
	    
