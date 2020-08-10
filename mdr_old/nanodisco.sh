#!/bin/bash

datadir=/uru/Data/Nanopore/projects/mdr/
samps="MDRstool_16 MDRstool_19"


if [ $1 == preprocess ] ; then
    for i in $samps ;
    do
	mkdir -p $datadir/$i/nanodisco/preprocess
	nanodisco preprocess \
		  -p 36 \
		  -f $datadir/$i/called/workspace/ \
		  -s $i \
		  -o $datadir/$i/nanodisco/preprocess \
		  -r $datadir/$i/metaflye/150m/assembly.fasta
    done
fi

	
if [ $1 == coverage ] ; then
    for i in $samps
    do
	mkdir -p $datadir/$i/nanodisco/coverage
	nanodisco coverage \
		  -b $datadir/$i/nanodisco/preprocess/$i.sorted.bam \
		  -r $datadir/$i/metaflye/150m/assembly.fasta \
		  -o $datadir/$i/nanodisco/coverage
    done
fi


if [ $1 == current_diff ] ; then
    for i in $samps
    do
	nanodisco chunk_info -r $datadir/$i/metaflye/150m/assembly.fasta

	nanodisco difference \
		  -nj 16 \
		  -nc 4 \
		  -p 2 \
		  -r $datadir/$i/metaflye/150m/assembly.fasta \
		  -i $datadir/$i/nanodisco/preprocess \
		  
