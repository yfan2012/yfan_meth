#!/bin/bash

datadir=/uru/Data/Nanopore/projects/mdr
samps="MDRstool_16 MDRstool_19"


if [ $1 == call ] ; then
datadir=~/data/mdr
    for i in $samps ;
    do
	guppy_basecaller \
	    -i $datadir/$i/multiraw \
	    -s $datadir/$i/called \
	    --flowcell FLO-MIN106 --kit SQK-LSK108 \
	    --device 'cuda:0' \
	    --fast5_out
    done
fi

if [ $1 == gather ] ; then
    for i in $samps ;
    do
	mkdir -p $datadir/$i/fastqs
	cat $datadir/$i/called/*fastq > $datadir/$i/fastqs/$i.fq
    done
fi

	     
if [ $1 == assemble ] ; then
    sizes="150m 100m 500m 10g"
    for i in $samps ;
    do
	for gsize in $sizes ;
	do
	    mkdir -p $datadir/$i/metaflye/$gsize
	    flye \
		--nano-raw $datadir/$i/fastqs/$i.fq \
		-o $datadir/$i/metaflye/$gsize \
		-t 36 \
		-g $gsize \
		--plasmids \
		--meta
	done
    done
fi

if [ $1 == recswab_assemble ] ; then

    gsize=1m
    mkdir -p $datadir/RECswab_1/metaflye/$gsize
    flye \
	--nano-raw $datadir/RECswab_1/fastqs/RECswab_1.fq \
	-o $datadir/RECswab_1/metaflye/$gsize \
	-t 36 \
	-g $gsize \
	--plasmids \
	--meta
fi

	
