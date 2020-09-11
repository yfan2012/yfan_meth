#!/bin/bash

rawdir=/dilithium/Data/Nanopore/oxford/181010_ecoli_methylation_standards/all_data
datadir=/mithril/Data/Nanopore/projects/methbin
allsamps='ecoli_CpG ecoli_CpGGpC ecoli_GpC ecoli_Unmethylated'

if [ $1 == untar ] ; then
    for i in $allsamps ;
    do
	##these tar files of isac are already in multi fast5 format
	tar -xzf $rawdir/$i.fast5.tgz -C $datadir/multiraw
    done
fi

if [ $1 == multiraw_sub ] ; then
    for i in $allsamps ;
    do
	mkdir -p ~/data/rerio/$i/multiraw_sub
	for j in {1..50} ;
	do
	    cp $datadir/multiraw/$i/${i}_$j.fast5 ~/data/rerio/$i/multiraw_sub/
	done
    done
fi
    
