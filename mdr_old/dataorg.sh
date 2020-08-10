#!/bin/bash

datadir=/uru/Data/Nanopore/projects/mdr
metaraw="MDRstool_16 MDRstool_19 20180524_2329_RECswab_1"
meta="MDRstool_16 MDRstool_19 RECswab_1"
iso="AB_531 Ecoli_526 KLPN_529 PRRE_530 PSAE_534"

if [ $1 == untar_meta ] ; then
    mkdir -p $datadir/raw/
    for i in $metaraw ;
    do
	tar -xzf $datadir/data/$i.tar.gz -C $datadir/raw/
    done
    mv $datadir/raw/20180524_2329_RECswab_1 $datadir/raw/RECswab_1
fi

if [ $1 == groupfast5 ] ; then
	mkdir -p ~/data/mdr/MDRstool_16/multiraw
	single_to_multi_fast5 \
	    -t 36 \
	    --recursive \
	    -i $datadir/raw/MDRstool_16/GA40000/reads \
	    -s ~/data/mdr/MDRstool_16/multiraw 

	mkdir -p ~/data/mdr/MDRstool_19/multiraw
	single_to_multi_fast5 \
	    -t 36 \
	    --recursive \
	    -i $datadir/raw/MDRstool_19/GA50000/reads \
	    -s ~/data/mdr/MDRstool_19/multiraw

	mkdir -p ~/data/mdr/RECswab_1/multiraw
	single_to_multi_fast5 \
	    -t 36 \
	    --recursive \
	    -i $datadir/raw/RECswab_1/fast5/pass \
	    -s ~/data/mdr/RECswab_1/multiraw
fi
