#!/bin/bash

datasource=/uru/Data/Nanopore/projects/methbin/multiraw
datadir=~/data/methbin
gdna="neb11 neb12 neb13 neb14 neb15 neb16 neb17 neb19 nebdcm"
ref=$datadir/reference/allsamps.fa

if [ $1 == get_test ] ; then
    mkdir -p $datadir/multiraw_test
    for i in $gdna ;
    do
	mkdir -p $datadir/multiraw_test/$i
	for j in {51..100} ;
	do
	    if [ -f $datasource/$i/${i}_$j.fast5 ] ; then
		cp $datasource/$i/${i}_$j.fast5 $datadir/multiraw_test/$i/
	    fi
	done
    done
fi

if [ $1 == 4mc_call ] ; then
    mkdir -p $datadir/calls
    mkdir -p $datadir/calls/neb12

    model=$datadir/train/neb12/training2/model_final.checkpoint
        
    for i in neb12 neb2 neb11 ;
    do
	mkdir -p $datadir/calls/neb12/$i
	megalodon $datadir/multiraw_test/$i \
		  --taiyaki-model-filename $model \
		  --reference $ref \
		  --devices 0 \
		  --outputs mod_basecalls mods \
		  --mod-motif X CCGG 0 \
		  --write-mods-text \
		  --processes 18 \
		  --overwrite \
		  --verbose-read-progress 3 \
		  --output-directory $datadir/calls/neb12/$i
    done
fi

if [ $1 == 6ma_call ] ; then
    for i in neb19 neb9 neb11 ;
    do
	mkdir -p $datadir/calls/neb12/$i
	megalodon $datadir/multiraw_test/$i \
		  --taiyaki-model-filename $model \
		  --reference $ref \
		  --devices 0 \
		  --outputs mod_basecalls mods \
		  --mod-motif X CCGG 0 \
		  --write-mods-text \
		  --processes 18 \
		  --overwrite \
		  --verbose-read-progress 3 \
		  --output-directory $datadir/calls/neb12/$i
    done


    for i in neb17 neb6 neb11 ;
    do
	mkdir -p $datadir/calls/neb12/$i
	megalodon $datadir/multiraw_test/$i \
		  --taiyaki-model-filename $model \
		  --reference $ref \
		  --devices 0 \
		  --outputs mod_basecalls mods \
		  --mod-motif X CCGG 0 \
		  --write-mods-text \
		  --processes 18 \
		  --overwrite \
		  --verbose-read-progress 3 \
		  --output-directory $datadir/calls/neb12/$i
    done
fi
