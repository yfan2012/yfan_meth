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
        
    ##for i in neb12 neb2 neb11 neb1 ;
    for i in neb2 ;
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
    for i in neb19 neb9 neb11 neb1 ;
    do
	model=$datadir/train/neb19/training2/model_final.checkpoint
	mkdir -p $datadir/calls/neb19/$i
	megalodon $datadir/multiraw_test/$i \
		  --taiyaki-model-filename $model \
		  --reference $ref \
		  --devices 0 \
		  --outputs mod_basecalls mods \
		  --mod-motif Z GATC 1 \
		  --write-mods-text \
		  --processes 18 \
		  --overwrite \
		  --verbose-read-progress 3 \
		  --output-directory $datadir/calls/neb19/$i
    done
    
    
    for i in neb17 neb6 neb11 neb1 ;
    do
	model=$datadir/train/neb17/training2/model_final.checkpoint
	mkdir -p $datadir/calls/neb17/$i
	megalodon $datadir/multiraw_test/$i \
		  --taiyaki-model-filename $model \
		  --reference $ref \
		  --devices 0 \
		  --outputs mod_basecalls mods \
		  --mod-motif Z GANTC 1 \
		  --write-mods-text \
		  --processes 18 \
		  --overwrite \
		  --verbose-read-progress 3 \
		  --output-directory $datadir/calls/neb17/$i
    done
fi


if [ $1 == 5mc_call ] ; then
    nebsamps=/uru/Data/Nanopore/projects/methbin/reference/nebsamps.tsv
    
    while read samp ; do
        name=`echo $samp | cut -d ' ' -f 1`
        motif=` echo $samp | cut -d ' ' -f 2 `
        modtype=` echo $samp | cut -d ' ' -f 3 `
        canonical=` echo $samp | cut -d ' ' -f 4 `
        mod=` echo $samp | cut -d ' ' -f 5 `
        pos=` echo $samp | cut -d ' ' -f 6 `
        dna=` echo $samp | cut -d ' ' -f 7 `
        seq=` echo $samp | cut -d ' ' -f 8 `
	pair=` echo $samp | cut -d ' ' -f 9 `
	model=$datadir/train/$name/training2/model_final.checkpoint

	if [ -f $model ] ; then

	    if [ $pair == none ] ; then
		nameset="$name neb11"
	    else
		nameset="$name $pair neb11 neb1"
	    fi
	    
	    for i in $nameset ;
	    do
		mkdir -p $datadir/calls/$name/$i
		megalodon $datadir/multiraw_test/$i \
			  --taiyaki-model-filename $model \
			  --reference $ref \
			  --devices 0 \
			  --outputs mod_basecalls mods \
			  --mod-motif $mod $motif $pos \
			  --write-mods-text \
			  --processes 36 \
			  --overwrite \
			  --verbose-read-progress 3 \
			  --output-directory $datadir/calls/$name/$i
	    done
	fi
    done < $nebsamps
fi
	    

if [ $1 == nebdcm ] ; then
    for i in nebdcm neb11 ;
    do
	model=$datadir/train/nebdcm/training2/model_final.checkpoint
	mkdir -p $datadir/calls/nebdcm/$i
	megalodon $datadir/multiraw_test/$i \
		  --taiyaki-model-filename $model \
		  --reference $ref \
		  --devices 0 \
		  --outputs mod_basecalls mods \
		  --mod-motif Y CCWGG 1 \
		  --write-mods-text \
		  --processes 18 \
		  --overwrite \
		  --verbose-read-progress 3 \
		  --output-directory $datadir/calls/nebdcm/$i
    done
fi


