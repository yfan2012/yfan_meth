#!/bin/bash

datasource=/uru/Data/Nanopore/projects/methbin/multiraw
datadir=/home/yfan/data/methbin
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


if [ $1 == cross_call ] ; then
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

	mkdir -p $datadir/calls/$name

	echo $mod $motif $pos
	
	if [ $name == neb14 ] ; then
	    for i in neb19 neb12 neb13 neb14 neb15 neb16 neb17 ;
	    do
		echo $name $i =======================================================================================
		if [ $i != $name ] ; then
		    rmdir $datadir/calls/$name/$i
		    megalodon $datadir/multiraw_test/$i \
			      --taiyaki-model-filename $model \
			      --reference $ref \
			      --devices 0 \
			      --outputs mod_basecalls mods \
			      --mod-motif Y GATC 3 \
			      --write-mods-text \
			      --processes 36 \
			      --verbose-read-progress 3 \
			      --output-directory $datadir/calls/$name/$i
		fi
	    done
	fi
    done < $nebsamps
fi


if [ $1 == neb14 ] ; then
    mkdir -p $datadir/calls
    mkdir -p $datadir/calls/neb14
    
    model=$datadir/train/neb14/training2/model_final.checkpoint
        
    ##for i in neb14 neb10 neb11 neb1 neb19 ;
    for i in neb10 ;
    do
	mkdir -p $datadir/calls/neb12/$i
	megalodon $datadir/multiraw_test/$i \
		  --taiyaki-model-filename $model \
		  --reference $ref \
		  --devices 0 \
		  --outputs mod_basecalls mods \
		  --mod-motif Y GATC 3 \
		  --write-mods-text \
		  --processes 18 \
		  --overwrite \
		  --verbose-read-progress 3 \
		  --output-directory $datadir/calls/neb14/$i
    done
fi

if [ $1 == neb15 ] ; then
    mkdir -p $datadir/calls
    mkdir -p $datadir/calls/neb15
    
    model=$datadir/train/neb15/training2/model_final.checkpoint
        
    for i in neb15 neb4 neb11 neb1 ;
    do
	mkdir -p $datadir/calls/neb15/$i
	megalodon $datadir/multiraw_test/$i \
		  --taiyaki-model-filename $model \
		  --reference $ref \
		  --devices 0 \
		  --outputs mod_basecalls mods \
		  --mod-motif Y GCNGC 1 \
		  --write-mods-text \
		  --processes 18 \
		  --overwrite \
		  --verbose-read-progress 3 \
		  --output-directory $datadir/calls/neb15/$i
    done
fi

if [ $1 == neb16 ] ; then
    mkdir -p $datadir/calls
    mkdir -p $datadir/calls/neb16
    
    model=$datadir/train/neb16/training2/model_final.checkpoint
        
    for i in neb16 neb5 neb11 neb1 ;
    do
	mkdir -p $datadir/calls/neb16/$i
	megalodon $datadir/multiraw_test/$i \
		  --taiyaki-model-filename $model \
		  --reference $ref \
		  --devices 0 \
		  --outputs mod_basecalls mods \
		  --mod-motif Y CCNGGC 5 \
		  --write-mods-text \
		  --processes 18 \
		  --overwrite \
		  --verbose-read-progress 3 \
		  --output-directory $datadir/calls/neb16/$i
    done
fi

if [ $1 == neb13 ] ; then
    mkdir -p $datadir/calls
    mkdir -p $datadir/calls/neb13
    
    model=$datadir/train/neb13/training2/model_final.checkpoint
        
    for i in neb13 neb3 neb11 neb1 ;
    do
	rmdir $datadir/calls/neb13/$i
	megalodon $datadir/multiraw_test/$i \
		  --taiyaki-model-filename $model \
		  --reference $ref \
		  --devices 0 \
		  --outputs mod_basecalls mods \
		  --mod-motif Y CTGCAG 3 \
		  --write-mods-text \
		  --processes 18 \
		  --verbose-read-progress 3 \
		  --output-directory $datadir/calls/neb13/$i
    done
fi


if [ $1 == neb14_cross ] ; then
    mkdir -p $datadir/calls
    mkdir -p $datadir/calls/neb14
    
    model=$datadir/train/neb14/training2/model_final.checkpoint
        
    ##for i in neb14 neb10 neb11 neb1 neb19 ;
    ##for i in neb10 ;
    for i in neb12 neb13 neb15 neb16 neb17 neb19 neb2 neb3 neb4 neb5 neb6 neb8 ;
    do
	megalodon $datadir/multiraw_test/$i \
		  --taiyaki-model-filename $model \
		  --reference $ref \
		  --devices 0 \
		  --outputs mod_basecalls mods \
		  --mod-motif Y GATC 3 \
		  --write-mods-text \
		  --processes 18 \
		  --verbose-read-progress 3 \
		  --output-directory $datadir/calls/neb14/$i
    done
fi
