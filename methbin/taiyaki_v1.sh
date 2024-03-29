#!/bin/bash

datadir=~/data/methbin
gdna="neb11 neb12 neb13 neb14 neb15 neb16 neb17 neb19 nebdcm"

if [ $1 == cp_sub ] ; then
    mkdir -p $datadir/multiraw_sub
    for i in $datadir/multiraw/* ;
    do
	prefix=`echo $i | rev | cut -d / -f 1 | rev`
	mkdir -p $datadir/multiraw_sub/$prefix

	for j in {1..50} ;
	do
	    if [ -f $i/${prefix}_$j.fast5 ] ; then
		cp $i/${prefix}_$j.fast5 $datadir/multiraw_sub/$prefix/
	    fi
	done
	
    done
fi


if [ $1 == gatherfq ] ; then
    mkdir -p $datadir/fastqs
    
    for i in $datadir/called/*;
    do
	prefix=`echo $i | rev | cut -d / -f 1 | rev`
	mkdir -p $datadir/fastqs/$prefix
	for j in {1..50} ;
	do
	    cat $datadir/called/$prefix/*_${j}_*fastq >> $datadir/fastqs/$prefix/${prefix}_100k.fq
	done
    done
fi


if [ $1 == create_ref_align ] ; then
    for i in $gdna;
    do
	mkdir -p $datadir/read_ref/$i
	minimap2 -I 16G -x map-ont -t 54 -a --secondary=no $datadir/reference/er2796.fa $datadir/fastqs/$i/${i}_100k.fq |\
	    samtools view -@ 36 -bST $datadir/reference/er2796.fa - > $datadir/read_ref/$i/${i}_100k.bam
    done
fi


if [ $1 == create_ref ] ; then
    ##get ground truth reads 
    for i in $gdna ;
    do
	get_refs_from_sam.py $datadir/reference/er2796.fa $datadir/read_ref/$i/${i}_100k.bam --min_coverage 0.8 > $datadir/read_ref/$i/${i}_100k.ref.fasta &
    done

fi


if [ $1 == mangle_ref ] ; then
    ##6ma will be Z, 5mc will be Y, 4mc will be X

    ##neb12
    sed -i -e 's/CCGG/XCGG/g' $datadir/read_ref/neb12/neb12_100k.ref.fasta

    ##neb13
    sed -i -e 's/CTGCAG/CTGYAG/g' $datadir/read_ref/neb13/neb13_100k.ref.fasta

    ##neb14
    sed -i -e 's/GATC/GATY/g' $datadir/read_ref/neb14/neb14_100k.ref.fasta

    ##neb15
    sed -i -e 's/GCAGC/GYAGC/g' $datadir/read_ref/neb15/neb15_100k.ref.fasta
    sed -i -e 's/GCCGC/GYCGC/g' $datadir/read_ref/neb15/neb15_100k.ref.fasta
    sed -i -e 's/GCGGC/GYGGC/g' $datadir/read_ref/neb15/neb15_100k.ref.fasta
    sed -i -e 's/GCTGC/GYTGC/g' $datadir/read_ref/neb15/neb15_100k.ref.fasta

    ##neb16
    sed -i -e 's/CCAGGC/CCAGGY/g' $datadir/read_ref/neb16/neb16_100k.ref.fasta
    sed -i -e 's/CCCGGC/CCCGGY/g' $datadir/read_ref/neb16/neb16_100k.ref.fasta
    sed -i -e 's/CCGGGC/CCGGGY/g' $datadir/read_ref/neb16/neb16_100k.ref.fasta
    sed -i -e 's/CCTGGC/CCTGGY/g' $datadir/read_ref/neb16/neb16_100k.ref.fasta

    ##neb17
    sed -i -e 's/GAATC/GZATC/g' $datadir/read_ref/neb17/neb17_100k.ref.fasta
    sed -i -e 's/GACTC/GZCTC/g' $datadir/read_ref/neb17/neb17_100k.ref.fasta
    sed -i -e 's/GAGTC/GZGTC/g' $datadir/read_ref/neb17/neb17_100k.ref.fasta
    sed -i -e 's/GATTC/GZTTC/g' $datadir/read_ref/neb17/neb17_100k.ref.fasta

    ##neb19
    sed -i -e 's/GATC/GZTC/g' $datadir/read_ref/neb19/neb19_100k.ref.fasta

    ##nebdcm
    sed -i -e 's/CCAGG/CYAGG/g' $datadir/read_ref/nebdcm/nebdcm_100k.ref.fasta
    sed -i -e 's/CCTGG/CYTGG/g' $datadir/read_ref/nebdcm/nebdcm_100k.ref.fasta
fi

if [ $1 == catref ] ; then
    ##model needs to recognize unmodified bases also, so each set of training reads needs 
    for i in $gdna ;
    do
	rm $datadir/read_ref/$i/${i}_100k_all.ref.fasta 
	head -n 200000 $datadir/read_ref/$i/${i}_100k.ref.fasta > $datadir/read_ref/$i/${i}_100k_all.ref.fasta
	head -n 200000 $datadir/read_ref/neb11/neb11_100k.ref.fasta >> $datadir/read_ref/$i/${i}_100k_all.ref.fasta 
    done
fi

if [ $1 == cp_raw ] ; then
    ##copy unmeth raw data to meth raw data
    ##not using for loop to save the time of copying neb11 to itself
    
    cp $datadir/multiraw_sub/neb11/*fast5 $datadir/multiraw_sub/neb12/

    ##5mc
    for i in neb13 neb14 neb15 neb16 nebdcm;
    do
	cp $datadir/multiraw_sub/neb11/*fast5 $datadir/multiraw_sub/$i/
    done

    ##6mA
    for i in neb17 neb19 ;
    do
	cp $datadir/multiraw_sub/neb11/*fast5 $datadir/multiraw_sub/$i/
    done
fi

if [ $1 == get_params ] ; then
    mkdir -p $datadir/train
    for i in $gdna ;
    do
	mkdir -p $datadir/train/$i
	generate_per_read_params.py --jobs 36 $datadir/multiraw_sub/$i > $datadir/train/$i/${i}_100k_modbase.tsv
    done
fi

if [ $1 == get_example ] ; then
    ##obtain pretrained checkpoint file 
    ##wget  https://s3-eu-west-1.amazonaws.com/ont-research/taiyaki_modbase.tar.gz
    cp ~/software/taiyaki/taiyaki_modbase/pretrained/r941_dna_minion.checkpoint $datadir/train/
fi

if [ $1 == 4mc_map_read_file ] ; then
    ##have to copy the pretrained model from the example data (see get_example)
    
    ##neb12 is the only 4mC (X to C)  
    prepare_mapped_reads.py \
	--overwrite \
	--jobs 36 \
	--mod X C neb12 \
	$datadir/multiraw_sub/neb12 \
	$datadir/train/neb12/neb12_100k_modbase.tsv \
	$datadir/train/neb12/neb12_100k_modbase.hdf5 \
	$datadir/train/r941_dna_minion.checkpoint \
	$datadir/read_ref/neb12/neb12_100k_all.ref.fasta
fi


if [ $1 == 5mc_map_read_file ] ; then
    ##5mC samples
    ##for i in neb13 neb14 neb15 neb16 nebdcm;
    for i in nebdcm ;
    do
	prepare_mapped_reads.py \
	    --overwrite \
	    --jobs 36 \
	    --mod Y C $i \
	    $datadir/multiraw_sub/$i \
	    $datadir/train/$i/${i}_100k_modbase.tsv \
	    $datadir/train/$i/${i}_100k_modbase.hdf5 \
	    $datadir/train/r941_dna_minion.checkpoint \
	    $datadir/read_ref/$i/${i}_100k_all.ref.fasta
    done
fi


if [ $1 == 6ma_map_read_file ] ; then
    ##6mA samples
    for i in neb17 neb19 ;
    do
	prepare_mapped_reads.py \
	    --overwrite \
	    --jobs 36 \
	    --mod Z A $i \
	    $datadir/multiraw_sub/$i \
	    $datadir/train/$i/${i}_100k_modbase.tsv \
	    $datadir/train/$i/${i}_100k_modbase.hdf5 \
	    $datadir/train/r941_dna_minion.checkpoint \
	    $datadir/read_ref/$i/${i}_100k_all.ref.fasta &> $datadir/train/$i/$i.mapreads.log
    done
fi


if [ $1 == 4mc_train ] ; then
    for i in neb12 ;
    do
	mkdir -p $datadir/train/$i/training
	##train_mod_flipflop.py --overwrite --device 0 --mod_factor 0.01 --outdir $datadir/train/$i/training ~/software/taiyaki/models/mGru_cat_mod_flipflop.py $datadir/train/$i/${i}_100k_modbase.hdf5

	mkdir -p $datadir/train/$i/training2
	train_mod_flipflop.py --overwrite --device 0 --mod_factor 1.0 --outdir $datadir/train/$i/training2 \
			      $datadir/train/$i/training/model_final.checkpoint \
			      $datadir/train/$i/${i}_100k_modbase.hdf5
    done
fi
