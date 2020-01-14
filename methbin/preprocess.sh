#!/bin/bash

datadir=/uru/Data/Nanopore/projects/methbin
gdna="neb11 neb12 neb13 neb14 neb15 neb16 neb17 neb19 nebdcm"


if [ $1 == call ] ; then
    ##using whack:guppy_bcall docker
    ##docker run --runtime=nvidia --name yfan -i -t -v /uru/Data/Nanopore/projects/methbin:/media yfan2012/whack:guppy_bcall /bin/bash
    mkdir -p /media/called
    for i in /media/multiraw/* ;
    do
	prefix=`echo $i | rev | cut -d / -f 1 | rev`
	mkdir /media/called/$prefix
	guppy_basecaller -i $i -s /media/called/$prefix --flowcell FLO-MIN106 --kit SQK-LSK108 -x "cuda:0"
    done
fi


if [ $1 == gatherfq ] ; then
    mkdir -p $datadir/fastqs
    for i in $datadir/called/*;
    do
	(prefix=`echo $i | rev | cut -d / -f 1 | rev`
	mkdir -p $datadir/fastqs/$prefix
	cat $datadir/called/$prefix/*fastq > $datadir/fastqs/$prefix/$prefix.fq ) &
    done
fi


if [ $1 == create_ref_align ] ; then
    mkdir -p $datadir/fastqs
    for i in $gdna;
    do
	mkdir -p $datadir/read_ref/$i
	minimap2 -I 16G -x map-ont -t 54 -a --secondary=no $datadir/reference/er2796.fa $datadir/fastqs/$i/$i.fq |\
	    samtools view -@ 36 -bST $datadir/reference/er2796.fa - > $datadir/read_ref/$i/$i.bam
    done
fi


if [ $1 == create_ref ] ; then
    ##get ground truth reads 
    for i in $gdna ;
    do
	get_refs_from_sam.py $datadir/reference/er2796.fa $datadir/read_ref/$i/$i.bam --min_coverage 0.8 > $datadir/read_ref/$i/$i.ref.fasta &
    done
fi

if [ $1 == mangle_ref ] ; then
    ##6ma will be Z, 5mc will be Y, 4mc will be X

    ##neb12
    sed -i -e 's/CCGG/XCGG/g' $datadir/read_ref/neb12/neb12.ref.fasta

    ##neb13
    sed -i -e 's/CTGCAG/CTGYAG/g' $datadir/read_ref/neb13/neb13.ref.fasta

    ##neb14
    sed -i -e 's/GATC/GATY/g' $datadir/read_ref/neb14/neb14.ref.fasta

    ##neb15
    sed -i -e 's/GCAGC/GYAGC/g' $datadir/read_ref/neb15/neb15.ref.fasta
    sed -i -e 's/GCCGC/GYCGC/g' $datadir/read_ref/neb15/neb15.ref.fasta
    sed -i -e 's/GCGGC/GYGGC/g' $datadir/read_ref/neb15/neb15.ref.fasta
    sed -i -e 's/GCTGC/GYTGC/g' $datadir/read_ref/neb15/neb15.ref.fasta

    ##neb16
    sed -i -e 's/CCAGGC/CCAGGY/g' $datadir/read_ref/neb16/neb16.ref.fasta
    sed -i -e 's/CCCGGC/CCCGGY/g' $datadir/read_ref/neb16/neb16.ref.fasta
    sed -i -e 's/CCGGGC/CCGGGY/g' $datadir/read_ref/neb16/neb16.ref.fasta
    sed -i -e 's/CCTGGC/CCTGGY/g' $datadir/read_ref/neb16/neb16.ref.fasta

    ##neb17
    sed -i -e 's/GAATC/GZATC/g' $datadir/read_ref/neb17/neb17.ref.fasta
    sed -i -e 's/GACTC/GZCTC/g' $datadir/read_ref/neb17/neb17.ref.fasta
    sed -i -e 's/GAGTC/GZGTC/g' $datadir/read_ref/neb17/neb17.ref.fasta
    sed -i -e 's/GATTC/GZTTC/g' $datadir/read_ref/neb17/neb17.ref.fasta

    ##neb19
    sed -i -e 's/GATC/GZTC/g' $datadir/read_ref/neb19/neb19.ref.fasta

    ##nebdcm
    sed -i -e 's/CCAGG/CYAGG/g' $datadir/read_ref/nebdcm/nebdcm.ref.fasta
    sed -i -e 's/CCTGG/CYTGG/g' $datadir/read_ref/nebdcm/nebdcm.ref.fasta
fi

if [ $1 == catref ] ; then
    ##model needs to recognize unmodified bases also, so each set of training reads needs 
    for i in $gdna ;
    do
	cat $datadir/read_ref/$i/$i.ref.fasta $datadir/read_ref/neb11/neb11.ref.fasta > $datadir/read_ref/$i/${i}_all.ref.fasta &
    done
fi

if [ $1 == catref_100k ] ; then
    ##get only 100k reads
    head -n 200000 $datadir/read_ref/neb11/neb11.ref.fasta > $datadir/read_ref/neb11/neb11_100k.ref.fasta 
    for i in $gdna ;
    do
	head -n 200000 $datadir/read_ref/$i/${i}.ref.fasta > $datadir/read_ref/$i/${i}_100k.ref.fasta
	cat $datadir/read_ref/$i/${i}_100k.ref.fasta $datadir/read_ref/neb11/neb11_100k.ref.fasta > $datadir/read_ref/$i/${i}_all200k.ref.fasta
    done
fi
	

if [ $1 == get_params ] ; then
    mkdir -p $datadir/train
    for i in $gdna ;
    do
	mkdir -p $datadir/train/$i
	generate_per_read_params.py --jobs 72 $datadir/multiraw/$i > $datadir/train/$i/${i}_modbase.tsv
    done
fi

if [ $1 == get_example ] ; then
    ##obtain pretrained checkpoint file 
    ##wget  https://s3-eu-west-1.amazonaws.com/ont-research/taiyaki_modbase.tar.gz
    cp ~/software/taiyaki/taiyaki_modbase/pretrained/r941_dna_minion.checkpoint $datadir/train/
fi

    
if [ $1 == map_read_file ] ; then
    ##have to copy the pretrained model from the example data

    ##neb12 is the only 4mC (X to C)  
    prepare_mapped_reads.py \
	--overwrite \
	--jobs 24 \
	--mod X C neb12 \
	$datadir/multiraw \
	$datadir/train/neb12/neb12_modbase.tsv \
	$datadir/train/neb12/neb12_modbase.hdf5 \
	$datadir/train/r941_dna_minion.checkpoint \
	$datadir/read_ref/neb12/neb12_all200k.ref.fasta

    ##5mC samples
    for i in neb13 neb14 neb15 neb16 ;
    do
	prepare_mapped_reads.py \
	    --overwrite \
	    --jobs 24 \
	    --mod Y C $i \
	    $datadir/multiraw \
	    $datadir/train/$i/${i}_modbase.tsv \
	    $datadir/train/$i/${i}_modbase.hdf5 \
	    $datadir/train/r941_dna_minion.checkpoint \
	    $datadir/read_ref/$i/${i}_all200k.ref.fasta
    done

    ##6mA samples
    for i in neb17 neb19 ;
    do
	prepare_mapped_reads.py \
	    --overwrite \
	    --jobs 24 \
	    --mod Z A $i \
	    $datadir/multiraw \
	    $datadir/train/$i/${i}_modbase.tsv \
	    $datadir/train/$i/${i}_modbase.hdf5 \
	    $datadir/train/r941_dna_minion.checkpoint \
	    $datadir/read_ref/$i/${i}_all200k.ref.fasta
    done
fi
