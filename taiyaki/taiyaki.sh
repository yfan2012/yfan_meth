#!/bin/bash

dat=/data/yfan
mod=$2
old_letter=$3
new_letter=$4
motif=$5
new_motif=$6

## bash taiyaki.sh make_multi dam A Z GATC GZTC

if [ $1 == make_multi ] ; then
    for i in $mod unmeth ;
    do
	mkdir -p $dat/$mod/$i/multiraw
	single_to_multi_fast5 -t 36 -i $dat/$mod/$i/raw -s $dat/$mod/$i/multiraw --recursive
    done
fi

##in docker image, called using guppy
##https://medium.com/@kepler_00/nanopore-gpu-basecalling-using-guppy-on-ubuntu-18-04-and-nvidia-docker-v2-with-a-rtx-2080-d875945e5c8d

if [ $1 == pick_50 ] ; then
    for i in $mod unmeth ;
    do
	for j in {0..49} ;
	do
	    mv $dat/$mod/$i/multiraw/batch_$j.fast5 $dat/$mod/$i/multiraw/${i}_batch_$j.fast5
	done
	rm $dat/$mod/$i/multiraw/batch*
    done
fi
	     

if [ $1 == call ] ; then
    for i in $mod unmeth ;
    do
	guppy_basecaller -i $dat/$mod/$i/multiraw -s $dat/$mod/$i/called --flowcell FLO-MIN106 --kit SQK-LSK108 -x 'cuda:1'
    done
fi


if [ $1 == create_ref_align ] ; then
    ##align to create ground truth reads based on alignment to reference genome
    for i in $mod unmeth ;
	     ##for i in unmeth ;
    do
	mkdir -p $dat/$mod/$i/read_ref
	cat $dat/$mod/$i/called/*fastq > $dat/$mod/$i/$i.fastq
	minimap2 -I 16G -x map-ont -t 48 -a --secondary=no $dat/$mod/*.fasta $dat/$mod/$i/$i.fastq | samtools view -bST $dat/$mod/er2796.fasta - > $dat/$mod/$i/read_ref/$i.bam
	samtools sort -@ 48 -o $dat/$mod/$i/read_ref/$i.sorted.bam $dat/$mod/$i/read_ref/$i.bam
	samtools index $dat/$mod/$i/read_ref/$i.sorted.bam
    done
fi



if [ $1 == create_ref ] ; then
    ##get ground truth reads
    for i in $mod unmeth ;
	     ##for i in unmeth ;
    do
	get_refs_from_sam.py $dat/$mod/*.fasta $dat/$mod/$i/read_ref/$i.bam --min_coverage 0.8 > $dat/$mod/$i/read_ref/$i.ref.fasta
    done
fi

if [ $1 == mangle_ref ] ; then
    ##sed -i -e 's/GATC/GZTC/g' ~/data/damlike/read_ref/damlike.ref.fasta
    sed -i -e 's/'$motif'/'$new_motif'/g' $dat/$mod/$mod/read_ref/$mod.ref.fasta
fi

if [ $1 == cat_ref ] ; then
    mkdir -p $dat/$mod/all_raw
    cat $dat/$mod/$mod/read_ref/$mod.ref.fasta $dat/$mod/unmeth/read_ref/unmeth.ref.fasta > $dat/$mod/all_raw/all.ref.fasta

    mkdir -p $dat/$mod/all_raw/multiraw
    cp $dat/$mod/$mod/multiraw/*fast5 $dat/$mod/all_raw/multiraw/
    cp $dat/$mod/unmeth/multiraw/*fast5 $dat/$mod/all_raw/multiraw/
fi

if [ $1 == get_params ] ; then
    mkdir -p $dat/$mod/train
    generate_per_read_params.py --jobs 48 $dat/$mod/all_raw/multiraw > $dat/$mod/train/modbase.tsv
fi

if [ $1 == map_read_file ] ; then
    ##have to copy the pretrained model from the example data
    prepare_mapped_reads.py --jobs 24 --mod $new_letter $old_letter $mod $dat/$mod/all_raw/multiraw $dat/$mod/train/modbase.tsv $dat/$mod/train/modbase.hdf5 $dat/pretrained/r941_dna_minion.checkpoint $dat/$mod/all_raw/all.ref.fasta
fi
