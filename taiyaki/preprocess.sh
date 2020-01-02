#!/bin/bash

if [ $1 == make_multi ] ; then
    ##for i in damlike unmeth ;
    for i in unmeth ;
    do
	mkdir -p ~/data/$i/multiraw
	single_to_multi_fast5 -t 36 -i ~/data/$i/raw -s ~/data/$i/multiraw --recursive
    done
fi

##in docker image, called using guppy
##https://medium.com/@kepler_00/nanopore-gpu-basecalling-using-guppy-on-ubuntu-18-04-and-nvidia-docker-v2-with-a-rtx-2080-d875945e5c8d


if [ $1 == create_ref_align ] ; then
    ##align to create ground truth reads based on alignment to reference genome
    for i in damlike unmeth ;
    ##for i in unmeth ;
    do
	mkdir -p ~/data/$i/read_ref
	cat ~/data/$i/called/*fastq > ~/data/$i/$i.fastq
	minimap2 -I 16G -x map-ont -t 72 -a --secondary=no ~/data/refs/er2796.fasta ~/data/$i/$i.fastq | samtools view -bST ~/data/refs/er2796.fasta - > ~/data/$i/read_ref/$i.bam
    done
fi

if [ $1 == create_ref ] ; then
    ##get ground truth reads 
    for i in damlike unmeth ;
    ##for i in unmeth ;
    do
	get_refs_from_sam.py ~/data/refs/er2796.fasta ~/data/$i/read_ref/$i.bam --min_coverage 0.8 > ~/data/$i/read_ref/$i.ref.fasta
    done
fi

if [ $1 == mangle_ref ] ; then
    sed -i -e 's/GATC/GZTC/g' ~/data/damlike/read_ref/damlike.ref.fasta
fi

if [ $1 == cat_ref ] ; then
    cat ~/data/damlike/read_ref/damlike.ref.fasta ~/data/unmeth/read_ref/unmeth.ref.fasta > ~/data/all_raw/all.ref.fasta
fi

if [ $1 == get_params ] ; then
    generate_per_read_params.py --jobs 72 ~/data/all_raw/multiraw > ~/data/train/modbase.tsv
fi

if [ $1 == map_read_file ] ; then
    ##have to copy the pretrained model from the example data
    prepare_mapped_reads.py --jobs 72 --mod Z A dam ~/data/all_raw/multiraw ~/data/train/modbase.tsv ~/data/train/modbase.hdf5 ~/data/train/pretrained/r941_dna_minion.checkpoint ~/data/all_raw/all.ref.fasta
fi