#!/bin/bash

if [ $1 == make_multi ] ; then
    for i in cpg unmeth ;
    do
	mkdir -p ~/cpg/$i/multiraw
	single_to_multi_fast5 -t 36 -i ~/cpg/$i/raw -s ~/cpg/$i/multiraw --recursive
    done
fi

##in docker image, called using guppy
##https://medium.com/@kepler_00/nanopore-gpu-basecalling-using-guppy-on-ubuntu-18-04-and-nvidia-docker-v2-with-a-rtx-2080-d875945e5c8d


if [ $1 == create_ref_align ] ; then
    ##align to create ground truth reads based on alignment to reference genome
    for i in cpg unmeth ;
    ##for i in unmeth ;
    do
	mkdir -p ~/cpg/$i/read_ref
	cat ~/cpg/$i/called/*fastq > ~/cpg/$i/$i.fastq
	minimap2 -I 16G -x map-ont -t 36 -a --secondary=no ~/cpg/refs/GRCH38.fa ~/cpg/$i/$i.fastq | samtools view -bST ~/cpg/refs/GRCH38.fa - > ~/cpg/$i/read_ref/$i.bam
    done
fi

if [ $1 == create_ref ] ; then
    ##get ground truth reads 
    for i in cpg unmeth ;
    ##for i in unmeth ;
    do
	get_refs_from_sam.py ~/cpg/refs/GRCH38.fa ~/cpg/$i/read_ref/$i.bam --min_coverage 0.8 > ~/cpg/$i/read_ref/$i.ref.fasta
    done
fi

if [ $1 == mangle_ref ] ; then
    sed -i -e 's/CG/YG/g' ~/cpg/cpg/read_ref/cpg.ref.fasta
fi

if [ $1 == cat_ref ] ; then
    cat ~/cpg/cpg/read_ref/cpg.ref.fasta ~/cpg/unmeth/read_ref/unmeth.ref.fasta > ~/cpg/all_raw/all.ref.fasta
fi

if [ $1 == get_params ] ; then
    generate_per_read_params.py --jobs 36 ~/cpg/all_raw/multiraw > ~/cpg/train/modbase.tsv
fi

if [ $1 == map_read_file ] ; then
    ##have to copy the pretrained model from the example cpg
    prepare_mapped_reads.py --jobs 72 --mod Y C CpG ~/cpg/all_raw/multiraw ~/cpg/train/modbase.tsv ~/cpg/train/modbase.hdf5 ~/cpg/train/pretrained/r941_dna_minion.checkpoint ~/cpg/all_raw/all.ref.fasta
fi
