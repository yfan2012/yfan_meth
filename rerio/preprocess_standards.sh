#!/bin/bash

rawdir=/dilithium/Data/Nanopore/oxford/181010_ecoli_methylation_standards/all_data
datadir=/mithril/Data/Nanopore/projects/methbin
allsamps='ecoli_CpG ecoli_CpGGpC ecoli_GpC ecoli_Unmethylated'

if [ $1 == call ] ; then
    for i in $allsamps ;
    do
	mkdir -p ~/data/rerio/$i/called
	
	guppy_basecaller \
	    -i ~/data/rerio/$i/multiraw_sub \
	    -s ~/data/rerio/$i/called \
	    -d ~/software/rerio/basecall_models/ \
	    -c res_dna_r941_min_modbases-all-context_v001.cfg \
	    -x "cuda:0" \
	    --fast5_out
	##if u give config, u can't give kit/fcell

	##move to mithril
	mkdir -p $datadir/rerio/$i
	mv ~/data/rerio/$i/called $datadir/rerio/$i/
    done
fi

ref=/dilithium/Data/Nanopore/Analysis/181010_ecoli_methylation_standards/Reference/ecoli_mg1655.fasta
if [ $1 == gather_align ] ; then
    for i in $allsamps ;
    do
	mkdir -p $datadir/fastqs/$i
	mkdir -p $datadir/align/$i
	
	cat $datadir/rerio/$i/called/*fastq > $datadir/fastqs/$i/${i}_rerio.fq

	minimap2 -t 36 -ax map-ont $ref $datadir/fastqs/$i/${i}_rerio.fq |\
	    samtools view -@ 36 -b |
	    samtools sort -@ 36 -o $datadir/align/$i/${i}_rerio.sorted.bam
	samtools index $datadir/align/$i/${i}_rerio.sorted.bam

	samtools calmd -@ 36 $datadir/align/$i/${i}_rerio.sorted.bam $ref |
	    samtools view -@ 36 -b |
	    samtools sort -@ 36 -o $datadir/align/$i/${i}_rerio.md.sorted.bam
	samtools index $datadir/align/$i/${i}_rerio.md.sorted.bam
    done
fi

	
if [ $1 == per_read_errors ] ; then
   for i in $allsamps ;
   do
       ~/Code/timp_nanopore/oxford/bam_extract.py -i $datadir/align/$i/${i}_rerio.md.sorted.bam &
   done
   gunzip $datadir/align/*/*_rerio.md.sorted.csv.gz
fi
   
if [ $1 == assess ] ; then

    for i in neb17 ;
    do
	mkdir -p $datadir/rerio/$i/assess
	python ~/Code/methylation/utils/methcall_check.py \
	       -r $ref \
	       -b $datadir/align/$i/${i}_rerio.md.sorted.bam \
	       -o ~/data/$i.GANTC.csv \
	       -f $datadir/rerio/$i/called/workspace \
	       -m GANTC \
	       -p 1 \
	       -s 0 \
	       -t 54
    done
fi

if [ $1 == assess_single ] ; then

    for i in neb17 ;
    do
	mkdir -p $datadir/rerio/$i/assess
	python ~/Code/methylation/utils/methcall_check_single.py \
	       -r $ref \
	       -b $datadir/align/$i/${i}_rerio.md.sorted.bam \
	       -o ~/data/$i.GANTC_singletest.csv \
	       -f $datadir/rerio/$i/called/workspace \
	       -m GANTC \
	       -p 1 \
	       -s 0

    done
fi

	       
    

