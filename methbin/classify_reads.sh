#!/bin/bash

##analysis for using eventaligh to classify reads independently of each other
##idea is to aggregate kmers over a read

datadir=/uru/Data/Nanopore/projects/methbin
gdna="neb11 neb12 neb13 neb14 neb15 neb16 neb17 neb19 nebdcm"
plasmids="neb1 neb2 neb3 neb4 neb5 neb6 neb8 neb9 neb10"

if [ $1 == gather40k ] ; then
    ##align 40k reads per sample
    for i in $datadir/fastqs/* ;
    do
	(prefix=`echo $i | rev | cut -d / -f 1 | rev `
	head -n 160000 $i/$prefix.fq > $i/${prefix}_40k.fq) &
    done
fi

if [ $1 == rename_ref ] ; then
    refdir=$datadir/reference
    cp $refdir/pRRSlac.fa $refdir/neb1.fa
    cp $refdir/pRRStetMPspJDRI.fa $refdir/neb2.fa
    cp $refdir/pRRStetMSin395ORF119.fa $refdir/neb3.fa
    cp $refdir/pBR322_MFnu4H.fa $refdir/neb4.fa
    cp $refdir/pLacZZMSdeAII.fa $refdir/neb5.fa
    cp $refdir/pUC19HinfI.fa $refdir/neb6.fa
    cp $refdir/XP12.fa $refdir/neb8.fa
    cp $refdir/pRRStetM3_1BstXII.fa $refdir/neb9.fa
    cp $refdir/pT7MSin395ORF667.fa $refdir/neb10.fa    
fi
    
if [ $1 == align ] ; then
    ##align 40k reads per sample
    mkdir -p $datadir/align
    for i in $gdna ;
    do
	ref=$datadir/reference/er2796.fa
	mkdir -p $datadir/align/$i
	
	minimap2 -t 36 -ax map-ont $ref $datadir/fastqs/$i/${i}_40k.fq |
	    samtools view -@ 36 -b |
	    samtools sort -@ 36 -o $datadir/align/$i/${i}_40k.sorted.bam -T $datadir/align/$i/$i.tmp
	samtools index $datadir/align/$i/${i}_40k.sorted.bam
    done
    
    for i in $plasmids ;
    do
	ref=$datadir/reference/$i.fa
	mkdir -p $datadir/align/$i
	
	minimap2 -t 36 -ax map-ont $ref $datadir/fastqs/$i/${i}_40k.fq |
	    samtools view -@ 36 -b |
	    samtools sort -@ 36 -o $datadir/align/$i/${i}_40k.sorted.bam -T $datadir/align/$i/$i.tmp
	samtools index $datadir/align/$i/${i}_40k.sorted.bam
    done
fi

if [ $1 == np_idx ] ; then
    for i in $gdna $plasmids ;
    do
	##make the index
	(nanopolish index \
		   -d $datadir/multiraw/$i \
		   -s $datadir/called/$i/sequencing_summary.txt \
		   $datadir/fastqs/$i/${i}_40k.fq) &
    done
fi

if [ $1 == eventalign ] ; then
    mkdir -p $datadir/eventalign

    for i in $gdna ;
    do
	mkdir -p $datadir/eventalign/$i
	ref=$datadir/reference/er2796.fa

	nanopolish eventalign \
		   --scale-events \
		   --samples \
		   -t 36 \
		   -r $datadir/fastqs/$i/${i}_40k.fq \
		   -g $ref \
		   -b $datadir/align/$i/${i}_40k.sorted.bam > $datadir/eventalign/$i/${i}_40kraw.eventalign.tsv
    done

    for i in $plasmids ;
    do
	mkdir -p $datadir/eventalign/$i
	ref=$datadir/reference/$i.fa

	nanopolish eventalign \
		   --scale-events \
		   --samples \
		   -t 36 \
		   -r $datadir/fastqs/$i/${i}_40k.fq \
		   -g $ref \
		   -b $datadir/align/$i/${i}_40k.sorted.bam > $datadir/eventalign/$i/${i}_40kraw.eventalign.tsv
    done
fi

if [ $1 == compore ] ; then
    mkdir -p $datadir/eventalign_collapsed

    for i in $gdna $plasmids ;
    ##for i in neb11 ;
    do
	mkdir -p $datadir/eventalign_collapsed/$i
	NanopolishComp Eventalign_collapse \
		       -i $datadir/eventalign/$i/${i}_40kraw.eventalign.tsv \
		       -o $datadir/eventalign_collapsed/$i \
		       -t 36 \
		       --max_reads 20 \
		       -p $i
    done
fi
