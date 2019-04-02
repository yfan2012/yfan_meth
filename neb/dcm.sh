#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/meth/180628_neb_dcm
srcdir=~/Code/utils/marcc
ref=/scratch/groups/mschatz1/cpowgs/meth/refs/er2796.fasta

if [ $1 == untar ] ; then
    mkdir -p $datadir/raw
    sbatch --output=$datadir/batch_logs/untar.out --job-name=ut_ncm $srcdir/untar.scr $datadir/180628_neb_dcm.tar.gz $datadir
fi


if [ $1 == call ] ; then
    mkdir -p $datadir/called
    mkdir -p $datadir/call_done
    mkdir -p $datadir/call_logs
    sbatch --array=0-290 --output=$datadir/call_logs/180628_neb_dcm.%A_%a.out --job-name=180628_neb_dcm $srcdir/call.scr $datadir
fi

if [ $1 == fastq ] ; then
    mkdir -p $datadir/fastqs
    cat $datadir/called/*/workspace/pass/*fastq > $datadir/fastqs/180628_neb_dcm.pass.fastq
    cat $datadir/called/*/workspace/fail/*fastq > $datadir/fastqs/180628_neb_dcm.fail.fastq
    cat $datadir/called/*/workspace/*/*fastq > $datadir/fastqs/180628_neb_dcm.all.fastq
fi

if [ $1 == align ] ; then
    mkdir -p $datadir/align
    sbatch --output=$datadir/batch_logs/align_pass.out --job-name=a_dcmp ./align.scr $datadir/fastqs/180628_neb_dcm.pass.fastq $ref
    sbatch --output=$datadir/batch_logs/align_fail.out --job-name=a_dcmf ./align.scr $datadir/fastqs/180628_neb_dcm.fail.fastq $ref
    sbatch --output=$datadir/batch_logs/align_all.out --job-name=a_dcma ./align.scr $datadir/fastqs/180628_neb_dcm.all.fastq $ref
fi

if [ $1 == npidx ] ; then
    mkdir -p $datadir/eventalign
    rm -f $datadir/eventalign/workspace.fofn
    touch $datadir/eventalign/workspace.fofn

    for i in $datadir/called/*/sequencing_summary.txt ;
    do
	echo $i >> $datadir/eventalign/workspace.fofn
    done

    seqtk seq -a $datadir/fastqs/180628_neb_dcm.pass.fastq > $datadir/fastqs/180628_neb_dcm.pass.fa
    nanopolish index -d $datadir/raw -f $datadir/eventalign/workspace.fofn $datadir/fastqs/180628_neb_dcm.pass.fa
fi


if [ $1 == eventalign ] ; then
    ml samtools
    samtools index $datadir/align/180628_neb_dcm.pass.sorted.bam
    nanopolish eventalign \
	       --scale-events \
	       -t 36 \
	       -r $datadir/fastqs/180628_neb_dcm.pass.fa \
	       -b $datadir/align/180628_neb_dcm.pass.sorted.bam \
	       -g $ref > $datadir/eventalign/180628_neb_dcm.eventalign.tsv
fi
    
if [ $1 == subalign ] ; then
    ml samtools

    head -n 200000 $datadir/fastqs/180628_neb_dcm.pass.fastq > $datadir/fastqs/180628_neb_dcm.sub50k.fastq
    sbatch --output=$datadir/batch_logs/align_sub50k.out --job-name=subdcmalign ./align.scr $datadir/fastqs/180628_neb_dcm.sub50k.fastq $ref
fi

if [ $1 == subeventalign ] ; then
    mkdir -p $datadir/eventalign
    rm -f $datadir/eventalign/workspace.fofn
    touch $datadir/eventalign/workspace.fofn

    for i in $datadir/called/*/sequencing_summary.txt ;
    do
	echo $i >> $datadir/eventalign/workspace.fofn
    done
    
    seqtk seq -a $datadir/fastqs/180628_neb_dcm.sub50k.fastq > $datadir/fastqs/180628_neb_dcm.sub50k.fa
    nanopolish index -d $datadir/raw -f $datadir/eventalign/workspace.fofn $datadir/fastqs/180628_neb_dcm.sub50k.fa

    nanopolish eventalign \
	       --scale-events \
	       -t 36 \
	       -r $datadir/fastqs/180628_neb_dcm.sub50k.fa \
	       -b $datadir/align/180628_neb_dcm.sub50k.sorted.bam \
	       -g $ref > $datadir/eventalign/180628_neb_dcm.sub50k.eventalign.tsv
fi

if [ $1 == shortref ] ; then
    ##use first 40k bases of genome and 50k reads
    sbatch --output=$datadir/batch_logs/align_sub50k_genome40k.out --job-name=subdcmalign ./align.scr $datadir/fastqs/180628_neb_dcm.sub50k.fastq $datadir/align/er2796_40k.fasta
fi

if [ $1 == subgenomeeventalign ] ; then
    mkdir -p $datadir/eventalign
    rm -f $datadir/eventalign/workspace.fofn
    touch $datadir/eventalign/workspace.fofn

    nanopolish eventalign \
	       --scale-events \
	       -t 36 \
	       -r $datadir/fastqs/180628_neb_dcm.sub50k.fa \
	       -b $datadir/align/180628_neb_dcm.sub50k.sorted.bam \
	       -g $datadir/align/er2796_40k.fasta > $datadir/eventalign/180628_neb_dcm.sub50k.genome40k.eventalign.tsv
fi
