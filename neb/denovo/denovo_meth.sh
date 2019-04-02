#!/bin/bash

##generalized dcm.sh for all other meth samps

srcdir=~/Code/utils/marcc
##gDNA="171020_neb11 171005_neb12 171012_neb13 170906_neb14 171020_neb15 171003_neb16 171019_neb17 171019_neb19 180628_neb_dcm"
gDNA='180628_neb_dcm'

if [ $1 == untar ] ; then
    mkdir -p $datadir/raw
    sbatch --output=$datadir/batch_logs/untar.out --job-name=ut_ncm $srcdir/untar.scr $datadir/$prefix.tar.gz $datadir
fi


if [ $1 == call ] ; then
    mkdir -p $datadir/called
    mkdir -p $datadir/call_done
    mkdir -p $datadir/call_logs
    sbatch --array=0-290 --output=$datadir/call_logs/$prefix.%A_%a.out --job-name=$prefix $srcdir/call.scr $datadir
fi

if [ $1 == fastq ] ; then
    for prefix in $gDNA ;
    do
	echo $prefix
	datadir=/scratch/groups/mschatz1/cpowgs/meth/$prefix
	mkdir -p $datadir/fastqs
	cat $datadir/called/*/workspace/pass/*fastq > $datadir/fastqs/$prefix.fq
	cat $datadir/called/*/workspace/fail/*fastq > $datadir/fastqs/$prefix.fail.fq
	cat $datadir/called/*/workspace/*/*fastq > $datadir/fastqs/$prefix.all.fq
    done
fi

if [ $1 == align ] ; then
    mkdir -p $datadir/align
    sbatch --output=$datadir/batch_logs/align_pass.out --job-name=a_dcmp ./align.scr $datadir/fastqs/$prefix.pass.fastq $ref
    sbatch --output=$datadir/batch_logs/align_fail.out --job-name=a_dcmf ./align.scr $datadir/fastqs/$prefix.fail.fastq $ref
    sbatch --output=$datadir/batch_logs/align_all.out --job-name=a_dcma ./align.scr $datadir/fastqs/$prefix.all.fastq $ref
fi

if [ $1 == npidx ] ; then
    for prefix in $gDNA ;
    do
	(datadir=/scratch/groups/mschatz1/cpowgs/meth/$prefix
	mkdir -p $datadir/fastqs
	rm -f $datadir/fastqs/workspace.fofn
	touch $datadir/fastqs/workspace.fofn
	
	for i in $datadir/called/*/sequencing_summary.txt ;
	do
	    echo $i >> $datadir/fastqs/workspace.fofn
	done

	seqtk seq -a $datadir/fastqs/$prefix.fq > $datadir/fastqs/$prefix.fa
	nanopolish index -d $datadir/raw -f $datadir/eventalign/workspace.fofn $datadir/fastqs/$prefix.fa ) &
    done
fi


if [ $1 == eventalign ] ; then

    basedir=/scratch/groups/mschatz1/cpowgs/meth
    storedir=/work-zfs/mschatz1/cpowgs/analysis/meth/neb/eventalign
    ref=$basedir/refs/er2796.fasta
    mkdir -p $storedir
    for prefix in $gDNA ;
    do
	datadir=/scratch/groups/mschatz1/cpowgs/meth/$prefix
	nanopolish eventalign \
		   --scale-events \
		   -t 72 \
		   -r $datadir/fastqs/$prefix.fa \
		   -b $datadir/align/$prefix.sorted.bam \
		   -g $ref > $storedir/$prefix.eventalign.tsv
    done
fi
    
if [ $1 == subalign ] ; then
    ml samtools
    head -n 200000 $datadir/fastqs/$prefix.fq > $datadir/fastqs/$prefix.sub50k.fastq
    sbatch --output=$datadir/batch_logs/align_sub50k.out --job-name=subdcmalign ./align.scr $datadir/fastqs/$prefix.sub50k.fastq $ref
fi

if [ $1 == subeventalign ] ; then
    mkdir -p $datadir/eventalign
    rm -f $datadir/eventalign/workspace.fofn
    touch $datadir/eventalign/workspace.fofn

    for i in $datadir/called/*/sequencing_summary.txt ;
    do
	echo $i >> $datadir/eventalign/workspace.fofn
    done
    
    seqtk seq -a $datadir/fastqs/$prefix.sub50k.fastq > $datadir/fastqs/$prefix.sub50k.fa
    nanopolish index -d $datadir/raw -f $datadir/eventalign/workspace.fofn $datadir/fastqs/$prefix.sub50k.fa

    nanopolish eventalign \
	       --scale-events \
	       -t 36 \
	       -r $datadir/fastqs/$prefix.sub50k.fa \
	       -b $datadir/align/$prefix.sub50k.sorted.bam \
	       -g $ref > $datadir/eventalign/$prefix.sub50k.eventalign.tsv
fi

if [ $1 == shortref ] ; then
    ##use first 40k bases of genome and 50k reads
    head -n 200000 $datadir/fastqs/$prefix.fq > $datadir/fastqs/$prefix.sub50k.fastq
    mv $datadir/bams $datadir/align
    sbatch --output=$datadir/batch_logs/align_sub50k_genome40k.out --job-name=subdcmalign ./align.scr $datadir/fastqs/$prefix.sub50k.fastq $datadir/../er2796_40k.fasta
fi

if [ $1 == subgenomeeventalign ] ; then
    mkdir -p $datadir/eventalign
    rm -f $datadir/eventalign/workspace.fofn
    touch $datadir/eventalign/workspace.fofn

    for i in $datadir/called/*/sequencing_summary.txt ;
    do
	echo $i >> $datadir/eventalign/workspace.fofn
    done

    seqtk seq -a $datadir/fastqs/$prefix.sub50k.fastq > $datadir/fastqs/$prefix.sub50k.fa
    nanopolish index -d $datadir/raw -f $datadir/eventalign/workspace.fofn $datadir/fastqs/$prefix.sub50k.fa
    
    
    nanopolish eventalign \
	       --scale-events \
	       -t 36 \
	       -r $datadir/fastqs/$prefix.sub50k.fa \
	       -b $datadir/align/$prefix.sub50k.sorted.bam \
	       -g $datadir/../er2796_40k.fasta > $datadir/eventalign/$prefix.sub50k.genome40k.eventalign.tsv
fi

if [ $1 == sub120kreads ] ; then
    mkdir -p $datadir/eventalign
    rm -f $datadir/eventalign/workspace.fofn
    touch $datadir/eventalign/workspace.fofn

    for i in $datadir/called/*/sequencing_summary.txt ;
    do
	echo $i >> $datadir/eventalign/workspace.fofn
    done

    head -n 480000 $datadir/fastqs/$prefix.fq > $datadir/fastqs/$prefix.sub120k.fastq
    seqtk seq -a $datadir/fastqs/$prefix.sub120k.fastq > $datadir/fastqs/$prefix.sub120k.fa
    nanopolish index -d $datadir/raw -f $datadir/eventalign/workspace.fofn $datadir/fastqs/$prefix.sub120k.fa
    
    ml samtools
    minimap2 -a -x map-ont -t 36 $datadir/../er2796.$genome.fasta $datadir/fastqs/$prefix.sub120k.fastq | samtools view -b | samtools sort -o $datadir/align/$prefix.sub120k.sorted.bam -T $datadir/align/$prefix.reads.tmp -
    samtools index $datadir/align/$prefix.sub120k.sorted.bam
    
    nanopolish eventalign \
	       --scale-events \
	       -t 36 \
	       -r $datadir/fastqs/$prefix.sub120k.fa \
	       -b $datadir/align/$prefix.sub120k.sorted.bam \
	       -g $datadir/../er2796_40k.fasta > $datadir/eventalign/$prefix.sub120k.$genome.eventalign.tsv
fi

if [ $1 == subgenome10k ] ; then
    mkdir -p $datadir/eventalign
    rm -f $datadir/eventalign/workspace.fofn
    touch $datadir/eventalign/workspace.fofn

    for i in $datadir/called/*/sequencing_summary.txt ;
    do
	echo $i >> $datadir/eventalign/workspace.fofn
    done

    seqtk seq -a $datadir/fastqs/$prefix.fq > $datadir/fastqs/$prefix.fa
    nanopolish index -d $datadir/raw -f $datadir/eventalign/workspace.fofn $datadir/fastqs/$prefix.fa

    for i in $datadir/../refs/er2796_2930000_2940000.fasta ;
    do
	refpre=`basename $i .fasta`
	ml samtools
	minimap2 -a -x map-ont -t 36 $i $datadir/fastqs/$prefix.fq | samtools view -b | samtools sort -o $datadir/align/$prefix.$refpre.sorted.bam -T $datadir/align/$prefix.$refpre.reads.tmp -
	samtools index $datadir/align/$prefix.$refpre.sorted.bam
	
	nanopolish eventalign \
		   --scale-events \
		   -t 36 \
		   -r $datadir/fastqs/$prefix.fa \
		   -b $datadir/align/$prefix.$refpre.sorted.bam \
		   -g $datadir/../er2796_40k.fasta > $datadir/eventalign/$prefix.$refpre.eventalign.tsv
    done
fi
