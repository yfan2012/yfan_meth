#!/bin/bash

srcdir=~/software/timp_nanopolish
npscript=~/Code/methylation-analysis
datadir=/scratch/groups/mschatz1/cpowgs
ref=$datadir/methref/ecoli_er2796.fasta



ml python/2.7.10

prefix=`echo $1 | rev | cut -d '/' -f 1 | rev`
enzyme=$2
echo $prefix
echo $enzyme

echo METH REF
##methylate the reference genome
refpre=`echo ${ref%.fasta}`
python $npscript/methylate_reference.py --recognition $enzyme $ref > $refpre.$enzyme.fasta


ml gcc/5.1.0
##setup
if [ "$2" != "unmeth" ] ; then
    echo SETUP
    for i in $srcdir/etc/r9-models/r9.4_450bps.nucleotide*6mer* ;
    do
	file=`echo $i | rev | cut -d '/' -f 1 | rev`
	python $npscript/expand_model_alphabet.py --alphabet $enzyme $i > $srcdir/etc/r9-models/input.$enzyme.$file
    done
    
    in_models=$datadir/$1/models/input_models.fofn
    ls $srcdir/etc/r9-models/input.$enzyme*.model | tr " " "\n" > $in_models
else
    $in_models
fi



##train
echo TRAINING
fq=$datadir/$1/fastqs/$1.all.fq
mkdir -p $datadir/$1/models

cd $datadir/$1/models
$srcdir/nanopolish methyltrain -v -t 48 \
		   --train-kmers all \
		   --out-fofn $datadir/$1/models/$1.$enzyme.fofn \
		   --out-suffix $1.$enzyme.model \
		   -m $in_models \
		   --reads $fq \
		   --bam $datadir/$1/bams/$1.$enzyme.all.sorted.bam \
		   --genome $refpre.$enzyme.fasta 

