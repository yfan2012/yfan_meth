#!/bin/bash

srcdir=~/Code/methylation/neb
datadir=/scratch/groups/mschatz1/cpowgs/meth

if [ $1 == npidx ] ; then
    for i in 170906_neb14 171003_neb16 171005_neb12 171012_neb13 171019_neb17 171019_neb19 171020_neb11 171020_neb15 ;
    do
	bash denovo_meth.sh npidx $i
    done
fi

if [ $1 == shortref ] ; then
    for i in 170906_neb14 171003_neb16 171005_neb12 171012_neb13 171019_neb17 171019_neb19 171020_neb11 171020_neb15 ;
    do
	bash denovo_meth.sh shortref $i
    done
fi

if [ $1 == subgenomeeventalign ] ; then
    for i in 170906_neb14 171003_neb16 171005_neb12 171012_neb13 171019_neb17 171019_neb19 171020_neb11 171020_neb15 ;
    do
	bash denovo_meth.sh subgenomeeventalign $i
    done
fi

if [ $1 == sub120kreads ] ; then
    for i in 170906_neb14 171003_neb16 171005_neb12 171012_neb13 171019_neb17 171019_neb19 171020_neb11 171020_neb15 ;
    do
	bash denovo_meth.sh sub120kreads $i
    done
fi


if [ $1 == subgenomeregions ] ; then
    for i in `seq 1 1 500` ;
    do
	end=$(($i * 10000))
	start=$(($end - 10000))
	python ~/Code/utils/fasta_utils.py -i $datadir/refs/er2796.fasta -s $start -e $end -o $datadir/refs/er2796_${start}_${end}.fasta
    done
fi
	 
if [ $1 == subgenome10k ] ; then
    for i in 170906_neb14 171003_neb16 171005_neb12 171012_neb13 171019_neb17 171019_neb19 171020_neb11 171020_neb15 ;
    do
	bash denovo_meth.sh subgenome10k $i
    done
fi
