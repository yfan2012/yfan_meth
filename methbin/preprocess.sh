#!/bin/bash

##datadir=/uru/Data/Nanopore/projects/methbin
datadir=~/data/methbin
gdna="neb11 neb12 neb13 neb14 neb15 neb16 neb17 neb19 nebdcm"


if [ $1 == call ] ; then
    ##using whack:guppy_bcall docker
    ##docker run --runtime=nvidia --name yfan -i -t -v ~/data/methbin:/media yfan2012/whack:guppy_bcall /bin/bash
    mkdir -p /media/called
    for i in /media/multiraw_sub/* ;
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


