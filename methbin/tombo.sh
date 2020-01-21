#!/bin/bash

rawdir=/dilithium/Data/Nanopore/oxford
datadir=/uru/Data/Nanopore/projects/methbin

if [ $1 == untar ] ; then
    mkdir -p $datadir/raw

    for i in $rawdir/*neb*/*neb*.t* ;
    do
	prefix=`basename $i | cut -d . -f 1 | cut -d _ -f 2`
	echo untaring $i
	
	mkdir -p $datadir/raw/$prefix
	tar -xzf $i -C $datadir/raw/$prefix
    done
fi

if [ $1 == rename ] ; then
    mv $datadir/raw/neb2mod $datadir/raw/neb2
    mv $datadir/raw/neb5mod $datadir/raw/neb5
    mv $datadir/raw/nebmod6 $datadir/raw/neb6
    mv $datadir/raw/nebmod9 $datadir/raw/neb9
    mv $datadir/raw/neb12mod10 $datadir/raw/neb12
    mv $datadir/raw/neb12mod $datadir/raw/neb14
    mv $datadir/raw/neb16mod $datadir/raw/neb16
    mv $datadir/raw/neb $datadir/raw/nebdcm
fi

if [ $1 == flatraw ] ; then
    for i in $datadir/raw/* ;
    do
	echo $i
	for j in $i/* ;
	do
	    echo $j
	    for k in $j/fast5/* ;
	    do
		mv $k/*.fast5 $i/
		rmdir $k
	    done
	    rmdir $j/fast5
	    rmdir $j
	done
    done
fi

if [ $1 == preprocess ] ; then
    for i in $datadir/raw/* ;
    do
	prefix=`echo $i | rev | cut -d / -f 1 | rev`
	echo $prefix
	tombo preprocess annotate_raw_with_fastqs \
	      --fast5-basedir $i \
	      --fastq-filenames $datadir/fastqs/$prefix/${prefix}_40k.fq
    done
fi
