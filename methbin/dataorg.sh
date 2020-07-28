#!/bin/bash

rawdir=/dilithium/Data/Nanopore/oxford
datadir=/mithril/Data/Nanopore/projects/methbin


if [ $1 == make_multi ] ; then
    ##untar the raw data, make the multifast5s, and delete the single fast5s
    mkdir -p $datadir/multiraw

    for i in $rawdir/*neb*/*neb*.t* ;
    do
	prefix=`basename $i | cut -d . -f 1 | cut -d _ -f 2`
	echo untaring $i
	
	mkdir -p $datadir/multiraw/$prefix
	mkdir -p $datadir/multiraw/$prefix/raw
	tar -xzf $i -C $datadir/multiraw/$prefix/raw

	echo making multi
	single_to_multi_fast5 -t 36 -i $datadir/multiraw/$prefix/raw -s $datadir/multiraw/$prefix -f $prefix --recursive

	echo deleting $datadir/multiraw/$prefix/raw
	rm -rf $datadir/multiraw/$prefix/raw
	
    done
fi

if [ $1 == renametest ] ; then
    mv $datadir/multiraw/neb2mod $datadir/multiraw/neb2
    rename 's/neb2mod/neb2/' $datadir/multiraw/neb2/*.fast5
fi

if [ $1 == rename ] ; then
    mv $datadir/multiraw/neb5mod $datadir/multiraw/neb5
    rename 's/neb5mod/neb5/' $datadir/multiraw/neb5/*.fast5
    
    mv $datadir/multiraw/nebmod6 $datadir/multiraw/neb6
    rename 's/nebmod6/neb6/' $datadir/multiraw/neb6/*.fast5

    mv $datadir/multiraw/nebmod9 $datadir/multiraw/neb9
    rename 's/nebmod9/neb9/' $datadir/multiraw/neb9/*.fast5
    
    mv $datadir/multiraw/neb12mod10 $datadir/multiraw/neb12
    rename 's/neb12mod10/neb12/' $datadir/multiraw/neb12/*.fast5
    
    mv $datadir/multiraw/neb12mod $datadir/multiraw/neb14
    rename 's/neb12mod/neb14/' $datadir/multiraw/neb14/*.fast5

    mv $datadir/multiraw/neb16mod $datadir/multiraw/neb16
    rename 's/neb16mod/neb16/' $datadir/multiraw/neb16/*.fast5
    
    mv $datadir/multiraw/neb $datadir/multiraw/nebdcm
    rename 's/neb/nebdcm/' $datadir/multiraw/nebdcm/*.fast5
fi

if [ $1 == renamedcm ] ; then
    rename 's/neb_/nebdcm_/' $datadir/multiraw/nebdcm/*.fast5
fi

