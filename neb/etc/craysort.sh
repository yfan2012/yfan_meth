#!/bin/bash

datadir=/work-zfs/mschatz1/cpowgs/analysis/meth/neb/eventalign
scratchdir=/scratch/groups/mschatz1/cpowgs/meth/eventalign

if [ $1 == divide ] ; then
    i=/work-zfs/mschatz1/cpowgs/analysis/meth/neb/eventalign/171019_neb19.eventalign.tsv
    ##i=/work-zfs/mschatz1/cpowgs/analysis/meth/neb/eventalign/test.txt
    parallel --tmpdir /scratch/groups/mschatz1/cpowgs/meth/tmp_cray -k -q --block 200M --pipepart -a $i awk '{ print $0 >> "/work-zfs/mschatz1/cpowgs/analysis/meth/neb/eventalign/frags/neb19_"$2".txt" 
        if(FNR % 10000000 ==0 ) {printf ("processed %d lines \n", FNR)  }}'
fi

if [ $1 == delete ] ; then
    for i in $datadir/frags/* ;
    do
	rm $i &
    done
fi

if [ $1 == gnusort ] ; then
    for i in $datadir/*.eventalign.tsv ;
    do
	prefix=`basename $i .eventalign.tsv | cut -d . -f 1`
	if [ ! -f $datadir/$prefix.sorted.eventalign.tsv ] ; then
	    echo $prefix
	    sort -S 100G -t $'\t' -k2,2 -k4,4 -T /scratch/groups/mschatz1/cpowgs/meth/tmp_gnu -o $datadir/$prefix.sorted.eventalign.tsv $i
	fi
    done
fi

if [ $1 == gnu_resort ] ; then
    for i in $datadir/*sorted.eventalign.tsv ;
    do
	prefix=`basename $i .eventalign.tsv | cut -d . -f 1`
	if [ ! -f $datadir/$prefix.eventalign.sorted.tsv ] ; then
	    echo $prefix
	    sort -S -n 100G -t $'\t' -k2,2 -k4,4 -T /scratch/groups/mschatz1/cpowgs/meth/tmp_gnu -o $datadir/$prefix.sorted.eventalign.tsv $i
	fi
    done
fi
