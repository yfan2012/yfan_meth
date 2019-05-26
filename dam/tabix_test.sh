#!/bin/bash

datadir=/work-zfs/mschatz1/cpowgs/analysis/meth/neb/eventalign

if [ $1 == bgzip ] ; then
    bgzip $datadir/171019_neb19.eventalign.sorted.tsv
fi

if [ $1 == tabix ] ; then
    tabix -p bed $datadir/171019_neb19.eventalign.sorted.tsv.gz
fi

if [ $1 == re-sort ] ; then
    sort -n -S 100G -t $'\t' -k2,2 -k4,4 -T /scratch/groups/mschatz1/cpowgs/meth/tmp_gnu -o $datadir/$prefix.eventalign.sorted.tsv $datadir/171019_neb19.sorted.eventalign.tsv
fi

if [ $1 == query ] ; then
    datadir=/scratch/groups/mschatz1/cpowgs/meth/eventalign
    tabix $datadir/171019_neb19.eventalign.sorted.tsv.gz 'gi|730582171|gb|CP009644.1|':4458547-4458548 > $datadir/test.txt
fi
