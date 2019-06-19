#!/bin/bash

datadir=/work-zfs/mschatz1/cpowgs/analysis/meth/neb/eventalign
##prefix=171019_neb19
prefix=171020_neb11

if [ $1 == bgzip ] ; then
    bgzip -@ 32 $datadir/$prefix.eventalign.sorted.bed
fi

if [ $1 == tabix ] ; then
    tabix -p bed -S 1 $datadir/$prefix.eventalign.sorted.bed.gz
fi

if [ $1 == bgzip_rsorted ] ; then
    bgzip -@ 32 $datadir/$prefix.eventalign.readsorted.tsv
fi

if [ $1 == tabix_rsorted ] ; then
    tabix -b 4 -e 4 -S 1 $datadir/$prefix.eventalign.readsorted.tsv.gz
fi

if [ $1 == re-sort ] ; then
    ##need to take col as int
    sort -n -S 100G -t $'\t' -k2,2 -k4,4 -T ~/work/methtmp -o $datadir/$prefix.eventalign.sorted.tsv $datadir/$prefix.sorted.eventalign.tsv
fi

if [ $1 == query ] ; then
    datadir=/scratch/groups/mschatz1/cpowgs/meth/eventalign
    tabix $datadir/$prefix.eventalign.sorted.tsv.gz 'gi|730582171|gb|CP009644.1|':4458547-4458548 > $datadir/test.txt
fi


if [ $1 == bedans ] ; then
    awk '{NF-=3}1' $datadir/$prefix.eventalign.sorted.tsv | awk '{$2=$2"\t"$2}1' OFS=$'\t' > $datadir/$prefix.eventalign.sorted.bed
fi

if [ $1 == readsort ] ; then
    ##stupidly got rid of the read sorted file
    mkdir -p /scratch/groups/mschatz1/cpowgs/meth/tmp_gnu
    sort -n -S 100G -t $'\t' -k4,4 -T ~/work/methtmp -o $datadir/$prefix.eventalign.readsorted.tsv $datadir/$prefix.eventalign.sorted.tsv
fi

