#!/bin/bash

datadir=/work-zfs/mschatz1/cpowgs/analysis/meth/neb/eventalign
ref=/scratch/groups/mschatz1/cpowgs/meth/refs/er2796.fasta

if [ $1 == getnorm ] ; then
    ##querying is different between readidx and posidx because I forced conversion to bed file on this one    
    posidx=$datadir/171019_neb19.eventalign.sorted.bed.gz
    readidx=$datadir/171019_neb19.eventalign.readsorted.tsv.gz


    ##find all positions that are 10 bases upstream of a gatc site
    dampos=$datadir/dampos.txt
    ustream=$datadir/dampos_ustream.txt
    python ~/Code/methylation/neb/etc/find_motif.py -r $ref -m GATC > $dampos
    awk -v s=10 '{print $1-s}' $dampos | head > $ustream



    
    norm=$datadir/dam_norm.csv
    echo read,avgdwell >> $norm
    while read i ; do
	tabpos=$datadir/tmp.$i.tabix.tsv
	end=$(($i + 1))

	##query each position and save out the tmp file
	tabix $posidx 'gi|730582171|gb|CP009644.1|':$i-$end &> $tabpos

	##check if each read is in the norm info doc
	cut -f 3 $tmppos | uniq > $datadir/uniq_reads.txt
	while read u ; do
	    if $u in
	
    done 
fi
