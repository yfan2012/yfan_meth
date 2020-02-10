#!/bin/bash

##analysis for using eventaligh to classify reads independently of each other
##idea is to aggregate kmers over a read

datadir=/uru/Data/Nanopore/projects/methbin
gdna="neb11 neb12 neb13 neb14 neb15 neb16 neb17 neb19 nebdcm"
plasmids="neb1 neb2 neb3 neb4 neb5 neb6 neb9 neb10"
virus="neb8"

if [ $1 == aggregate_test ] ; then

    for i in neb12;
    do
	python3 ~/Code/methylation/methbin/aggregate_events.py \
		-r $datadir/reference/er2796.fa \
		-c $datadir/eventalign_collapsed/neb11/neb11_eventalign_collapse.tsv \
		-d ~/software/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model \
		-m CCGG \
		-o $datadir/eventalign_collapsed/$i/neb11.positionpvals.tsv \
		-p $datadir/eventalign_collapsed/$i/neb11.readpvals.tsv \
		-t 12 >& $datadir/eventalign_collapsed/$i/weirdreads_neb11.txt
	python3 ~/Code/methylation/methbin/aggregate_events.py \
		-r $datadir/reference/er2796.fa \
		-c $datadir/eventalign_collapsed/$i/${i}_eventalign_collapse.tsv \
		-d ~/software/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model \
		-m CCGG \
		-o $datadir/eventalign_collapsed/$i/$i.positionpvals.tsv \
		-p $datadir/eventalign_collapsed/$i/$i.readpvals.tsv \
		-t 12
    done
fi

if [ $1 == aggregate_test2 ] ; then
    for i in neb2 ;
    do
	python3 ~/Code/methylation/methbin/aggregate_events.py \
		-r $datadir/reference/allsamps.fa \
		-c $datadir/eventalign_collapsed/$i/${i}_eventalign_collapse.tsv \
		-d ~/software/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model \
		-m CCGG \
		-o $datadir/eventalign_collapsed/$i/$i.positionpvals.tsv \
		-p $datadir/eventalign_collapsed/$i/$i.readpvals.tsv \
		-t 36
	
	python3 ~/Code/methylation/methbin/aggregate_events.py \
		-r $datadir/reference/allsamps.fa \
		-c $datadir/eventalign_collapsed/neb1/neb1_eventalign_collapse.tsv \
		-d ~/software/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model \
		-m CCGG \
		-o $datadir/eventalign_collapsed/$i/neb1.positionpvals.tsv \
		-p $datadir/eventalign_collapsed/$i/neb1.readpvals.tsv \
		-t 36
    done
fi


if [ $1 == aggregate_all ] ; then
    nebsamps=$datadir/reference/nebsamps.tsv
    
    while read samp ; do
	name=`echo $samp | cut -d ' ' -f 1`
	motif=` echo $samp | cut -d ' ' -f 2 `
	modtype=` echo $samp | cut -d ' ' -f 3 `
	canonical=` echo $samp | cut -d ' ' -f 4 `
	mod=` echo $samp | cut -d ' ' -f 5 `
	pos=` echo $samp | cut -d ' ' -f 6 `
	dna=` echo $samp | cut -d ' ' -f 7 `
	
	if [ $motif != 'none' ] && [ ! -f $datadir/eventalign_collapsed/$name/$name.readpvals.tsv ] ; then
	    python3 ~/Code/methylation/methbin/aggregate_events.py \
		    -r $datadir/reference/allsamps.fa \
		    -c $datadir/eventalign_collapsed/$name/${name}_eventalign_collapse.tsv \
		    -d ~/software/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model \
		    -m $motif \
		    -o $datadir/eventalign_collapsed/$name/$name.positionpvals.tsv \
		    -p $datadir/eventalign_collapsed/$name/$name.readpvals.tsv \
		    -t 12
	    if [ $dna == gDNA ] ; then
		python3 ~/Code/methylation/methbin/aggregate_events.py \
			-r $datadir/reference/allsamps.fa \
			-c $datadir/eventalign_collapsed/neb1/neb11_eventalign_collapse.tsv \
			-d ~/software/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model \
			-m $motif \
			-o $datadir/eventalign_collapsed/$name/neb11.positionpvals.tsv \
			-p $datadir/eventalign_collapsed/$name/neb11.readpvals.tsv \
			-t 12
	    elif [ $dna == plas ] ; then
		python3 ~/Code/methylation/methbin/aggregate_events.py \
			-r $datadir/reference/allsamps.fa \
			-c $datadir/eventalign_collapsed/neb1/neb1_eventalign_collapse.tsv \
			-d ~/software/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model \
			-m $motif \
			-o $datadir/eventalign_collapsed/$name/neb1.positionpvals.tsv \
			-p $datadir/eventalign_collapsed/$name/neb1.readpvals.tsv \
			-t 12
	    else
		echo skipping $name
	    fi
	fi
    done < $nebsamps
fi

    
