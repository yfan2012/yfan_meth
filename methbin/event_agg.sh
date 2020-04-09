#!/bin/bash

##analysis for using eventaligh to classify reads independently of each other
##idea is to aggregate kmers over a read

datadir=/uru/Data/Nanopore/projects/methbin
gdna="neb11 neb12 neb13 neb14 neb15 neb16 neb17 neb19 nebdcm"
plasmids="neb1 neb2 neb3 neb4 neb5 neb6 neb9 neb10"
virus="neb8"

if [ $1 == aggregate_test ] ; then

    for i in neb19;
    do
	python3 ~/Code/methylation/methbin/aggregate_events.py \
		-r $datadir/reference/allsamps.fa \
		-c $datadir/eventalign_collapsed/$i/${i}_eventalign_collapse.tsv \
		-d ~/software/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model \
		-m GATC \
		-o $datadir/eventalign_collapsed/$i/$i.positionpvals.tsv \
		-p $datadir/eventalign_collapsed/$i/$i.readpvals.tsv \
		-l 3000 \
		-t 36
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
	seq=` echo $samp | cut -d ' ' -f 8 `

	if [ $motif != 'none' ] ; then
	##if [ $name == neb13 ] ; then
	    echo $name $name
	    python3 ~/Code/methylation/methbin/aggregate_events.py \
		    -r $datadir/reference/allsamps.fa \
		    -c $datadir/eventalign_collapsed/$name/${name}_eventalign_collapse.tsv \
		    -d ~/software/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model \
		    -m $motif \
		    -o $datadir/eventalign_collapsed/$name/$name.positionpvals.tsv \
		    -p $datadir/eventalign_collapsed/$name/$name.readpvals.tsv \
		    -l 3000 \
       		    -t 36
	    if [ $dna == gDNA ] ; then
		echo $name neb11
		python3 ~/Code/methylation/methbin/aggregate_events.py \
			-r $datadir/reference/allsamps.fa \
			-c $datadir/eventalign_collapsed/neb11/neb11_eventalign_collapse.tsv \
			-d ~/software/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model \
			-m $motif \
			-o $datadir/eventalign_collapsed/$name/neb11.positionpvals.tsv \
			-p $datadir/eventalign_collapsed/$name/neb11.readpvals.tsv \
			-l 3000 \
			-t 36
	    elif [ $dna == plas ] ; then
		echo $name neb1
		python3 ~/Code/methylation/methbin/aggregate_events.py \
			-r $datadir/reference/allsamps.fa \
			-c $datadir/eventalign_collapsed/neb1/neb1_eventalign_collapse.tsv \
			-d ~/software/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model \
			-m $motif \
			-o $datadir/eventalign_collapsed/$name/neb1.positionpvals.tsv \
			-p $datadir/eventalign_collapsed/$name/neb1.readpvals.tsv \
			-l 3000 \
			-t 36
	    fi
	fi
    done < $nebsamps
fi

    
if [ $1 == testneb13 ] ; then
    for i in neb13 neb11;
    do
	echo neb13 $i
	python3 ~/Code/methylation/methbin/aggregate_events.py \
		-r $datadir/reference/allsamps.fa \
		-c $datadir/eventalign_collapsed/$i/${i}_eventalign_collapse.tsv \
		-d ~/software/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model \
		-m CTGCAG \
		-o $datadir/eventalign_collapsed/neb13/$i.positionpvals.tsv \
		-p $datadir/eventalign_collapsed/neb13/$i.readpvals.tsv \
		-t 36
    done

    for i in neb3 neb1 ;
    do
	echo neb3 $i
	python3 ~/Code/methylation/methbin/aggregate_events.py \
		-r $datadir/reference/allsamps.fa \
		-c $datadir/eventalign_collapsed/$i/${i}_eventalign_collapse.tsv \
		-d ~/software/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model \
		-m CTGCAG \
		-o $datadir/eventalign_collapsed/neb3/$i.positionpvals.tsv \
		-p $datadir/eventalign_collapsed/neb3/$i.readpvals.tsv \
		-t 36
    done
fi   
