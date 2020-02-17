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
	seq=` echo $samp | cut -d ' ' -f 8 `
		
	##if [ $name == neb13 ] && [ ! -f $datadir/eventalign_collapsed/$name/neb11.readpvals.tsv ] ; then
	if [ $name == neb10 ] ; then
	    echo $name $name
	    #python3 ~/Code/methylation/methbin/aggregate_events.py \
	#	    -r $datadir/reference/allsamps.fa \
	#	    -c $datadir/eventalign_collapsed/$name/${name}_eventalign_collapse.tsv \
	#	    -d ~/software/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model \
	#	    -m $motif \
	#	    -s $seq \
	#	    -o $datadir/eventalign_collapsed/$name/$name.positionpvals.tsv \
	#	    -p $datadir/eventalign_collapsed/$name/$name.readpvals.tsv \
       	#	    -t 36
	    if [ $dna == sfiojeofj ] ; then
		echo $name neb11
		python3 ~/Code/methylation/methbin/aggregate_events.py \
			-r $datadir/reference/allsamps.fa \
			-c $datadir/eventalign_collapsed/neb11/neb11_eventalign_collapse.tsv \
			-d ~/software/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model \
			-m $motif \
			-s "gi|730582171|gb|CP009644.1|" \
			-o $datadir/eventalign_collapsed/$name/neb11.positionpvals.tsv \
			-p $datadir/eventalign_collapsed/$name/neb11.readpvals.tsv \
			-t 36
	    elif [ $name == neb10 ] && [ ! -f $datadir/eventalign_collapsed/$name/neb1.readpvals.tsv ] ; then
		echo $name neb1
		python3 ~/Code/methylation/methbin/aggregate_events.py \
			-r $datadir/reference/allsamps.fa \
			-c $datadir/eventalign_collapsed/neb1/neb1_eventalign_collapse.tsv \
			-d ~/software/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model \
			-m $motif \
			-s pT7MSin395ORF667 \
			-o $datadir/eventalign_collapsed/$name/neb1.positionpvals.tsv \
			-p $datadir/eventalign_collapsed/$name/neb1.readpvals.tsv \
			-t 36
	    else
		echo skipping $name
	    fi
	fi
    done < $nebsamps
fi

    
if [ $1 == testneb10 ] ; then
    for i in neb10 ;
    do
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
