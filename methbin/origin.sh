#!/bin/bash

datadir=/uru/Data/Nanopore/projects/methbin
##just looking at the biggest shifts from the roc on the poster
gdna="neb12 neb14 neb15 neb19"


if [ $1 == find_oriC ] ; then
    ##already checked that there's only one line with 'origin' in it
    grep origin $datadir/reference/ecoli_k12.gff > $datadir/reference/ecoli_k12_origin.gff
    bedtools getfasta \
	     -fi $datadir/reference/ecoli_k12.fa \
	     -bed $datadir/reference/ecoli_k12_origin.gff \
	     -fo $datadir/reference/ecoli_k12_origin.fa
    blastn -query $datadir/reference/ecoli_k12_origin.fa \
	   -subject $datadir/reference/er2796.fa \
	   -outfmt 7 \
	   -out $datadir/reference/er2796_origin.tsv
fi

if [ $1 == find_allpvals ] ; then
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

	if [[ $gdna == *$name* ]] && [ $name != neb1 ] ; then
	    echo $name $name
	    mkdir -p $datadir/origin/$name
            python3 ~/Code/methylation/methbin/aggregate_events.py \
                    -r $datadir/reference/allsamps.fa \
                    -c $datadir/eventalign_collapsed/$name/${name}_eventalign_collapse.tsv \
                    -d ~/software/nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model \
                    -m $motif \
                    -o $datadir/origin/$name/$name.positionpvals.tsv \
                    -p $datadir/origin/$name/$name.readpvals.tsv \
                    -t 54
	fi
    done < $nebsamps
fi
