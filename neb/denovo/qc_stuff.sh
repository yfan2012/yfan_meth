#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/meth
gDNA="171020_neb11 171005_neb12 171012_neb13 170906_neb14 171020_neb15 171003_neb16 171019_neb17 171019_neb19 180628_neb_dcm"
ref=$datadir/refs/er2796.fasta

##this is forward looking to make it easy to switch between gDNA and plasmids using $2
samps=$gDNA

if [ $1 == cat_workspace ] ; then
    for i in $samps ;
    do
	touch $datadir/$i/seq_sum.txt
	rm $datadir/$i/seq_sum.txt
	touch $datadir/$i/seq_sum.txt
	cat $datadir/$i/called/*/sequencing_summary.txt >> $datadir/$i/seq_sum.txt
	sed -i '/^filename/ d' $datadir/$i/seq_sum.txt
    done
fi


      
	     
   
