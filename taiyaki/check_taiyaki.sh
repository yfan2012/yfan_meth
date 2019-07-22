#!/bin/bash

##on whack
if [ $1 == write_probs ] ; then
   ##write cond probs to txt file for roc plotting in r
   for i in cpg dam ;
   do
       if [ $i == cpg ] ; then
	   motif=CG
	   pos=0
       elif [ $i == dam ] ; then
	   motif=GATC
	   pos=1
       fi
       python ~/Code/methylation/taiyaki/write_probs.py -i ~/data/test_$i/$i/basecall/basecalls.hdf5 -o ~/data/test_$i/${i}_condprobs.txt -m $motif -p $pos -f ~/data/test_$i/$i/basecall/basecalls.fa
       python ~/Code/methylation/taiyaki/write_probs.py -i ~/data/test_$i/unmeth/basecall/basecalls.hdf5 -o ~/data/test_$i/unmeth_condprobs.txt -m $motif -p $pos -f ~/data/test_$i/unmeth/basecall/basecalls.fa
   done
fi
   
