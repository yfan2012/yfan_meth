#!/bin/bash

##on whack

if [ $1 == write_probs ] ; then
   ##write cond probs to txt file for roc plotting in r
   for i in cpg dam ;
   do
       (
       python ~/Code/methylation/taiyaki/write_probs.py -i ~/data/test_$i/$i/basecall/basecalls.hdf5 -o ~/data/test_$i/${i}_condprobs.txt &
       python ~/Code/methylation/taiyaki/write_probs.py -i ~/data/test_$i/unmeth/basecall/basecalls.hdf5 -o ~/data/test_$i/unmeth_condprobs.txt
       )&
   done
fi
   
