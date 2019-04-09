#!/bin/bash

datadir=/dilithium/Data/Nanopore/oxford
for i in neb1 neb2mod neb3 neb4 neb5mod nebmod6 neb8 nebmod9 neb11 neb12mod10 neb13 neb12mod neb15 neb16mod neb17 neb19 neb10 ;
do
    aws s3 cp $datadir/*$i/* s3://neb-meth/
done
