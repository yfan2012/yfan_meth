#!/bin/bash

datadir=/scratch/groups/mschatz1/cpowgs/meth

if [ $1 == get_fq ] ; then
    for i in $datadir/*neb* ;
    do
	prefix=`echo $i | rev | cut -d / -f 1 | rev`
	scp $i/fastqs/$prefix.fq smaug:/atium/Data/Nanopore/methyl/neb/fastqs
    done
fi
