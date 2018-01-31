#!/bin/bash

##run like this: bash call_cre.sh /work-zfs/mschatz1/cpowgs/171011_170_KLPN strip 

##Set up the run name and input options
if [ -d "$1" ]; then
    prefix=` echo $1 | rev |  cut -d '/' -f1 | rev `
    if [ -z $prefix ]; then
	prefix=` echo $1 | rev |  cut -d '/' -f2 | rev `
    fi
else
    echo Enter sample directory
    exit
fi


##set up dirs
raw=$1/raw
srcpath=~/Code/methylation/scripts
outdir=/scratch/groups/mschatz1/cpowgs/$prefix/call_logs
mkdir -p  $outdir


##figure out array size
numdirs=`find $raw/* -maxdepth 0 -type d | wc -l `
dummy=1
maxdir=`expr $numdirs - $dummy`
echo maxdir is $maxdir
echo $prefix


sbatch --array=0-$maxdir --job-name=$prefix --output=$outdir/$prefix.%A_%a.out $srcpath/call.scr $1

