#!/bin/bash

##wrap for easy running of megalodon with how I've set up the neb data

samp=$1 ##name of neb samp
mod=$2 ##symbol of mod
motif=$3 ##motif
pos=$4 ##mod location

ref=/dilithium/Data/Nanopore/Analysis/181010_ecoli_methylation_standards/Reference/ecoli_mg1655.fasta

if [ $5 == rerio_allmod ] ; then
	mkdir -p ~/data/rerio/$samp/megalodon/${samp}_${motif}_$mod

	megalodon \
	    ~/data/rerio/$samp/multiraw_sub/ \
	    --overwrite \
	    --guppy-server-path "/usr/bin/guppy_basecall_server" \
	    --guppy-params "-d /home/yfan/software/rerio/basecall_models/ --num_callers 5 --ipc_threads 6" \
	    --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg \
	    --reference $ref \
	    --mod-motif $mod $motif $pos \
	    --outputs basecalls mod_basecalls per_read_mods mods \
	    --output-directory ~/data/rerio/$samp/megalodon/${samp}_${motif}_$mod \
	    --write-mods-text \
	    --devices "cuda:0" \
	    --processes 36
fi

if [ $5 == megalodon_vanilla ] ; then
	mkdir -p ~/data/rerio/$samp/megalodon/${samp}_vanilla_${motif}_$mod

	megalodon \
	    ~/data/rerio/$samp/multiraw_sub/ \
	    --overwrite \
	    --guppy-server-path "/usr/bin/guppy_basecall_server" \
	    --reference $ref \
	    --mod-motif $mod $motif $pos \
	    --outputs basecalls mod_basecalls per_read_mods mods \
	    --output-directory ~/data/rerio/$samp/megalodon/${samp}_vanilla_${motif}_$mod \
	    --write-mods-text \
	    --devices "cuda:0" \
	    --processes 36

fi


if [ $5 == test_model ] ; then
    mkdir -p ~/data/rerio/$samp/megalodon/${samp}_vanilla_${motif}_$mod
    
    megalodon \
	~/data/rerio/$samp/multiraw_sub/ \
	--guppy-config dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg \
	--guppy-server-path "/usr/bin/guppy_basecall_server" \
	--reference $ref \
	--mod-motif $mod $motif $pos \
	--outputs basecalls mod_basecalls per_read_mods mods \
	--output-directory ~/data/rerio/$samp/megalodon/${samp}_vanilla_${motif}_$mod \
	--write-mods-text \
	--devices "cuda:0" \
	--processes 36 \
	--overwrite
fi

    

