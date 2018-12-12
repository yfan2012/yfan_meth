#!/bin/bash

datadir=/scratch/groups/mschatz/cpowgs/meth
samp=171019_neb19

##try 5mC and 6mA detection with model
tombo detect_modifications alternative_model --fast5-basedirs $datadir/$samp \
      --per-read-statistics-basename --alternate-bases 5mC 6mA --statistics-file-basename sample_alt_model
