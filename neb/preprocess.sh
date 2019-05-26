#!/bin/bash

##Run through training pipeline
##structured analyis dir is first arg

npdir=~/software/timp_nanopolish
srcdir=~/Code/utils/marcc
datadir=/scratch/groups/mschatz1/cpowgs/meth
ref=${datadir}ref/ecoli_er2796.fasta 

gDNA="171020_neb11 171005_neb12 171012_neb13 170906_neb14 171020_neb15 171003_neb16 171019_neb17 171019_neb19 180628_neb_dcm"
plasmids="180104_neb1 170922_neb2 180104_neb3 180104_neb4 170922_neb5 171122_neb6 180104_neb8 171122_neb9"

if [ $2 == gDNA ] ; then
    samps=$gDNA
elif [ $2 == plasmid ] ; then
    samps=$plasmids
else
    echo gDNA or plasmids?
    exit 1
fi  


if [ $1 == clean ] ; then
    for i in $samps ;
    do
	(
	    if [ -d $datadir/$i ] ; then
	       rm -r $datadir/$i
	    fi
	)&
    done
fi

	     

##untar
if [ $1 == untar ] ; then
    for i in $samps ;
    do
	mkdir -p $datadir/$i
	mkdir -p $datadir/$i/raw
	mkdir -p $datadir/$i/batch_logs
	sbatch --output=$datadir/$i/batch_logs/$i.untar.out $srcdir/untar.scr $datadir/tarballs/$i.tgz $datadir/$i
    done
fi


##basecall
if [ $1 == call ] ; then
    for i in $samps ;
    do
	mkdir -p $datadir/$i/called
	bash $srcdir/call.sh $datadir/$i
    done
fi 


##gather fastqs
if [ $1 == fqs ] ; then
    for i in $samps ;
    do
	mkdir -p $datadir/$i/fastqs
	sbatch --output=$datadir/$i/batch_logs/$i.fqs.out $srcdir/fqs.scr $datadir/$i
    done
fi


##align
if [ $1 == align ] ; then
    for i in $samps ;
    do
	mkdir -p $datadir/$i/bams
	sbatch --output=$datadir/$i/batch_logs/$i.align.out $srcdir/align.scr $datadir/$i $ref 
    done
fi


##nanopolish index
if [ $1 == npindex ] ; then
    for i in $samps ;
    do 
	mkdir -p $datadir/$i/np_index
n	bash $srcdir/np_index.sh $datadir/$i
    done
fi


##cat index
if [ $1 == npindex_cat ] ; then
    for i in $samps ;
    do
	sbatch --output=$datadir/$i/batch_logs/$i.npindex_cat.out $srcdir/np_index_cat.scr $datadir/$i
    done
fi

	
##make alphabet/model files
if [ $1 == alphabet ] ; then
    sbatch --output=$srcdir/alphabet.%A_%a.out $srcdir/alphabet.scr
fi


##train
if [ $1 == train ] ; then
    ##neb19 was used to test this code, so excluding it from the list below
    ##unmeth is just too hard to deal with. Doing that one on it's own

    for i in 171020_neb11 171019_neb17 ;
    do
	if [ $i == 171005_neb12 ] ; then
	    enzyme=pspjdri
	elif [ $i == 171012_neb13 ] ; then
	    enzyme=sin395
	elif [ $i == 170906_neb14 ] ; then
	    enzyme=sin395
	elif [ $i == 171020_neb15 ] ; then
	    enzyme=fnu4h
	elif [ $i == 171003_neb16 ] ; then
	    enzyme=sdeaII
	elif [ $i == 171019_neb17 ] ; then
	    enzyme=hinfI
	elif [ $i == 171019_neb19 ] ; then
	    enzyme=dam
	elif [ $i == 171020_neb11 ] ; then
	    enzyme=unmeth
	fi
	
	##set up model files
	mkdir -p $datadir/$i/models
	in_models=$datadir/$i/models/input_models.fofn
	if [ $enzyme != unmeth ] ; then
	    ls $npdir/etc/r9-models/input.$enzyme*.model | tr " " "\n" > $in_models
	else
	    ls $npdir/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model | tr " " "\n" > $in_models
	fi

	
	##sbatch --output=$datadir/$i/batch_logs/$i.train.out $srcdir/meth_train.scr $datadir/$i $enzyme $ref
	bash $srcdir/meth_train.scr $datadir/$i $enzyme $ref
	
    done
fi

if [ $1 == assemble ] ; then
    ##if assembly doesn't exist, submit canu
    ml gnuplot
    for i in $samps ;
    do
	mkdir -p $datadir/$i/canu_assembly
	if [ ! -f $datadir/$i/canu_assembly/$i.contigs.fasta ] ; then
	    canu \
		-p $i -d $datadir/$i/canu_assembly \
		-gridOptions="--time=22:00:00 --account=mschatz1 --partition=parallel" \
		genomeSize=4.6m \
		-nanopore-raw $datadir/$i/fastqs/$i.fq
	fi
    done
fi
	
if [ $1 == fix_canu ] ; then
    ml gnuplot
    ##171012_neb13 ran out of time
    canu \
	-p 171012_neb13 -d $datadir/171012_neb13/canu_assembly \
	-gridOptions="--time=48:00:00 --account=mschatz1 --partition=parallel" \
	genomeSize=4.6m \
	-nanopore-raw $datadir/171012_neb13/fastqs/171012_neb13.fq

    ##171020_neb15 ran out of time
    canu \
	-p 171020_neb15 -d $datadir/171020_neb15/canu_assembly \
	-gridOptions="--time=48:00:00 --account=mschatz1 --partition=parallel" \
	genomeSize=4.6m \
	-nanopore-raw $datadir/171020_neb15/fastqs/171020_neb15.fq

    ##180628_neb_dcm had some weird mem error at unitigging. changed the unitigger.sh to submit
    canu \
	-p 180628_neb_dcm -d $datadir/180628_neb_dcm/canu_assembly \
	-gridOptions="--time=22:00:00 --account=mschatz1 --partition=parallel" \
	genomeSize=4.6m \
	-nanopore-raw $datadir/180628_neb_dcm/fastqs/180628_neb_dcm.fq

    ##171019_neb17 had some weird mem error at unitigging. changed the unitigger.sh to submit
    canu \
	-p 171019_neb17 -d $datadir/171019_neb17/canu_assembly \
	-gridOptions="--time=22:00:00 --account=mschatz1 --partition=parallel" \
	genomeSize=4.6m \
	-nanopore-raw $datadir/171019_neb17/fastqs/171019_neb17.fq
fi    
    
if [ $1 == fix_neb13 ] ; then
   ml gnuplot

   ##neb13 is still bitchy, so isolate reads higher than 2kb and redo
   rm -r $datadir/171012_neb13/canu_assembly
   python ~/Code/utils/fastq_long.py -i $datadir/171012_neb13/fastqs/171012_neb13.fq -o $datadir/171012_neb13/fastqs/171012_neb13_5kb.fq -l 5000

   canu \
       -p 171012_neb13 -d $datadir/171012_neb13/canu_assembly \
       -gridOptions="--time=48:00:00 --account=mschatz1 --partition=parallel" \
       genomeSize=4.6m \
       -nanopore-raw $datadir/171012_neb13/fastqs/171012_neb13_5kb.fq
fi
