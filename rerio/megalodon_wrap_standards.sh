#!/bin/bash

if [ $1 == meth ] ; then
    ##testing true motif, all 6mA 

    ##CpG
    bash ~/Code/yfan_meth/rerio/megalodon_standards.sh ecoli_CpG Z CG 0 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon_standards.sh ecoli_CpG Z CG 0 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon_standards.sh ecoli_CpG Z GC 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon_standards.sh ecoli_CpG Z GC 1 megalodon_vanilla
    
    ##GpC
    bash ~/Code/yfan_meth/rerio/megalodon_standards.sh ecoli_GpC Z GC 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon_standards.sh ecoli_GpC Z GC 1 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon_standards.sh ecoli_GpC Z CG 0 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon_standards.sh ecoli_GpC Z CG 0 megalodon_vanilla
    
    
fi

if [ $1 == control ] ; then
    ##all motifs for control samp
    bash ~/Code/yfan_meth/rerio/megalodon_standards.sh ecoli_Unmethylated Z GC 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon_standards.sh ecoli_Unmethylated Z GC 1 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon_standards.sh ecoli_Unmethylated Z CG 0 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon_standards.sh ecoli_Unmethylated Z CG 0 megalodon_vanilla
fi
