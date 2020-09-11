#!/bin/bash


if [ $1 == test ] ; then
    bash ~/Code/yfan_meth/rerio/megalodon.sh nebdcm Z CCWGG 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh nebdcm Z CCWGG 1 megalodon_vanilla
fi

if [ $1 == 6mA ] ; then
    ##testing true motif, all 6mA 
    
    ##neb17
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb17 Y GANTC 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb17 Y GATC 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb17 Z CCWGG 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb17 Z GATC 3 rerio_allmod
    
    
    ##neb19
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb19 Y GATC 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb19 Y GANTC 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb19 Z CCWGG 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb19 Z GATC 3 rerio_allmod

fi

if [ $1 == 5mC_all ] ; then
    ##testing true motif, all 5mC

    ##neb14
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb14 Y GATC 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb14 Z GATC 3 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb14 Z CCWGG 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb14 Z GCNGC 1 rerio_allmod

    ##neb15
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb15 Y GATC 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb15 Z GATC 3 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb15 Z CCWGG 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb15 Z GCNGC 1 rerio_allmod

    ##nebdcm
    bash ~/Code/yfan_meth/rerio/megalodon.sh nebdcm Y GATC 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh nebdcm Z GATC 3 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh nebdcm Z CCWGG 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh nebdcm Z GCNGC 1 rerio_allmod

    ##vanilla
    ##neb14
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb14 Y GATC 1 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb14 Z GATC 3 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb14 Z CCWGG 1 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb14 Z GCNGC 1 megalodon_vanilla

    ##neb15
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb15 Y GATC 1 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb15 Z GATC 3 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb15 Z CCWGG 1 megalodon_vanilla 
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb15 Z GCNGC 1 megalodon_vanilla

    ##nebdcm
    bash ~/Code/yfan_meth/rerio/megalodon.sh nebdcm Y GATC 1 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon.sh nebdcm Z GATC 3 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon.sh nebdcm Z CCWGG 1 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon.sh nebdcm Z GCNGC 1 megalodon_vanilla
fi


if [ $1 == 6mA_vanilla ] ; then
    ##testing true motif, all 6mA 
    
    ##neb17
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb17 Y GANTC 1 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb17 Y GATC 1 megalodon_vanilla 
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb17 Z CCWGG 1 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb17 Z GATC 3 megalodon_vanilla
    
    
    ##neb19
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb19 Y GATC 1 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb19 Y GANTC 1 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb19 Z CCWGG 1 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb19 Z GATC 3 megalodon_vanilla

fi



if [ $1 == control ] ; then
    ##all motifs for control samp

    ##4mc
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb11 X CCGG 0 rerio_allmod

    ##5mc
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb11 Z GATC 3 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb11 Z GCNGC 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb11 Z CCNGGC 5 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb11 Z CCWGG 1 rerio_allmod

    ##6ma
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb11 Y GANTC 1 rerio_allmod
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb11 Y GATC 1 rerio_allmod


    ##4mc
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb11 X CCGG 0 megalodon_vanilla

    ##5mc
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb11 Z GATC 3 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb11 Z GCNGC 1 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb11 Z CCNGGC 5 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb11 Z CCWGG 1 megalodon_vanilla

    ##6ma
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb11 Y GANTC 1 megalodon_vanilla
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb11 Y GATC 1 megalodon_vanilla

fi
