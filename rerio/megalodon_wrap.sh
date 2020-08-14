#!/bin/bash


if [ $1 == test ] ; then
    bash ~/Code/yfan_meth/rerio/megalodon.sh neb17 Y GANTC 1 rerio_allmod
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
fi
