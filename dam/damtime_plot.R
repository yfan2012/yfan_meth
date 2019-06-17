library(ggplot2)

timefile='/work-zfs/mschatz1/cpowgs/analysis/meth/neb/eventalign/171019_neb19.dwells.csv'

time=read.table(file=timefile, header=F, sep=',', stringsAsFactors=FALSE)
