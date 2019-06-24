library(ggplot2)

timefile='/work-zfs/mschatz1/cpowgs/analysis/meth/neb/eventalign/171019_neb19.dwells.csv'

time=read.table(file=timefile, header=F, sep=',', stringsAsFactors=FALSE)
colnames(time)=c('chrom', 'read', 'time')


outfile='/home-4/yfan7@jhu.edu/scratch/plots/meth/damtime.pdf'
pdf(outfile
