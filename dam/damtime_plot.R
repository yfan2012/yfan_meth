library(ggplot2)

meth='~/data/dwelltime/171019_neb19.dwells.csv'
umeth='~/data/dwelltime/171020_neb11.dwells.csv'

methtime=read.table(file=meth, header=F, sep=',', stringsAsFactors=FALSE)
colnames(methtime)=c('chrom', 'pos', 'time')
umethtime=read.table(file=umeth, header=F, sep=',', stringsAsFactors=FALSE)
colnames(umethtime)=c('chrom', 'pos', 'time')


outfile='~/Dropbox/yfan/methylation/neb_meth/time/dam_offset.pdf'
pdf(outfile, width=11, height=8.5)

##plot a few dwells, paired by position
##sample from umeth, since that's not been processed as long
locs=sample(umethtime$pos, 5)
locsmeth=methtime[methtime$pos %in% locs,]
locsmeth$status='meth'
locsumeth=umethtime[umethtime$pos %in% locs,]
locsumeth$status='umeth'
pospaired=rbind(locsmeth, locsumeth)
pospaired$pos=as.character(pospaired$pos)

plot=ggplot(pospaired, aes(x=pos, y=time, fill=status)) +
    geom_violin(scale='width') +
    ylim(0,5) +
    xlab('Genome Position') +
    ylab('Normalized Time (0 to 5)') +
    ggtitle('Dwell Times by Position') +
    theme_bw()
print(plot)


dev.off()
