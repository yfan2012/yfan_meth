library(ggplot2)
library(plyr)
library(tidyverse)
library(doParallel)
cl=makeCluster(10)
registerDoParallel(cl, cores=10)

datadir='/dilithium/Data/Nanopore/dam/'
outdir='~/Dropbox/yfan/methylation/neb_meth/time/'

damraw=paste0(datadir,'ecolidam.dam.eventalign.tsv')
umraw=paste0(datadir, 'ecoliUM.dam.eventalign.tsv')

dwell_times <- function(tsvpath) {
    ##from isac tsv, collapse events into single dwell times for each kmer in each read
    cols=c('contig', 'position', 'reference_kmer', 'read_index', 'strand', 'event_idx', 'eventmean', 'event_std', 'length', 'model_kmer', 'model_mean', 'model_std', 'std_level')
    raw=read_tsv(tsvpath, col_names=cols)
    dwell=ddply(raw, .(position, read_index), .parallel=TRUE, function(x) {
        duration=sum(x$length)
        return(data.frame(refmer=x$reference_kmer[1], time=duration))
    })
    return(dwell)
}

position_times <- function(dwell) {
    ##from collapsed dwell times, get mean/sd dwell times for each position
    posdwell=ddply(dwell, .(position), .parallel=TRUE, function(x) {
        duration_mean=mean(x$time)
        duration_std=sd(x$time)
        return(data.frame(refmer=x$refmer[1], meantime=duration_mean, sdtime=duration_std))
    })
    return(posdwell)
}

kmer_times <- function(dwell) {
    kmerdwell=ddply(dwell, .(refmer), .parallel=TRUE, function(x) {
        duration_mean=mean(x$time)
        duration_std=sd(x$time)
        return(data.frame(refmer=x$refmer[1], meantime=duration_mean, sdtime=duration_std))
    })
    return(kmerdwell)  
}


damdwell=dwell_times(damraw)
dampos=position_times(damdwell)
damkmer=kmer_times(damdwell)

umdwell=dwell_times(umraw)
umpos=position_times(umdwell)
umkmer=kmer_times(umdwell)

damdwell$samp='dam'
umdwell$samp='um'
alldwell=rbind(damdwell, umdwell)

comparepos=foreach(i=1:dim(umpos)[1], .combine='rbind') %dopar% {
    pos=umpos$position[i]
    damline=dampos[dampos$position==pos,]
    diff=umpos$meantime[i]-damline$meantime
    if (dim(damline)[1]>0) {
        return(data.frame(pos=pos, refmer=umpos$refmer[i], dammean=damline$meantime[1], damsd=damline$sdtime[1], ummean=umpos$meantime[i], umsd=umpos$sdtime[i], diff=diff))
    }
}
comparepos$absdiff=abs(comparepos$diff)
ranked=head(arrange(comparepos, -absdiff))

comparekmer=foreach(i=1:dim(umkmer)[1], .combine='rbind') %dopar% {
    mer=umkmer$refmer[i]
    damline=damkmer[damkmer$refmer==mer,]
    diff=umkmer$meantime[i]-damline$meantime
    if (dim(damline)[1]>0) {
        return(data.frame(mer=mer, dammean=damline$meantime[1], damsd=damline$sdtime[1], ummean=umkmer$meantime[i], umsd=umkmer$sdtime[i], diff=diff))
    }
}
comparekmer$absdiff=abs(comparekmer$diff)
rankedmer=head(arrange(comparekmer, -absdiff))


pdf(paste0(outdir, 'dam_ecolistandard_aggregate.pdf'), width=11, height=8.5)
##dam dwell time dist and um dwell time dist
plot(ggplot(alldwell, aes(x=time)) +
     geom_density(aes(fill=samp, alpha=.4)) +
     xlab('Dwell Time (s)') +
     ggtitle('dam vs normal (all)') +
     theme_bw())
##hist of mean difference in dwell time by position
plot(ggplot(comparepos, aes(x=diff)) +
     geom_density(fill='lightblue', aes(alpha=.4)) +
     ggtitle('Difference in mean dwell time aggregated over position') +
     xlim(-.00075, .00075) +
     xlab('Difference') +
     theme_bw())
##hist of mean difference in dwell time by kmer
plot(ggplot(comparekmer, aes(x=diff)) +
     geom_density(fill='lightblue', aes(alpha=.4))+
     ggtitle('Difference in mean dwell time aggregated over kmer') +
     xlim(-.00075, .00075) +
     xlab('Difference') +
     theme_bw())
dev.off()


##plot dwell distributions for each largest difference kmer
pdf(paste0(outdir, 'dam_largest_diffs_grid.pdf'), width=11, height=8.5)
merdam=damdwell[damdwell$refmer %in% rankedmer$mer,]
merum=umdwell[umdwell$refmer %in% rankedmer$mer,]
merall=rbind(merdam, merum)
plot(ggplot(merall, aes(x=time)) +
     geom_density(aes(fill=samp, alpha=.4)) +
     facet_grid(rows=vars(refmer)) +
     ggtitle('Dwell times for kmers ') +
     xlim(0,.025) +
     xlab('time') +
     theme_bw())

dev.off()
    
     
     
