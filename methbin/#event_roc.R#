library(tidyverse)
library(doParallel)
library(foreach)

cl=makeCluster(10)
registerDoParallel(cl, cores=10)

datadir='/uru/Data/Nanopore/projects/methbin/eventalign_collapsed/'
samps=tibble(gDNA=c('neb12', 'neb13', 'neb14', 'neb15', 'neb16', 'neb17', 'neb19'),
             plas=c('neb2', 'neb3', 'neb10', 'neb4', 'neb5', 'neb6', 'neb9'))


getconfusion <- function(methfile, unmethfile) {
    meth=read_tsv(methfile, col_names=c('read', 'ratio')) %>%
        mutate(meth=TRUE)
    unmeth=read_tsv(unmethfile, col_names=c('read', 'ratio')) %>%
        mutate(meth=FALSE)

    all=rbind(meth, unmeth)

    thresholds=sort(all$ratio)[seq(1, length(all$ratio), 20)]
    pos=sum(all$meth)
    neg=sum(!all$meth)

    confusions=foreach(i=1:length(thresholds), .combine=rbind) %dopar% {
        thresh=thresholds[i]

        call=all$ratio > thresh
        tp=call & all$meth
        fp=call & !all$meth

        tpr=sum(tp)/pos
        fpr=sum(fp)/neg

        conf=data.frame(thresh=thresh, tpr=tpr, fpr=fpr)
        return(conf)
    }
}




for (i in 1:dim(samps)[1]) {
    modfile=paste0(datadir, samps$gDNA[i], '/', samps$gDNA[i], '.readpvals.tsv')
    umodfile=paste0(datadir, samps$gDNA[i], '/neb11.readpvals.tsv')

    plsfile=paste0(datadir, samps$plas[i], '/', samps$plas[i], '.readpvals.tsv')
    uplsfile=paste0(datadir, samps$plas[i], '/neb1.readpvals.tsv')
    
    
    plasconf=getconfusion(plsfile, uplsfile)
    plasconf$samp='plasmid'


    modconf=getconfusion(modfile, umodfile)
    modconf$samp='gDNA'

    plotconf=rbind(modconf, plasconf)
    plotconf$mtase=samps$gDNA[i]

    plotfile=paste0('~/Dropbox/yfan/methylation/methbin/event/roc_', samps$gDNA[i], '.pdf')
    pdf(plotfile)
    ggplot(plotconf, aes(x=fpr, y=tpr, colour=mtase)) +
        geom_line(aes(linetype=samp)) +
        geom_abline(slope=1, intercept=0) +
        xlim(0,1) +
        ylim(0,1) +
        xlab('False Positive Rate') +
        ylab('True Positive Rate') +
        theme_bw()
    dev.off()
}

