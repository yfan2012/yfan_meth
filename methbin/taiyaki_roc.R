library(tidyverse)
library(doParallel)
library(foreach)

cl=makeCluster(8)
registerDoParallel(cl, cores=8)

datadir='~/data/methbin/calls/'
samps=tibble(gDNA=c('neb12', 'neb13', 'neb14', 'neb15', 'neb16', 'neb17', 'neb19'), 
             plas=c('neb2', 'neb3', 'neb10', 'neb4', 'neb5', 'neb6', 'neb9'))


getconfusion <- function(methfile, unmethfile) {
    ##computes tpr and fpr from megalodon output file
    meth=read_tsv(methfile) %>%
        select(read_id, pos, mod_log_prob, can_log_prob) %>%
        mutate(ratio=can_log_prob/mod_log_prob) %>%
        mutate(meth=TRUE)
    unmeth=read_tsv(unmethfile) %>%
        select(read_id, pos, mod_log_prob, can_log_prob) %>%
        mutate(ratio=can_log_prob/mod_log_prob) %>%
        mutate(meth=FALSE)
    len=min(c(dim(meth)[1], dim(unmeth)[1]))
    all=rbind(meth[1:len,], unmeth[1:len,])

    thresholds=sort(all$ratio)[seq(1, length(all$ratio), 2000)]
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


classread <- function(methfile, unmethfile) {
    meth=read_tsv(methfile) %>%
        select(read_id, mod_log_prob, can_log_prob) %>%
        group_by(read_id) %>%
        summarise(mod=sum(mod_log_prob), can=sum(can_log_prob)) %>%
        mutate(ratio=can/mod) %>%
        mutate(meth=TRUE)
    unmeth=read_tsv(unmethfile) %>%
        select(read_id, mod_log_prob, can_log_prob) %>%
        group_by(read_id) %>%
        summarise(mod=sum(mod_log_prob), can=sum(can_log_prob)) %>%
        mutate(ratio=can/mod) %>%
        mutate(meth=FALSE)
    all=rbind(meth,unmeth)
    
    thresholds=sort(all$ratio)[seq(1, length(all$ratio), 100)]
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


##gets per base roc
allposconf=data.frame(matrix(ncol=3, nrow=0))
for (i in 1:dim(samps)[1]) {
    modfile=paste0(datadir, samps$gDNA[i], '/', samps$gDNA[i], '/per_read_modified_base_calls.txt')
    plsfile=paste0(datadir, samps$gDNA[i], '/', samps$plas[i], '/per_read_modified_base_calls.txt')
    umodfile=paste0(datadir, samps$gDNA[i], '/neb11/per_read_modified_base_calls.txt')
    uplsfile=paste0(datadir, samps$gDNA[i], '/neb1/per_read_modified_base_calls.txt')

    plasconf=getconfusion(plsfile, uplsfile)
    plasconf$samp='plasmid'
    

    modconf=getconfusion(modfile, umodfile)
    modconf$samp='gDNA'
    
    plotconf=rbind(modconf, plasconf)
    plotconf$mtase=samps$gDNA[i]

    plotfile=paste0('~/Dropbox/yfan/methylation/methbin/taiyaki/roc_', samps$gDNA[i], '.pdf')
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
    allposconf=rbind(allposconf, plotconf)
}    
pdf('~/Dropbox/yfan/methylation/methbin/taiyaki/roc_all.pdf')
ggplot(allposconf, aes(x=fpr, y=tpr, colour=mtase)) +
    geom_line(aes(linetype=samp)) +
    geom_abline(slope=1, intercept=0) +
    xlim(0,1) +
    ylim(0,1) +
    xlab('False Positive Rate') +
    ylab('True Positive Rate') +
    theme_bw()
dev.off()

   
##gets per read roc
allreadconf=data.frame(matrix(ncol=3, nrow=0))
for (i in 1:dim(samps)[1]) {
    modfile=paste0(datadir, samps$gDNA[i], '/', samps$gDNA[i], '/per_read_modified_base_calls.txt')
    plsfile=paste0(datadir, samps$gDNA[i], '/', samps$plas[i], '/per_read_modified_base_calls.txt')
    umodfile=paste0(datadir, samps$gDNA[i], '/neb11/per_read_modified_base_calls.txt')
    uplsfile=paste0(datadir, samps$gDNA[i], '/neb1/per_read_modified_base_calls.txt')

    plasconf=classread(plsfile, uplsfile)
    plasconf$samp='plasmid'
    

    modconf=classread(modfile, umodfile)
    modconf$samp='gDNA'

    plotconf=rbind(modconf, plasconf)
    plotconf$mtase=samps$gDNA[i]

    plotfile=paste0('~/Dropbox/yfan/methylation/methbin/taiyaki/readroc_', samps$gDNA[i], '.pdf')
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
    allreadconf=rbind(allreadconf, plotconf)
}    
pdf('~/Dropbox/yfan/methylation/methbin/taiyaki/readroc_all.pdf')
ggplot(allreadconf, aes(x=fpr, y=tpr, colour=mtase)) +
    geom_line(aes(linetype=samp)) +
    geom_abline(slope=1, intercept=0) +
    xlim(0,1) +
    ylim(0,1) +
    xlab('False Positive Rate') +
    ylab('True Positive Rate') +
    theme_bw()
dev.off()



allsamps=samps$gDNA[c(1:2, 4:7)]
allconf=data.frame(matrix(ncol=5, nrow=0))
for (i in 1:length(allsamps)) {
    unmodfile=paste0(datadir, 'neb14/', allsamps[i], '/per_read_modified_base_calls.txt')
    modfile=paste0(datadir, 'neb14/neb14/per_read_modified_base_calls.txt')

    conf=classread(modfile, unmodfile)
    conf$mtase=allsamps[i]
    conf$samp='gDNA'

    allconf=rbind(allconf,conf)
}

pdf('~/Dropbox/yfan/methylation/methbin/taiyaki/crossgatc_roc.pdf')
ggplot(allconf, aes(x=fpr, y=tpr, colour=mtase)) +
    geom_line(aes(linetype=samp)) +
    geom_abline(slope=1, intercept=0) +
    xlim(0,1) +
    ylim(0,1) +
    xlab('False Positive Rate') +
    ylab('True Positive Rate') +
    theme_bw()
dev.off()

