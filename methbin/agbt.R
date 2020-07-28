library(tidyverse)
library(doParallel)
library(foreach)
library(scattermore)
library(ggpubr)
library(purrr)
library(R.utils)

cl=makeCluster(10)
registerDoParallel(cl, cores=10)

datadir='/uru/Data/Nanopore/projects/methbin/eventalign_collapsed/'
samps=tibble(gDNA=c('neb12', 'neb13', 'neb14', 'neb15', 'neb16', 'neb17', 'neb19', 'nebdcm'),
             plas=c('neb2', 'neb3', 'neb10', 'neb4', 'neb5', 'neb6', 'neb9', 'none'))


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

getlengths  <- function(idxfile, pvalfile) {
    pvals=read_tsv(pvalfile, col_names=c('readname', 'ratio')) %>%
        mutate(read=substring(readname, 2))
    idxcols=c('refid', 'refstart', 'refend', 'read', 'kmers', 'dwell', 'N', 'mismatch', 'missing', 'offset', 'length')
    idxinfo=read_tsv(idxfile, col_names=idxcols) %>%
        filter(read %in% pvals$read)
    lenratio=inner_join(pvals, idxinfo, by=c('read'))

    lenratio=as_tibble(transform(lenratio, kmers=as.numeric(kmers)))
    return(lenratio)
}
    


##get some length/pval ratio correlations
plotfile=paste0('~/Dropbox/yfan/methylation/methbin/event/lenratio_corr2.pdf')
pdf(plotfile)
allratios=data.frame(matrix(ncol=4, nrow=0))
plots=vector(mode='list')
for (i in 1:dim(samps)[1]) {
    pvalfile=paste0(datadir, samps$gDNA[i], '/', samps$gDNA[i], '.readpvals.tsv')
    idxfile=paste0(datadir, samps$gDNA[i], '/', samps$gDNA[i], '_eventalign_collapse.tsv.idx')
    lenratio=getlengths(idxfile, pvalfile)
    lenratio$mtase=samps$gDNA[i]

    npvalfile=paste0(datadir, samps$gDNA[i], '/neb11.readpvals.tsv')
    nidxfile=paste0(datadir,'neb11/neb11_eventalign_collapse.tsv.idx')
    nlenratio=getlengths(nidxfile, npvalfile)
    nlenratio$mtase='neb11'

    ##numpoints=min(dim(lenratio)[1], dim(nlenratio)[1])
    numpoints=500
    alenratios=rbind(lenratio[1:numpoints,], nlenratio[1:numpoints,])
    
    plot=ggplot(alenratios, aes(x=kmers, y=ratio)) +
        geom_point(aes(colour=mtase, alpha=.2)) +
        ggtitle(samps$gDNA[i]) +
        xlab('Read Length') +
        ylab('ratio_logpval') +
        theme_bw()
    plots[[i]]=plot
    allratios=rbind(allratios, lenratio)
}
print(ggplot(allratios, aes(x=kmers, y=ratio)) +
      geom_point(aes(colour=mtase, alpha=.2)) +
      ggtitle('All') +
      xlab('Read Length') +
      ylab('ratio_logpval') +
      theme_bw())
dev.off()

plotfile=paste0('~/Dropbox/yfan/methylation/methbin/event/lenratio_correlations_all2.pdf')
pdf(plotfile, width=15, height=5)
ggarrange(plotlist=plots, ncol=4, nrow=2, align='v')
dev.off()



##get pval distributions
plots=vector(mode='list')
allratios=data.frame(matrix(ncol=4, nrow=0))
for (i in 1:dim(samps)[1]) {
    pvalfile=paste0(datadir, samps$gDNA[i], '/', samps$gDNA[i], '.readpvals.tsv')
    pvals=read_tsv(pvalfile, col_names=c('readname', 'ratio')) %>%
        mutate(mtase=samps$gDNA[i])
    
    npvalfile=paste0(datadir, samps$gDNA[i], '/neb11.readpvals.tsv')
    npvals=read_tsv(npvalfile, col_names=c('readname', 'ratio')) %>%
        mutate(mtase='neb11')

    numpoints=min(dim(pvals)[1], dim(npvals)[1])
    
    allpvals=rbind(pvals[1:numpoints,], npvals[1:numpoints,])
    
    plot=ggplot(allpvals, aes(x=ratio)) +
        geom_density(aes(colour=mtase, fill=mtase, alpha=.5)) +
        ggtitle(samps$gDNA[i]) +
        xlab('pval ratio') +
        xlim(0, 10) +
        theme_bw()
    plots[[i]]=plot
    allratios=rbind(allratios, lenratio)
}
plotfile=paste0('~/Dropbox/yfan/methylation/methbin/event/pvaldensities_group.pdf')
pdf(plotfile, width=18, height=5)
ggarrange(plotlist=plots, ncol=4, nrow=2, align="hv")
dev.off()

 
plots=vector(mode='list')
for (i in 1:(dim(samps)[1]-1)) {
    pvalfile=paste0(datadir, samps$plas[i], '/', samps$plas[i], '.readpvals.tsv')
    pvals=read_tsv(pvalfile, col_names=c('readname', 'ratio')) %>%
        mutate(mtase=samps$plas[i])
    
    npvalfile=paste0(datadir, samps$plas[i], '/neb1.readpvals.tsv')
    npvals=read_tsv(npvalfile, col_names=c('readname', 'ratio')) %>%
        mutate(mtase='neb1')

    numpoints=min(dim(pvals)[1], dim(npvals)[1], 2000)
    
    allpvals=rbind(pvals[1:numpoints,], npvals[1:numpoints,])
    
    plot=ggplot(allpvals, aes(x=ratio)) +
        geom_density(aes(colour=mtase, fill=mtase, alpha=.5)) +
        ggtitle(samps$plas[i]) +
        xlab('pval ratio') +
        xlim(0, 10) +
        theme_bw()
    plots[[i]]=plot
    allratios=rbind(allratios, lenratio)
}
plotfile=paste0('~/Dropbox/yfan/methylation/methbin/event/pvaldensities_group_plas.pdf')
pdf(plotfile, width=18, height=5)
ggarrange(plotlist=plots, ncol=4, nrow=2, align="hv")
dev.off()




samps$locs=c(8824, 0, 2018, 858, 0, 5728, 2018, 0)

pospvalfile=paste0('~/Dropbox/yfan/methylation/methbin/event/pospvals.pdf')
pdf(pospvalfile, width=15, height=4)
for (i in 1:dim(samps)[1]) {
    ##grab examples per genomic position
    name=samps$gDNA[i]
    loc=samps$locs[i]
    if ( loc != 0 ) {
        pospvalfile=paste0(datadir, name, '/', name, '.positionpvals.tsv')
        npospvalfile=paste0(datadir, name,'/neb11.positionpvals.tsv')
        pospval=read_tsv(pospvalfile, col_names=c('read', 'chr', 'pos', 'kmer', 'pval')) %>%
            filter( pos > loc-6 & pos < loc+6) %>%
            mutate(mtase=name) %>%
            mutate(meth=TRUE)
        npospval=read_tsv(npospvalfile, col_names=c('read', 'chr', 'pos', 'kmer', 'pval')) %>%
            filter( pos > loc-6 & pos < loc+6) %>%
            mutate(mtase=name) %>%
            mutate(meth=FALSE)
        alllocs=rbind(pospval, npospval) %>%
            mutate(label=paste0(as.character(pos), '\n', kmer))
        print(ggplot(alllocs, aes(x=label, y=pval, colour=meth, fill=meth)) +
              ##geom_violin(data=alllocs[alllocs$meth==TRUE,], alpha=.1) +
              ##geom_violin(data=alllocs[alllocs$meth==FALSE,], alpha=.1) +
              ##geom_violin(alpha=.25) +
              ##geom_dotplot(aes(x=label, group=meth), binaxis='y', stackdir='center',  dotsize=.5, alpha=.2) +
              geom_jitter(aes(x=meth), height=0, width=.1) +
              geom_violin(aes(x=meth), alpha=.2) +
              facet_grid(. ~ label) + 
              ggtitle(name) +
              xlab('Position') +
              ylab('pval') +
              theme_bw())
    }
}
dev.off()


##show gatc on one plot
gatcpvalfile=paste0('~/Dropbox/yfan/methylation/methbin/event/pospvals_gatc.pdf')
pdf(gatcpvalfile, width=15, height=4, useDingbats=FALSE)
loc=2018
mcfile=paste0(datadir, 'neb14/neb14.positionpvals.tsv')
mc=read_tsv(mcfile, col_names=c('read', 'chr', 'pos', 'kmer', 'pval')) %>%
    filter( pos > loc-6 & pos < loc+6) %>%
    mutate(mtase='5mC') %>%
    mutate(meth=TRUE)
mafile=paste0(datadir, 'neb19/neb19.positionpvals.tsv')
ma=read_tsv(mafile, col_names=c('read', 'chr', 'pos', 'kmer', 'pval')) %>%
    filter( pos > loc-6 & pos < loc+6) %>%
    mutate(mtase='6mA') %>%
    mutate(meth=TRUE)
nofile=paste0(datadir, 'neb14/neb11.positionpvals.tsv')
no=read_tsv(nofile, col_names=c('read', 'chr', 'pos', 'kmer', 'pval')) %>%
    filter( pos > loc-6 & pos < loc+6) %>%
    mutate(mtase='none') %>%
    mutate(meth=TRUE)
all=rbind(mc, ma, no) %>%
    mutate(label=paste0(as.character(pos), '\n', kmer))
colours=c('#F8766D', '#619CFF', '#00BA38')
print(ggplot(all, aes(x=label, y=pval, colour=mtase, fill=mtase)) +
      geom_jitter(aes(x=mtase), height=0, width=.1) +
      geom_violin(aes(x=mtase), alpha=.2) +
      scale_x_discrete(limits=c('none', '5mC', '6mA')) +
      scale_fill_manual(values=c('#1CBDC2','#2AB34B','#F3766E')) +
      scale_colour_manual(values=c('#1CBDC2','#2AB34B','#F3766E')) +
      facet_grid(. ~ label) + 
      ggtitle('GATC') +
      xlab('Position') +
      ylab('pval') +
      theme_bw())
dev.off()
