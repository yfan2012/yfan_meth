library(ggplot2)
library(tidyverse)
library(doParallel)
library(foreach)
cl=makeCluster(10)
registerDoParallel(cl, cores=10)

plotdir='~/Dropbox/timplab_data/dnamod/plots/'

datadir='/dilithium/Data/Nanopore/Analysis/171120_ecoliStandards/methcall/'
dammeth=paste0(datadir,'dam/subset/171122_ecolidam.0.meth.dam.tsv')
damunmeth=paste0(datadir,'dam/subset/171122_ecoliUM.0.meth.dam.tsv')
cgmeth=paste0(datadir, 'cpg/subset/171124_ecoliCpG.0.meth.cpg.tsv')
cgunmeth=paste0(datadir,'cpg/subset/171122_ecoliUM.0.meth.cpg.tsv')
gcmeth=paste0(datadir,'gpc/subset/171120_ecoliGpC.0.meth.gpc.tsv.subset')
gcunmeth=paste0(datadir,'gpc/subset/171122_ecoliUM.0.meth.gpc.tsv.subset')

getconfusion <- function(methfile, unmethfile) {
    meth=read_tsv(methfile) %>%
        mutate(meth=TRUE)
    unmeth=read_tsv(unmethfile) %>%
        mutate(meth=FALSE)
    all=rbind(meth, unmeth)
    
    thresholds=seq(min(all$log_lik_ratio), max(all$log_lik_ratio), .1)
    pos=sum(all$meth)
    neg=sum(!all$meth)
    
    confusions=foreach(i=1:length(thresholds), .combine=rbind) %dopar% {
        thresh=thresholds[i]
        
        call=all$log_lik_ratio > thresh
        tp=call & all$meth
        fp=call & !all$meth
    
        tpr=sum(tp)/pos
        fpr=sum(fp)/neg
        
        conf=data.frame(thresh=thresh, tpr=tpr, fpr=fpr)
        return(conf)
    }
}

damconf=getconfusion(dammeth, damunmeth)
damconf$motif='dam'
cpgconf=getconfusion(cgmeth, cgunmeth)
cpgconf$motif='cpg'
gpcconf=getconfusion(gcmeth, gcunmeth)
gpcconf$motif='gpc'

allconf=rbind(damconf,cpgconf, gpcconf)

pdf(paste0(plotdir,'rocs.pdf'), width=11, height=8.5)
ggplot(allconf, aes(x=fpr, y=tpr, colour=motif)) +
    geom_line() +
    ggtitle('ROC meth calling') +
    xlab('False Positive Rate') +
    ylab('True Positive Rate') +
    theme_bw()
dev.off()
