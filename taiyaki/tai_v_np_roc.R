library(ggplot2)
library(tidyverse)
library(doParallel)
library(foreach)
cl=makeCluster(10)
registerDoParallel(cl, cores=10)


##read in stuff
taidir='/kyber/Data/Nanopore/projects/taiyaki_yfan/'
cpgfile=paste0(taidir, 'test_cpg/cpg_condprobs.txt')
ucpgfile=paste0(taidir, 'test_cpg/unmeth_condprobs.txt')
damfile=paste0(taidir, 'test_dam/dam_condprobs.txt')
udamfile=paste0(taidir, 'test_dam/unmeth_condprobs.txt')

cpg=read.table(cpgfile, sep=',')
colnames(cpg)=c('prob', 'context')

##getting confusion thing for taiyaki
getconfusion_tai <- function(methfile, unmethfile) {
    meth=read_csv(methfile, col_names=c('prob','context')) %>%
        mutate(meth=TRUE)
    unmeth=read_csv(unmethfile, col_names=c('prob', 'context')) %>%
        mutate(meth=FALSE)
    all=rbind(meth, unmeth)

    thresholds=seq(min(all$prob), max(all$prob), .1)
    pos=sum(all$meth)
    neg=sum(!all$meth)

    confusions=foreach(i=1:length(thresholds), .combine=rbind) %dopar% {
        thresh=thresholds[i]

        call=all$prob > thresh
        tp=call & all$meth
        fp=call & !all$meth

        tpr=sum(tp)/pos
        fpr=sum(fp)/neg

        conf=data.frame(thresh=thresh, tpr=tpr, fpr=fpr)
        return(conf)
    }
}

cpgconf=getconfusion_tai(cpgfile, ucpgfile)
cpgconf$motif='cpg'
damconf=getconfusion_tai(damfile, udamfile)
damconf$motif='dam'


taiconf=rbind(cpgconf, damconf)
taiconf$caller='taiyaki'


datadir='/dilithium/Data/Nanopore/Analysis/171120_ecoliStandards/methcall/'
dammeth=paste0(datadir,'dam/subset/171122_ecolidam.0.meth.dam.tsv')
damunmeth=paste0(datadir,'dam/subset/171122_ecoliUM.0.meth.dam.tsv')
cgmeth=paste0(datadir, 'cpg/subset/171124_ecoliCpG.0.meth.cpg.tsv')
cgunmeth=paste0(datadir,'cpg/subset/171122_ecoliUM.0.meth.cpg.tsv')

getconfusion_np <- function(methfile, unmethfile) {
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

damconf=getconfusion_np(dammeth, damunmeth)
damconf$motif='dam'
cpgconf=getconfusion_np(cgmeth, cgunmeth)
cpgconf$motif='cpg'

npconf=rbind(cpgconf, damconf)
npconf$caller='np'

allconf=rbind(npconf, taiconf)

rocfile='~/Dropbox/yfan/methylation/taiyaki/tai_v_np_roc.pdf'
pdf(rocfile, height=8.5, width=11)
ggplot(allconf, aes(x=fpr, y=tpr, colour=motif, linetype=caller)) +
    geom_line() +
    geom_abline(slope=1, intercept=0) +
    xlim(0,1) +
    ylim(0,1) +
    ggtitle('Taiyaki vs Nanopolish ROC') +
    xlab('False Positive Rate') +
    ylab('True Positive Rate') +
    theme_bw()
dev.off()


##plot distributions
getdists_np <- function(methfile, unmethfile) {
    meth=read_tsv(methfile) %>%
        mutate(meth=TRUE)
    unmeth=read_tsv(unmethfile) %>%
        mutate(meth=FALSE)
    all=rbind(meth, unmeth)
    dists=data.frame(meth=all$meth, loglik=all$log_lik_ratio)
    return(dists)
}


damnpdists=getdists_np(dammeth, damunmeth)
damnpdists$motif='dam'
damnpdists$caller='np'
cgnpdists=getdists_np(cgmeth, cgunmeth)
cgnpdists$motif='cpg'
cgnpdists$caller='np'

all=rbind(damnpdists, cgnpdists)

getdists_tai <- function(methfile, unmethfile) {
    meth=read_csv(methfile, col_names=c('prob','context')) %>%
        mutate(meth=TRUE)
    unmeth=read_csv(unmethfile, col_names=c('prob', 'context')) %>%
        mutate(meth=FALSE)
    all=rbind(meth, unmeth)
    return(all)
}    

damtaidists=getdists_tai(damfile, udamfile)
damtaidists$motif='dam'
damtaidists$caller='taiyaki'
damtaidists$probmeth=exp(damtaidists$prob)
cgtaidists=getdists_tai(cpgfile, ucpgfile)
cgtaidists$motif='cpg'
cgtaidists$caller='taiyaki'
cgtaidists$probmeth=exp(cgtaidists$prob)


dists='~/Dropbox/yfan/methylation/taiyaki/probdists.pdf'
pdf(dists, height=8.5, width=11)
npdam=ggplot(damnpdists, aes(x=loglik)) +
    geom_histogram(data=subset(damnpdists, meth==TRUE), fill='red', colour='red', binwidth=1, alpha=.2) +
    geom_histogram(data=subset(damnpdists, meth==FALSE), fill='blue', colour='blue', binwidth=1, alpha=.2) +
    ggtitle('Nanopolish dam Log Likelihood Ratios') +
    xlab('Log Lik Ratio') +
    ylab('Density') +
    theme_bw()
npcg=ggplot(cgnpdists, aes(x=loglik)) +
    geom_histogram(data=subset(cgnpdists, meth==TRUE), fill='red', colour='red', binwidth=1, alpha=.2) +
    geom_histogram(data=subset(cgnpdists, meth==FALSE), fill='blue', colour='blue', binwidth=1, alpha=.2) +
    ggtitle('Nanopolish cpg Log Likelihood Ratios') +
    xlab('Log Lik Ratio') +
    ylab('Density') +
    theme_bw()
taidam=ggplot(damtaidists, aes(x=probmeth)) +
    geom_histogram(data=subset(damtaidists, meth==TRUE), fill='red', colour='red', binwidth=.01, alpha=.2) +
    geom_histogram(data=subset(damtaidists, meth==FALSE), fill='blue', colour='blue', binwidth=.01, alpha=.2) +
    ggtitle('Taiyaki dam Conditional Probabilities') +
    xlab('Cond Prob') +
    ylab('Density') +
    theme_bw()
taicg=ggplot(cgtaidists, aes(x=probmeth)) +
    geom_histogram(data=subset(cgtaidists, meth==TRUE), fill='red', colour='red', binwidth=.01, alpha=.2) +
    geom_histogram(data=subset(cgtaidists, meth==FALSE), fill='blue', colour='blue', binwidth=.01, alpha=.2) +
    ##scale_y_log10() +
    ggtitle('Taiyaki cpg Conditional Probabilities') +
    xlab('Cond Prob') +
    ylab('Density') +
    theme_bw()
library(gridExtra)
grid.arrange(npdam, npcg, taidam, taicg, ncol=2, nrow=2)
dev.off()



##calculate auc
library(plyr)
auc=ddply(allconf, .(caller), function(x) {
    fpr=seq(0,1,.01)
    dam=x[x$motif=='dam',]

    damval=0
    for (i in fpr) {
        damval=damval+dam$tpr[which.min(abs(dam$fpr-i))]
    }
    damauc=damval/101

    cg=x[x$motif=='cpg',]
    cgval=0
    for (i in fpr) {
        cgval=cgval+cg$tpr[which.min(abs(cg$fpr-i))]
    }
    cgauc=cgval/101
    return(data.frame(dam=damauc, cpg=cgauc))
})

table='~/Dropbox/yfan/methylation/taiyaki/auc.pdf'
pdf(table)
grid.table(auc)
dev.off()
