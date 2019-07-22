library(ggplot2)
library(tidyverse)
library(doParallel)
library(foreach)
cl=makeCluster(10)
registerDoParallel(cl, cores=10)


##read in stuff
cpgfile='~/data/test_cpg/cpg_condprobs.txt'
ucpgfile='~/data/test_cpg/unmeth_condprobs.txt'
damfile='~/data/test_dam/dam_condprobs.txt'
udamfile='~/data/test_dam/unmeth_condprobs.txt'

cpg=read.table(cpgfile, sep=',')
colnames(cpg)=c('prob', 'context')


getconfusion <- function(methfile, unmethfile) {
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

cpgconf=getconfusion(cpgfile, ucpgfile)
cpgconf$motif='cpg'
damconf=getconfusion(damfile, udamfile)
damconf$motif='dam'

allconf=rbind(cpgconf, damconf)


rocfile='~/Dropbox/yfan/methylation/taiyaki/methcall_roc.pdf'

pdf(rocfile, height=8.5, width=11)
ggplot(allconf, aes(x=fpr, y=tpr, colour=motif)) +
    geom_line() +
    geom_abline(slope=1, intercept=0) +
    xlim(0,1) +
    ylim(0,1) +
    ggtitle('Taiyaki ROC') +
    xlab('False Positive Rate') +
    ylab('True Positive Rate') +
    theme_bw()
dev.off()
   
