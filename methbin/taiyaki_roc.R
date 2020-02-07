library(tidyverse)
library(doParallel)
library(foreach)

cl=makeCluster(10)
registerDoParallel(cl, cores=10)

datadir='~/data/methbin/calls/'
samps=tibble(gDNA=c('neb12', 'neb13', 'neb14', 'neb15', 'neb16', 'neb17', 'neb19'), 
             plas=c('neb2', 'neb3', 'neb10', 'neb4', 'neb5', 'neb6', 'neb9'))




getconfusion <- function(methfile, unmethfile) {
    meth=read_tsv(methfile) %>%
        select(read_id, pos, mod_log_prob, can_log_prob) %>%
        mutate(ratio=can_log_prob/mod_log_prob) %>%
        mutate(meth=TRUE)
    unmeth=read_tsv(unmethfile) %>%
        select(read_id, pos, mod_log_prob, can_log_prob) %>%
        mutate(ratio=can_log_prob/mod_log_prob) %>%
        mutate(meth=FALSE)
    all=rbind(meth, unmeth)

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




for (i in 1:dim(samps)[1]) {
    modfile=paste0(datadir, samps$gDNA[i], '/', samps$gDNA[i], '/per_read_modified_base_calls.txt')
    plsfile=paste0(datadir, samps$gDNA[i], '/', samps$plas[i], '/per_read_modified_base_calls.txt')
    ctrfile=paste0(datadir, samps$gDNA[i], '/neb11/per_read_modified_base_calls.txt')

    conf=getconfusion(modfile, ctrfile)

    

   



