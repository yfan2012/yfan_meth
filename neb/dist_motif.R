library(ggplot2)
library(tidyverse)
library(gridExtra)
library(doParallel)
library(foreach)
library(reshape2)
cl=makeCluster(10)
registerDoParallel(cl, cores=10)

plot='~/scratch/plots/neb/'

pvals <- function(kmerevent) {
    ##take in model mean, std, and set of event means. Sample from normal and get t test, ks, man whitney u, kldiv, etc
    ttest.pval=-log10(t.test(kmerevent$event_level_mean, mu=kmerevent$model_mean[1])$p.value)
    kstest.pval=-log10(ks.test(kmerevent$event_level_mean, 'pnorm', kmerevent$model_mean[1], kmerevent$model_stdv[1])$p.value)
    wilcox.pval=-log10(wilcox.test(kmerevent$event_level_mean, kmerevent$model_mean[1])$p.value)
    zval=mean(kmerevent$event_level_mean-kmerevent$model_mean[1])/kmerevent$event_stdv[1]
    return(c(ttest.pval, kstest.pval, wilcox.pval, zval))
}

allsamps=c('170906_neb14', '171003_neb16', '171005_neb12', '171012_neb13', '171019_neb17', '171019_neb19', '171020_neb11', '171020_neb15', '180628_neb_dcm')
motifkey=read_csv('~/Code/methylation/neb/meth_key.csv')
motifs=unique(motifkey$motif)

##for (prefix in prefixes) {
allmotifdists=foreach(prefix=allsamps, .combine=rbind) %dopar% { 
    library(tidyverse)
    library(foreach)
    eventtsv=paste0('/scratch/groups/mschatz1/cpowgs/meth/', prefix, '/eventalign/', prefix, '.er2796_2930000_2940000.eventalign.tsv')
    eventalign=read_tsv(eventtsv) %>%
        filter(model_kmer != 'NNNNNN') %>%
        select(-contig) %>%
        select(-read_index) %>%
        select(-event_index) %>%
        select(-strand) %>%
        select(-standardized_level) %>%
        mutate(diff=abs(model_mean-event_level_mean))

    ##calculate distance to nearest motif
    methdist=foreach(seq=1:length(motifs), .combine=rbind) %do% {
        library(tidyverse)
        ##for (seq in 1:length(motifs)){
        motifdist=eventalign %>%
            mutate(motif=grepl(motifs[seq], reference_kmer, fixed=TRUE))
        
        mpos=unique(motifdist$position[motifdist$motif==TRUE])
        
        dist=motifdist$position
        for (i in 1:length(dist)) {
            diffs=dist[i]-mpos
            dist[i]=diffs[which.min(abs(diffs))]
        }
        
        motifdist$dist=dist
        
        motifdist=motifdist %>%
            filter(dist > -20 & dist < 20) %>%
            mutate(motif=motifs[seq]) %>%
            mutate(sample=prefix)
        
        return(motifdist)
    }

    return(methdist)
}
stopCluster(cl)

##figure out a pval for each motif/sample combination at each position
for (i in motifs) {
    cl=makeCluster(10)
    registerDoParallel(cl, cores=10)

    prefixes=c(unique(motifkey$sample[motifkey$motif==i]), '171020_neb11')
    motifsamps=foreach(samp=1:length(prefixes), .combine=rbind) %dopar% {
        library(tidyverse)
        library(foreach)
        ##for (samp in prefixes) {
        ##get relevant data (within 20bp of one motif in one sample, template only)
        sampmotif=allmotifdists %>%
            filter(motif==i) %>%
            filter(sample==prefixes[samp]) %>%
            filter(reference_kmer==model_kmer)
        
        positions=unique(sampmotif$position)
        kmers=unique(sampmotif$reference_kmer)
        
        ##figure out pvals for each position 
        motifvals=foreach(pos=1:length(positions), .combine=rbind) %do% {
            thispos=sampmotif %>%
                filter(position==positions[pos]) %>%
                filter(reference_kmer==model_kmer) %>%
                mutate(diff=event_level_mean-model_mean)
            
            if (dim(thispos)[1]>1) {
                ##data frame where kmers and dists should match
                diffvals=pvals(thispos)
            }
            
            mval=data.frame(mer=thispos$reference_kmer[1], dist=thispos$dist[1], samp=thispos$sample[1], cov=dim(thispos)[1], diffmean=mean(thispos$diff), diffstd=sd(thispos$diff), eventmeans=mean(thispos$event_level_mean), eventstd=sd(thispos$event_level_mean), zval=diffvals[4], tval=diffvals[1], ks=diffvals[2], wilcox=diffvals[3])
            return(mval)
        }
        
        return(motifvals)
    }
    stopCluster(cl)

    pdf(paste0(plot, i, '_distance.pdf'))    
    ##Coverage
    print(ggplot(motifsamps, aes(x=dist, y=cov, colour=samp)) +
          geom_point(alpha=.5) +
          ggtitle(paste0('Coverage within 20 bp of ', i)) +
          xlab(paste0('Distance to ', i)) +
          ylab('Event Coverge') +
          theme_bw())

    ##mean difference
    print(ggplot(motifsamps, aes(x=dist, y=diffmean, colour=samp)) +
          geom_point(alpha=.5) +
          ggtitle(paste0('Mean event level difference  within 20 bp of ', i)) +
          xlab(paste0('Distance to ', i)) +
          ylab('Mean Difference') +
          theme_bw())

    ##zvalue
    print(ggplot(motifsamps, aes(x=dist, y=zval, colour=samp)) +
          geom_point(alpha=.5) +
          ggtitle(paste0('Z value within 20 bp of ', i)) +
          xlab(paste0('Distance to ', i)) +
          ylab('Z value') +
          theme_bw())

    ##T test
    print(ggplot(motifsamps, aes(x=dist, y=tval, colour=samp)) +
          geom_point(alpha=.5) +
          ggtitle(paste0('T-test p value within 20 bp of ', i)) +
          xlab(paste0('Distance to ', i)) +
          ylab('P value') +
          theme_bw())

    ##KS test
    print(ggplot(motifsamps, aes(x=dist, y=ks, colour=samp)) +
          geom_point(alpha=.5) +
          ggtitle(paste0('KS test within 20 bp of ', i)) +
          xlab(paste0('Distance to ', i)) +
          ylab('P value') +
          theme_bw())

    ##wilcox
    print(ggplot(motifsamps, aes(x=dist, y=wilcox, colour=samp)) +
          geom_point(alpha=.5) +
          ggtitle(paste0('Wilcox Rank Sum within 20 bp of ', i)) +
          xlab(paste0('Distance to ', i)) +
          ylab('P value') +
          theme_bw())
    dev.off()
}
    
