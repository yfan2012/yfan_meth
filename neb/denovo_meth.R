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
    zval=(kmerevent$event_level_mean-kmerevent$model_mean[1])/kmerevent$event_stdv[1]
    return(c(ttest.pval, kstest.pval, wilcox.pval, zval))
}

prefixes=c('170906_neb14', '171003_neb16', '171005_neb12', '171012_neb13', '171019_neb17', '171019_neb19', '171020_neb11', '171020_neb15', '180628_neb_dcm')
##prefixes=c('180628_neb_dcm')
for (prefix in prefixes) {
    eventtsv=paste0('/scratch/groups/mschatz1/cpowgs/meth/', prefix, '/eventalign/', prefix, '.er2796_2930000_2940000.eventalign.tsv')
    eventalign=read_tsv(eventtsv) %>%
        filter(model_kmer != 'NNNNNN') %>%
        select(-contig) %>%
        select(-read_index) %>%
        select(-event_index) %>%
        select(-strand) %>%
        select(-standardized_level) %>%
        mutate(diff=abs(model_mean-event_level_mean))

    
    ##group by ref position
    pos=unique(eventalign$position)

    posinfo=foreach(i=1:length(pos), .combine=rbind) %dopar% {
    ##for (i in 1:length(pos)) {
        library(tidyverse)
        posevent=eventalign %>%
            filter(position == pos[i])
        posevent_temp=posevent %>%
            filter(reference_kmer == model_kmer)
        ##check if there are enough events to add to merinfo
        if (dim(posevent_temp)[1] > 1 ) {
            pvaltemp=pvals(posevent_temp)
            posevent_comp=posevent %>%
                filter(reference_kmer!= model_kmer)
            if (dim(posevent_comp)[1] > 1) {
                pvalcomp=pvals(posevent_comp)
            } else {
                pvalcomp=pvaltemp
            }
        
            merinfo=data.frame(pos=pos[i], cov_temp=dim(posevent_temp)[1],cov_comp=dim(posevent_comp)[1], diffmean=mean(posevent$diff), diffstd=sd(posevent$diff), eventmeans=mean(posevent$event_level_mean), eventstd=sd(posevent$event_level_mean),compz=pvaltemp[4], tempz=pvaltemp[4], tempt=pvaltemp[1], tempks=pvaltemp[2], tempwil=pvaltemp[3], compt=pvalcomp[1], compks=pvalcomp[2], compwil=pvalcomp[3])

            return(merinfo)
        }
        
    }
    
    pdf(paste0(plot, prefix, '_genome_pvals.pdf'), width=21, height=7)
    print(ggplot(posinfo, aes(x=pos, y=cov_temp)) +
          geom_point(aes(alpha=.5)) +
          ggtitle('Events per position: Template Strand') +
          xlab('K-mer Coordinate') +
          ylab('Coverage') +
          theme_bw())
    
    print(ggplot(posinfo, aes(x=pos, y=cov_comp)) +
          geom_point(aes(alpha=.5)) +
          ggtitle('Events per position: Complement Strand') +
          xlab('K-mer Coordinate') +
          ylab('Coverage') +
          theme_bw())
    
    print(ggplot(posinfo, aes(x=pos, y=diffmean)) +
          geom_point(aes(alpha=.5)) +
          ggtitle('Mean Current Diff (model-event)') +
          xlab('K-mer Coordinate') +
          ylab('Event level difference') +
          theme_bw())
    
    zval=tibble(template=posinfo$tempz, complement=posinfo$compz, pos=posinfo$pos) %>%
        gather(strand, zval, -pos)
    print(ggplot(zval, aes(x=pos, y=zval, colour=strand)) +
          geom_point(aes(alpha=.5)) +
          ggtitle('Z values') +
          xlab('K-mer Coordinate') +
          ylab('z value') +
          theme_bw())
    
    ttest=tibble(template=posinfo$tempt, complement=posinfo$compt, pos=posinfo$pos) %>%
        gather(strand, pval, -pos)
    print(ggplot(ttest, aes(x=pos, y=pval, colour=strand)) +
          geom_point(aes(alpha=.5)) +
          ggtitle('T test') +
          xlab('K-mer Coordinate') +
          ylab('-log P value') +
          theme_bw())
    
    kstest=tibble(template=posinfo$tempks, complement=posinfo$compks, pos=posinfo$pos) %>%
        gather(strand, pval, -pos)
    print(ggplot(kstest, aes(x=pos, y=pval, colour=strand)) +
          geom_point(aes(alpha=.5)) +
          ggtitle('KS test') +
          xlab('K-mer Coordinate') +
          ylab('-log P value') +
          theme_bw())

    
    wiltest=tibble(template=posinfo$tempwil, complement=posinfo$compwil, pos=posinfo$pos) %>%
        gather(strand, pval, -pos)
    print(ggplot(wiltest, aes(x=pos, y=pval, colour=strand)) +
        geom_point(aes(alpha=.5)) +
        ggtitle('Wilcox rank sum test') +
        xlab('K-mer Coordinate') +
        ylab('-log P value') +
        theme_bw())

    dev.off()
}
