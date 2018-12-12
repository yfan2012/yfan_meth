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


    kmers=unique(eventalign$reference_kmer)

    
    ##group by events
    merinfo=foreach(i=1:length(kmers), .combine=rbind) %dopar% {
    ##for (i in 1:length(kmers)) {
        library(tidyverse)
        kmerevent=eventalign %>%
            filter(reference_kmer == kmers[i])
        kmerevent_temp=kmerevent %>%
            filter(reference_kmer == model_kmer)
        pvaltemp=pvals(kmerevent_temp)
        kmerevent_comp=kmerevent %>%
            filter(reference_kmer != model_kmer)
        if (dim(kmerevent_comp)[1] >= 1) {
            pvalcomp=pvals(kmerevent_comp)
        } else {
            pvalcomp=pvaltemp
        }
        
        merinfo=data.frame(mer=kmers[i], diffmean=mean(kmerevent$diff), diffstd=sd(kmerevent$diff), eventmeans=mean(kmerevent$event_level_mean), eventstd=sd(kmerevent$event_level_mean), tempz=mean(kmerevent_temp$event_level_mean)/kmerevent_temp$model_stdv[1], compz=mean(kmerevent_comp$event_level_mean)/kmerevent_temp$model_stdv[1], tempt=pvaltemp[1], tempks=pvaltemp[2], tempwil=pvaltemp[3], compt=pvalcomp[1], compks=pvalcomp[2], compwil=pvalcomp[3])
        return(merinfo)
    }

    
    largediff=as_tibble(merinfo) %>%
        arrange(desc(diffmean))

    largetempz=as_tibble(merinfo) %>%
        arrange(desc(tempz))

    largecompz=as_tibble(merinfo) %>%
        arrange(desc(compz))
    
    largetempt=as_tibble(merinfo) %>%
        arrange(desc(tempt))

    largetempks=as_tibble(merinfo) %>%
        arrange(desc(tempks))

    largetempwil=as_tibble(merinfo) %>%
        arrange(desc(tempwil))

    largecompt=as_tibble(merinfo) %>%
        arrange(desc(compt))

    largecompks=as_tibble(merinfo) %>%
        arrange(desc(compks))

    largecompwil=as_tibble(merinfo) %>%
        arrange(desc(compwil))
    
    view=data.frame(diffmer=largediff$mer[1:10], diff=largediff$diffmean[1:10],
                    tempzmer=largetempz$mer[1:10], temp_z=largetempz$tempz[1:10],
                    compzmer=largecompz$mer[1:10], comp_z=largecompz$compz[1:10],
                    temptmer=largetempt$mer[1:10], temp_t=largetempt$tempt[1:10],
                    comptmer=largecompt$mer[1:10], comp_t=largecompt$compt[1:10],
                    tempksmer=largetempks$mer[1:10], temp_ks=largetempks$tempks[1:10],
                    compksmer=largecompks$mer[1:10], comp_ks=largecompks$compks[1:10],
                    tempwilmer=largetempwil$mer[1:10], temp_wil=largetempwil$tempwil[1:10],
                    compwilmer=largecompwil$mer[1:10], comp_wil=largecompwil$compwil[1:10])


    pdf(paste0(plot, prefix, '_rankmer.pdf'), width=30, height=5)
    grid.table(view, rows=NULL)
    dev.off()

}
