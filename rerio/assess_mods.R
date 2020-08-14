library(tidyverse)
library(cowplot)

datadir='/mithril/Data/Nanopore/projects/methbin/rerio/'
dbxdir='~/Dropbox/yfan/methylation/methbin/rerio/'

##Just doing modA comparison
samps=c('neb17', 'neb19')
amods=c('GANTC', 'GATC')
cmods=c('CCWGG', 'GATC')

plot_modprobs <- function(samp, motifs, base) {
    ###takes the samp and motifs that were extracted, and plots
    cols=c('read_name', 'read_index', 'read_base', 'ref_name','ref_index', 'ref_base', 'pAmod', 'pCmod', 'motifmatch', 'strand', 'readmotif')
    sampmodinfo=tibble(
        read_name=as.character(),
        read_index=as.numeric(),
        read_base=as.character(),
        ref_name=as.character(),
        ref_index=as.character(),
        ref_base=as.numeric(),
        pAmod=as.numeric(),
        pCmod=as.numeric(), 
        motifmatch=as.character(),
        strand=as.character(), 
        readmotif=as.character(),
        motif=as.character(),
        samp=as.character()
    )
    
    for (motif in motifs) {
        modfile=paste0(datadir, samp, '/assess/', samp, '.', motif, '.csv')
        modinfo=read_csv(modfile, col_names=cols) %>%
            mutate(motif=motif) %>%
            mutate(samp=samp) %>%
            mutate(pAmod=as.numeric(pAmod)) %>%
            mutate(pCmod=as.numeric(pCmod)) %>%
            filter(motifmatch==TRUE)
        sampmodinfo=rbind(sampmodinfo, modinfo)
    }
    
    ##basic density of pA
    if (base=='A') {
        plot=ggplot(sampmodinfo, aes(x=pAmod, colour=motif, fill=motif, alpha=.3)) +
            geom_density() +
            scale_y_continuous(trans='pseudo_log') +
            ggtitle(samp) +
            xlab('Prob mod (out of 255)') +
            theme_bw()
        return(plot)
    }
    if (base=='C') {
        plot=ggplot(sampmodinfo, aes(x=probCmod, colour=motif, fill=motif, alpha=.3)) +
            geom_density() +
            scale_y_continuous(trans='pseudo_log') +
            ggtitle(samp) +
            xlab('Prob mod (out of 255)') +
            theme_bw()
        return(plot)
    }
}

neb17a=plot_modprobs('neb17', amods, 'A')
neb19a=plot_modprobs('neb19', amods, 'A')


outfile=paste0(dbxdir, 'rerio_Amod.pdf')
pdf(outfile, h=10, w=10)
plot_grid(neb17a, neb19a, ncol=1, align='v')
dev.off()


neb17c=plot_modprobs('neb17', amods, 'C')
neb19c=plot_modprobs('neb19', amods, 'C')

outfile=paste0(dbxdir, 'rerio_Amod_Ctest.pdf')
pdf(outfile, h=10, w=10)
plot_grid(neb17c, neb19c, ncol=1, align='v')
dev.off()


plot_moderrors <- function(samp, motifs, base) {
    ###takes the samp and motifs that were extracted, and plots
    cols=c('read_name', 'read_index', 'read_base', 'ref_name','ref_index', 'ref_base', 'pAmod', 'pCmod', 'motifmatch', 'strand', 'readmotif')
    errmodinfo=tibble(
        read_name=as.character(),
        numerr=as.numeric(),
        errfrac=as.numeric(),
        nummotifs=as.numeric(),
        motif=as.character(),
        samp=as.character()
    )
    
    for (motif in motifs) {
        modfile=paste0(datadir, samp, '/assess/', samp, '.', motif, '.csv')
        errinfo=read_csv(modfile, col_names=cols) %>%
            group_by(read_name) %>%
            summarise(numerr=sum(motifmatch), errfrac=sum(motifmatch)/length(read_name), nummotifs=length(read_name)) %>%
            mutate(motif=motif) %>%
            mutate(samp=samp) %>%
            filter(numerr>=5)
            
        errmodinfo=rbind(errmodinfo, errinfo)
    }
    return(errmodinfo)
}

neb17err=plot_moderrors('neb17', amods, 'A')
neb19err=plot_moderrors('neb19', amods, 'A')
allerr=rbind(neb17err, neb19err)

outfile=paste0(dbxdir, 'rerio_Amod_errors.pdf')

pdf(outfile, h=13, w=10)
plot=ggplot(allerr, aes(x=samp, y=errfrac, colour=motif, fill=motif, alpha=.3)) +
    geom_violin(scale='width') +
    ggtitle('Errors at meth motifs') +
    xlab('Sample') +
    ylab('Fraction of motifs that are correct per read (with >5 motifs)') +
    theme_bw()
print(plot)
dev.off()
