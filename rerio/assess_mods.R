library(tidyverse)
library(cowplot)

datadir='/mithril/Data/Nanopore/projects/methbin/rerio/'
dbxdir='~/Dropbox/yfan/methylation/methbin/rerio/'

##Just doing modA comparison
samps=c('neb17', 'neb19')
motifs=c('GANTC', 'GATC')

plot_modprobs <- function(samp) {
    cols=c('readid', 'readpos', 'refpos', 'base', 'chr', 'pA', 'pAmod', 'pC', 'pCmod')
    sampmodinfo=tibble(
        readid=as.character(),
        readpos=as.numeric(),
        refpos=as.numeric(),
        base=as.character(),
        chr=as.character(),
        pA=as.numeric(),
        pAmod=as.numeric(),
        pC=as.numeric(),
        pCmod=as.numeric(),
        samp=as.character(),
        motif=as.character()
    )
    
    for (motif in motifs) {
        modfile=paste0(datadir, samp, '/assess/', samp, '.', motif, '.csv')
        modinfo=read_csv(modfile, col_names=cols) %>%
            mutate(motif=motif) %>%
            mutate(samp=samp) %>%
            filter(pAmod!='None') %>%
            mutate(probAmod=as.numeric(pAmod))
        sampmodinfo=rbind(sampmodinfo, modinfo)
    }
    
    ##basic density of pA 
    plot=ggplot(sampmodinfo, aes(x=probAmod, colour=motif, fill=motif, alpha=.3)) +
        geom_density() +
        scale_y_continuous(trans='pseudo_log') +
        ggtitle(samp) +
        xlab('Prob mod (out of 255)') +
        theme_bw()
    return(plot)
}

neb17=plot_modprobs('neb17')
neb19=plot_modprobs('neb19')

outfile=paste0(dbxdir, 'rerio_Amod.pdf')
pdf(outfile, h=10, w=10)
plot_grid(neb17, neb19, ncol=1, align='v')
dev.off()



