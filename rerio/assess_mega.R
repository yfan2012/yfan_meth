library(tidyverse)
library(cowplot)

datadir='/mithril/Data/Nanopore/projects/methbin/rerio/'
dbxdir='~/Dropbox/yfan/methylation/methbin/rerio/'

amods=c('GATC', 'GANTC') ##Z
cmods=c('CCWGG', 'GATC', 'GCNGC') ##Y
                 

plot_mega <- function(samp, motifs, mod, mode) {
    allmegainfo=tibble(
        read_id=as.character(),
        chrm=as.character(),
        strand=as.numeric(),
        pos=as.numeric(),
        mod_log_prob=as.numeric(),
        can_log_prob=as.numeric(),
        mod_base=as.character(),
        motif=as.character(),
        modratio=as.numeric(),
        samp=as.character())
    
    for (motif in motifs) {
        if (mode=='vanilla') {
            megafile=paste0(datadir, '/', samp, '/megalodon/', samp, '_vanilla_', motif, '_', mod, '/per_read_modified_base_calls.txt')
            title=paste0(samp, ' old model')
        }else{
            megafile=paste0(datadir, '/', samp, '/megalodon/', samp, '_', motif, '_', mod, '/per_read_modified_base_calls.txt')
            title=paste0(samp, ' rerio')
        }
        megainfo=read_tsv(megafile) %>%
            mutate(modratio=can_log_prob/mod_log_prob) %>%
            mutate(samp=samp)
        allmegainfo=rbind(allmegainfo, megainfo)
    }

    plot=ggplot(allmegainfo, aes(x=modratio, colour=motif, fill=motif, alpha=.3)) +
        geom_density() +
        ggtitle(title) +
        xlab('modprob ratio') +
        theme_bw()
    
    return(plot)
}


##amods
neb19=plot_mega('neb19', amods, 'Y', '')
neb19v=plot_mega('neb19', amods, 'Y', 'vanilla')
neb17=plot_mega('neb17', amods, 'Y', '')
neb17v=plot_mega('neb17', amods, 'Y', 'vanilla')
neb11=plot_mega('neb11', amods, 'Y', '')
neb11v=plot_mega('neb11', amods, 'Y', 'vanilla')

##cmods
cneb11=plot_mega('neb11', cmods, 'Z', '')
cneb11v=plot_mega('neb11', cmods, 'Z', 'vanilla')
cneb14=plot_mega('neb14', cmods, 'Z', '')
cneb14v=plot_mega('neb14', cmods, 'Z', 'vanilla')
cneb15=plot_mega('neb15', cmods, 'Z', '')
cneb15v=plot_mega('neb15', cmods, 'Z', 'vanilla')
cnebdcm=plot_mega('nebdcm', cmods, 'Z', '')
cnebdcmv=plot_mega('nebdcm', cmods, 'Z', 'vanilla')

##plot nebX and neb11, regular and vanilla, windowpane style
outfile=paste0(dbxdir, 'rerio_mega_Amod.pdf')
pdf(outfile, h=7, w=10)
plot_grid(neb19, neb11, neb19v, neb11v, ncol=2, align='v')
plot_grid(neb17, neb11, neb17v, neb11v, ncol=2, align='v')
dev.off()

outfile=paste0(dbxdir, 'rerio_mega_Cmod.pdf')
pdf(outfile, h=7, w=13)
plot_grid(cneb14, cneb11, cneb14v, cneb11v, ncol=2, align='v')
plot_grid(cneb15, cneb11, cneb15v, cneb11v, ncol=2, align='v')
plot_grid(cnebdcm, cneb11, cnebdcmv, cneb11v, ncol=2, align='v')
dev.off()


