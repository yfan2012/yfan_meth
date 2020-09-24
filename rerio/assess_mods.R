library(tidyverse)
library(cowplot)
library(doParallel)
library(foreach)

datadir='/mithril/Data/Nanopore/projects/methbin/rerio/'
dbxdir='~/Dropbox/yfan/methylation/methbin/rerio/'

##Just doing modA comparison
samps=c('neb17', 'neb19', 'neb11')
amods=data.frame(motifs=c('GANTC', 'GATC'), positions=c(1,1))
cmods=data.frame(motifs=c('CCWGG', 'GATC', 'GCNGC'), positions=c(1,3,1))

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
    
    for (i in 1:dim(motifs)[1]) {
        motif=motifs$motifs[i]
        pos=as.character(motifs$positions[i])
        modfile=paste0(datadir, samp, '/assess/', samp, '.', motif, '.', pos, '.csv')
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
        plot=ggplot(sampmodinfo, aes(x=pCmod, colour=motif, fill=motif, alpha=.3)) +
            geom_density() +
            scale_y_continuous(trans='pseudo_log') +
            ggtitle(samp) +
            xlab('Prob mod (out of 255)') +
            theme_bw()
        return(plot)
    }
}

rerio_roc  <- function(samp, motif, basepos, base) {
    modfile=paste0(datadir, samp, '/assess/', samp, '.', motif, '.', basepos, '.csv')
    unmodfile=paste0(datadir, 'neb11/assess/neb11', '.', motif, '.', basepos, '.csv')
    ###takes the samp and motifs that were extracted, and makes an roc
    cols=c('read_name', 'read_index', 'read_base', 'ref_name','ref_index', 'ref_base', 'pAmod', 'pCmod', 'motifmatch', 'strand', 'readmotif')

    mod=read_csv(modfile, col_names=cols) %>%
        mutate(meth=TRUE)
    unmod=read_csv(unmodfile, col_names=cols) %>%
        mutate(meth=FALSE)
    all=rbind(mod, unmod)

    thresholds=seq(0,255,1)
    pos=sum(all$meth)
    neg=sum(!all$meth)


    confusion=foreach(i=1:length(thresholds), .combine=rbind) %dopar% {
        thresh=thresholds[i]

        if (base=='A') {
            call=all$pAmod > thresh
        }else if (base=='C') {
            call=all$pCmod > thresh
        }

        
        tp=call & all$meth
        fp=call & !all$meth

        tpr=sum(tp)/pos
        fpr=sum(fp)/neg
        conf=data.frame(thresh=thresh, tpr=tpr, fpr=fpr)
        return(conf)
    }
}

        

neb17a=plot_modprobs('neb17', amods, 'A')
neb19a=plot_modprobs('neb19', amods, 'A')
neb11a=plot_modprobs('neb11', amods, 'A')

outfile=paste0(dbxdir, 'rerio_Amod.pdf')
pdf(outfile, h=10, w=7)
plot_grid(neb11a, neb17a, neb19a, ncol=1, align='v')
dev.off()

neb17c=plot_modprobs('neb17', amods, 'C')
neb19c=plot_modprobs('neb19', amods, 'C')

neb11c=plot_modprobs('neb11', cmods, 'C')
neb14c=plot_modprobs('neb14', cmods, 'C')
neb15c=plot_modprobs('neb15', cmods, 'C')
nebdcmc=plot_modprobs('nebdcm', cmods, 'C')

outfile=paste0(dbxdir, 'rerio_Cmod.pdf')
pdf(outfile, h=10, w=18)
plot_grid(neb11c, neb14c, neb15c, nebdcmc, ncol=2, align='v')
dev.off()


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
    
    for (i in 1:dim(motifs)[1]) {
        motif=motifs$motif[i]
        pos=as.character(motifs$positions[i])
        modfile=paste0(datadir, samp, '/assess/', samp, '.', motif, '.', pos, '.csv')
        errinfo=read_csv(modfile, col_names=cols) %>%
            group_by(read_name) %>%
            summarise(numerr=sum(motifmatch), errfrac=sum(motifmatch)/length(read_name), nummotifs=length(read_name)) %>%
            mutate(motif=motif) %>%
            mutate(samp=samp) %>%
            filter(nummotifs>=5)
            
        errmodinfo=rbind(errmodinfo, errinfo)
    }
    return(errmodinfo)
}

neb17err=plot_moderrors('neb17', amods, 'A')
neb19err=plot_moderrors('neb19', amods, 'A')
neb11err=plot_moderrors('neb11', amods, 'A')
allerr=rbind(neb11err, neb17err, neb19err)

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

neb11err=plot_moderrors('neb11', cmods, 'C')
neb14err=plot_moderrors('neb14', cmods, 'C')
neb15err=plot_moderrors('neb15', cmods, 'C')
nebdcmerr=plot_moderrors('nebdcm', cmods, 'C')
allCerr=rbind(neb11err, neb14err, neb15err, nebdcmerr)

outfile=paste0(dbxdir, 'rerio_Cmod_errors.pdf')

pdf(outfile, h=7, w=13)
plot=ggplot(allCerr, aes(x=samp, y=errfrac, colour=motif, fill=motif, alpha=.3)) +
    geom_violin(scale='width') +
    ggtitle('Errors at meth motifs') +
    xlab('Sample') +
    ylab('Fraction of motifs that are correct per read (with >5 motifs)') +
    theme_bw()
print(plot)
dev.off()


##ROC
sampinfo=tibble(
    samps=c('neb14', 'neb15', 'neb17', 'neb19', 'nebdcm'),
    motifs=c('GATC', 'GCNGC', 'GANTC', 'GATC', 'CCWGG'),
    basepos=c('3', '1', '1', '1', '1'),
    base=c('C', 'C', 'A', 'A', 'C')
    )
    
rocinfo=tibble(
    thresh=as.numeric(),
    tpr=as.numeric(),
    fpr=as.numeric(),
    samp=as.character(),
    model=as.character(),
    base=as.character())
