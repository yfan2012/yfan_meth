library(tidyverse)
library(cowplot)
library(foreach)
library(doParallel)
cl=makeCluster(10)
registerDoParallel(cl, cores=10)


datadir='/mithril/Data/Nanopore/projects/methbin/rerio/'
##datadir='~/data/rerio/'
dbxdir='~/Dropbox/yfan/methylation/methbin/rerio/'

amods=c('GATC', 'GANTC') ##Z
cmods=c('CCWGG', 'GATC', 'GCNGC') ##Y
cgmods=c('CG', 'GC') ##Y                  


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
        megainfo=read_tsv(megafile, n_max=5000000) %>%
            mutate(modratio=log(can_log_prob/mod_log_prob)) %>%
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


plot_mega_control <- function(samp, motif, mod, mode) {
    if (mode=='vanilla') {
        megafile=paste0(datadir, '/', samp, '/megalodon/', samp, '_vanilla_', motif, '_', mod, '/per_read_modified_base_calls.txt')
        controlfile=paste0(datadir, '/neb11/megalodon/neb11_vanilla_', motif, '_', mod, '/per_read_modified_base_calls.txt')
        title=paste0(samp, ' old model')
    }else{
        megafile=paste0(datadir, '/', samp, '/megalodon/', samp, '_', motif, '_', mod, '/per_read_modified_base_calls.txt')
        controlfile=paste0(datadir, '/neb11/megalodon/neb11_', motif, '_', mod, '/per_read_modified_base_calls.txt')
        title=paste0(samp, ' rerio')
    }
    megainfo=read_tsv(megafile, n_max=5000000) %>%
        mutate(modratio=log(can_log_prob/mod_log_prob)) %>%
        mutate(samp=samp)
    controlinfo=read_tsv(controlfile, n_max=5000000) %>%
        mutate(modratio=log(can_log_prob/mod_log_prob)) %>%
        mutate(samp='neb11')
    allmegainfo=rbind(megainfo, controlinfo)
    

    plot=ggplot(allmegainfo, aes(x=modratio, colour=samp, fill=samp, alpha=.3)) +
        geom_density() +
        ggtitle(title) +
        xlab('modprob ratio') +
        theme_bw()
    
    return(plot)
}


mega_roc <- function(modfile, unmodfile) {
    mod=read_tsv(modfile, n_max=15000000) %>%
        mutate(modratio=log(can_log_prob/mod_log_prob)) %>%
        mutate(meth=TRUE)
    unmod=read_tsv(unmodfile, n_max=15000000) %>%
        mutate(modratio=log(can_log_prob/mod_log_prob)) %>%
        mutate(meth=FALSE)
    
    all=rbind(mod, unmod)
    thresholds=seq(min(all$modratio)-.05, max(all$modratio)+.05, .05)
    
    pos=sum(all$meth)
    neg=sum(!all$meth)

    
    confusion=foreach(i=1:length(thresholds), .combine=rbind) %dopar% {
        thresh=thresholds[i]

        call=all$modratio > thresh
        tp=call & all$meth
        fp=call & !all$meth

        tpr=sum(tp)/pos
        fpr=sum(fp)/neg
        conf=data.frame(thresh=thresh, tpr=tpr, fpr=fpr)
        return(conf)
    }

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


##standards
ecolicg=plot_mega('ecoli_CpG', cgmods, 'Z', '')
ecolicgv=plot_mega('ecoli_CpG', cgmods, 'Z', 'vanilla')

outfile=paste0(dbxdir, 'rerio_mega_ecolistandards.pdf')
pdf(outfile, h=5, w=13)
plot_grid(ecolicg, ecolicgv, ncol=2, align='v')
dev.off()

ecoligc=plot_mega('ecoli_GpC', cgmods, 'Z', '')
ecolicgv=plot_mega('ecoli_GpC', cgmods, 'Z', 'vanilla')


##amods vs neb11, true motifs
neb19=plot_mega_control('neb19', 'GATC', 'Y', '')
neb19v=plot_mega_control('neb19', 'GATC', 'Y', 'vanilla')
neb17=plot_mega_control('neb17', 'GANTC', 'Y', '')
neb17v=plot_mega_control('neb17', 'GANTC', 'Y', 'vanilla')
outfile=paste0(dbxdir, 'rerio_mega_Amod_v_control.pdf')
pdf(outfile, h=6, w=13)
plot_grid(neb19, neb19v, neb17, neb17v, ncol=2, align='v')
dev.off()


##cmods vs neb11, true motifs
neb14=plot_mega_control('neb14', 'GATC', 'Z', '')
neb14v=plot_mega_control('neb14', 'GATC', 'Z', 'vanilla')
neb15=plot_mega_control('neb15', 'GCNGC', 'Z', '')
neb15v=plot_mega_control('neb15', 'GCNGC', 'Z', 'vanilla')
nebdcm=plot_mega_control('nebdcm','CCWGG' , 'Z', '')
nebdcmv=plot_mega_control('nebdcm', 'CCWGG', 'Z', 'vanilla')
outfile=paste0(dbxdir, 'rerio_mega_Cmod_v_control.pdf')
pdf(outfile, h=9, w=13)
plot_grid(neb14, neb14v, neb15, neb15v, nebdcm, nebdcmv, ncol=2, align='v')
dev.off()


##ROC
samptable=tibble(samp=c('neb11', 'neb14', 'neb15', 'neb17', 'neb19', 'nebdcm')) %>%
    mutate(motif=c('', 'GATC', 'GCNGC', 'GANTC', 'GATC', 'CCWGG')) %>%
    mutate(mod=c('', 'Z', 'Z', 'Y', 'Y', 'Z')) %>%
    mutate(vmodfile=paste0(datadir, samp, '/megalodon/', samp, '_vanilla_', motif, '_', mod, '/per_read_modified_base_calls.txt')) %>%
    mutate(vunmodfile=paste0(datadir, 'neb11/megalodon/neb11_vanilla_', motif, '_', mod, '/per_read_modified_base_calls.txt')) %>%
    mutate(modfile=paste0(datadir, samp, '/megalodon/', samp, '_', motif, '_', mod, '/per_read_modified_base_calls.txt')) %>%
    mutate(unmodfile=paste0(datadir, 'neb11/megalodon/neb11_', motif, '_', mod, '/per_read_modified_base_calls.txt'))

rocinfo=tibble(
    thresh=as.numeric(),
    tpr=as.numeric(),
    fpr=as.numeric(),
    samp=as.character(),
    model=as.character(),
    base=as.character())

for (i in 2:dim(samptable)[1]) {
    neb=samptable[i,]
    samp=neb$samp
    
    nebmega=as_tibble(mega_roc(neb$modfile, neb$unmodfile)) %>%
        mutate(samp=samp) %>%
        mutate(model='rerio') %>%
        mutate(base=neb$mod)
    nebmegav=as_tibble(mega_roc(neb$vmodfile, neb$vunmodfile)) %>%
        mutate(samp=samp) %>%
        mutate(model='vanilla') %>%
        mutate(base=neb$mod)

    rocinfo=rbind(rocinfo, nebmega, nebmegav)
}

outfile=paste0(dbxdir, 'mega_roc.pdf')
pdf(outfile, w=13, h=9)
plot=ggplot(rocinfo, aes(x=fpr, y=tpr, colour=samp, linetype=model)) +
    geom_step() +
    ggtitle('Megalodon ROC') +
    theme_bw()
print(plot)
dev.off()

##estimate auc
auc=rocinfo %>%
    group_by(samp, model) %>%
    mutate(differ=c(0,diff(fpr))) %>%
    mutate(prod=abs(differ)*tpr) %>%
    summarise(auc=sum(prod))


