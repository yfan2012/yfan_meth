library(ggplot2)
library(tidyverse)
library(argparse)
library(foreach)

datadir='/scratch/groups/mschatz1/cpowgs/meth/'
prefixes=c('170906_neb14', '171003_neb16', '171005_neb12', '171012_neb13', '171019_neb17', '171019_neb19', '171020_neb11', '171020_neb15', '180628_neb_dcm')


colnames=c('filename','read_id','run_id', 'channel', 'start_time','duration','num_events','passes_filtering','template_start','num_events_template','template_duration','num_called_template','sequence_length_template','mean_qscore_template','strand_score_template','calibration_strand_genome_template','calibration_strand_identity_template','calibration_strand_accuracy_template','aligned_speed_bps_template')



sumall=foreach(i=prefixes, .combine=rbind) %dopar% {
    sumfile=paste0(datadir,i,'/seq_sum.txt')
    seqsums=read_tsv(sumfile, col_names=colnames)


    iname=tail(strsplit(i, '_', fixed=TRUE)[[1]],1)
    lensum=sum(seqsums$sequence_length_template)
    namesum=paste0(iname, ' ', as.character(round(lensum/1000000000, digits=2)), 'Gb')

    seqsums=seqsums %>%
        mutate(samp=namesum)
    
    lengthinfo=seqsums[,c('samp','sequence_length_template')]
    return(lengthinfo)
}

pdf('~/scratch/plots/meth/read_length_hist.pdf', width=20, height=8.5)
ggplot(sumall, aes(x=samp, y=sequence_length_template, fill=samp)) +
    geom_violin() +
    xlab('Sample') +
    ylab('Length') +
    ylim(0, 50000) +
    ggtitle('Read Length') +
    guides(fill=FALSE) +
    theme_bw()
dev.off()
