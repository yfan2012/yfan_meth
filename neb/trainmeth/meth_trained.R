library(tidyverse)
library(ggridges)
library(ggjoy)
library(cowplot)

##Compare current distros of all relevant kmers
datadir='/scratch/groups/mschatz1/cpowgs/meth/'
plotdir='~/plots/meth'

##excluding hinfI 
enzymes=data.frame(enzymes=c('pspjdri', 'fnu4h', 'sin395', 'sdeaII' , 'dam'), motif=c('CCGG', 'GCNGC', 'GATC', 'CCNGGC', 'GATC'), meth=c('MCGG', 'GMNGC', 'GATM', 'CCNGGM', 'GMTC'), data=c('171005_neb12', '171020_neb15', '170906_neb14', '171003_neb16', '171019_neb19'))
enzymes$tsv=paste0(datadir,enzymes$data, '/models/methyltrain' , enzymes$enzymes, '.model.round4.events.tsv')
enzymes$type=c('4mc', '5mC', '5mC', '5mC', '6mA')


##Get data for meth
allmeth=read_tsv(enzymes$tsv[1], n_max=5000000)
add_column(meth, type=enzymes$type[1])
for (i in 2:length(enzymes$tsv)) {
    meth=read_tsv(enzymes$tsv[i], n_nmax=5000000) 
    add_column(meth, type=enzymes$type[i])

    allmeth=rbind(allmeth, meth)
}


##all=all %>%
##select(model_kmer, level_mean, type)


##Get model for later
modelpath='~/software/timp_nanopolish/etc/r9-models/r9.4_450bps.nucleotide.6mer.template.model'
model=read_tsv(modelpath, comment='#')




##Make ultimate key for all enzyme motifs
motifkey=tibble(enzyme=character(), motif=character(), methmotif=character())
for (i in c(1:6)) {
    if (grepl('N', enzymes$motif[i])) {
        for (j in c('A','T','C', 'G')) {
            expanded_motif=gsub('N',j , enzymes$motif[i])
            expanded_methmotif=gsub('N',j, enzymes$meth[i])
            enzymekey=tibble(enzyme=enzymes$enzymes[i], motif=expanded_motif, meth=expanded_methmotif)
            motifkey=rbind(motifkey, enzymekey)
        }        
    }else{
        enzymekey=tibble(enzyme=enzymes$enzymes[i], motif=enzymes$motif[i], meth=enzymes$meth[i])
        motifkey=rbind(motifkey, enzymekey)
    } 
}    



all=all %>%
    mutate(enzyme='unknown')
all=all %>%
    mutate(mstatus='unknown')
all=all %>%
    mutate(group='unknown')

for (i in 1:dim(motifkey)[1]) {
    unmeth=grepl(as.character(motifkey$motif[i]), all$model_kmer)
    meth=grepl(as.character(motifkey$meth[i]), all$model_kmer)
    all$enzyme[unmeth]=as.character(motifkey$enzyme[i])
    all$enzyme[meth]=as.character(motifkey$enzyme[i])
    all$mstatus[unmeth]='unmeth'
    all$mstatus[meth]='meth'
    all$group[unmeth]=as.character(motifkey$meth[i])
    all$group[meth]=as.character(motifkey$meth[i])
}

all=all %>%
    filter(enzyme!='unknown')

ameth=all %>%
    filter(enzyme=='hinfI' | enzyme=='dam')

four=all %>%
    filter(enzyme=='pspjdri')

five=all %>%
    filter(enzyme=='sin395' | enzyme=='fnu4h' | enzyme=='sdeaII')



##not really a thing
pdf(paste0(plotdir, '/distros.pdf'), width=11, height=8.5)
print(ggplot(all, aes(x=level_mean, y=model_kmer, fill=type))+geom_joy(alpha=.5))
dev.off()
pdf(paste0(plotdir, '/grouped_distros.pdf'), width=11, height=8.5)
print(ggplot(all, aes(x=level_mean, y=group, fill=type))+geom_joy(alpha=.4))
dev.off()


sub=all %>%
    filter(group=='MCGG')

pdf(paste0(plotdir, '/MCGG.pdf'), width=11, height=8.5)
print(ggplot(sub, aes(x=level_mean, y=group, fill=type))+geom_joy(alpha=.4))
dev.off()




four=all %>%
    filter(group=='MCGG')
methfour=all %>%
    filter(type=='pspjdri') %>%
    filter(grepl('M', model_kmer)) %>%
    mutate(mpos=regexpr('M', model_kmer))

four=all %>%
    mutate(kmer='unknown')
targets=c('ACCGGA', 'TCCGGA', 'TCCGGG','ATTCCG', 'TAGCCG', 'CCGGAG', 'CCGGCA', 'ACTCCG')
methtar=c('AMCGGA', 'TMCGGA', 'TMCGGG','ATTMCG', 'TAGMCG', 'MCGGAG', 'MCGGCA', 'ACTMCG')
four=four %>%
    filter(model_kmer %in% targets | model_kmer %in% methtar)


model_current=as.numeric()
for (i in 1:dim(methfour)[1]) {
    model_current[i]=model$level_mean[model$kmer==gsub('M', 'C', methfour$model_kmer[i])]
}
methfour=methfour %>%
    mutate(model_level=model_current) %>%
    mutate(diff=level_mean-model_level)

pdf(paste0(plotdir, '/4mc_diff.pdf'), width=11, height=8.5)
print(ggplot(methfour, aes(x=mpos, y=diff)) +
      geom_boxplot(aes(group=mpos)) +
      theme_bw())
dev.off()


for (i in 1:8) {
    four$kmer[four$model_kmer==targets[i]]=targets[i]
    four$kmer[four$model_kmer==methtar[i]]=targets[i]
}

foursub=four %>%
    filter(type=='pspjdri' | type=='unmeth' | type=='dam' | type=='sdeaII')
    
pdf(paste0(plotdir, '/4mc.pdf'), width=11, height=8.5)
print(ggplot(four[four$type!='unmeth',], aes(x=level_mean, y=kmer, fill=type))+geom_joy(alpha=.4))
print(ggplot(foursub, aes(x=level_mean, y=kmer, fill=type))+geom_joy(alpha=.4))
print(ggplot(foursub, aes(x=level_mean, y=kmer, fill=type))+geom_density_ridges(stat="binline", scale=0.95, draw_baseline=FALSE))
dev.off()



fivesites=c('GATM', 'GMAGC', 'GMCGC', 'GMGGC', 'GMTGC','CCAGGC', 'CCCGGC','CCGGGC', 'CCTGGC')
five=all %>%
    filter(group %in% fivesites)

methfive=all %>%
    filter(type=='sin395' | type=='fnu4h' | type=='sdeaII') %>%
    filter(grepl('M', model_kmer)) %>%
    mutate(mpos=regexpr('M', model_kmer))

five=all %>%
    mutate(kmer='unknown')
#targets=c('CCTGGC', 'CCAGGC', 'CCGGGC','GATCTT', 'GGATCG', 'TGCAGC','GCAGCG', 'AGCAGC')
#methtar=c('CCTGGM', 'CCAGGM', 'CCGGGM','GATMTT', 'GGATMG', 'TGMAGC','GMAGCG', 'AGMAGC')

targets=c('CCTGGC', 'CCAGGC', 'CCGGGC', 'CCCGGC', 'CAGGCT', 'CTGGCA', 'CGGGCT', 'TGGCAT')
methtar=c('CCTGGM', 'CCAGGM', 'CCGGGM', 'CCCGGM', 'CAGGMT', 'CTGGMA', 'CGGGMT', 'TGGMAT')
five=five %>%
    filter(model_kmer %in% targets | model_kmer %in% methtar)


model_current=as.numeric()
for (i in 1:dim(methfive)[1]) {
    model_current[i]=model$level_mean[model$kmer==gsub('M', 'C', methfive$model_kmer[i])]
}
methfive=methfive %>%
    mutate(model_level=model_current) %>%
    mutate(diff=level_mean-model_level)

pdf(paste0(plotdir, '/5mc_diff.pdf'), width=11, height=8.5)
print(ggplot(methfive, aes(x=mpos, y=diff)) +
      geom_boxplot(aes(group=mpos)) +
      theme_bw())
dev.off()


for (i in 1:8) {
    five$kmer[five$model_kmer==targets[i]]=targets[i]
    five$kmer[five$model_kmer==methtar[i]]=targets[i]
}

fivesub=five %>%
    filter(type=='unmeth' | type=='sdeaII' | type=='dam' | type=='pspjdri')
pdf(paste0(plotdir, '/5mc.pdf'), width=11, height=8.5)
print(ggplot(five[five$type!='unmeth',], aes(x=level_mean, y=kmer, fill=type))+geom_joy(alpha=.4))
print(ggplot(fivesub, aes(x=level_mean, y=kmer, fill=type))+geom_joy(alpha=.4))
print(ggplot(fivesub, aes(x=level_mean, y=kmer, fill=type))+geom_density_ridges(stat="binline", scale=0.95, draw_baseline=FALSE))
dev.off()



sixsites=c('GAATC', 'GACTC', 'GAGTC', 'GATTC', 'GMTC')
six=all %>%
    filter(group %in% sixsites)

methsix=all %>%
    filter(grepl('M', model_kmer)) %>%
    filter(type=='dam' | type=='hinfI') %>%
    mutate(mpos=regexpr('M', model_kmer))
    n
six=all %>%
    mutate(kmer='unknown')
targets=c('TGATTC', 'GATCAT', 'AGATCA','CGGATC', 'AGAGTC', 'GATCAA','AGTGAG', 'CCGGAA')
methtar=c('TGMTTC', 'GMTCAT', 'AGMTCA','CGGMTC', 'AGMGTC', 'GMTCAA','AGTGMG', 'CCGGMA')
six=six %>%
    filter(model_kmer %in% targets | model_kmer %in% methtar)


model_current=as.numeric()
for (i in 1:dim(methsix)[1]) {
    model_current[i]=model$level_mean[model$kmer==gsub('M', 'A', methsix$model_kmer[i])]
}
methsix=methsix %>%
    mutate(model_level=model_current) %>%
    mutate(diff=level_mean-model_level)

pdf(paste0(plotdir, '/6ma_diff.pdf'), width=11, height=8.5)
print(ggplot(methsix, aes(x=mpos, y=diff)) +
      geom_boxplot(aes(group=mpos)) +
      theme_bw())
dev.off()


for (i in 1:8) {
    six$kmer[six$model_kmer==targets[i]]=targets[i]
    six$kmer[six$model_kmer==methtar[i]]=targets[i]
}



pdf(paste0(plotdir, '/6ma.pdf'), width=11, height=8.5)
print(ggplot(six, aes(x=level_mean, y=kmer, fill=type))+geom_joy(alpha=.4))
dev.off()


##all_diffs by position
methfour=methfour %>%
    filter(type=='pspjdri') %>%
    mutate(tag="4mc")
methfive=methfive %>%
    filter(type=='sin395' |  type=='fnu4h' | type=='sdeaII') %>%
    mutate(tag="5mc")
methsix=methsix %>%
    filter(type=='dam' | type=='hinfI') %>%
    mutate(tag="6ma")

all_diff=rbind(methfour, methfive, methsix)
all_diff=all_diff %>%
    mutate(pos_label=as.character(mpos))

pdf(paste0(plotdir, '/all_diffs.pdf'), width=11, height=8.5)
print(ggplot(all_diff, aes(x=pos_label, y=diff, colour=tag)) +
      geom_boxplot() +
      theme_bw())
print(ggplot(all_diff, aes(x=diff, y=tag ,fill=pos_label)) +
      geom_joy(alpha=.4) +
      xlim(-20,20) +
      theme_bw())
print(ggplot(all_diff, aes(x=diff, y=pos_label,fill=tag)) +
      geom_joy(alpha=.4) +
      xlim(-20,20) +
      theme_bw())
dev.off()


system("scp ~/plots/meth/*.pdf  yfan@smaug.timplab.jhu.edu:/home/yfan/Dropbox/Lab/talks/2017_ncm/")
