library(gridExtra)
library(tidyverse)

datadir='~/Dropbox/yfan/methylation/neb_meth/asm_diffs/counts/'
##prefixes=c('170906_neb14', '171003_neb16', '171005_neb12', '171012_neb13', '171019_neb17', '171019_neb19', '171020_neb11', '171020_neb15', '180628_neb_dcm')
prefixes=c('170906_neb14', '171003_neb16', '171005_neb12', '171019_neb17', '171019_neb19', '171020_neb11', '171020_neb15', '180628_neb_dcm')

for (i in prefixes) {
    print(i)
    enrichfile=paste0(datadir, i, '.ref.6mer.hpfilt.sorted.csv')


    enrich_cols=c('motif', '#in_errors', '#in_genome', 'proportion')
    enrich=read_csv(enrichfile, enrich_cols)

    outfile=paste0(datadir, i, '.ref.top_6mers.pdf')
    pdf(outfile)
    grid.table(enrich[1:5,], rows=NULL)
    dev.off()
}
