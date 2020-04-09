library(tidyverse)
library(doParallel)
library(foreach)

cl=makeCluster(12)
registerDoParallel(cl, cores=12)

datadir='/uru/Data/Nanopore/projects/methbin/origin/'
##samps=c('neb12', 'neb14', 'neb15', 'neb19')
samps=c('neb14')
columns=c('read', 'seqname', 'position', 'kmer', 'pval')
outfile='~/Dropbox/yfan/methylation/methbin/genome_pvals.pdf'

pdf(outfile, height=5, width=16)
for (samp in samps) {
    file=paste0(datadir, samp, '/', samp, '.positionpvals.tsv')

    pvals=read_tsv(file, col_names=columns)
    lims=seq(1000, max(pvals$position), 1000)


    distpvals=foreach(i=lims, .combine=rbind) %dopar% {
        library(tidyverse)
        min=i-1000
        max=i
        inpos=pvals %>%
            filter(position<max & position>min & seqname=='gi|730582171|gb|CP009644.1|') %>%
            mutate(logpval=log(pval, base=10)) %>%
            mutate(logrecip=log(1-pval, base=10))
        if (dim(inpos)[1] > 500) {
            summary=tibble(pos=abs(3897400-max), ratio=sum(inpos$logpval)/sum(inpos$logrecip))
            return(summary)
        }
    }
        
    print(ggplot(distpvals, aes(x=pos, y=ratio, alpha=.2)) +
          geom_bin2d(bins=c(100, 20)) +
          ggtitle(samp) +
          xlab('distance from oriC') +
          ylab('pval ratio') +
          theme_bw()
          )
}
dev.off()
