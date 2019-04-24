##Try in R
library(seqLogo)
library(foreach)
library(doParallel)
cl=makeCluster(10)
registerDoParallel(cl, cores=10)

findmotif <- function(file_name, W, iterations=100){
    seqs=as.vector(as.matrix(read.table(file_name, sep=",", stringsAsFactors=FALSE)))
    bases=do.call(rbind, strsplit(seqs, split=""))

    ##Initialize Pck, Bc, prob of each nucleotide at each position in the motif and background
    Pck=matrix(.25, nrow=4, ncol=W)
    rownames(Pck)=c('A', 'T', 'C', 'G')
    Bc=matrix(.25, nrow=4, ncol=1)
    rownames(Bc)=c('A', 'T', 'C', 'G')
    Mc=matrix(c(sum(bases=='A'), sum(bases=='T'), sum(bases=='C'), sum(bases=='G')), nrow=4, ncol=1)
    rownames(Mc)=c('A', 'T', 'C', 'G')
    
    L=dim(bases)[2]

    iter=0
    loglike=c()
    eps=1
    while ((eps>.001) && (iter<iterations)){

        ##Estep
        pregam=foreach(i=1:dim(bases)[1], .combine=rbind) %dopar% {
            ##initialize the row to be returned
            pregam_i=matrix(0, nrow=1, ncol=L-W+1)
            for (j in 1:(L-W+1)) {
                if (j==1) {
                    beginprob=1
                }else{
                    beginseq=bases[i,1:j-1]
                    startprobs=Bc[beginseq,]
                    beginprob=prod(startprobs)
                }
                if (j==(L-W+1)) {
                    endprob=1
                }else{
                    endseq=bases[i,(j+W):L]
                    lastprobs=Bc[endseq,]
                    endprob=prod(lastprobs)
                }
                motifseq=bases[i,j:(j+W-1)]
                motifprob=prod(Pck[motifseq,])
                
                pregam_i[j]=beginprob*motifprob*endprob
            }
            return(pregam_i)
        }
        gamma=pregam/rowSums(pregam)

        ##Mstep======================================

        ##Get Nck from gamma and bases
        Nck=foreach(i=c('A', 'T', 'C', 'G'), .combine=rbind) %dopar% {
            Nck_i=matrix(0, nrow=1, ncol=W)

            for (k in 1:W) {
                pos=which(bases[,k:(L-W+1)]==i, arr.ind=TRUE)
                Nck_i[k]=sum(gamma[pos])
            }
            return(Nck_i)
        }
        rownames(Nck)=c('A', 'T', 'C', 'G')

        ##Get Pck from Nck
        Pck=sweep((Nck+1),2,colSums(Nck+1), '/')

        ##Get Gc from Mc and Nck
        Gc=Mc-rowSums(Nck)

        Bc=(Gc+1)/sum(Gc+1)

        ##Log like
        loglike=c(loglike,sum(log(rowSums(pregam))))
        if (length(loglike)>1){
            eps=abs(loglike[length(loglike)]-loglike[length(loglike)-1])
        }

        iter=iter+1
        print(iter)
    }
    
    return(Pck)
}



datadir='/scratch/groups/mschatz1/cpowgs/meth/'
prefixes=c('170906_neb14', '171003_neb16', '171005_neb12', '171012_neb13', '171019_neb17', '171019_neb19', '171020_neb11', '171020_neb15', '180628_neb_dcm')

for (i in prefixes) {
    seqlist=paste0(datadir,i,'/asm_diffs/',i,'.12merlist.txt')
    if (file.exists(seqlist)) {
        pck6=findmotif(seqlist, 6, 100)
        pck5=findmotif(seqlist, 5, 100)
        pck4=findmotif(seqlist, 4, 100)

        pdf(paste0(datadir,i,'/asm_diffs/',i,'.motifs.pdf'))
        seqLogo(pck6[c('A', 'C', 'G', 'T'),], ic.scale=FALSE)
        seqLogo(pck5[c('A', 'C', 'G', 'T'),], ic.scale=FALSE)
        seqLogo(pck4[c('A', 'C', 'G', 'T'),], ic.scale=FALSE)
        dev.off()
    }
}
