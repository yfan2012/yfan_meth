import csv
import argparse
import sys
import multiprocessing as mp
import argparse
import re
import scipy
import math
from scipy import stats
import pysam
import time

parser=argparse.ArgumentParser(description='calculate event pvals')
parser.add_argument('-r', '--ref', type=str, required=True, help='reference aligned to for nanopolish')
parser.add_argument('-c', '--col', type=str, required=True, help='collapsed eventalign file')
parser.add_argument('-d', '--dists', type=str, required=True, help='model file')
parser.add_argument('-m', '--motif', type=str, required=True, help='mod motif. standard iupac nucleotide codes')
parser.add_argument('-o', '--out', type=str, required=True, help='output file for kmer pvals')
parser.add_argument('-p', '--pout', type=str, required=True, help='output file for pvals aggregated across each read')
parser.add_argument('-t', '--threads', type=int, required=True, help='number of threads to use')
args=parser.parse_args()

iupacnt={
    'N': ['A', 'C', 'G', 'T'],
    'W': ['A', 'T']
}

def expand_motifs(motifs):
    '''
    Take list of motifs
    return list of motifs with explicit nts
    '''
    ##check for any motifs
    ambig=0
    for i in iupacnt:
        if i in ''.join(motifs):
            ambig+=1
    if ambig==0 :
        return motifs
    else:
        newmotifs=[]
        for i in motifs:
            for j in iupacnt:
                if j in i:
                    newmotifs.extend([i.replace(j, k, 1) for k in iupacnt[j]])
        return expand_motifs(newmotifs)

    
def find_motifs(reffile, motifs, motiflen):
    '''
    need a fastafile and a list of motifs
    give dict with 'chrname':[list of positions where any motif occurs with flanking region]
    '''
    fa=pysam.FastaFile(reffile)
    tigs=fa.references
    fastadict={ x:fa.fetch(x) for x in tigs }
    motifpos={}
    for seq in fastadict:
        motifpos[seq]=[]
        for motif in motifs:
            for i in re.finditer(motif, fastadict[seq]):
                motifpos[seq].extend(range(i.start()-motiflen,i.start()+motiflen)) 
            ##motifpos[seq]+=[i.start() for i in re.finditer(motif, fastadict[seq])]
        for i in motifpos:
            motifpos[i]=list(set(motifpos[i]))
    return motifpos


def read_model(modelfile):
    '''
    grab the model file
    kmer:[mean,std]
    '''
    with open(modelfile) as f:
        content=f.read().split('\n')
    model={}
    for i in content:
        if len(i)>0:
            if i[0] != '#':
                model[i.split('\t')[0]]=[i.split('\t')[1], i.split('\t')[2]]
    return model


def read_index(indexfile):
    '''
    grab the index file
    read:[numkmers, byteoffset, bytelen]
    '''
    with open(indexfile) as f:
        content=f.read().split('\n')
    readidx={}
    for i in content[1:]:
        if len(i)>0:
            readidx[i.split('\t')[3]]=[int(i.split('\t')[4]), int(i.split('\t')[9]), int(i.split('\t')[10]), str(i.split('\t')[0])]
    return readidx


def per_read_pvals(collapsefile, readchrom, byteoffset, bytelen, motifpos, model, q, r):
    '''
    get collapsed file and index info
    return [readnum, chrom, refpos, kmer, pval] and [read_id, loglik]
    '''
    start_time=time.time()
    with open(collapsefile, 'r') as f:
        f.seek(byteoffset, 0)
        readcontent=f.read(bytelen).split('\n')
        f.close()
    kmervals=[]
    readnum=readcontent[0].split('\t')[0]
    readchr=readcontent[0].split('\t')[1]
    agg=0
    revagg=0
    strkmervals=''
    for i in readcontent[2:]:
        ##if the position is relevant, get info
        if int(i.split('\t')[0]) in motifpos[readchrom]:
            eventmean=float(i.split('\t')[6])
            zval=(eventmean-float(model[i.split('\t')[1]][0]))/float(model[i.split('\t')[1]][1])
            pval=scipy.stats.norm.sf(abs(zval))*2
            kmerinfo=[readnum, readchr, i.split('\t')[0], i.split('\t')[1], str(pval)]
            if float(pval) > 0:
                agg+=math.log(float(pval),10)
                revagg+=math.log(1-float(pval),10)
                strkmervals+='\t'.join(kmerinfo)+'\n'
            elif float(pval) == 0:
                agg+=-500
                revagg+=0
                strkmervals+='\t'.join(kmerinfo)+'\n'
            else:
                print('why is the pval ' + str(pval) + ' in read ' + readnum )
                strkmervals+='\t'.join(kmerinfo)+'\n'
    try:
        aggratio=agg/revagg
    except ZeroDivisionError:
        aggratio=agg/(revagg+.0001)
        print(str(readnum)+ ': sum of pvals is 0?')
    print(readnum+ ' took ' + str(time.time()-start_time) + ' seconds')
    q.put(strkmervals)
    r.put('\t'.join([str(readnum), str(aggratio)])+'\n')

    
def listener(q, outfile):
    '''
    writes from q, a manager.Queue()
    '''
    with open(outfile, 'w') as f:
        while True:
            m=q.get()
            if  m=='Done now, ty 4 ur service':
                break
            f.write(m)
            f.flush()

            
#https://stackoverflow.com/questions/13446445/python-multiprocessing-safely-writing-to-a-file
def main(reffile, motifs, collapsefile, modelfile, koutfile, poutfile):
    model=read_model(modelfile)
    motiflen=6
    allmotifs=expand_motifs(motifs)
    motifpos=find_motifs(reffile, allmotifs, motiflen)
    indexfile=collapsefile+'.idx'
    readidx=read_index(indexfile)

    manager=mp.Manager()
    q=manager.Queue()
    r=manager.Queue()
    pool=mp.Pool(args.threads)

    ##listener is like a capaciter for stuff that needs to be written to a file
    qwatcher=pool.apply_async(listener, (q,koutfile))
    rwatcher=pool.apply_async(listener, (r,poutfile))
    
    ##start jobs whose results will accumulate in the watcher
    jobs=[]
    for i in readidx:
        if readidx[i][0] > 2000:
            job=pool.apply_async(per_read_pvals, (collapsefile, readidx[i][3], readidx[i][1], readidx[i][2], motifpos, model, q, r))
            jobs.append(job)

    ##I guess this bit actually runs the jobs? idk
    for job in jobs:
        job.get()

    ##signal listeners to die
    q.put('Done now, ty 4 ur service')
    r.put('Done now, ty 4 ur service')
    pool.close()
    pool.join()


if __name__ == "__main__":
    main(args.ref, args.motif, args.col, args.dists, args.out, args.pout)
  
 
    
    
    
