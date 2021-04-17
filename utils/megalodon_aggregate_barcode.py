import argparse
import multiprocessing as mp
import pysam
import math
import numpy as np
import itertools
from megalodon_barcode import *
import re


class modAgg:
    '''
    aggregate mod pos info for a read
    '''
    def __init__(self, read, bcpos):
        self.readname=read.readname
        self.atot=self.get_a_totals(read, bcpos)
        self.ctot=self.get_c_totals(read, bcpos)
        self.combined=self.get_combined()
    def get_a_totals(self, read, bcpos):
        '''
        take read, and get A averages for each barcode
        '''
        atot={}
        for i in bcpos:
            fil=np.array([x in bcpos[i] for x in read.apos])
            aprobs=list(itertools.compress(read.apmod, fil))
            if len(aprobs)>0:
                atot[i]=sum(aprobs)/float(len(aprobs))
            else:
                atot[i]='none'
        return atot
    def get_c_totals(self, read, bcpos):
        '''
        take read, and get C averages for each barcode
        '''
        ctot={}
        for i in bcpos:
            fil=np.array([x in bcpos[i] for x in read.cpos])
            cprobs=list(itertools.compress(read.cpmod, fil))
            if len(cprobs)>0:
                ctot[i]=sum(cprobs)/float(len(cprobs))
            else:
                ctot[i]='none'
        return ctot
    def get_combined(self):
        '''
        get combined averages
        '''
        combined={}
        for i in self.atot:
            if self.ctot[i]!='none' and self.atot[i]!='none':
                combined[i]=self.atot[i]+self.ctot[i]/2
            else:
                combined[i]='none'
        return combined

        
def find_motifs(ref, barcodes):
    '''
    for each barcode get positions as list of lists
    barcode:[[refchr, refpos], [refchr, refpos]]
    '''
    barcode_locations={}
    for i in barcodes:
        motiflen=len(i)
        barcode_locations[i]=[]
        for j in barcodes[i]:
            for chrom in ref:
                for pos in re.finditer(j, ref[chrom]):
                    barcode_locations[i].extend([[chrom, x] for x in range(pos.start(), pos.start()+motiflen)])
    final_locations={}
    for i in barcode_locations:
        final_locations[i]=[]
        for loc in barcode_locations[i]:
            if ref[loc[0]][loc[1]] == 'A' or ref[loc[0]][loc[1]] == 'C':
                final_locations[i].append(loc)
    return final_locations
            

def agg_per_read(readinfo, modfile, ref, bcpos, q):
    '''
    for a given read
    get the barcode scores
    '''
    read=grab_read(modfile, readinfo[0], readinfo[1], readinfo[2])
    agg=modAgg(read, bcpos)
    callinfo=[read.readname]
    for i in bcpos:
        callinfo.append(str(agg.combined[i]))
    q.put('\t'.join(callinfo)+'\n')
        

def main(reffile, modfile, idxfile, barcodefile, outfile, threads, numreads):
    '''
    a lot of this is the same in megalodon_barcode
    '''
    ref=fasta_dict(reffile)
    barcodes=expand_barcodes(barcodefile)
    bcpos=find_motifs(ref, barcodes)
    readidx=read_megalodon_index(idxfile)
    
    manager=mp.Manager()
    q=manager.Queue()
    pool=mp.Pool(threads)

    watcher=pool.apply_async(listener, (q, outfile))

    jobs=[]
    readcounts=0
    for i in readidx:
        if readcounts < numreads:
            readcounts+=1
            job=pool.apply_async(agg_per_read, (i, modfile, ref, bcpos, q))
            jobs.append(job)

    for job in jobs:
        job.get()

    q.put('Done now, ty 4 ur service')
    pool.close()
    pool.join()


    
if __name__=='__main__':
    args=parseArgs()
    main(args.reffile, args.modfile, args.idxfile, args.barcodefile, args.outfile, args.threads, args.numreads)


