import itertools
import argparse
import multiprocessing as mp
import pysam
import math
import numpy as np
from itertools import repeat
import math
from megalodon_barcode_functions import *


def parseArgs():
    parser=argparse.ArgumentParser(description='get a meth barcode per read for megalodon')
    parser.add_argument('-m', '--modfile', type=str, required=True, help='megalodon output file with mod motif probs')
    parser.add_argument('-i', '--idxfile', type=str, required=True, help='index of megalodon outputfile')
    parser.add_argument('-r', '--reffile', type=str, required=True, help='reference genome used in megaldon')
    parser.add_argument('-b', '--barcodefile', type=str, required=True, help='motif list file')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='output file that lists each read and each barcode number')
    parser.add_argument('-a', '--abound', type=float, required=False, help='A threshold')
    parser.add_argument('-c', '--cbound', type=float, required=False, help='C threshold')
    parser.add_argument('-n', '--numreads', type=int, required=False, help='how many reads to consider')
    parser.add_argument('-t', '--threads', type=int, required=True, help='number of threads to use')
    args=parser.parse_args()
    return args

class logmodCalls:
    '''
    package the mod calls
    '''
    def __init__(self, readname, chrname, strand, calledmotifs, barcodes, motifcounts):
        self.readname=readname
        self.chrname=chrname
        self.strand=strand
        self.calledmotifs=calledmotifs
        self.motifcounts=motifcounts
        self.bc_counts=self.assign_barcode(barcodes, strand)
        self.bc_norm=self.norm_barcode(barcodes)
    def assign_barcode(self, barcodes, strand):
        '''
        take called motifs
        count how many have each motif
        '''
        bc_counts={}
        for i in barcodes:
            bc_counts[i]=0
        for i in self.calledmotifs:
            for j in barcodes:
                for k in barcodes[j]:
                    if strand=='-':
                        seq=revcomp(k)
                    else:
                        seq=k
                    if seq in i:
                        bc_counts[j]+=1
                        break
        return bc_counts
    def norm_barcode(self, barcodes):
        '''
        normalize motif counts in called motifs
        returns counts as a proportion of the number you'd expect randomly 
        '''
        bc_norm={}
        if len(self.calledmotifs) > 0:
            numcalled=len(self.calledmotifs)
            lencalled=len(self.calledmotifs[0])
            for i in self.bc_counts:
                nummotifs=len(barcodes[i])
                pseudocount=float((self.bc_counts[i]+1) / (1+(lencalled-len(i)+1) * (4**(lencalled-len(i))) * numcalled * nummotifs / float(4**lencalled)))
                normconst=math.log10(pseudocount)
                bc_norm[i]=normconst
        else:
            for i in self.bc_counts:
                bc_norm[i]=0
        for i in self.bc_counts:
            if self.motifcounts[i]==0:
                bc_norm[i]=None
        return bc_norm


def log_call_read(read, thresh, ref, barcodes, k):
    '''
    take in a readMod object
    return a modCalls object
    '''
    ccallpos=np.array(read.cpmod) > thresh[0]
    acallpos=np.array(read.apmod) > thresh[1]
    poscalled=list(itertools.compress(read.cpos, ccallpos.tolist()))
    poscalled.extend(list(itertools.compress(read.apos, acallpos.tolist())))
    calledmotifs=[]
    for i in poscalled:
        calledmotifs.append(ref[i[0]][i[1]-k:i[1]+k])
    motifcounts=get_num_motifs(read.chrname, read.minpos, read.maxpos, ref, barcodes, read.strand)
    calls=logmodCalls(read.readname, read.chrname, read.strand, calledmotifs, barcodes, motifcounts)
    return calls

def per_read_log(idxchunk, modfile, thresh, ref, barcodes, k, L):
    '''
    for a given read
    run through the pipe, but with log scores
    '''
    holdreads=[]
    for readinfo in idxchunk:
        read=grab_read(modfile, readinfo[0], readinfo[1], readinfo[2], readinfo[3], readinfo[4])
        calls=log_call_read(read, thresh, ref, barcodes, k)
        callinfo=[calls.readname, calls.chrname]
        for i in barcodes:
            callinfo.append(str(calls.bc_norm[i]))
        holdreads.append('\t'.join(callinfo)+'\n')
    L+=holdreads 
    
def main(reffile, modfile, idxfile, barcodefile, outfile, abound, cbound, threads):
    ref=fasta_dict(reffile)
    readidx=read_megalodon_index(idxfile)

    ##everything keeps in the order of the barcodes
    barcodes=expand_barcodes(barcodefile)
    if abound is not None and cbound is not None:
        thresh=[cbound, abound]
    else:
        thresh=find_thresh()
    k=4

    ##split idx into chunks
    idxchunks=[]
    chunksize=math.ceil(len(readidx)/threads)
    for i in range(0,threads):
        start=i*chunksize
        end=(i+1)*chunksize
        if len(readidx)<end:
            idxchunks.append(readidx[start:])
        else:
            idxchunks.append(readidx[start:end])
    
    manager=mp.Manager()
    L=manager.list()
    pool=mp.Pool(threads)
    pool.starmap(per_read_log, zip(idxchunks, repeat(modfile), repeat(thresh), repeat(ref), repeat(barcodes), repeat(k), repeat(L)))
    with open (outfile, 'w') as f:
        for i in L:
            f.write(i)
    pool.close()
    pool.join()

if __name__ == "__main__":
    args=parseArgs()
    main(args.reffile, args.modfile, args.idxfile, args.barcodefile, args.outfile, args.abound, args.cbound, args.threads)
