import itertools
import argparse
import multiprocessing as mp
import pysam
import math
import sys
import numpy as np
from itertools import repeat
from megalodon_barcode_functions import *
import time

def parseArgs():
    parser=argparse.ArgumentParser(description='get a meth barcode per read for megalodon')
    parser.add_argument('-m', '--modfile', type=str, required=True, help='megalodon output file with mod motif probs')
    parser.add_argument('-i', '--idxfile', type=str, required=True, help='index of megalodon outputfile')
    parser.add_argument('-r', '--reffile', type=str, required=True, help='reference genome used in megaldon')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='output file that lists each read and each barcode number')
    parser.add_argument('-t', '--threads', type=int, required=True, help='number of threads to use')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose outputs')
    args=parser.parse_args()
    return args


def get_empty_ref(ref):
    '''
    take ref dict
    return ref dict with np array full of zeros
    '''
    aggcalls={}
    for i in ref:
        aggcalls[i]=np.zeros(len(ref[i]))
    return aggcalls


def get_reads(readinfo, modfile):
    '''
    take read info from read index
    return batch of reads
    '''
    try:
        startbyte=readinfo[0][2]
    except IndexError:
        print(readinfo)
    endbyte=readinfo[-1][2]+readinfo[-1][3]
    bytelen=endbyte-startbyte
    readbatch=[]
    with open(modfile, 'r') as f:
        f.seek(startbyte,0)
        readcontent=f.read(bytelen).split('\n')
        f.close()
    for i in readcontent:
        readbatch.append(i.split('\t'))
    return readbatch
    
        
def aggregate_reads(idxchunk, ref, modfile, thresh, L):
    '''
    take a batch of reads (in idxchunk)
    return aggregated read methylation calls
    '''
    aggmeth=get_empty_ref(ref)
    aggunmeth=get_empty_ref(ref)
    readbatch=get_reads(idxchunk, modfile)
    for i in readbatch[0:-1]:
        ratio=math.log10(float(i[5])/float(i[4]))
        pos=float(i[3])
        base=ref[i[1]][int(i[3])]
        if base=='C' or base=='G':
            if ratio > thresh[0]:
                aggmeth[i[1]][int(i[3])]+=1
            else:
                aggunmeth[i[1]][int(i[3])]+=1
        if base=='A' or base=='T':
            if ratio > thresh[1]:
                aggmeth[i[1]][int(i[3])]+=1
            else:
                aggunmeth[i[1]][int(i[3])]+=1
    L.append([aggmeth, aggunmeth])

    
def pre_sum_refs(aglist, ref, L):
    '''
    take list of aggregated meth calls
    combine them into one final thing
    '''
    fullmeth=get_empty_ref(ref)
    for i in aglist:
        for j in i:
            fullmeth[j]=np.add(i[j], fullmeth[j])
    L.append(fullmeth)

    
def sum_refs(aglist, ref):
    '''
    take list of aggregated meth calls
    combine them into one final thing
    '''
    fullmeth=get_empty_ref(ref)
    for i in aglist:
        for j in i:
            fullmeth[j]=np.add(i[j], fullmeth[j])
    return fullmeth


def main(reffile, modfile, idxfile, outfile, threads, verbose):
    start_time=time.time()
    if verbose:
        print('reading reference and index files, splitting into jobs')
    ref=fasta_dict(reffile)
    readidx=read_megalodon_index(idxfile)
    ##potentially add custom thresholds later
    thresh=find_thresh()

    ##split idx into chunks
    idxchunks=[]
    chunksize=3000
    numchunks=math.ceil(len(readidx)/3000)
    for i in range(numchunks):
        start=i*chunksize
        end=(i+1)*chunksize
        if len(readidx)<end:
            if len(readidx[start:]) > 0:
                idxchunks.append(readidx[start:])
        else:
            if len(readidx[start:end]) > 0:
                idxchunks.append(readidx[start:end])

    if verbose:
        read_duration=round(time.time()-start_time)
        print('finished reading in %d seconds' % (read_duration))
        print('starting parallel')
        par_start=time.time()

    manager=mp.Manager()
    L=manager.list()
    pool=mp.Pool(threads)
    pool.starmap(aggregate_reads, zip(idxchunks, repeat(ref), repeat(modfile), repeat(thresh), repeat(L)))
    pool.close()

    '''
    results=[]
    args=list(zip(idxchunks, repeat(ref), repeat(modfile), repeat(thresh)))
    Parallel(n_jobs=threads)(delayed(aggregate_reads)(idxchunk, ref, modfile, thresh) for idxchunk, ref, modfile, thresh in args)
    '''
            
    if verbose:
        par_duration=round(time.time()-par_start)
        print('finished parallel in %d seconds' % (par_duration))
        print('starting aggregation')
        agg_start=time.time()


    am=[item[0] for item in L]
    au=[item[1] for item in L]

    ##split so aggregation is parallel
    amsplit=[]
    ausplit=[]
    splitsize=math.ceil(len(am)/threads)
    for i in range(threads):
        start=i*splitsize
        end=(i+1)*splitsize
        if len(am)<end:
            if len(am[start:]) > 0:
                amsplit.append(am[start:])
                ausplit.append(au[start:])
        else:
                amsplit.append(am[start:end])
                ausplit.append(au[start:end])
                
    manager=mp.Manager()
    mL=manager.list()
    uL=manager.list()
    pool=mp.Pool(threads)
    pool.starmap(pre_sum_refs, zip(amsplit, repeat(ref), repeat(mL)))
    pool.close()
    
    pool=mp.Pool(threads)
    pool.starmap(pre_sum_refs, zip(ausplit, repeat(ref), repeat(uL)))
    pool.close()
    
    fullmeth=sum_refs(mL, ref)
    fullunmeth=sum_refs(uL, ref)

    if verbose:
        agg_duration=round(time.time()-agg_start)
        print('finished aggregation in %d seconds' % (agg_duration))
        print('starting writing')
        write_start=time.time()
        

    with open (outfile, 'w') as f:
        for i in fullmeth:
            for pos in range(0, len(fullmeth[i])):
                towrite=[i, str(pos), str(fullmeth[i][pos]), str(fullunmeth[i][pos])]
                f.write('\t'.join(towrite)+'\n')

    if verbose:
        write_duration=round(time.time()-write_start)
        print('finished writing in %d seconds' % (write_duration))

if __name__ == "__main__":
    args=parseArgs()
    main(args.reffile, args.modfile, args.idxfile, args.outfile, args.threads, args.verbose)
