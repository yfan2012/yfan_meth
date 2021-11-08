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


def get_strand(ref):
    '''
    take the refernce 
    return strand info
    '''
    strands=get_empty_ref(ref)
    for i in ref:
        for j in range(len(ref[i])):
            if ref[i][j]=='A' or ref[i][j]=='C':
                strands[i][j]+=1
    return strands

def get_reads(idxchunks, q, threads):
    '''
    take read info from read index
    return batch of reads
    '''
    for readinfo in idxchunks:
        startbyte=readinfo[0][2]
        endbyte=readinfo[-1][2]+readinfo[-1][3]
        bytelen=endbyte-startbyte
        readbatch=[]
        with open(modfile, 'r') as f:
            f.seek(startbyte,0)
            readcontent=f.read(bytelen).split('\n')
            f.close()
        for i in readcontent:
            readbatch.append(i.split('\t'))
        q.put(readbatch)
    for i in range(threads):
        q.put('done now')
    
        
def aggregate_reads(q, r):
    '''
    take a batch of reads (in idxchunk)
    return aggregated read methylation calls
    '''
    while True:
        readbatch=q.get()
        if readbatch=='done now':
            r.put('done now')
            break
        
        aggmeth=get_empty_ref(ref)
        aggunmeth=get_empty_ref(ref)

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
        r.put([aggmeth, aggunmeth])
    r.put('done now')


def gather_results(r, aggmeth, aggunmeth, threads):
    '''
    get results from queue and put them in lists
    '''
    dones=0
    while True:
        res=r.get()
        if res == 'done now':
            dones+=1
            if dones==threads:
                break
        else:
            aggmeth.append(res[0])
            aggunmeth.append(res[1])
    return aggmeth, aggunmeth


def sum_refs(agg):
    '''
    take list of aggregated meth calls
    combine them into one final thing
    '''
    fullmeth=get_empty_ref(ref)
    for i in agg:
        for chrom in fullmeth:
                fullmeth[chrom]=np.add(i[chrom], fullmeth[chrom])
    return fullmeth



if __name__ == "__main__":
    args=parseArgs()
    reffile=args.reffile
    modfile=args.modfile
    idxfile=args.idxfile
    threads=args.threads
    verbose=args.verbose

    start_time=time.time()
    if verbose:
        print('reading reference and index files, splitting into jobs')
    ref=fasta_dict(reffile)
    strands=get_strand(ref)
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
    q=manager.Queue()
    r=manager.Queue()
    methagg=[]
    unmethagg=[]

    processes=[]
    processes.append(mp.Process(target=get_reads, args=(idxchunks, q, threads)))
    for i in range(threads):
        processes.append(mp.Process(target=aggregate_reads, args=(q, r)))
    processes.append(mp.Process(target=gather_results, args=(r, methagg, unmethagg, threads)))

    for process in processes:
        process.start()

    for process in processes:
        process.join()

    q.close()
    r.close()


        
    if verbose:
        par_duration=round(time.time()-par_start)
        print('finished parallel in %d seconds' % (par_duration))
        print('starting aggregation')
        agg_start=time.time()

    fullmeth=sum_refs(methagg)
    fullunmeth=sum_refs(unmethagg)

    if verbose:
        agg_duration=round(time.time()-agg_start)
        print('finished aggregation in %d seconds' % (agg_duration))
        print('starting writing')
        write_start=time.time()
        

    with open (outfile, 'w') as f:
        for i in fullmeth:
            for pos in range(0, len(fullmeth[i])):
                towrite=[i, str(pos), str(strands[i][j]), str(fullmeth[i][pos]), str(fullunmeth[i][pos])]
                f.write('\t'.join(towrite)+'\n')

    if verbose:
        write_duration=round(time.time()-write_start)
        print('finished writing in %d seconds' % (write_duration))
