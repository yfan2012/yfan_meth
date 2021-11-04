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
    parser.add_argument('-o', '--outfile', type=str, required=True, help='output file that lists each read and each barcode number')
    parser.add_argument('-t', '--threads', type=int, required=True, help='number of threads to use')
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
    return readbatch
    
    
        
def aggregate_reads(idxchunk, ref, modfile, thresh, am, au):
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
    am.append(aggmeth)
    au.append(aggunmeth)


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


def main(reffile, modfile, idxfile, outfile, threads):
    ref=fasta_dict(reffile)
    readidx=read_megalodon_index(idxfile)

    ##potentially add custom thresholds later
    thresh=find_thresh()


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
    am=manager.list()
    au=manager.list()
    pool=mp.Pool(threads)
    pool.starmap(aggregate_reads, zip(idxchunks, repeat(ref), repeat(modfile), repeat(thresh), repeat(am), repeat(au)))

    fullmeth=sum_refs(am, ref)
    fullunmeth=sum_refs(au, ref)

    with open (outfile, 'w') as f:
        for i in fullmeth:
            chrom=i
            for pos in fullmeth[i]:
                towrite=[i, str(pos), str(fullmeth[i][pos]), str(fullunmeth[i][pos])]
                f.write('\t'.join(towrite)+'\n')

    pool.close()
    pool.join()

if __name__ == "__main__":
    args=parseArgs()
    main(args.reffile, args.modfile, args.idxfile, args.outfile, args.threads)
