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
    parser=argparse.ArgumentParser(description='pull out relevant positions from cx style report based on barcode')
    parser.add_argument('-r', '--reffile', type=str, required=True, help='reference genome used in megaldon')
    parser.add_argument('-b', '--barcodefile', type=str, required=True, help='list of motifs in the barcode')
    parser.add_argument('-m', '--modfile', type=str, required=True, help='per read mod calls form megalodon')
    parser.add_argument('-i', '--idxfile', type=str, required=True, help='index of per read mod calls from megalodon')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='output file')
    parser.add_argument('-t', '--threads', type=int, required=True, help='number of threads to use')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose outputs')
    args=parser.parse_args()
    return args


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
    

##idxfile='/mithril/Data/Nanopore/projects/methbin/zymo/megalodon/20190809_zymo_control/per_read_modified_base_calls.txt.small.idx'
##barcodefile='/home/yfan/Code/yfan_nanopore/mdr/zymo/barcodes_zymo_curated.txt'
##modfile='/mithril/Data/Nanopore/projects/methbin/zymo/megalodon/20190809_zymo_control/per_read_modified_base_calls.txt'
##reffile='/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa'


def main(reffile, barcodefile, modfile, idxfile, outfile, threads, verbose):
    start_time=time.time()

        
    ref=fasta_dict(reffile)
    barcodes=expand_barcodes(barcodefile)
    readidx=read_megalodon_index(idxfile)
    

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


    manager=mp.Manager()
    L=manager.list()
    pool=mp.Pool(threads)
    pool.starmap(aggregate_reads, zip(idxchunks, repeat(ref), repeat(modfile), repeat(thresh), repeat(L)))
    pool.close()

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
                towrite=[i, str(pos), str(strands[i][pos]), str(fullmeth[i][pos]), str(fullunmeth[i][pos])]
                f.write('\t'.join(towrite)+'\n')

    if verbose:
        write_duration=round(time.time()-write_start)
        print('finished writing in %d seconds' % (write_duration))

if __name__ == "__main__":
    args=parseArgs()
    main(args.reffile, args.modfile, args.idxfile, args.outfile, args.threads, args.verbose)
