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
    pool.starmap(per_read_calls, zip(idxchunks, repeat(modfile), repeat(thresh), repeat(ref), repeat(barcodes), repeat(k), repeat(L)))
    with open (outfile, 'w') as f:
        for i in L:
            f.write(i)
    pool.close()
    pool.join()

if __name__ == "__main__":
    args=parseArgs()
    main(args.reffile, args.modfile, args.idxfile, args.barcodefile, args.outfile, args.abound, args.cbound, args.threads)
