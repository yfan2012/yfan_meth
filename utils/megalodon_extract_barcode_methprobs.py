import itertools
import argparse
import multiprocessing as mp
import pysam
import math
import sys
import numpy as np
from itertools import repeat
from megalodon_barcode_functions import *
from megalodon_agg import get_reads
import time
import re


def parseArgs():
    parser=argparse.ArgumentParser(description='pull out relevant positions from cx style report based on barcode')
    parser.add_argument('-r', '--reffile', type=str, required=True, help='reference genome used in megalodon')
    parser.add_argument('-b', '--barcodefile', type=str, required=True, help='list of motifs in the barcode')
    parser.add_argument('-m', '--modfile', type=str, required=True, help='per read mod calls form megalodon')
    parser.add_argument('-i', '--idxfile', type=str, required=True, help='index of per read mod calls from megalodon')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='output file')
    parser.add_argument('-t', '--threads', type=int, required=True, help='number of threads to use')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose outputs')
    args=parser.parse_args()
    return args


def find_motif_positions(barcodes, ref):
    '''
    take motifs in barocdes and ref
    find all motif related positions
    {motif:{chr:[positions in chr]}}
    '''
    motif_positions={}
    for chrom in ref:
        motif_positions[chrom]={}
        for motif in barcodes:
            positions=[]
            for i in barcodes[motif]:
                mlen=len(i)
                ranges=[list(range(m.start(), m.start()+mlen)) for m in re.finditer(i, ref[chrom])]
                for j in ranges:
                    positions.extend(j)
            motif_positions[chrom][motif]=set(positions)
    return(motif_positions)


def check_position(chrom, pos, motif_positions):
    '''
    given an chrom and a pos, check if it's associated with any of the meth motifs
    return the motif(s) associated
    returns blank list if there are no motifs assocated with the position
    '''
    motifhits=[]
    for i in motif_positions[chrom]:
        if pos in motif_positions[chrom][i]:
            motifhits.append(i)
    return(motifhits)

    
def label_positions(p, motif_positions, modfile, ref, q):
    '''
    read the appropriate piece of the modfile to get the read info
    for each meth call in the read info, check if it's part of a relevant meth motif
    save out the positions with relevant motifs
    '''
    while True:
        idxchunk=p.get()
        if idxchunk!='done':
            readchunk=get_reads(idxchunk, modfile)
            meth_overlaps=[]
            for i in readchunk[:-1]:
                chrom=i[1]
                if chrom in list(motif_positions):
                    pos=int(i[3])
                    motif=check_position(chrom, pos, motif_positions)
                    base=ref[chrom][pos]
                    if len(motif)>0:
                        for j in motif:
                            meth_overlaps.append([chrom, pos, i[2], float(i[4]), j, base])
            q.put(meth_overlaps)
        else:
            q.put('done')
            break


def get_idxchunks(readidx, chunksize, p, threads):
    '''
    split readidx into chunks for processing
    put in to p (an input queue)
    '''
    numchunks=math.ceil(len(readidx)/3000)
    for i in range(numchunks):
        start=i*chunksize
        end=(i+1)*chunksize
        if len(readidx)<end:
            if len(readidx[start:]) > 0:
                p.put(readidx[start:])
        else:
            if len(readidx[start:end]) > 0:
                p.put(readidx[start:end])
    for i in range(threads):
        p.put('done')

        
def writer(q, outfile, threads):
    '''
    writes stuff from the out q to the file
    '''
    donecount=0
    with open(outfile, 'w') as f:
        while donecount<threads:
            m=q.get()
            if m!='done':
                for i in m:
                    f.write(','.join([str(x) for x in i])+'\n')
                    f.flush()
            else:
                donecount+=1
                
##idxfile='/mithril/Data/Nanopore/projects/methbin/zymo/megalodon/20190809_zymo_control/per_read_modified_base_calls.txt.small.idx'
##barcodefile='/home/yfan/Code/yfan_nanopore/mdr/zymo/barcodes_zymo_curated.txt'
##modfile='/mithril/Data/Nanopore/projects/methbin/zymo/megalodon/20190809_zymo_control/per_read_modified_base_calls.txt'
##reffile='/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa'

def main(reffile, barcodefile, modfile, idxfile, outfile, threads, verbose):
    ref=fasta_dict(reffile)
    ##for our purposes, don't bother with the yeast genome
    ##this could become a prob later for bigger metagenomes though
    
    barcodes=expand_barcodes(barcodefile)
    motif_positions=find_motif_positions(barcodes, ref)
    
    readidx=read_megalodon_index(idxfile)
    chunksize=3000

    manager=mp.Manager()
    p=manager.Queue()
    q=manager.Queue()

    ps=[]
    ps.append(mp.Process(target=get_idxchunks, args=(readidx, chunksize, p, threads)))
    for i in range(threads):
        ps.append(mp.Process(target=label_positions, args=(p, motif_positions, modfile, ref, q)))
    ps.append(mp.Process(target=writer, args=(q, outfile, threads)))

    for i in ps:
        i.start()
    for i in ps:
        i.join()
        


if __name__ == "__main__":
    args=parseArgs()
    main(args.reffile, args.barcodefile, args.modfile, args.idxfile, args.outfile, args.threads, args.verbose)
