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
    parser.add_argument('-c', '--cxfile', type=str, required=True, help='bisulfite cytosine report file (*CX_report.txt)')
    parser.add_argument('-r', '--reffile', type=str, required=True, help='reference genome used in bismark analysis')
    parser.add_argument('-b', '--barcodefile', type=str, required=True, help='motif list file')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='output file that lists each read and each barcode number')
    parser.add_argument('-p', '--percentcall', type=int, required=True, help='number of threshold to consider a position methylated')
    parser.add_argument('-l', '--lencontext', type=int, required=True, help='number of threshold to consider a position methylated')
    parser.add_argument('-t', '--threads', type=int, required=True, help='number of threads to use')
    args=parser.parse_args()
    return args


def read_meth_file(cxfile, percentcall):
    '''
    read in cytosine report file
    return positions with meth count greater than certain threshold
    '''
    methcalled=[]
    with open(cxfile, 'r') as f:
        content=f.read().split('\n')
    for i in content:
        if len(i)>0:
            readinfo=i.split('\t')
            meth=float(readinfo[3])
            unmeth=float(readinfo[4])
            fracmeth=meth/(unmeth+meth)
            if fracmeth>percentcall:
                methcalled.append(i)
    return methcalled



def get_contexts(methcalled, ref, lencontext, chroms):
    '''
    take positions called as methylated
    return list of contexts
    '''
    methseqs=[]
    for i in methcalled:
        chrom=i[0]
        if chrom in chroms:
            pos=int(i[1])
            strand=i[2]
            seq=ref[chrom][pos-lencontext:pos+lencontext]
            if strand=='-':
                seq=revcomp(seq)
            methseqs.append(seq)
    return methseqs


def count_motif_occurences_meth(methseqs, barcodes):
    '''
    take list of meth seq contexts
    return barcode counts
    '''
    motifcounts={}
    for i in barcodes:
        motifcounts[i]=0
        for motif in barcodes[i]:
            for seq in methseqs:
                if motif in seq:
                    motifcounts[i]+=1
    return motifcounts


def count_motif_occurences_genome(ref, chroms, barcodes):
    '''
    take ref
    return counts of each barcode occurence
    '''
    refmotifs={}
    for i in barcodes:
        refmotifs[i]=0
        for motif in bracodes[i]:
            for chrom in chroms:
                refmotifs[i]+=ref[chrom].count(motif)
    return refmotifs
        

def main(reffile, cxfile, barcodefile, outfile, threads, percentcall, lencontext, choosechrom):
    ref=fasta_dict(reffile)
    barcodes=expand_barcodes(barcodefile)
    
    methcalled=read_meth_file(cxfile, percentcall)
    if choosechrom: 
        chroms=find_chroms(methcalled)
    else:
        chroms=[ref.keys()]
    methseqs=get_contexts(methcalled, ref, lencontext)
    motifcounts=count_motif_occurences_meth(methseqs, chroms, barcodes)
    refmotifs=count_motif_occurencees_genome(ref, chroms, barcodes)
        
    with open (outfile, 'w') as f:
        for i in L:
            f.write(i)


if __name__ == "__main__":
    args=parseArgs()
    main(args.reffile, args.modfile, args.idxfile, args.barcodefile, args.outfile, args.abound, args.cbound, args.threads)
