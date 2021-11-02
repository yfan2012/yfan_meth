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
    parser.add_argument('-p', '--percentcall', type=float, required=True, help='number of threshold to consider a position methylated')
    parser.add_argument('-l', '--lencontext', type=int, required=True, help='length of the context considered')
    parser.add_argument('-m', '--mincov', type=int, required=True, help='minimum coverage required to be considered')
    parser.add_argument('-s', '--seqname', type=str, required=False, help='specify chromosome if u want')
    args=parser.parse_args()
    return args


def read_meth_file(cxfile, percentcall, mincov, seqname):
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
            chrname=readinfo[0]
            if chrname==seqname:
                meth=float(readinfo[3])
                unmeth=float(readinfo[4])
                if meth+unmeth>=mincov:
                    fracmeth=meth/(unmeth+meth)
                    if fracmeth>percentcall:
                        methcalled.append(i.split('\t'))
    return methcalled


def get_contexts(methcalled, ref, lencontext, seqname):
    '''
    take positions called as methylated
    return list of contexts
    '''
    methseqs=[]
    for i in methcalled:
        pos=int(i[1])
        strand=i[2]
        seq=ref[seqname][pos-lencontext:pos+lencontext]
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


def count_motif_occurences_genome(ref, seqname, barcodes):
    '''
    take ref
    return counts of each barcode occurence
    '''
    refmotifs={}
    for i in barcodes:
        refmotifs[i]=0
        for motif in barcodes[i]:
            refmotifs[i]+=ref[seqname].count(motif)
    return refmotifs


def calc_normalized_occurences(motifcounts, refmotifs):
    '''
    take motifcounts and num of motifs in the reference (refmotifs)
    return motifcounts normalized by refmotifs
    '''
    normcounts={}
    for i in motifcounts:
        if refmotifs[i]==0:
            normcounts[i]=None
        else:
            normcounts[i]=motifcounts[i]/refmotifs[i]
    return normcounts

##for testing:
##reffile='/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa'
##cxfile='/mithril/Data/Nanopore/projects/methbin/zymo/truth/bisulfite/bismark/ecoli/ecoli_1_bismark_bt2_pe.CX_report.txt'
##barcodefile='/home/yfan/Code/yfan_nanopore/mdr/rebase/barcodes20.txt'
##seqname='Escherichia_coli_chromosome'

def main(reffile, cxfile, barcodefile, percentcall, lencontext, mincov, seqname):
    ref=fasta_dict(reffile)
    barcodes=expand_barcodes(barcodefile)
    
    methcalled=read_meth_file(cxfile, percentcall, mincov, seqname)
    methseqs=get_contexts(methcalled, ref, lencontext, seqname)
    motifcounts=count_motif_occurences_meth(methseqs, barcodes)
    refmotifs=count_motif_occurences_genome(ref, seqname, barcodes)

    normcounts=calc_normalized_occurences(motifcounts, refmotifs)
    
    normcountlist=[seqname]
    for i in normcounts:
        normcountlist.append(str(normcounts[i]))
    print(','.join(normcountlist))
    


if __name__ == "__main__":
    args=parseArgs()
    main(args.reffile, args.cxfile, args.barcodefile, args.percentcall, args.lencontext, args.mincov, args.seqname)
