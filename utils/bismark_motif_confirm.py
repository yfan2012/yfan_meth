import itertools
import argparse
import multiprocessing as mp
import pysam
import math
import numpy as np
from itertools import repeat
import math
from megalodon_barcode_functions import *

'''
cxfile='/mithril/Data/Nanopore/projects/methbin/zymo/truth/bisulfite/bismark/ecoli/ecoli_1_bismark_bt2_pe.CX_report.txt'
reffile='/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa'
'''

def parseArgs():
    parser=argparse.ArgumentParser(description='report how many methylated genomic loci are not accounted for by the motif list')
    parser.add_argument('-c', '--cxfile', type=str, required=True, help='cx report file from bismark')
    parser.add_argument('-r', '--reffile', type=str, required=True, help='reference genome used in bismark analysis')
    parser.add_argument('-m', '--motifs', type=str, required=True, help='comma separated list of motifs')
    parser.add_argument('-p', '--positions', type=str, required=True, help='comma separated list of positions in the order the motifs are listed')
    args=parser.parse_args()
    return args


def read_cx(cxfile):
    '''
    takes in a cxfile
    returns same info in list of lists 
    '''
    with open(cxfile, 'r') as f:
        contents=f.read().split('\n')
    cx=[]
    for i in contents:
        if len(i)>0:
            info=i.split('\t')
            meth=int(info[3])
            unmeth=int(info[4])
            if meth+unmeth > 10 :
                if meth/(meth+unmeth) > .9 :
                    cx.append([info[0], int(info[1]), info[2], meth, unmeth])
    return cx


def main(cxfile, reffile, motifs, positions):
    ref=fasta_dict(reffile)
    cx=read_cx(cxfile)
    motifpos=dict(zip(motifs.split(','), [int(x) for x in positions.split(',')]))

    counted=0
    total=len(cx)
    for i in cx:
        chrom=i[0]
        pos=i[1]
        strand=i[2]
        for motif in motifpos:
            motifs=expand_motif([motif])
            if strand=='+' :
                start=pos-motifpos[motif]-1
                end=start+len(motif)
                context=ref[chrom][start:end]
            else:
                start=pos-len(motif)+motifpos[motif]
                end=start+len(motif)
                context=revcomp(ref[chrom][start:end])
            if context in motifs:
                counted+=1

    print(' '.join([str(counted), str(total)]))

    
if __name__ == "__main__":
    args=parseArgs()
    main(args.cxfile, args.reffile, args.motifs, args.positions)
