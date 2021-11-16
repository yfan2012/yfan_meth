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
chrfile='/home/yfan/Code/yfan_nanopore/mdr/zymo/truth/chrlist_withlabels.txt'
motiffile='/mithril/Data/Nanopore/projects/methbin/zymo/truth/bisulfite/zymo_cmeth.csv'
cxfiledir='/mithril/Data/Nanopore/projects/methbin/zymo/truth/bisulfite/bismark'
'''

def parseArgs():
    parser=argparse.ArgumentParser(description='check if cmeth motif set is complete')
    parser.add_argument('-c', '--cxfiledir', type=str, required=True, help='directory containing cxfiles')
    parser.add_argument('-r', '--reffile', type=str, required=True, help='reference genome used in bismark analysis')
    parser.add_argument('-s', '--chrfile', type=str, required=True, help='chromosome key file')
    parser.add_argument('-m', '--motiffile', type=str, required=True, help='motif list for each chromosome')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='motif list for each chromosome')
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
            cx.append(i.split('\t'))
    return cx


def read_chroms(chrfile):
    '''
    takes in a chr key file
    returns label:[chr1, plas1, plas2, ...]
    '''
    with open(chrfile, 'r') as f:
        contents=f.read().split('\n')
    chroms={}
    for i in contents:
        if len(i)>0:
            chromsinfo=i.split(' ')
            if chromsinfo[1] in chroms:
                chroms[chromsinfo[1]].append(chromsinfo[0])
            else:
                chroms[chromsinfo[1]]=[chromsinfo[0]]
    return chroms


def read_motif(motiffile):
    '''
    takes in a motif info file
    returns label:[motif, pos]
    '''
    with open(motiffile, 'r') as f:
        contents=f.read().split('\n')
    motifs={}
    for i in contents:
        if len(i)>0:
            motifinfo=i.split(',')
            if motifinfo[0] in motifs:
                motifs[motifinfo[0]].append([motifinfo[1], motifinfo[2]])
            else:
                motifs[motifinfo[0]]=[[motifinfo[1], motifinfo[2]]]
    return motifs


def main(cxfiledir, reffile, chrfile, motiffile, outfile):
    ref=fasta_dict(reffile)
    chroms=read_chroms(chrfile)
    motifs=read_motif(motiffile)
    newchroms={}
    for species in chroms:
        if species in motifs:
            newchroms[species]=chroms[species]
    chroms=newchroms

    meth=[]
    
    for species in chroms:
        cxfile=cxfiledir+'/'+species+'/'+species+'_1_bismark_bt2_pe.CX_report.txt'
        cx=read_cx(cxfile)
        seqs=chroms[species]
        numcalled=0
        numcounted=0
        for i in cx:
            seq=i[0]
            pos=int(i[1])
            strand=i[2]
            called=float(i[3])
            uncalled=float(i[4])
            if called+uncalled>15:
                if seq in seqs and called/(called+uncalled)>.5:
                    numcalled+=1
                    for j in motifs[species]:
                        start=pos-int(j[1])+1
                        end=start+len(j[0])
                        motiflist=expand_motif([j[0]])
                        if ref[seq][start:end] in motiflist:
                            numcounted+=1

        meth.append([species, str(numcalled), str(numcounted)])

    with open(outfile, 'w') as f:
        for i in meth:
            f.write(','.join(i))

            
if __name__ == "__main__":
    args=parseArgs()
    main(args.cxfiledir, args.reffile, args.chrfile, args.motiffile, args.outfile)
