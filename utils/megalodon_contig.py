import itertools
from bisect import bisect_left
import argparse
import multiprocessing as mp
import pysam
import re
import math
import numpy as np
from itertools import repeat
from megalodon_barcode import fasta_dict, expand_barcodes, read_megalodon_index

def parseArgs():
    parser=argparse.ArgumentParser(description='get a meth barcode per read for megalodon')
    parser.add_argument('-m', '--modfile', type=str, required=True, help='megalodon output file with mod motif probs')
    parser.add_argument('-i', '--idxfile', type=str, required=True, help='index of megalodon outputfile')
    parser.add_argument('-r', '--reffile', type=str, required=True, help='reference genome used in megaldon')
    parser.add_argument('-b', '--barcodefile', type=str, required=True, help='motif list file')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='output file that lists each ref seqname and its motif scores')
    parser.add_argument('-v', '--covfile', type=str, required=True, help='output file that lists each ref seqname and motif coverage')
    parser.add_argument('-a', '--abound', type=float, required=False, help='A threshold')
    parser.add_argument('-c', '--cbound', type=float, required=False, help='C threshold')
    parser.add_argument('-t', '--threads', type=int, required=True, help='number of threads to use')
    args=parser.parse_args()
    return args


iupacnt={
    'N': ['A', 'C', 'G', 'T'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'S': ['G', 'C'],
    'M': ['A', 'C'],
    'V': ['G', 'C', 'A'],
    'Y': ['C', 'T'], 
    'R': ['A', 'G']
}

'''
for testing
modfile='/mithril/Data/Nanopore/projects/methbin/zymo/megalodon/20190809_zymo_control/per_read_modified_base_calls.txt'
idxfile='/mithril/Data/Nanopore/projects/methbin/zymo/megalodon/20190809_zymo_control/per_read_modified_base_calls.txt.idx'
reffile='/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa'
barcodefile='/home/yfan/Code/yfan_nanopore/mdr/rebase/barcodes15.txt'
'''

def find_motifs(ref, barcodes):
    '''
    take in reference dictionary and barcode info
    return nested dict:
    {seqname1: {motif1:[start1, start2...], motif2:[start1, start2...]}, 
    seqname2: {motif1:[start1, start2...], motif2:[start1, start2...]},
    seqname3: ...}
    '''
    motifstarts={}
    for i in ref:
        motifstarts[i]={}
        for motif in barcodes:
            motifstarts[i][motif]=[]
            for j in barcodes[motif]:
                motifstarts[i][motif].extend([x.start() for x in re.finditer(j, ref[i])])
    return motifstarts


def check_position(starts, position):
    '''
    See if the position is a relevant motif
    returns list of motifs that are relevant
    considers 3 bases before and after the full motif
    '''
    meth=[]
    for i in starts:
        mlen=len(i)
        rank=bisect_left(starts[i], position)
        if rank==0:
            nearest=0
        elif rank==len(starts[i]):
            nearest=rank-1
        else:
            abovediff=abs(starts[i][rank]-position)
            belowdiff=abs(starts[i][rank-1]-position)
            if abovediff < belowdiff:
                nearest=rank
            else:
                nearest=rank-1
        pos=starts[i][nearest]
        if position >= pos-1 and position <= pos+mlen+1:
            meth.append(i)
    return meth
            
    
def add_read_list(read_index, modfile, athresh, cthresh, motifstarts, C, H, headers):
    '''
    take each read
    return motifcounts info 
    {motif1:[meth, unmeth], motif2:[meth, unmeth], ...}
    '''
    seqname=read_index[1]
    byteoffset=read_index[2]
    bytelen=read_index[3]
    with open(modfile, 'r') as f:
        f.seek(byteoffset,0)
        readcontent=f.read(bytelen).split('\n')
        f.close()
    headerlen=len(headers)
    meth=[seqname]+[0]*(headerlen-1) #query C or H is slow?
    cov=[seqname]+[0]*(headerlen-1) #query C or H is slow?
    for i in readcontent:
        if len(i)>0:
            readinfo=i.split('\t')
            ratio=math.log10(float(readinfo[5])/float(readinfo[4]))
            modtype=readinfo[6]
            starts=motifstarts[seqname]
            methinfo=check_position(starts, int(readinfo[3]))
            if len(methinfo)>0:
                for j in methinfo:
                    pos=headers.index(j)
                    cov[pos]+=1
                    if modtype=='Y' and ratio>athresh:
                        meth[pos]+=1
                    if modtype=='Z' and ratio>cthresh:
                        meth[pos]+=1
    C.append(cov)
    H.append(meth)
            


def main(reffile, modfile, idxfile, barcodefile, outfile, covfile, abound, cbound, threads):
    '''
    main function
    '''
    print('reading stuff')
    ref=fasta_dict(reffile)
    barcodes=expand_barcodes(barcodefile)
    motifstarts=find_motifs(ref, barcodes)
    
    readidx=read_megalodon_index(idxfile)

    headers=['chrname']
    headers.extend(barcodes.keys())
    
    ##set up parallel dict so the read processing func can write to it in parallel
    manager=mp.Manager()
    pool=mp.Pool(threads)
    
    C=manager.list()
    C.append(headers)
    H=manager.list()
    H.append(headers)
    pool.starmap(add_read_list, zip(readidx, repeat(modfile), repeat(abound), repeat(cbound), repeat(motifstarts), repeat(C), repeat(H), repeat(headers)))

    with open (outfile, 'w') as f:
        for i in H:
            f.write('\t'.join([str(x) for x in i])+'\n')
        f.close()
    with open(covfile, 'w') as f:
        for i in C:
            f.write('\t'.join([str(x) for x in i])+'\n')
        f.close()


        
if __name__ == "__main__":
    args=parseArgs()
    main(args.reffile, args.modfile, args.idxfile, args.barcodefile, args.outfile, args.covfile, args.abound, args.cbound, args.threads)
