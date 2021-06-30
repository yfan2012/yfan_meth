import itertools
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
    parser.add_argument('-c', '--covfile', type=str, required=True, help='output file that lists each ref seqname and motif coverage')
    parser.add_argument('-a', '--abound', type=float, required=False, help='A threshold')
    parser.add_argument('-c', '--cbound', type=float, required=False, help='C threshold')
    parser.add_argument('-n', '--numreads', type=int, required=False, help='how many reads to consider')
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
        nearest=bisect_left(starts[i], position)
        pos=starts[i][nearest]
        if position >= pos-1 and position <= pos+mlen+1:
            meth.append(i)
    return meth
            
    
def add_read(read_index, modfile, athresh, cthresh, motifcounts):
    '''
    take each read
    return motifcounts info 
    {motif1:[meth, unmeth], motif2:[meth, unmeth], ...}
    '''
    byteoffset=read_index[2]
    bytelen=read_index[3]
    with open(modfile, 'r') as f:
        f.seek(byteoffset,0)
        readcontent=f.read(bytelen).split('\n')
        f.close()
    for i in readcontent:
        if len(i)>0:
            readinfo=i.split('\t')
            ratio=math.log10(float(readinfo[5])/float(readinfo[4]))
            modtype=readinfo[6]
            seqname=readinfo[1]
            starts=motifstarts[seqname]
            methinfo=check_position(starts, int(readinfo[3]))
            if len(methinfo)>0:
                for j in methinfo:
                    if modtype=='Y':
                        if ratio>athresh:
                            motifcounts[seqname][j][0]+=1
                        else:
                            motifcounts[seqname][j][1]+=1
                    if modtype=='Z':
                        if ratio>cthresh:
                            motifcounts[seqname][j][0]+=1
                        else:
                            motifcounts[seqname][j][1]+=1
                        
            


def main(reffile, modfile, idxfile, barcodefile, outfile, covfile, abound, cbound, threads):
    '''
    main function
    '''
    ref=fasta_dict(reffile)
    barcodes=expand_barcodes(barcodefile)
    motifstarts=find_motifs(ref, barcodes)
    
    readidx=read_megalodon_index(idxfile)

    ##set up parallel dict so the read processing func can write to it in parallel
    manager=mp.Manager()
    pool=mp.Pool(threads)
    motifcounts=manager.dict()
    for i in ref:
        motifcounts[i]=manager.dict()
        for j in barcodes:
            motifcounts[i][j]=manager.list([0,0])

    pool.starmap(add_read, zip(readidx, repeat(modfile), repeat(abound), repeat(cbound), repeat(motifcounts)))

    motiftotals=[]
    headers=['chrname']
    first=motifcounts.keys()[0]
    headers.extend(motifcounts[first].keys())

    allinfo=[]
    allsums=[]
    for i in motifcounts:
        chrinfo=[i]
        chrsum=[i]
        for j in headers[1:]:
            counts=motifcounts[i][j]
            tot=sum(counts)
            ratio=counts[0]/tot
            chrinfo.append(ratio)
            chrsum.append(tot)
        allinfo.append(chrinfo)
        allsums.append(chrsum)
        
    with open (outfile, 'w') as f:
        f.write('\t'.join(headers))
        for i in allinfo:
            f.write('\t'.join(i))
        f.close()
    with open(covfile, 'w') as f:
        f.write('\t'.join(headers))
        for i in allsums:
            f.write('\t'.join(i))
        f.close()
        

if __name__ == "__main__":
    args=parseArgs()
    main(args.reffile, args.modfile, args.idxfile, args.barcodefile, args.outfile, args.covfile, args.abound, args.cbound, args.threads):
