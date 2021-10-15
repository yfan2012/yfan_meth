import itertools
import argparse
import multiprocessing as mp
import pysam
import math
import numpy as np
from itertools import repeat

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
for testing:
reffile='/uru/Data/Nanopore/projects/read_class/zymo/ref/zymo_all.fa'
idxfile='/mithril/Data/Nanopore/projects/methbin/zymo/megalodon/20190809_zymo_control/per_read_modified_base_calls.txt.small.idx'
modfile='/mithril/Data/Nanopore/projects/methbin/zymo/megalodon/20190809_zymo_control/per_read_modified_base_calls.txt'
barcodefile='/home/yfan/Code/yfan_nanopore/mdr/rebase/barcodes20.txt'
'''


class readMods:
    '''
    initial thing that packages all the modinfo
    '''
    def __init__(self, readname, chrname, strand, positions, pmodratio, modtype, minpos, maxpos):
        self.readname=readname
        self.chrname=chrname
        self.strand=strand
        self.minpos=minpos
        self.maxpos=maxpos
        self.afil=self.get_afil(modtype)
        self.cfil=self.get_cfil(modtype)
        self.apos=self.get_apos(positions)
        self.cpos=self.get_cpos(positions)
        self.apmod=self.get_apmod(pmodratio)
        self.cpmod=self.get_cpmod(pmodratio)
    def get_afil(self, modtype):
        afil=np.array(modtype)=='Y'
        return afil
    def get_cfil(self, modtype):
        cfil=np.array(modtype)=='Z'
        return cfil
    def get_apos(self, positions):
        apos=list(itertools.compress(positions, self.afil.tolist()))
        return apos
    def get_cpos(self, positions):
        cpos=list(itertools.compress(positions, self.cfil.tolist()))
        return cpos
    def get_apmod(self, pmodratio):
        apmod=list(itertools.compress(pmodratio, self.afil.tolist()))
        return apmod
    def get_cpmod(self, pmodratio):
        cpmod=list(itertools.compress(pmodratio, self.cfil.tolist()))
        return cpmod

    
class modCalls:
    '''
    package the mod calls
    '''
    def __init__(self, readname, chrname, strand, calledmotifs, barcodes, motifcounts):
        self.readname=readname
        self.chrname=chrname
        self.strand=strand
        self.calledmotifs=calledmotifs
        self.motifcounts=motifcounts
        self.bc_counts=self.assign_barcode(barcodes, strand)
        self.bc_norm=self.norm_barcode(barcodes)
    def assign_barcode(self, barcodes, strand):
        '''
        take called motifs
        count how many have each motif
        '''
        bc_counts={}
        for i in barcodes:
            bc_counts[i]=0
        for i in self.calledmotifs:
            for j in barcodes:
                for k in barcodes[j]:
                    if strand=='-':
                        seq=revcomp(k)
                    else:
                        seq=k
                    if seq in i:
                        bc_counts[j]+=1
                        break
        return bc_counts
    def norm_barcode(self, barcodes):
        '''
        normalize motif counts in called motifs
        returns counts as a proportion of the number you'd expect randomly 
        '''
        bc_norm={}
        if len(self.calledmotifs) > 0: 
            numcalled=len(self.calledmotifs)
            lencalled=len(self.calledmotifs[0])
            for i in self.bc_counts:
                nummotifs=len(barcodes[i])
                normconst=self.bc_counts[i]/ ( (lencalled-len(i)+1) * (4**(lencalled-len(i))) * numcalled * nummotifs /float(4**lencalled))
                bc_norm[i]=normconst
        else:
            for i in self.bc_counts:
                bc_norm[i]=0
        for i in self.bc_counts:
            if self.motifcounts[i]==0:
                bc_norm[i]=None
        return bc_norm
    

def fasta_dict(reffile):
    '''
    take reference filepath
    read into dict
    '''
    fa=pysam.FastaFile(reffile)
    tigs=fa.references
    fastadict={x:fa.fetch(x) for x in tigs}
    return fastadict


def revcomp(seq):
    '''
    reverse complement a sequence
    '''
    key={'A': 'T', 'T':'A', 'G':'C',  'C':'G'}
    newseq=''
    for i in seq:
        newseq+=key[i]
    return newseq[::-1]

    
def find_thresh():
    '''
    determine what logpmod ratio to use as a threshold
    [cmod_thresh, amod_thresh]
    for now just take from the rerio roc u made before
    mean of the cmods and mean of the amods
    cmods are nebdcm and neb15 (not sure why i never did roc for neb16? and neb14 specificity is uncertain)
    '''
    cmods=[1.0416759935039135, 0.9416759935039138]
    amods=[0.011380116536200191, -0.03861988346379963]
    cthresh=sum(cmods)/2
    athresh=sum(amods)/2
    thresh=[cthresh,athresh]
    return thresh


##for testing: idxfile='/mithril/Data/Nanopore/projects/methbin/zymo/megalodon/20190809_zymo_control/per_read_modified_base_calls.txt.idx'
def read_megalodon_index(idxfile):
    '''
    read the mod indexfile into a dictionary
    [readname, chrname, byteoffset, bytelen]
    '''
    readidx=[]
    with open(idxfile, 'r') as f:
        content=f.read().split('\n')
    for i in content[1:]:
        if len(i)>0:
            readinfo=i.split('\t')
            readidx.append([readinfo[0], readinfo[1], int(readinfo[2]), int(readinfo[3]), readinfo[4]])
    return readidx


def get_num_motifs(chrname, minpos, maxpos, ref, barcodes, strand):
    '''
    count how many of each motif a read covers
    '''
    motifcounts={}
    for bc in barcodes:
        counts=0
        for motif in barcodes[bc]:
            if strand=='-':
                revmotif=revcomp(motif)
                counts+=ref[chrname][minpos:maxpos].count(revmotif)
            else:
                counts+=ref[chrname][minpos:maxpos].count(motif)
        motifcounts[bc]=counts
    return motifcounts
            

def grab_read(modfile, readname, chrname, byteoffset, bytelen, strand):
    '''
    get modinfo for the read
    '''
    with open(modfile, 'r') as f:
        f.seek(byteoffset,0)
        readcontent=f.read(bytelen).split('\n')
        f.close()
    positions=[]
    pmodratio=[]
    modtype=[]
    for i in readcontent:
        if len(i)>0:
            readinfo=i.split('\t')
            positions.append([readinfo[1], int(readinfo[3])])
            pmodratio.append(math.log10(float(readinfo[5])/float(readinfo[4])))
            modtype.append(readinfo[6])
    minpos=min(positions[0][1], positions[-1][1])
    maxpos=max(positions[0][1], positions[-1][1])
    read=readMods(readname, chrname, strand, positions, pmodratio, modtype, minpos, maxpos)
    return read


def expand_motif(motifs):
    '''
    take in motif with possible ambiguous nts
    return list of non-ambigusous nts
    '''
    ambig=0
    for i in iupacnt:
        if i in ''.join(motifs):
            ambig+=1
    if ambig==0:
        return motifs
    else:
        newmotifs=[]
        for i in motifs:
            for j in iupacnt:
                if j in i:
                    newmotifs.extend([i.replace(j,k,1) for k in iupacnt[j]])
        return expand_motif(newmotifs)


##for testing: barcodefile='/home/yfan/Code/yfan_nanopore/mdr/rebase/barcodes20.txt'
def expand_barcodes(barcodefile):
    '''
    take barcodefile
    return dict original:expanded
    '''
    with open(barcodefile, 'r') as f:
        bcinfo=f.read().split('\n')
    barcodes={}
    for i in bcinfo:
        if len(i)>0:
            barcodes[i]=expand_motif([i])
    return barcodes


def call_read(read, thresh, ref, barcodes, k):    
    '''
    take in a readMod object
    return a modCalls object
    '''
    ccallpos=np.array(read.cpmod) > thresh[0]
    acallpos=np.array(read.apmod) > thresh[1]
    poscalled=list(itertools.compress(read.cpos, ccallpos.tolist()))
    poscalled.extend(list(itertools.compress(read.apos, acallpos.tolist())))
    calledmotifs=[]
    for i in poscalled:
        calledmotifs.append(ref[i[0]][i[1]-k:i[1]+k])
    motifcounts=get_num_motifs(read.chrname, read.minpos, read.maxpos, ref, barcodes, read.strand)
    calls=modCalls(read.readname, read.chrname, read.strand, calledmotifs, barcodes, motifcounts)
    return calls


def per_read_calls(idxchunk, modfile, thresh, ref, barcodes, k, L):
    '''
    for a given read
    run through the pipe
    this is the function to parallel later
    '''
    holdreads=[]
    for readinfo in idxchunk:
        read=grab_read(modfile, readinfo[0], readinfo[1], readinfo[2], readinfo[3], readinfo[4])
        calls=call_read(read, thresh, ref, barcodes, k)
        callinfo=[calls.readname, calls.chrname]
        for i in barcodes:
            callinfo.append(str(calls.bc_norm[i]))
        holdreads.append('\t'.join(callinfo)+'\n')
    L+=holdreads    

    
