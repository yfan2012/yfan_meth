import argparse
import multiprocessing
import pysam
import math
import numpy as np
import itertools

def parseArgs():
    parser=argparse.ArgumentParser(description='get a meth barcode per read for megalodon')
    parser.add_argument('-i', '--input', type=str, required=True, help='megalodon output file with mod motif probs')
    parser.add_argument('-r', '--ref', type=str, required=True, help='reference genome used in megaldon')
    parser.add_argument('-m', '--motif', type=str, required=True, help='motif list file')
    parser.add_argument('-o', '--out', type=str, required=True, help='output file that lists each read and each barcode number')
    parser.add_argument('-t', '--threads', type=int, required=True, help='number of threads to use')
    parser.add_argument('-l', '--numreads', type=int, required=False, help='how many reads to consider')
    args=parser.parse_args()
    return args

class readMods:
    '''
    initial thing that packages all the modinfo
    '''
    def __init__(self, readname, positions, pmodratio, modtype):
        self.readname=readname
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
        apos=list(itertools.compress(positions, afil.tolist()))
        return apos
    def get_cpos(self, positions):
        cpos=list(itertools.compress(positions, cfil.tolist()))
        return cpos
    def get_apmod(self, pmodratio):
        apmod=list(itertools.compress(pmodratio, afil.tolist()))
        return apmod
    def get_cpmod(self, pmodratio):
        cpmod=list(itertools.compress(pmodratio, cfil.tolist()))
        return cpmod

    
class modCalls:
    '''
    package the mod calls
    '''
    def __init__(self, readname, motifs):
        self.readname=readname
        self.motifs=motifs

        
def fasta_dict(reffile):
    '''
    take reference filepath
    read into dict
    '''
    fa=pysam.FastaFile(reffile)
    tigs=fa.references
    fastadict={x:fa.fetch(x) for x in tigs}
    return fastadict


def find_thresh(pmodratio):
    '''
    determine what logpmod ratio to use as a threshold
    [cmod_thresh, amod_thresh]
    '''
    #for now just take from the rerio roc u made before
    #mean of the cmods and mean of the amods
    #cmods are nebdcm and neb15 (not sure why i never did roc for neb16? and neb14 specificity is uncertain)
    cmods=[1.0416759935039135, 0.9416759935039138]
    amods=[0.011380116536200191, -0.03861988346379963]
    athresh=sum(amods)/2
    cthresh=sum(cmods)/2
    thresh=[cthresh,athresh]
    return thresh

    
def grab_read(readname, byteoffset, bytelen):
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
    read=readMods(readname, positions, pmodratio, modtype):
    return read


def call_read(read, thresh, ref, k):    
    '''
    take in a readMod object
    return a list of calledmotifs
    '''
    ccallpos=np.array(read.cpmod) > thresh[0]
    acallpos=np.array(read.apmod) > thresh[1]
    poscalled=list(itertools.compress(read.cpos, ccallpos.tolist()))
    poscalled.extend(list(itertools.compress(read.apos, acallpos.tolist())))
    calledmotifs=[]
    for i in poscalled:
        calledmotifs.append(ref[i[0]][i[1]-k:i[1]+k])
    calls=modCalls(readname, calledmotifs)
    return calls



def per_read_calls(readinfo):
    '''
    for a given read get
    '''
    #load data from mod file

    #figure out which ones you consider to be meth

    #figure out which motifs your meth matches to
    

def read_megalodon_index(idxfile):
    '''
    read the mod indexfile into a dictionary
    readname:[byteoffset, bytelen]
    '''
    readidx={}
    with open(idxfile, 'r') as f:
        content=f.read().split('\n')
    for i in content[1:]:
        if len(i)>0:
            readinfo=i.split('\t')
            readidx[readinfo[0]]=[readinfo[1], readinfo[2]]
    return readidx




