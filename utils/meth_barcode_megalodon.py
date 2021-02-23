import argparse
import multiprocessing
import pysam
import numpy as np

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

def grab_read(readname, byteoffset, bytelen):
    '''
    get modinfo for the read
    a numpy array of [chr]
    '''
    
with open(modfile, 'r') as f:
    f.seek(byteoffset,0)
    readcontent=f.read(bytelen).split('\n')
    f.close()

class readMods:
    def __init__(self, readname, readcontent):
        self.readname=readname
        self.positions=get_positions(self)
    
    
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




