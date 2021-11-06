from functools import partial
import sys
import itertools
import argparse
import multiprocessing as mp
import pysam
import math
import numpy as np
from itertools import repeat
from megalodon_barcode_functions import *
import time
from joblib import Parallel, delayed 

def parseArgs():
    parser=argparse.ArgumentParser(description='get a meth barcode per read for megalodon')
    parser.add_argument('-m', '--modfile', type=str, required=True, help='megalodon output file with mod motif probs')
    parser.add_argument('-i', '--idxfile', type=str, required=True, help='index of megalodon outputfile')
    parser.add_argument('-r', '--reffile', type=str, required=True, help='reference genome used in megaldon')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='output file that lists each read and each barcode number')
    parser.add_argument('-t', '--threads', type=int, required=True, help='number of threads to use')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose outputs')
    args=parser.parse_args()
    return args

##try nanocompore parallel method

class methylation_summary ():
    from megalodon_barcode_functions import fasta_dict, revcomp, find_thresh
    import time
    def __init__(self,
                 modfile:str,
                 idxfile:str,
                 reffile:str,
                 outfile:str,
                 threads:int=12,
                 verbose:bool=False):
        self.__modfile=modfile
        self.__idxfile=idxfile
        self.__reffile=reffile
        self.__outfile=outfile
        self.__threads=threads
        self.__verbose=verbose
        self.__ref=fasta_dict(reffile)

        
    def __call__(self):
        in_q=mp.Queue() ##full of read batches
        out_q=mp.Queue() ##full of (aggmeth, aggunmeth)
        err_q=mp.Queue()
        
        
        processes=[]
        processes.append(mp.Process(target=self.__get_reads, args=(in_q, err_q)))
        for i in range(self.__threads):
            processes.append(mp.Process(target=self.__aggregate_reads, args=(in_q, out_q)))
        processes.append(mp.Process(target=self.__write_reads, args=(out_q, err_q)))

        for p in processes:
            p.start()

        for p in processes:
            p.join()
        
        in_q.close()
        out_q.close()

        
    def __get_reads(self, in_q, err_q):
        '''
        take chunks of read info from readidx
        return batch of reads
        '''
        threads=self.__threads
        ref=self.__ref
        readidx=read_megalodon_index(self.__idxfile)

        ##separate reads into batches
        idxchunks=[]
        chunksize=math.ceil(len(readidx)/threads/100)
        for i in range(0,threads*100):
            start=i*chunksize
            end=(i+1)*chunksize
            if len(readidx)<end:
                idxchunks.append(readidx[start:])
            else:
                idxchunks.append(readidx[start:end])

        ##read in the read batch
        for readinfo in idxchunks:
            if len(readinfo)>0:
                startbyte=readinfo[0][2]
                endbyte=readinfo[-1][2]+readinfo[-1][3]
                bytelen=endbyte-startbyte
                readbatch=[]
                with open(self.__modfile, 'r') as f:
                    f.seek(startbyte,0)
                    readcontent=f.read(bytelen).split('\n')
                    f.close()
                for i in readcontent:
                    readbatch.append(i.split('\t'))
                in_q.put(readbatch)
        for i in range(self.__threads):
            in_q.put(None)

            
    def __aggregate_reads(self, in_q, out_q):
        '''
        take a batch of reads (in idxchunk)
        return aggregated read methylation calls
        '''
        for readbatch in iter(in_q.get, None):
            '''
            print('starting batch of length %d ' % len(readbatch))
            sys.stdout.flush()
            tstart=time.time()
            '''
            aggmeth=self.__get_empty_ref()
            aggunmeth=self.__get_empty_ref()
            thresh=find_thresh()
            
            for i in readbatch[0:-1]:
                ratio=math.log10(float(i[5])/float(i[4]))
                pos=float(i[3])
                base=self.__ref[i[1]][int(i[3])]
                if base=='C' or base=='G':
                    if ratio > thresh[0]:
                        aggmeth[i[1]][int(i[3])]+=1
                    else:
                        aggunmeth[i[1]][int(i[3])]+=1
                if base=='A' or base=='T':
                    if ratio > thresh[1]:
                        aggmeth[i[1]][int(i[3])]+=1
                    else:
                        aggunmeth[i[1]][int(i[3])]+=1
            out_q.put((aggmeth, aggunmeth))
        out_q.put(None)
        
    def __get_empty_ref(self):
        '''
        take ref dict
        return ref dict with np array full of zeros
        '''
        aggcalls={}
        for i in self.__ref:
            aggcalls[i]=np.zeros(len(self.__ref[i]))
        return aggcalls

    

    def __write_reads(self, out_q, err_q):
        '''
        write from the out queue
        '''
        
        methagg=[]
        unmethagg=[]
        for _ in range(self.__threads):
            for am, au in iter(out_q.get, None):
                methagg.append(am)
                unmethagg.append(au)
                
        with open (self.__outfile, 'w') as f:
                fullmeth=self.__sum_refs(methagg, self.__ref)
                fullunmeth=self.__sum_refs(unmethagg, self.__ref)
                for i in fullmeth:
                    for pos in range(0, len(fullmeth[i])):
                        towrite=[i, str(pos), str(fullmeth[i][pos]), str(fullunmeth[i][pos])]
                        f.write('\t'.join(towrite)+'\n')


    def __sum_refs(aglist):
        '''
        take list of aggregated meth calls
        combine them into one final thing
        '''
        fullmeth=get_empty_ref()
        for i in aglist:
            for j in i:
                fullmeth[j]=np.add(i[j], fullmeth[j])
        return fullmeth
                 


def main(reffile, modfile, idxfile, outfile, threads, verbose):
    m=methylation_summary(modfile=modfile,
                        reffile=reffile,
                        idxfile=idxfile,
                        outfile=outfile,
                        threads=threads,
                        verbose=verbose)
    m()

if __name__ == "__main__":
    args=parseArgs()
    main(args.reffile, args.modfile, args.idxfile, args.outfile, args.threads, args.verbose)
