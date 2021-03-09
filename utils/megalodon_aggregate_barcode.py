import argparse
import multiprocessing as mp
import pysam
import math
import numpy as np
import itertools
from megalodon_barcode import *




def agg_per_read(readinfo, modfile, ref, barcodes, q):
    '''
    for a given read
    get the barcode scores
    '''
    read=grab_read(modfile, readinfo[0], readinfo[1], readinfo[2])


def main(reffile, modfile, idxfile, barcodefile, outfile, threads):
    '''
    a lot of this is the same in megalodon_barcode
    '''
    ref=fast_dict(reffile)
    readidx=read_megalodon_index(idxfile)
    barcodes=expandbarcodes(barcodefile)

    manager=mp.Manager()
    q=manager.Queue()
    pool=mp.Pool(threads)

    watcher=pool.apply_async(listener, (q, outfile))

    jobs=[]
    for i in readidx:
        job=pool.apply_async(agg_per_read, (i, modfile, ref, barcodes, q))
        jobs.append(job)

    for job in jobs:
        job.get()

    q.put('Done now, ty 4 ur service')
    pool.close()
    pool.join()


    
if __name__=='__main__':
    args=parseArgs()
    main(args.reffile, args.modfile, args.idxfile, args.barcodefile, args.outfile, args.threads)


