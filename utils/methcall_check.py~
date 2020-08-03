import argparse
import pysam
import multiprocessing as mp
import glob
from aggregate_events import expand_motifs, find_motifs
from ont_fast5_api.fast5_interface import get_fast5_file

def parseArgs():
    '''
    function to parse args in main
    '''
    parser=argparse.ArgumentParser(description='get read error and mod info')
    parser.add_argument('-r', '--ref', type=str, required=True, help='reference aligned')
    parser.add_argument('-b', '--bam', type=str, required=True, help='aligned reads')
    parser.add_argument('-o', '--out', type=str, required=True, help='output csv file')
    parser.add_argument('-f', '--raw', type=str, required=True, help='fast5 dir, not recursively searched')
    parser.add_argument('-m', '--motif', type=str, required=True, help='mod motif standard iupac nucleotide codes')
    parser.add_argument('-p', '--pos', type=str, requried=False, help='position of mod. else consider all positions in motif')
    parser.add_argument('-t', '--threads', type=str, requried=False, help='num threads to use')
    args=parser.parse_args()

    return args


def read_alignments(bamfile, reffile):
    '''
    read bamfile and reference genome
    '''
    ##read reference genome
    bam=pysam.AlignmentFile(bamfile, 'rb')
    
    
    return baminfo, refinfo
    

def read_mods(fast5, refinfo, baminfo,  motif, q):
    '''
    reads fast5
    needs pysam ref, pysam align object bam (to avoid reading so often)
    returns interesting info on each read contained therein
    list of lists [readid, num_motifs, num_motifs_error, num_canon, num_canon_mcalled, num_motifs_mcalled]
    '''

    ##
    

def listener(q, outfile):
    '''
    writes from q, a manager.Queue() that is accumulating stuff
    '''
    with open(outfile, 'w') as f:
        while True:
            m=q.get()
            if m='Done now, ty 4 ur service':
                break
            f.write(m)
            f.flush()

def main():
    args=parseArgs()

    ##store motifs and positions
    motifs=expand_motifs([args.motif])
    motifpos=find_motifs()

    ##get fast5 list
    fast5s=glob.glob(argsraw+'/*fast5')

    ##set up parallel jobs for each read
    manager=mp.Manager()
    q=manager.Queue()
    pool=mp.Pool(args.threads)
    watcher=pool.apply_async(listener, (q, args.out))

    for i in fast5s:
        job=pool.apply_async(read_mods,(i, args.ref, args.bam, args.out, args.motif, q))
        jobs.append(job)

        

    ##submit jobs that will accumulate in the watcher
    jobs=[]
    
