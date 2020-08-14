import argparse
import pysam
import numpy as np
import multiprocessing as mp
import glob
import re
import time
from aggregate_events import expand_motifs, find_motifs
from ont_fast5_api.fast5_interface import get_fast5_file

def parseArgs():
    '''
    function to parse args in main
    '''
    parser=argparse.ArgumentParser(description='get read error and mod info')
    parser.add_argument('-r', '--ref', type=str, required=True,
                        help='reference aligned to')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='aligned reads')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output csv file')
    parser.add_argument('-f', '--raw', type=str, required=True,
                        help='fast5 dir, not recursively searched')
    parser.add_argument('-m', '--motif', type=str, required=True,
                        help='mod motif standard iupac nucleotide codes')
    parser.add_argument('-p', '--pos', type=int, required=False,
                        help='position of mod')
    parser.add_argument('-v', '--verbose', action='store_true', required=False, default=False,
                        help='list reads/alignments without md tags')
    parser.add_argument('-t', '--threads', type=int, required=False,
                        help='num threads to use')
    args=parser.parse_args()
    return args


def fasta_dict(reffile):
    '''
    from fastafile, get a dict of seqname:seq
    '''
    fa=pysam.FastaFile(reffile)
    tigs=fa.references
    fastadict={ x:fa.fetch(x) for x in tigs }
    return fastadict

def rev_comp(seq):
    key={'A': 'T', 'T':'A', 'G':'C',  'C':'G'}
    newseq=''
    for i in seq:
        newseq+=key[i]
    return newseq[::-1]

def motif_position(fastadict, motifs, modpos):
    '''
    need fastadict, list of motifs, position of modified base, and how many bases surrounding
    get dict of seqname:sorted np.array of motif positions of the putative mod
    '''
    motifpos={}
    for seq in fastadict:
        motifpos[seq]=[]
        for motif in motifs:
            for i in re.finditer(motif, fastadict[seq]):
                    motifpos[seq].extend(range(i.start()+modpos, i.start()+modpos+1))
    for i in motifpos:
        motifpos[i]=np.sort(np.array(list(set(motifpos[i]))))
    return motifpos


def read_mods(fast5, bamfile, motifpos, motifs, modpos, q):
    '''
    reads fast5file, indexed pysam align object bam, motifpos dict
    returns interesting info on each read contained therein
    '''
    start_time=time.time()

    ##store bamfiles
    bam=pysam.AlignmentFile(bamfile, 'rb')
    baminfo=pysam.IndexedReads(bam)
    baminfo.build()

    motiflen=len(motifs[0])

    with get_fast5_file(fast5, mode="r") as f5:
        for read_id in f5.get_read_ids():
            readff=f5.get_read(read_id)
            latest_basecall=readff.get_latest_analysis('Basecall_1D')
            mod_base_table=readff.get_analysis_dataset(latest_basecall, 'BaseCalled_template/ModBaseProbs')
            readrecord=baminfo.find(read_id)
            
            for read in readrecord:
                ##make sure the record is not unmapped or secondary (calmd doesn't give these md tags)
                if not read.is_unmapped and not read.is_secondary:
                    align=np.array(read.get_aligned_pairs(with_seq=True))
                    chrom=read.reference_name
                    if read.is_reverse:
                        seq=rev_comp(read.get_forward_sequence())
                        strand='-'
                    else:
                        seq=read.get_forward_sequence()
                        strand='+'
                    minrefpos=read.reference_start
                    maxrefpos=read.reference_end
                    firstmotif=sum(motifpos[chrom]<minrefpos)+1
                    lastmotif=sum(motifpos[chrom]<maxrefpos)
                    
                    ###go through each motif position in the alignemnt range (so you're not scanning thru the whole genome)
                    for pos in motifpos[chrom][firstmotif:lastmotif]:
                        alignidx=np.where(align[:,1]==pos)[0][0]
                        alignpos=align[alignidx]

                        ##if the mod position is a deletion, all bets are off
                        if alignpos[0] == None:
                            q.put(','.join(map(lambda x: str(x), [read_id, 'del', '-', chrom, alignpos[1], alignpos[2], 'del', 'del', 'False', strand, 'del']))+'\n')

                        else:
                            ##check for recognizable motif on the read, assume uniform length motifs for now
                            startidx=alignpos[0]-modpos
                            endidx=alignpos[0]-modpos+motiflen
                            readmotif=seq[startidx:endidx]
                            if readmotif in motifs:
                                motif_correct='True'
                            else:
                                motif_correct='False'
                                                
                            if read.is_reverse:
                                modprobs=mod_base_table[len(seq)-alignpos[0]-motiflen+modpos+modpos]
                            else:
                                modprobs=mod_base_table[alignpos[0]]
                            q.put(','.join(map(lambda x: str(x), [read_id, alignpos[0], seq[alignpos[0]], chrom, alignpos[1], alignpos[2], modprobs[1], modprobs[3], motif_correct, strand, readmotif]))+'\n')

                        
    print(fast5+' took '+str(time.time()-start_time) + ' seconds')

    

def listener(q, outfile):
    '''
    writes from q, a manager.Queue() that is accumulating stuff
    '''
    with open(outfile, 'w') as f:
        while True:
            m=q.get()
            if m=='Done now, ty 4 ur service':
                break
            f.write(m)
            f.flush()

            
def main():
    args=parseArgs()

    ##prep motifs and ref info
    motifs=expand_motifs([args.motif])
    fastadict=fasta_dict(args.ref)
    motifpos=motif_position(fastadict, motifs, args.pos)

    
    ##get fast5 list
    fast5s=glob.glob(args.raw+'/*fast5')
    
    ##set up parallel jobs for each read
    manager=mp.Manager()
    q=manager.Queue()
    pool=mp.Pool(args.threads)
    watcher=pool.apply_async(listener, (q, args.out))


    ##submit jobs that will accumulate in the watcher
    jobs=[]
    for i in fast5s:
        job=pool.apply_async(read_mods, (i, args.bam, motifpos, motifs, args.pos, q))
        jobs.append(job)

        
    for job in jobs:
        job.get()
    q.put('Done now, ty 4 ur service')
    pool.close()
    pool.join()

    
if __name__ == "__main__":
    main()



    
