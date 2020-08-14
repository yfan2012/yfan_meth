import argparse
import pysam
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
    parser.add_argument('-s', '--offset', type=int, required=False,
                        help='range +/-[int] to look at from the mod position')
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


def motif_position(fastadict, motifs, modpos, offset):
    '''
    need fastadict, list of motifs, position of modified base
    get dict of seqname:list of motif positions +/-4 from the putative mod
    '''
    motifpos={}
    for seq in fastadict:
        motifpos[seq]=[]
        for motif in motifs:
            for i in re.finditer(motif, fastadict[seq]):
                if offset==0:
                    motifpos[seq].extend(range(i.start()+modpos, i.start()+modpos+1))
                else:
                    motifpos[seq].extend(range(i.start()+modpos-offset, i.start()+modpos+offset))
            for i in motifpos:
                motifpos[i]=list(set(motifpos[i]))
    return motifpos


def read_mods(fast5, bamfile, motifpos, q):
##def read_mods(fast5, bamfile, motifpos):
    '''
    reads fast5file, indexed pysam align object bam, motifpos dict
    returns interesting info on each read contained therein
    '''
    start_time=time.time()

    print('Reading bam file')
    ##store bamfiles
    bam=pysam.AlignmentFile(bamfile, 'rb')
    baminfo=pysam.IndexedReads(bam)
    baminfo.build()

    ###motifinfo will be a list of lists [readname, index_on_read, read_base,refname,  index_on_ref, ref_base, pAmod, pCmod]
    motifinfo=''

    with get_fast5_file(fast5, mode="r") as f5:
        for read_id in f5.get_read_ids():
            readff=f5.get_read(read_id)
            latest_basecall=readff.get_latest_analysis('Basecall_1D')
            mod_base_table=readff.get_analysis_dataset(latest_basecall, 'BaseCalled_template/ModBaseProbs')

            readrecord=baminfo.find(read_id)
            for read in readrecord:
                ##make sure the record is not unmapped or secondary (calmd doesn't give these md tags)
                if not read.is_unmapped and not read.is_secondary:
                    align=read.get_aligned_pairs(with_seq=True)
                    chrom=read.reference_name
                    seq=read.get_forward_sequence()


                    ###go through each position in the read
                    for pos in align:
                        
                        ##store info if the reference position is listed 
                        if pos[1] in motifpos[chrom]:
                            ##if the interesting position is a deletion on the read
                            if pos[0]==None:
                                motifinfo+=','.join(map(lambda x: str(x), [read_id, 'del', 'del', chrom, pos[1], pos[2], 'del', 'del']))+'\n'
                            else:
                                modprobs=mod_base_table[pos[0]]
                                readbase=seq[pos[0]]
                                motifinfo+=','.join(map(lambda x: str(x), [read_id, pos[0], readbase, chrom, pos[1], pos[2], modprobs[1], modprobs[3]]))+'\n'

                                
                        ##append positioal info if ref position is None (insertion in read) and it's preceeded or proceeded by an interesting position
                        elif pos[1]==None and pos[0]!=None:
                            modprobs=mod_base_table[pos[0]]
                            readbase=seq[pos[0]]
                            prelist=list(filter(lambda x: pos[0]-1==x[0], align))
                            prolist=list(filter(lambda x: pos[0]+1==x[0], align))
                            if len(prelist)>0:
                                if prelist[0][1] in motifpos[chrom]:
                                    motifinfo+=','.join(map(lambda x: str(x), [read_id, pos[0], readbase, chrom, pos[1], pos[2], modprobs[1], modprobs[3]]))+'\n'
                            elif len(prolist)>0:
                                if prolist[0][1] in motifpos[chrom]:
                                    motifinfo+=','.join(map(lambda x: str(x), [read_id, pos[0], readbase, chrom, pos[1], pos[2], modprobs[1], modprobs[3]]))+'\n'
                                    

    q.put(motifinfo)
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
    motifpos=motif_position(fastadict, motifs, args.pos, args.offset)

    
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
        job=pool.apply_async(read_mods, (i, args.bam, motifpos, q))
        jobs.append(job)

        
    for job in jobs:
        job.get()
    q.put('Done now, ty 4 ur service')
    pool.close()
    pool.join()

    
if __name__ == "__main__":
    main()



    
