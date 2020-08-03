import argparse
import pysam
import multiprocessing as mp
import glob
import re
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
    parser.add_argument('-p', '--pos', type=int, required=False, help='position of mod')
    parser.add_argument('-s', '--offset', type=int, required=False, help='range +/-[int] to look at from the mod position')
    parser.add_argument('-t', '--threads', type=int, required=False, help='num threads to use')
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


def read_mods(fast5, baminfo, motifpos, q):
    '''
    reads fast5
    needs pysam ref, indexed pysam align object bam (to avoid reading in so often), motifpos dict
    returns interesting info on each read contained therein
    list of lists [readid, num_motifs, num_motifs_error, num_canon, num_canon_mcalled, num_motifs_mcalled]
    '''
    with get_fast5_file(fast5, mode="r") as f5:
        for read_id in f5.get_read_ids():
            readff=f5.get_read(read_id)
            latest_basecall=readff.get_latest_analysis('Basecall_1D')
            mod_base_table=readff.get_analysis_dataset(latest_basecall, 'BaseCalled_template/ModBaseProbs')

            ###extract read from alignment file
            readrecord=baminfo.find(read_id)

            ###motifinfo will be a list of tuples (index_on_read, index_on_ref, ref_base, ref_seqname)
            ###really just expecting on read in readrecord, but just in case
            motifinfo=[]
            for read in readrecord:
                align=read.get_aligned_pairs(with_seq=True)
                seq=read.reference_name
                
                ###go through each position in the read
                motifseq=[]
                for pos in align:
                    ##append the tuple if ref position is listed as interesting
                    if pos[1] in motifpos[seq]:
                        motifseq.append(pos+(seq,))
                    ##append if ref position is None (insertion in read) and it's preceeded by an interesting position
                    elif pos[1]==None and len(motifseq)>0:
                        motifseq.append(pos+(seq,))
                    ##flush if the position is not interesting and something has accumulated in motifseq
                    elif pos[1]!=None and pos[1] not in motifpos[seq] and len(motifseq)>0:
                        motifinfo.extend(motifseq)
                        motifseq=[]

            ###add modcall info for each and write (put in write queue)
            ###modinfo will be list of lists in str form  [read_id, index_read, index_ref, ref_base, ref_seqname, p(A), p(6mA), p(C), p(5mC)]
            modinfo=''
            for i in motifinfo:
                if i[0]==None:
                    modinfo+=','.join([read_id, 'del', str(i[1]), i[2], i[3], 'None', 'None', 'None', 'None'])+'\n'
                else:
                    readmodinfo=mod_base_table[i[0]]
                    modinfo+=','.join([read_id, str(i[0]), str(i[1]), i[2], i[3], str(readmodinfo[0]), str(readmodinfo[1]), str(readmodinfo[2]), str(readmodinfo[3])])+'\n'
            q.put(modinfo)

    

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

    
    ##store bamfiles
    bam=pysam.AlignmentFile(args.bam, 'rb')
    baminfo=pysam.IndexedReads(bam)
    baminfo.build()
    
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
        job=pool.apply_async(read_mods,(i, baminfo, motifpos, q))
        jobs.append(job)

    for job in jobs:
        job.get()
    q.put('Done now, ty 4 ur service')
    pool.close()
    pool.join()




    
if __name__ == "__main__":
    main()



    
