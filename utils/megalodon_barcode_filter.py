import argparse
from megalodon_barcode import fasta_dict, expand_motif, expand_barcodes
import time

def parseArgs():
    parser=argparse.ArgumentParser(description='filter barcoded reads based on alignment')
    parser.add_argument('-a', '--alignfile', type=str, required=True, help='alignment file, type assumed based on extention')
    parser.add_argument('-r', '--reffile', type=str, required=True, help='reference genome used in megaldon')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='output file for filtered barcoded reads')
    parser.add_argument('-m', '--megafile', type=str, required=False, help='file with barcodes')
    parser.add_argument('-b', '--barcodefile', type=str, required=False, help='motif list file')
    parser.add_argument('-q', '--mapq', type=int, required=False, help='min mapq')
    parser.add_argument('-l', '--minlen', type=int, required=False, help='motif list file')
    parser.add_argument('-n', '--nummotifs', type=int, required=False, help='motif list file')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose outputs')
    args=parser.parse_args()
    return args

iupacnt={
    'N': ['A', 'C', 'G', 'T'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'V': ['G', 'C', 'A']
}

##for testing
##alignfile='/mithril/Data/Nanopore/projects/methbin/align/neb11/neb11_sub.paf'
##reffile='/mithril/Data/Nanopore/projects/methbin/reference/allsamps.fa'
##megafile='/mithril/Data/Nanopore/projects/methbin/barcode/qc/neb11_barcodes.txt'
##barcodefile='/home/yfan/Code/yfan_nanopore/mdr/qc/barcodes.txt'


class alignment:
    '''
    stores alignment info needed for potentially filtering
    '''
    def __init__(self, readname, refname, refstart, refend, mapq):
        self.readname=readname
        self.refname=refname
        self.refstart=int(refstart)
        self.refend=int(refend)
        self.mapq=int(mapq)
        self.barcode=None
    def getbarcode(self, barcode):
        self.barcode=barcode

    
def get_align_ranges_paf(paffile):
    '''
    take in paf alignment file
    return dict with readname:[]
    '''
    with open(paffile, 'r') as f:
        content=f.read().split('\n')
    pafaligns=[]
    for i in content:
        if len(i)>0:
            r=i.split('\t')
            pafaligns.append(alignment(r[0], r[5], r[7], r[8], r[11]))
    return pafaligns


def filter_length(aligns, length):
    '''
    filter for alignment length
    takes in list of reads as alignment objects
    returns list of reads as alignment objects
    '''
    readlist=[]
    for i in aligns:
        if i.refend-i.refstart > length:
            readlist.append(i)
    return readlist


def filter_mapq(aligns, minmapq):
    '''
    filter for mapq
    takes in list of reads as alignment objects
    returns list of reads as alignment objects
    '''
    readlist=[]
    for i in aligns:
        if i.mapq > minmapq:
            readlist.append(i)
    return readlist


def filter_nummotifs(aligns, ref, barcodes, minmotifs, verbose):
    '''
    filter for number of motifs that need to be present
    takes in list of reads as alignment objects
    returns list of reads as alignment objects
    '''
    start_time=time.time()
    readlist=[]
    count=0
    for i in aligns:
        count+=1
        if verbose and count % 50000==0:
            elapsed=round(time.time()-start_time)
            start_time=time.time()
            print('%d reads in %d seconds' % (count, elapsed))
        refseq=ref[i.refname][i.refstart:i.refend]
        motifcounts={}
        for motif in barcodes:
            motifcounts[motif]=0
            for j in barcodes[motif]:
                motifcounts[motif]+=refseq.count(j)
        if all(x>minmotifs for x in motifcounts.values()):
            readlist.append(i)
    return readlist


def barcode_info(aligns, megafile, verbose):
    '''
    add barcode info to alignments
    takes in list of alignment objects wihtout barcode info
    returns list of alignment objects with barcode info
    '''
    start_time=time.time()
    barcodeinfo={}
    with open(megafile, 'r') as f:
        content=f.read().split('\n')
    for i in content:
        readinfo=i.split('\t')
        barcodeinfo[readinfo[0]]=readinfo[1:]
    count=0
    for i in aligns:
        count+=1
        if verbose and count % 50000==0:
            elapsed=round(time.time()-start_time)
            start_time=time.time()
            print('%d reads in %d seconds' % (count, elapsed))
        if i.readname in barcodeinfo: #aligned read might have been filtered out already
            i.getbarcode(barcodeinfo[i.readname])
    return aligns


def main(alignfile, reffile, megafile, barcodefile, outfile, minmapq, minlen, minmotifs, verbose):
    ref=fasta_dict(reffile)
    barcodes=expand_barcodes(barcodefile)
    
    if alignfile.split('.')[-1]=='paf':
        if verbose:
            print('Reading paf file')
        aligns=get_align_ranges_paf(alignfile)

    if minmapq is not None:
        if verbose:
            print('Filtering for mapq')
        aligns=filter_mapq(aligns, minmapq)
        numpass=len(aligns)
        if verbose:
            print('Number of reads passing mapq filter: %d' % (numpass))
    if minlen is not None:
        print('Filtering for length')
        aligns=filter_length(aligns, minlen)
        numpass=len(aligns)
        if verbose:
            print('Number of reads passing length filter: %d' % (numpass))
    if minmotifs is not None:
        print('Filtering for number of motifs')
        aligns=filter_nummotifs(aligns, ref, barcodes, minmotifs, verbose)
        numpass=len(aligns)
        if verbose:
            print('Number of reads passing min motif filter: %d' % (numpass))

    if verbose:
        print('Matching barcode info')
    aligned_barcodes=barcode_info(aligns, megafile, verbose)
        
    with open(outfile, 'w') as f:
        for i in aligned_barcodes:
            if i.barcode is not None:
                barcodes='\t'.join(i.barcode)
                f.write('\t'.join([i.readname, i.refname, str(i.refstart), str(i.refend), str(i.mapq), barcodes])+'\n')
            
        
                
if __name__ == "__main__":
    args=parseArgs()
    main(args.alignfile, args.reffile, args.megafile, args.barcodefile, args.outfile, args.mapq, args.minlen, args.nummotifs, args.verbose)
