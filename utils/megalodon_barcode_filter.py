import argparse
from megalodon_barcode import fasta_dict, expand_motif, expand_barcodes

def parseArgs():
    parser=argparse.ArgumentParser(description='filter barcoded reads based on alignment')
    parser.add_argument('-a', '--alignfile', type=str, required=True, help='alignment file, type assumed based on extention')
    parser.add_argument('-r', '--reffile', type=str, required=True, help='reference genome used in megaldon')
    parser.add_argument('-m', '--megafile', type=str, required=True, help='file with barcodes')
    parser.add_argument('-b', '--barcodefile', type=str, required=True, help='motif list file')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='output file for filtered barcoded reads')
    args=parser.parse_args()
    return args


iupacnt={
    'N': ['A', 'C', 'G', 'T'],
    'W': ['A', 'T']
}


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
        if i.end-i.start > length:
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


def filter_nummotifs(aligns, ref, barcodes, minmotifs):
    '''
    filter for number of motifs that need to be present
    takes in list of reads as alignment objects
    returns list of reads as alignment objects
    '''
    readlist=[]
    for i in aligns:
        refseq=ref[i.refname][i.refstart:i.refend]
        motifcounts={}
        for motif in barcodes:
            motifcounts[motif]=0
            for j in barcodes[motif]:
                motifcounts[motif]+=refseq.count(j)
        if all(x>minmotifs for x in your_dict.values()):
            readlist.append(i)

    
def main(alignfile, reffile, megafile, barcodefile, outfile):
    ref=fasta_dict(reffile)
    barcodes=expand_barcodes(barcodefile)
    minmapq=30
    minlen=5000
    minmotifs=15
    
    if alignfile.split('.')[-1]=='paf':
        pafaligns=get_align_ranges_paf(alignfile)
    
    
                
if __name__ == "__main__":
    args=parseArgs()
    main(args.alignfile, args.reffile, args.megafile, args.barcodefile, args.outfile)
