def mumseqs(snps, ref, seqlen):
    '''get list of sequences of length 2*seqlen. don't need seq identifiers'''
    import sys
    sys.path.insert(0, '/home-4/yfan7@jhu.edu/Code/utils')
    from fasta_utils import fasta_dict

    with open(snps, 'r') as f:
        content=f.read().splitlines()
    pos=[]
    for i in content:
        pos.append([i.split('\t')[10], int(i.split('\t')[0])])

    refseqs=fasta_dict(ref)
    ##get rid of extras in chrom name
    for i in refseqs:
        newkey=i.split(' ')[0]
        refseqs[newkey]=refseqs.pop(i)

    ##list the seqs
    seqs=[]
    for i in pos:
        seqs.append(refseqs['>'+i[0]][i[1]-seqlen:i[1]+seqlen])
    return seqs


def main(snps, ref, seqlen, outfile):
    seqs=mumseqs(snps, ref, seqlen)
    with open(outfile, 'w') as f:
        for x in seqs:
            if 'AAAA' not in x and 'TTTT' not in x and 'GGGG' not in x and 'CCCC' not in x:
                f.write(x+'\n')

        
if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser(description = 'make a file thats just a list of sequences where a sindel has occured')
    parser.add_argument('--snps', '-s', type=str,  help='input mummer snps file path')
    parser.add_argument('--ref', '-r', type=str, help='path the reference used to make the mummer snp file')
    parser.add_argument('--seqlen', '-l', type=int, help='half the sequence length you want')
    parser.add_argument('--outfile', '-o', type=str, required=True,  help='output path, txt file')
    args=parser.parse_args()
    main(args.snps, args.ref, args.seqlen, args.outfile)

