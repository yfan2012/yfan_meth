def find_motif(fasta_dictionary, motif):
    import re
    import sys
    positions={}
    for i in fasta_dictionary:
        positions[i]=[m.start() for m in re.finditer(motif, fasta_dictionary[i])]
        sys.stdout.write(i)
        for j in positions[i]:
            sys.stdout.write(str(j)+'\n')

def main(fastafile, motif):
    import sys
    sys.path.insert(0, '/home-4/yfan7@jhu.edu/Code/utils')
    from fasta_utils import fasta_dict
    fasta_dictionary=fasta_dict(fastafile)
    find_motif(fasta_dictionary, motif)
    

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='list out each motif in a ref')
    parser.add_argument('--ref','-r',  help='input ref path', type=str, required=True)
    parser.add_argument('--motif', '-m', help='motif', type=str, required=True)
    args = parser.parse_args()

    main(args.ref, args.motif)
