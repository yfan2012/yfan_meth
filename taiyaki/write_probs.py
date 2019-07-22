import h5py
import numpy as np
import pysam

def get_cond_probs(bcfile, fafile, motif, methposition):
    '''
    take the basecalls.hdf5 file and return a list of lists
    [condprob, seq context]
    '''
    hf=h5py.File(bcfile, 'r')
    hfreads=hf.get('Reads')
    rnames=np.array(hfreads)
    fa=pysam.FastaFile(fafile)
    allprobs=[]
    numreads=0
    ##===================limit to 10k reads rn because it won't write more than 30k for some reason??? ===============
    for i in rnames[0:1000]:
        if len(rnames)>0:
            condprobs=hf['Reads/'+i]
        rseq=fa.fetch(i)
        whereprobs=np.ndarray.tolist(np.where(~np.isnan(condprobs))[0])
        methprobs=condprobs[whereprobs]
        for x in whereprobs:
            if rseq[x-methposition:x+len(motif)-methposition] == motif:
                allprobs+=[[str(condprobs[x][0]),rseq[x-methposition:x+len(motif)-methposition]]]
        numreads+=1
        if numreads % 10000 == 0:
            print(str(numreads))
    return allprobs


def main(bcfile, fafile, motif, methposition, outfile):
    '''
    write out a list of meth cond probs
    '''
    probs=get_cond_probs(bcfile, fafile, motif, methposition)
    with open(outfile, 'w') as f:
        for i in probs:
            f.write(','.join(i)+'\n')

            
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='write out a list of log condprobs from the basecall hdf file')
    parser.add_argument('--bcfile','-i',  help='input hdf path', type=str, required=True)
    parser.add_argument('--outfile','-o',  help='output text path', type=str, required=True)
    parser.add_argument('--fafile', '-f', help='basecalled fa file', type=str, required=True)
    parser.add_argument('--motif', '-m', help='motif', type=str, required=True)
    parser.add_argument('--pos', '-p', help='methylated position in motif, 0 based indexing', type=int, required=True)    
    args = parser.parse_args()
    main(args.bcfile, args.fafile, args.motif, args.pos, args.outfile)
    
