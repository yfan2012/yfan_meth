import h5py
import numpy as np

def get_cond_probs(bcfile):
    '''
    take the basecalls.hdf5 file and return a list of lists
    [meth_status, condprob]
    '''
    hf=h5py.File(bcfile, 'r')
    hfreads=hf.get('Reads')
    rnames=np.array(hfreads)
    allprobs=[]
    numreads=0
    ##===================limit to 10k reads rn because it won't write more than 30k for some reason??? ===============
    for i in rnames:
        if len(rnames)>0:
            condprobs=hf['Reads/'+i]
        methprobs=condprobs[np.ndarray.tolist(np.where(~np.isnan(condprobs))[0])]
        allprobs+=[str(x[0]) for x in np.ndarray.tolist(methprobs)]
        numreads+=1
        if numreads % 10000 == 0:
            print(str(numreads))
    return allprobs


def main(bcfile, outfile):
    '''
    write out a list of meth cond probs
    '''
    probs=get_cond_probs(bcfile)
    open(outfile, 'w').write('\n'.join(probs))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='write out a list of log condprobs from the basecall hdf file')
    parser.add_argument('--bcfile','-i',  help='input hdf path', type=str, required=True)
    parser.add_argument('--outfile','-o',  help='output text path', type=str, required=True)
    args = parser.parse_args()
    main(args.bcfile, args.outfile)
    
