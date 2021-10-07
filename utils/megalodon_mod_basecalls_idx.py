import argparse

def parseArgs():
    parser=argparse.ArgumentParser(description='index megalodon mod_basecalls text file so you can parallel later on')
    parser.add_argument('-i', '--inmods', type=str, required=True, help='megalodon mod_basecalls file to index')
    parser.add_argument('-o', '--outidx', type=str, required=True, help='index file to write')
    args=parser.parse_args()
    return args

##inspired by nanocompore indexting
##https://github.com/tleonardi/nanocompore
def make_idx(modfile, idxfile):
    '''
    index the mod_basecalls
    adapted from nanocompore
    note that each character is 1 byte
    '''
    with open(idxfile, 'w') as f:
        f.write('\t'.join(['readname', 'chrom', 'byte_offset', 'byte_len', 'strand'])+'\n')
        modinfo=open(modfile, 'r')
        byteoff=0
        bytelen=0
        chrname=''
        readname=''
        strand=''
        for line in modinfo:
            if byteoff!=0:
                ##if still on same read/chr, add to line
                if readname==line.split('\t')[0] and chrname==line.split('\t')[1]:
                    bytelen+=len(line)
                ##if first of file
                elif readname=='':
                    readname=line.split('\t')[0]
                    chrname=line.split('\t')[1]
                    strand=line.split('\t')[2]
                    bytelen+=len(line)
                ##if read ended
                else:
                    f.write('\t'.join([readname, chrname, str(byteoff), str(bytelen), strand])+'\n')
                    byteoff+=bytelen
                    readname=line.split('\t')[0]
                    chrname=line.split('\t')[1]
                    strand=line.split('\t')[2]
                    bytelen=len(line)
            else:
                ##take header into account
                byteoff+=len(line)
        ##write last read
        f.write('\t'.join([readname, chrname, str(byteoff), str(bytelen), strand]))
        f.close()

def main(modfile, idxfile):
    make_idx(modfile, idxfile)
    
if __name__ == "__main__":
    args=parseArgs()
    main(args.inmods, args.outidx)
            
        
                        
    
