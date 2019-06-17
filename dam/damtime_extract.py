import os, sys, re
import numpy as np
import pandas as pd

def find_motif(ref, motif):
    '''
    take reference and motif
    return list of positions of the motif in the reference
    assumes 1 chr reference
    '''
    sys.path.insert(0, '/home-4/yfan7@jhu.edu/Code/utils')
    from fasta_utils import fasta_dict
    fasta_dictionary=fasta_dict(ref)
    positions={}
    for i in fasta_dictionary:
        positions[i]=[m.start() for m in re.finditer(motif, fasta_dictionary[i])]
    return positions


def collapse_dwells(position_events, whichmode):
    '''
    take a list of lists of eventalign
    return list of dwell times (where events at the same position are combined)
    this only really applies for tabix files where there's only one read
    '''
    ##header if statement is just for dam because the read tabix and the position tabix aren't consistent
    ##make the two consistent in future
    if len(position_events[0])==13:
        headers=['chr', 'position', 'refmer', 'read', 'strand', 'event', 'current_mean', 'current_std', 'dwell', 'model_kmer', 'model_mean', 'model_std', 'level']
    else:
        headers=['chr', 'position', 'refmer', 'read', 'strand', 'event', 'current_mean', 'current_std', 'dwell', 'model_kmer']
    df=pd.DataFrame(position_events, columns=headers)
    df['dwell']=df.dwell.astype('float')
    if whichmode == 'read':
        dwell_sum=df.groupby(['position'])['dwell'].agg('sum')
    else:
        dwell_sum=df.groupby(['read'])['dwell'].agg('sum')
    dwell_sum_list=dwell_sum.values.tolist()
    return dwell_sum_list


def get_position_events(tabpos, locus):
    '''
    call system to query tabix and save results in a tmpfile
    return the data in a list of lists
    '''
    tmpname=tabpos.split('.')[0]+'.'+str(locus)+'.tmp.tsv'
    print 'Creating '+tmpname
    os.system('tabix '+ tabpos + ' "gi|730582171|gb|CP009644.1|":' + str(locus) + '-' + str(locus+1) + ' &> ' + tmpname)
    with open(tmpname) as f:
        content=f.read().split('\n')
        f.close()
    position_events=[]
    for i in content:
        a=i.split('\t')
        ##trailing newline in the tsvs sometimes, so check that there's actually stuff
        if len(i)>0:
            a.pop(2)
            position_events.append(a)
    print 'Removing '+tmpname
    os.system('rm ' + tmpname)
    ##position_dwells=collapse_dwells(position_events)
    ##return position_events
    return position_events


def get_read_norm(read_number, tabread):
    '''
    call system to query tabix for read number
    return a list [read_num, median_dwell]
    '''
    tmpread=tabread.split('.')[0]+'.read'+str(read_number)+'.tmp.tsv'
    os.system('tabix ' + tabread + ' "gi|730582171|gb|CP009644.1|":' + str(read_number) + '-' + str(read_number) + ' &> ' + tmpread)
    with open(tmpread) as f:
        content=f.read().split('\n')
        f.close()
    read_events=[]
    for i in content:
        if len(i)>0:
            read_events.append(i.split('\t'))
    print 'Removing '+tmpread
    os.system('rm ' + tmpread)
    collapsed=collapse_dwells(read_events, 'read')
    med=np.median(collapsed)
    return [read_number, med]


def main(ref, motif, tabpos, tabread, normout, dwellout):
    '''
    take reference and list of positions and corresponding meth data
    give back norm data, dwell times for each position queried (and relevant kmer)
    '''
    positions=find_motif(ref, motif)
    ##load norm data as dictionary of read_number:median_dwell
    with open(normout) as f:
        content=f.read().split('\n')
        f.close()
    norm={}
    for i in content:
        norm[i.split(',')[0]]=float(i.split(',')[1])
    ##position analysis
    for chrom in positions:
        for i in positions[chrom]:
            posdwells=[]
            locus=int(i)-10
            position_events=get_position_events(tabpos, locus)
            for j in position_events:
                ##there's something weird about calling tabix for read 0 and read 1 that I don't quite understand
                if j[3] != '0' and j[3] != '1':
                    ##get the normalization time for the read at this event if not already there
                    if j[3] not in norm:
                        newnorm=get_read_norm(j[3], tabread)
                        norm[newnorm[0]]=newnorm[1]
                        ##periodically save norm info 
                        if len(norm)%100 == 0:
                            with open(normout, 'w') as f:
                                for i in norm:
                                    f.write(','.join([i, str(norm[i])])+'\n')
                                f.close()
                    ##divide the dwell time by the normazlization
                    ##whoa, python will actually change the thing in place. cool.
                    j[8]=str(float(j[8])/norm[j[3]])
            position_dwells=collapse_dwells(position_events, 'position')
            for j in position_dwells:
                posdwells.append([str(chrom.split(' ')[0]), str(locus), str(j)])
            ##append this position's info
            with open(dwellout, 'a') as f:
                for i in posdwells:
                    f.write(','.join(i)+'\n')
                f.close()

                
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='find reads that aligned exactly one time')
    parser.add_argument('--ref','-r',  help='reference genome', type=str, required=True)
    parser.add_argument('--motif','-m',  help='motif', type=str, required=True)
    parser.add_argument('--norm','-n',  help='normalization csv file', type=str, required=True)
    parser.add_argument('--tabpos','-p',  help='tabix by position', type=str, required=True)
    parser.add_argument('--tabread','-s',  help='tabix by read', type=str, required=True)
    parser.add_argument('--outdwells','-o',  help='output of dwell times', type=str, required=True)
    args = parser.parse_args()
        
    main(args.ref, args.motif, args.norm, args.tabpos, args.tabread, args.outdwells)
    
