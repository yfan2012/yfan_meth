import argparse
import multiprocessing
import pysam

def parseArgs():
    parser=argparse.ArgumentParser(description='get a meth barcode per read for megalodon')
    parser.add_argument('-i', '--input', type=str, required=True, help='megalodon output file with mod motif probs')
    parser.add_argument('-r', '--ref', type=str, required=True, help='reference genome used in megaldon')
    parser.add_argument('-m', '--motif', type=str, required=True, help='motif list file')
    parser.add_argument('-o', '--out', type=str, required=True, help='output file that lists each read and each barcode number')
    parser.add_argument('-t', '--threads', type=int, required=True, help='number of threads to use')
    parser.add_argument('-l', '--numreads', type=int, required=False, help='how many reads to consider')
    args=parser.parse_args()
    return args


