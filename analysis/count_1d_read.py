import sys
import os
import pyBigWig
import pickle 
import numpy as np

def bigwig2count(input_bw):
    bw = pyBigWig.open(input_bw)
    chroms = bw.chroms()
    signal_dict = {}
    total_chrom = 0
    total_read = 0
    for chrom in chroms:
        print("Processing", chrom)
        chrom_size = chroms[chrom]
        print("Chrom size:", chrom_size)
        signal = bw.stats(chrom)
        total_read += np.sum(signal)*chrom_size
        total_chrom += chrom_size
    bw.close()
    print("Total reads:", total_read)
    print("Total chrom size:", total_chrom)
    print("Average read count per base:", total_read/total_chrom)
    return total_read/total_chrom
"""
This script is to calculate the average read count per base and total read count in the bigwig file.
```
python3 bigwig2count.py [input_bw]
```
[input_bw]: the input bigwig file. <br>

"""
if __name__ == '__main__':
    if len(sys.argv)!=2:
        print("Usage: python bigwig2count.py [input_bw]")
        print("input_bw: the input bigwig file")
        sys.exit(1)
    input_bw = sys.argv[1]
    coverage=bigwig2count(input_bw)



