import sys
import os
import pyBigWig
import pickle 
import numpy as np

def bigwig2coverage(input_bw):
    bw = pyBigWig.open(input_bw)
    chroms = bw.chroms()
    count_dict = {}
    for chrom in chroms:
        print("Processing", chrom)
        chrom_size = chroms[chrom]
        print("Chrom size:", chrom_size)
        #stats bigwig
        avg_cov=bw.stats(chrom,0,chrom_size, type="mean")

        coverage =  avg_cov*chrom_size
        count_dict[chrom] = coverage
        print(f"Finished processing {chrom}: {coverage} coverage")
    total_coverage = sum(count_dict.values())
    return total_coverage
"""
This script calculates the total coverage of a bigwig file.
```
python3 bigwig_coverage.py [input_bw]
```
[input_bw]: the input bigwig file. <br>

"""
if __name__ == '__main__':
    if len(sys.argv)!=2:
        print("Usage: python bigwig_coverage.py [input_bw]")
        print("input_bw: the input bigwig file")
        sys.exit(1)

    input_bw = os.path.abspath(sys.argv[1])
    cov=bigwig2coverage(input_bw)
    print("Total coverage:", cov)