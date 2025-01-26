import sys
import os
import pyBigWig
import pickle 
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
def read_bed(input_bed):
    locus_dict=defaultdict(list)
    with open(input_bed) as f:
        for line in f:
            line=line.strip()
            if line=="":
                continue
            lst=line.split()
            chrom=lst[0]
            start=int(lst[1])
            end=int(lst[2])
            locus_dict[chrom].append((start, end))
    return locus_dict


def plot_distribution(input_bw, input_peak, output_fig):
    #read peak file
    locus_dict=read_bed(input_peak)
    bw = pyBigWig.open(input_bw)
    chroms = bw.chroms()
    value_list = []
    for chrom in chroms:
        print("Processing", chrom)
        chrom_size = chroms[chrom]
        print("Chrom size:", chrom_size)
        signal = bw.values(chrom, 0, chrom_size, numpy=True)
        current_locus = locus_dict[chrom]
        if len(current_locus)==0:
            continue
        for locus in current_locus:
            start, end = locus
            signal_locus = signal[start:end]
            signal_locus = np.nan_to_num(signal_locus)
            signal_average = np.mean(signal_locus)
            value_list.extend(signal_average)
    bw.close()
    plt.hist(value_list, bins=100)
    plt.xlabel("Peak strength per base")
    plt.ylabel("Frequency")
    plt.title("Peak distribution")
    #plt.yscale('log') #very necessary because many 0 values
    plt.savefig(output_fig,dpi=600)
"""
This script plots the peak distribution of the bigwig file according to the peak region specified in .bed file.
```
python3 bigwig_peak_distribution.py [input.bw] [input.peak] [output_fig]
```
[input.bw]: the input bigwig file. <br>
[input.peak]: the input peak file,specify the peak region. <br>
[output_fig]: the output figure path to show the peak distribution. <br>


"""


if __name__ == '__main__':
    if len(sys.argv)!=4:
        print("Usage: python3 bigwig_peak_distribution.py [input.bw] [input.peak] [output_fig]")
        print("input.bw: the input bigwig file")
        print("input.peak: the input peak file,specify the peak region")
        print("output_fig: the output figure path to show the peak distribution")
        sys.exit(1)
    input_bw = os.path.abspath(sys.argv[1])
    input_peak = os.path.abspath(sys.argv[2])
    output_fig = os.path.abspath(sys.argv[3])
    output_dir = os.path.dirname(output_fig)
    os.makedirs(output_dir, exist_ok=True)
    plot_distribution(input_bw, input_peak, output_fig)