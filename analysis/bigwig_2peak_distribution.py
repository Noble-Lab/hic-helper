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


def plot_distribution(input_bw, positive_peak,negative_peak, output_fig,mode):
    #read peak file
    positive_dict =read_bed(positive_peak)
    negative_dict =read_bed(negative_peak)
    bw = pyBigWig.open(input_bw)
    chroms = bw.chroms()
    value_list = []
    neg_value_list = []
    for chrom in chroms:
        print("Processing", chrom)
        chrom_size = chroms[chrom]
        print("Chrom size:", chrom_size)
        signal = bw.values(chrom, 0, chrom_size, numpy=True)
        current_locus = positive_dict[chrom]
        if len(current_locus)==0:
            continue
        for locus in current_locus:
            start, end = locus
            signal_locus = signal[start:end]
            signal_locus = np.nan_to_num(signal_locus)
            signal_average = np.mean(signal_locus)
            value_list.append(signal_average)
        current_locus = negative_dict[chrom]
        if len(current_locus)==0:
            continue
        for locus in current_locus:
            start, end = locus
            signal_locus = signal[start:end]
            signal_locus = np.nan_to_num(signal_locus)
            signal_average = np.mean(signal_locus)
            neg_value_list.append(signal_average)
    bw.close()
    plt.hist(value_list, bins=100,label="Positive",alpha=0.5)
    plt.hist(neg_value_list, bins=100,label="Negative",alpha=0.5)
    plt.legend(fontsize=18)
    plt.xlabel("Peak strength per base")
    plt.title("Peak distribution")
    if mode==1:
        plt.yscale('log')
        
    plt.ylabel("Frequency")
    plt.savefig(output_fig,dpi=600)
"""
This script compares the positive and negative peak distribution of the bigwig file according to the peak region specified in .bed files.
```
python3 bigwig_2peak_distribution.py [input.bw] [positive.bed] [negative.bed] [output_fig] [mode]
```
[input.bw]: the input bigwig file. <br>
[positive.bed]: the input positive peak file,specify the positive peak region. <br>
[negative.bed]: the input negative peak file,specify the negative peak region. <br>
[output_fig]: the output figure path to show the peak distribution. <br>
[mode]: 0:raw_value, 1:log10_value. <br>


"""


if __name__ == '__main__':
    if len(sys.argv)!=6:
        print("Usage: python3 bigwig_2peak_distribution.py [input.bw] [positive.peak] [negative.peak] [output_fig] [mode]")
        print("input.bw: the input bigwig file")
        print("positive.peak: the input positive peak file,specify the positive peak region")
        print("negative.peak: the input negative peak file,specify the negative peak region")
        print("output_fig: the output figure path to show the peak distribution")
        print("mode:  0:raw_value, 1:log10_value")
        sys.exit(1)
    input_bw = os.path.abspath(sys.argv[1])
    positive_peak = os.path.abspath(sys.argv[2])
    negative_peak = os.path.abspath(sys.argv[3])
    output_fig = os.path.abspath(sys.argv[4])
    mode = int(sys.argv[5])
    output_dir = os.path.dirname(output_fig)
    os.makedirs(output_dir, exist_ok=True)
    plot_distribution(input_bw, positive_peak,negative_peak, output_fig,mode)