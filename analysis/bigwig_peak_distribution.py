import sys
import os
import pyBigWig
import pickle 
import numpy as np
import matplotlib.pyplot as plt
def plot_distribution(input_bw, window_size, stride, output_fig,plot_mode):
    bw = pyBigWig.open(input_bw)
    chroms = bw.chroms()
    value_list = []
    for chrom in chroms:
        print("Processing", chrom)
        chrom_size = chroms[chrom]
        print("Chrom size:", chrom_size)
        signal = bw.values(chrom, 0, chrom_size, numpy=True)
        signal = np.nan_to_num(signal)
        for i in range(0, chrom_size-window_size, stride):
            value = np.mean(signal[i:i+window_size])
            if plot_mode==1:
                value = np.log10(value+1)
            value_list.append(value)
        print("Accumulated value list length:", len(value_list))
        print("Signal stats: mean ",np.mean(value_list), "std ", np.std(value_list), "max ", np.max(value_list), "min ", np.min(value_list))
    bw.close()
    plt.hist(value_list, bins=100)
    plt.xlabel("Peak Strength")
    plt.ylabel("Frequency")
    plt.title("Peak Distribution")
    plt.yscale('log') #very necessary because many 0 values
    plt.savefig(output_fig,dpi=600)
"""
This script plots the peak distribution of the bigwig file.
```
python3 bigwig_peak_distribution.py [input.bw] [window_size] [stride] [output_fig] [mode]
```
[input.bw]: the input bigwig file. <br>
[window_size]: the window size for analyzing peak distribution. <br>
[stride]: the stride for analyzing peak distribution. <br>
[output_fig]: the output figure path to show the peak distribution. <br>
[mode]: 0:raw_value, 1:log10_value. <br>
"""
if __name__ == '__main__':
    if len(sys.argv)!=6:
        print("Usage: python3 bigwig_peak_distribution.py [input.bw] [window_size] [stride] [output_fig] [mode]")
        print("input.bw: the input bigwig file")
        print("window_size: the window size for analyzing peak distribution")
        print("stride: the stride for analyzing peak distribution")
        print("output_fig: the output figure path to show the peak distribution")
        print("mode: 0:raw_value, 1:log10_value")
        sys.exit(1)
    input_bw = os.path.abspath(sys.argv[1])
    window_size = int(sys.argv[2])
    stride = int(sys.argv[3])
    output_fig = os.path.abspath(sys.argv[4])
    mode=int(sys.argv[5])
    output_dir = os.path.dirname(output_fig)
    os.makedirs(output_dir, exist_ok=True)
    plot_distribution(input_bw, window_size, stride, output_fig, mode)

