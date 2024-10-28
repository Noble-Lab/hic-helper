import sys
import os
import pyBigWig

import numpy as np
from collections import defaultdict
def read_peak_info(input_bed):
    peak_dict=defaultdict(list)
    with open(input_bed,'r') as f:
        for line in f:
            line = line.strip()
            if len(line)==0:
                continue
            arr = line.split('\t')
            chrom = arr[0]
            start = int(arr[1])
            end = int(arr[2])
            peak_dict[chrom].append((start,end))
    return peak_dict


def report_peak_strength(input_bw, input_bed, output_bed):
    #get the peak information in a dict
    peak_dict = read_peak_info(input_bed)
    bw = pyBigWig.open(input_bw)
    chroms = bw.chroms()
    return_peak_dict={}
    for chrom in peak_dict:
        chrom_size = chroms[chrom]
        print("Chrom size:", chrom_size)
        signal = bw.values(chrom, 0, chrom_size, numpy=True)
        signal = np.nan_to_num(signal)
        current_peak_list = peak_dict[chrom]
        for peak in current_peak_list:
            start = peak[0]
            end = peak[1]
            value = np.mean(signal[start:end])
            return_peak_dict[chrom].append((start,end,value))
    bw.close()

    with open(output_bed,'w') as f:
        for chrom in return_peak_dict:
            for peak in return_peak_dict[chrom]:
                f.write(f"{chrom}\t{peak[0]}\t{peak[1]}\t{peak[2]}\n")
"""
This script is used to report the peak strength of the input bed file.
```
python3 report_peak_strength.py [input.bw] [input.bed] [output.bed]
```
[input.bw]: the input bigwig file. <br>
[input.bed]: the input bed file. <br>
[output.bed]: the output bed file, with last column represents the peak strength. <br>

"""


if __name__ == '__main__':
    if len(sys.argv)!=4:
        print("Usage: python3 report_peak_strength.py [input.bw] [input.bed] [output.bed]")
        print("This script is used to report the peak strength of the input bed file.")
        print("[input.bw]: the input bigwig file.")
        print("[input.bed]: the input bed file.")
        print("[output.bed]: the output bed file, with last column represents the peak strength.")
        sys.exit(1)

    input_bw = os.path.abspath(sys.argv[1])
    input_bed = os.path.abspath(sys.argv[2])
    output_bed = os.path.abspath(sys.argv[3])
    output_dir = os.path.dirname(output_bed)
    os.makedirs(output_dir, exist_ok=True)
    report_peak_strength(input_bw, input_bed, output_bed)
    print("Peak strength reported successfully.", output_bed)