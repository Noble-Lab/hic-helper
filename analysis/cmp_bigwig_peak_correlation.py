
import os
import sys 
import pyBigWig
import numpy as np
from scipy import stats
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

def caclulate_correlation(input_bigwig1, input_bigwig2, input_bed):
    locus_dict=read_bed(input_bed)
    bw1 = pyBigWig.open(input_bigwig1)
    bw2 = pyBigWig.open(input_bigwig2)
    report_dict=defaultdict(list)
    chroms = bw1.chroms()
    for chrom in chroms:
        current_locus = locus_dict[chrom]
        if len(current_locus)==0:
            continue
        chrom_size = chroms[chrom]
        print("chrom",chrom, " Chrom size:", chrom_size)
        signal1 = bw1.values(chrom, 0, chrom_size, numpy=True)
        signal1 = np.nan_to_num(signal1)
        signal2 = bw2.values(chrom, 0, chrom_size, numpy=True)
        signal2 = np.nan_to_num(signal2)
        for locus in current_locus:
            start, end = locus
            signal1_locus = signal1[start:end]
            signal2_locus = signal2[start:end]
            if len(signal1_locus)==0 or len(signal2_locus)==0:
                continue
            if np.sum(signal1_locus)==0 or np.sum(signal2_locus)==0:
                report_dict['pearson'].append(0)
                report_dict['spearman'].append(0)
                report_dict['cosine'].append(0)
                continue
            #calculate pearson correlation
            res = stats.pearsonr(signal1_locus, signal2_locus)
            report_dict['pearson'].append(res[0])
            #add spearman correlation
            res = stats.spearmanr(signal1_locus, signal2_locus)
            report_dict['spearman'].append(res[0])
            #add cosine similarity
            res = np.dot(signal1_locus, signal2_locus)/(np.linalg.norm(signal1_locus)*np.linalg.norm(signal2_locus))
            report_dict['cosine'].append(res)
        print("finished chrom", chrom)
        print("Pearson correlation:", np.mean(report_dict['pearson']))
        print("Spearman correlation:", np.mean(report_dict['spearman']))
        print("Cosine similarity:", np.mean(report_dict['cosine']))
    bw1.close()
    bw2.close()
    print("Pearson correlation:", np.mean(report_dict['pearson']))
    print("Spearman correlation:", np.mean(report_dict['spearman']))
    print("Cosine similarity:", np.mean(report_dict['cosine']))
"""
This script is used to calculate the correlation between two bigwig files on the locus specified in .bed file.
```
python3 cmp_bigwig_peak_correlation.py [input1.bigWig] [input2.bigWig] [reference.bed]
```
[input1.bigWig]: the first bigwig file. <br>
[input2.bigWig]: the second bigwig file. <br>
[reference.bed]: the bed file containing the locus to calculate the correlation. <br>
This script will output pearson correlation, spearman correlation, and cosine similarity between the two bigwig files. <br>

"""

if __name__ == '__main__':
    if len(sys.argv)!=4:
        print("Usage: python3 cmp_bigwig_peak_correlation.py [input1.bigWig] [input2.bigWig] [reference.bed]")
        print("This script is used to calculate the correlation between two bigwig files on the locus specified in .bed file.")
        print("[input1.bigWig]: the first bigwig file.")
        print("[input2.bigWig]: the second bigwig file.")
        print("[reference.bed]: the bed file containing the locus to calculate the correlation.")
        sys.exit(1)

    input_bigwig1 = os.path.abspath(sys.argv[1])
    input_bigwig2 = os.path.abspath(sys.argv[2])
    input_bed = os.path.abspath(sys.argv[3])
    caclulate_correlation(input_bigwig1, input_bigwig2, input_bed)