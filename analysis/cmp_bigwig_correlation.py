import os
import sys 
import pyBigWig
import numpy as np
from scipy import stats
from collections import defaultdict
def calculate_correlation(input1, input2, resolution):
    bw1 = pyBigWig.open(input1)
    bw2 = pyBigWig.open(input2)
    chroms1 = bw1.chroms()
    chroms2 = bw2.chroms()
    report_dict=defaultdict(list)
    shared_chroms = set(chroms1.keys()) & set(chroms2.keys())
    for chrom in shared_chroms:
        print("Processing", chrom)
        chrom_size = chroms1[chrom]
        print("Chrom size:", chrom_size)
        chrom_size2 = chroms2[chrom]
        print("chrom",chrom, " Chrom size2:", chrom_size2)
        chrom_size = min(chrom_size, chrom_size2)
        signal1 = bw1.values(chrom, 0, chrom_size, numpy=True)
        signal2 = bw2.values(chrom, 0, chrom_size, numpy=True)
        
        #average the signal according to the resolution, and remove nan values
        signal1 = signal1[~np.isnan(signal1)]
        signal2 = signal2[~np.isnan(signal2)]
        if len(signal1)==0 or len(signal2)==0:
            continue
        #average the signal according to the resolution, every resolution window will be converted to one value
        signal1 = np.mean(signal1[:len(signal1)//resolution*resolution].reshape(-1, resolution), axis=1)
        signal2 = np.mean(signal2[:len(signal2)//resolution*resolution].reshape(-1, resolution), axis=1)
        #calculate pearson correlation
        res = stats.pearsonr(signal1, signal2)
        report_dict['pearson'].append(res[0])
        #add spearman correlation
        res = stats.spearmanr(signal1, signal2)
        report_dict['spearman'].append(res[0])
        #add cosine similarity
        res = np.dot(signal1, signal2)/(np.linalg.norm(signal1)*np.linalg.norm(signal2))
        report_dict['cosine'].append(res)
    bw1.close()
    bw2.close()
    print("Pearson correlation:", np.mean(report_dict['pearson']))
    print("Spearman correlation:", np.mean(report_dict['spearman']))
    print("Cosine similarity:", np.mean(report_dict['cosine']))
"""
This script is used to calculate the correlation between two bigwig files.
```
python3 cmp_bigwig_correlation.py [input1.bigWig] [input2.bigWig] [resolution]
```
[input1.bigWig]: the first bigwig file. <br>
[input2.bigWig]: the second bigwig file. <br>
[resolution]: the resolution to calculate the correlation. <br>
This script will output pearson correlation, spearman correlation, and cosine similarity between the two bigwig files. <br>
""" 



if __name__ == '__main__':
    if len(sys.argv)!=4:
        print("Usage: python3 cmp_bigwig_correlation.py [input1.bigWig] [input2.bigWig] [resolution]")
        print("This script is used to calculate the correlation between two bigwig files.")
        print("[input1.bigWig]: the first bigwig file.")
        print("[input2.bigWig]: the second bigwig file.")
        print("[resolution]: the resolution to calculate the correlation.")
        sys.exit(1)

    input1 = os.path.abspath(sys.argv[1])
    input2 = os.path.abspath(sys.argv[2])
    resolution = int(sys.argv[3])
    calculate_correlation(input1, input2, resolution)