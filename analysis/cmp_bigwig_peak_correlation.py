
import os
import sys 
import pyBigWig
import numpy as np
from scipy import stats
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
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

def caclulate_correlation(input_bigwig1, input_bigwig2, input_bed,output_png):
    locus_dict=read_bed(input_bed)
    bw1 = pyBigWig.open(input_bigwig1)
    bw2 = pyBigWig.open(input_bigwig2)
    report_dict=defaultdict(list)
    chroms = bw1.chroms()
    chroms2 = bw2.chroms()
    chroms_list = set(chroms.keys()) & set(chroms2.keys())
    count_cmp_dict={"x":[], "y":[]}
    for chrom in chroms_list:
        current_locus = locus_dict[chrom]
        if len(current_locus)==0:
            continue
        chrom_size = chroms[chrom]
        print("chrom",chrom, " Chrom size:", chrom_size)
        chrom_size2 = chroms2[chrom]
        print("chrom",chrom, " Chrom size2:", chrom_size2)
        chrom_size = min(chrom_size, chrom_size2)
        
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
            #if std is 0, then the correlation is 0
            if np.std(signal1_locus)==0 or np.std(signal2_locus)==0:
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
            count_cmp_dict["x"].append(np.sum(signal1_locus))
            count_cmp_dict["y"].append(np.sum(signal2_locus))


        print("finished chrom", chrom)
        print("Pearson correlation:", np.mean(report_dict['pearson']))
        print("Spearman correlation:", np.mean(report_dict['spearman']))
        print("Cosine similarity:", np.mean(report_dict['cosine']))
    bw1.close()
    bw2.close()
    print("Profile correlation:")
    print("Pearson correlation:", np.mean(report_dict['pearson']))
    print("Spearman correlation:", np.mean(report_dict['spearman']))
    print("Cosine similarity:", np.mean(report_dict['cosine']))
    print("*"*50)
    print("Peak total reads comparison:")
    pearon_count = stats.pearsonr(count_cmp_dict["x"], count_cmp_dict["y"])[0]
    spearman_count = stats.spearmanr(count_cmp_dict["x"], count_cmp_dict["y"])[0]
    print("Pearson correlation:", pearon_count)
    print("Spearman correlation:", spearman_count)
    plt.figure(figsize=(8, 8))
    #plt.hist2d(count_cmp_dict["x"], count_cmp_dict["y"], bins=100, cmap=plt.cm.jet)
    count_cmp_df = pd.DataFrame(count_cmp_dict)
    sns.kdeplot(
    data=count_cmp_df, x="x", y="y", fill=True,
    )
    #set x scale, y scale to 95% percentile xmax,ymax
    xmax = np.percentile(count_cmp_dict["x"], 95)
    ymax = np.percentile(count_cmp_dict["y"], 95)
    plt.xlim(0, xmax)
    plt.ylim(0, ymax)
    plt.xlabel("Peak reads in bigwig1",fontsize=18)
    plt.ylabel("Peak reads in bigwig2",fontsize=18)
    plt.title("Pearson: %.2f"%pearon_count,fontsize=18)
    plt.savefig(output_png,dpi=300)

"""
This script is used to calculate the correlation between two bigwig files on the locus specified in .bed file.
```
python3 cmp_bigwig_peak_correlation.py [input1.bigWig] [input2.bigWig] [reference.bed] [output_png]
```
[input1.bigWig]: the first bigwig file. <br>
[input2.bigWig]: the second bigwig file. <br>
[reference.bed]: the bed file containing the locus to calculate the correlation. <br>
[output_png]: the output png file to save the peak total reads comparison. <br>
This script will output pearson correlation, spearman correlation, and cosine similarity between the two bigwig files. <br>

"""

if __name__ == '__main__':
    if len(sys.argv)!=5:
        print("Usage: python3 cmp_bigwig_peak_correlation.py [input1.bigWig] [input2.bigWig] [reference.bed] [output_png]")
        print("This script is used to calculate the correlation between two bigwig files on the locus specified in .bed file.")
        print("[input1.bigWig]: the first bigwig file.")
        print("[input2.bigWig]: the second bigwig file.")
        print("[reference.bed]: the bed file containing the locus to calculate the correlation.")
        print("[output_png]: the output png file to save the peak total reads comparison.")
        sys.exit(1)

    input_bigwig1 = os.path.abspath(sys.argv[1])
    input_bigwig2 = os.path.abspath(sys.argv[2])
    input_bed = os.path.abspath(sys.argv[3])
    output_png = os.path.abspath(sys.argv[4])
    output_dir = os.path.dirname(output_png)
    os.makedirs(output_dir, exist_ok=True)
    caclulate_correlation(input_bigwig1, input_bigwig2, input_bed,output_png)