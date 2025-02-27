
import sys
import numpy as np
import os 
from collections import defaultdict



def extract_loop_loc(pred_detect_path):
    overall_dict = defaultdict(list)
    with open(pred_detect_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip().split()
            chrom1 = line[0]
            try:
                start1 = int(float(line[1]))
                end1 = int(float(line[2]))
            except:
                continue
            chrom2 = line[3]
            try:
                start2 = int(float(line[4]))
                end2 = int(float(line[5]))
            except:
                continue
            if "chr" not in chrom1:
                chrom1 = "chr"+chrom1
            overall_dict[chrom1].append([start1, start2])
            
    return overall_dict


def extract_peak_bed(pred_detect_path):
    overall_dict = defaultdict(list)
    with open(pred_detect_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip().split()
            chrom1 = line[0]
            try:
                start1 = int(float(line[1]))
                end1 = int(float(line[2]))
            except:
                continue
            
            if "chr" not in chrom1:
                chrom1 = "chr"+chrom1
            overall_dict[chrom1].append([start1, end1])
            
    return overall_dict

def check_overlap(start1, end1, peak_list):
    """
    start1: the start of the loop
    end1: the end of the loop
    peak_list: the list of peaks
    """
    for peak in peak_list:
        start2 = peak[0]
        end2 = peak[1]
        if start1 <= start2 and end1 >= start2:
            return True
        if start2 <= start1 and end2 >= start1:
            return True
    return False

def calculate_loop_ctcf_ratio(loop_dict, peak_dict, resolution):
    loop_ctcf_ratio={}
    loop_singlectcf_ratio={}
    for chrom in loop_dict:
        cur_loop_list= loop_dict[chrom]
        cur_peak_list = peak_dict[chrom]
        if len(cur_peak_list)==0:
            continue
        count_match= 0
        count_singlematch = 0
        for loop_loc in cur_loop_list:
            start1 = loop_loc[0]
            start2 = loop_loc[1]
            start1 = int(start1/resolution)
            start2 = int(start2/resolution)
            start1 = start1*resolution
            start2 = start2*resolution
            check_start1 = start1 - resolution #3 pixel to check if the loop overlap with CTCF peak.
            check_end1 = start1 + resolution*2 

            check_start2 = start2 - resolution
            check_end2 = start2 + resolution*2

            match_flag1 = check_overlap(check_start1, check_end1, cur_peak_list)
            match_flag2 = check_overlap(check_start2, check_end2, cur_peak_list)
            if match_flag1 and match_flag2:
                count_match += 1
            if match_flag1 or match_flag2:
                count_singlematch += 1
        loop_ctcf_ratio[chrom] = count_match / len(cur_loop_list)
        loop_singlectcf_ratio[chrom] = count_singlematch / len(cur_loop_list)
        print(chrom, "loop_ctcf_ratio:", loop_ctcf_ratio[chrom])
        print(chrom, "loop_singlectcf_ratio:", loop_singlectcf_ratio[chrom])
    avg_loop_ctcf_ratio = np.mean(list(loop_ctcf_ratio.values()))
    avg_loop_singlectcf_ratio = np.mean(list(loop_singlectcf_ratio.values()))
    print("Average loop_ctcf_ratio:", avg_loop_ctcf_ratio)
    print("Average loop_singlectcf_ratio:", avg_loop_singlectcf_ratio)

"""
This script is used to calculate the ratio of chromatin loops that overlap with CTCF ChIP peaks.
```
python3 loop_ctcfpeak_ratio.py [loop.bed] [ctcf_peak.bed] [resolution]
```
[loop.bed]: the chromatin loop coordinate, in bed format <br>
[ctcf_peak.bed]: the CTCF peak coordinate, in bed format <br>
[resolution]: the resolution of the Hi-C data <br>
"""

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python3 loop_ctcfpeak_ratio.py [loop.bed] [ctcf_peak.bed] [resolution]")
        print("loop.bed: the chromatin loop coordinate, in bed format")
        print("ctcf_peak.bed: the CTCF peak coordinate, in bed format")
        print("resolution: the resolution of the Hi-C data")
        sys.exit(1)
    loop_bed = os.path.abspath(sys.argv[1])
    ctcf_peak_bed = os.path.abspath(sys.argv[2])
    resolution = int(sys.argv[3])

    loop_dict = extract_loop_loc(loop_bed)
    peak_dict = extract_peak_bed(ctcf_peak_bed)
    calculate_loop_ctcf_ratio(loop_dict, peak_dict, resolution)






