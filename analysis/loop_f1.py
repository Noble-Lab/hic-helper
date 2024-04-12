
import sys
import numpy as np
import os 
from collections import defaultdict

def extract_loc(pred_detect_path):
    overall_dict = defaultdict(list)
    with open(pred_detect_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip().split()
            chrom1 = line[0]
            try:
                start1 = int(line[1])
                end1 = int(line[2])
            except:
                continue
            chrom2 = line[3]
            try:
                start2 = int(line[4])
                end2 = int(line[5])
            except:
                continue
            overall_dict[chrom1].append([start1, start2])
            
    return overall_dict
def calculate_chromosome_loop_f1(predict_loc, gt_loc,max_scope=5):
    """
    predict_loc: nd array N*3
    gt_loc: nd array M*3
    max_scope: max distance to match
    """
    recall = np.zeros(len(gt_loc))
    precision = np.zeros(len(predict_loc))
    for i,pred_peak in enumerate(predict_loc):
        for j,tgt_peak in enumerate(gt_loc):
            dist=np.linalg.norm(pred_peak-tgt_peak)
            if dist<=max_scope:
                recall[j]=1
                precision[i]=1
    recall = np.mean(recall)
    precision = np.mean(precision)
    f1score = 2*precision*recall/(precision+recall)
    return f1score
def calculate_loop_f1(true_dict, pred_dict, resolution,max_scope=5):
    max_scope = max_scope*resolution
    loop_f1 = {}
    for chrom in true_dict:
        if chrom not in pred_dict:
            loop_f1[chrom] = 0
            continue
        true_loc = np.array(true_dict[chrom])
        pred_loc = np.array(pred_dict[chrom])
        loop_f1[chrom] = calculate_chromosome_loop_f1(pred_loc, true_loc,max_scope)
    return loop_f1


#input two bed files, one with the true peaks and one with the predicted peaks
#output the f1 score of the predicted peaks
"""
```
python3 loop_f1.py [true.bed] [pred.bed] [resolution]
```
[true.bed]: the true peaks, in bed format <br>
[pred.bed]: the predicted peaks, in bed format <br>
[resolution]: the resolution of the Hi-C data <br>

"""
if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python3 loop_f1.py true.bed pred.bed resolution")
        print("ture.bed: the true peaks, in bed format")
        print("pred.bed: the predicted peaks, in bed format")
        print("resolution: the resolution of the Hi-C data")
        sys.exit(1)
    true_bed = sys.argv[1]
    pred_bed = sys.argv[2]
    resolution = int(sys.argv[3])
    true_dict = extract_loc(true_bed)
    pred_dict = extract_loc(pred_bed)
    loop_f1=calculate_loop_f1(true_dict, pred_dict, resolution)
    print("Chromosome-wise F1: ",loop_f1)
    #average across different chromosomes
    print("Overall F1:",np.mean(list(loop_f1.values() ) ) )

