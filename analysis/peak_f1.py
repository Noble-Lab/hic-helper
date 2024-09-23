
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
                start1 = int(float(line[1]))
                end1 = int(float(line[2]))
            except:
                continue
            
            if "chr" not in chrom1:
                chrom1 = "chr"+chrom1
            overall_dict[chrom1].append([start1, end1])
            
    return overall_dict

def check_close_dist(true_loc,pred_loc):
    start1 = true_loc[0]
    end1 = true_loc[1]
    start2 = pred_loc[0]
    end2 = pred_loc[1]
    dist1 = abs(start1 - start2)
    dist2 = abs(end1 - end2)
    dist3 = abs(start1 - end2)
    dist4 = abs(end1 - start2)
    return min(dist1, dist2, dist3, dist4)
def calculate_peak_f1(true_dict, pred_dict, max_dist):
    true_peaks = 0
    pred_peaks = 0
    matched_peaks = 0
    for chrom in true_dict:
        if chrom not in pred_dict:
            print("Chromosome ",chrom," not in prediction")
            continue
        true_loc = np.array(true_dict[chrom])
        pred_loc = np.array(pred_dict[chrom])
        true_peaks += len(true_loc)
        pred_peaks += len(pred_loc)
        for i in range(len(true_loc)):
            for j in range(len(pred_loc)):
                if check_close_dist(true_loc[i],pred_loc[j]) <= max_dist:
                    matched_peaks += 1
                    break
    recall = matched_peaks / true_peaks
    precision = matched_peaks / pred_peaks
    f1score = 2 * precision * recall / (precision + recall)
    print("True peaks: ", true_peaks)
    print("Predicted peaks: ", pred_peaks)
    print("Matched peaks: ", matched_peaks)
    return f1score


"""
```
python3 peak_f1.py [true.bed] [pred.bed] [max_dist]
```
[true.bed]: the true peaks, in bed format <br>
[pred.bed]: the predicted peaks, in bed format <br>
[max_dist]: the maximum distance to match the peaks <br>

"""
if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python3 peak_f1.py [true.bed] [pred.bed] [max_dist]")
        print("ture.bed: the true peaks, in bed format")
        print("pred.bed: the predicted peaks, in bed format")
        print("max_dist: the maximum distance to match the peaks")
        sys.exit(1)
    
    true_bed = sys.argv[1]
    pred_bed = sys.argv[2]
    max_dist = int(sys.argv[3])
    true_dict = extract_loc(true_bed)
    pred_dict = extract_loc(pred_bed)

    peak_f1 = calculate_peak_f1(true_dict, pred_dict, max_dist)
    print("Peak F1 score: ", peak_f1)