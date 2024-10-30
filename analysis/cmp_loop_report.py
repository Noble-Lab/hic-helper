
import os
import sys 
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd
import seaborn as sns 
import numpy as np
#import cdist
from scipy.spatial.distance import cdist
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
def compare_loops(control_loc, input_loc,max_scope):
    count_loop_lost = 0
    count_loop_gain = 0
    count_loop_same = 0
    dist_matrix = cdist(control_loc, input_loc)
    #other_merge_count_input=0
    visit_dict={}
    for i in range(len(control_loc)):
        min_dist = np.min(dist_matrix[i])
        min_idx = int(np.argmin(dist_matrix[i]))
        if min_dist<=max_scope:
            all_min_idx = np.argwhere(dist_matrix[i]<=max_scope).flatten()
            cur_visit = False
            for idx in all_min_idx:
                idx = int(idx)
                if idx in visit_dict:
                    cur_visit = True
                    continue
                visit_dict[idx] = 1
            if cur_visit is False:
                count_loop_same += 1
            # if min_idx in visit_dict:
            #     continue 
            # count_loop_same += 1
            # visit_dict[min_idx] = 1
        else:
            count_loop_lost += 1
    count_loop_gain = len(input_loc) - len(visit_dict)
    return count_loop_lost, count_loop_gain, count_loop_same
def cmp_loop_report(control_bed, input_bed, output_pdf,max_scope=5):
    control_dict = extract_loc(control_bed)
    input_dict = extract_loc(input_bed)
    count_loop_lost = 0
    count_loop_gain = 0
    count_loop_same = 0
    for chrom in control_dict:
        if chrom not in input_dict:
            print("Chromosome ",chrom," not in prediction")
            continue
        control_loc = np.array(control_dict[chrom])
        input_loc = np.array(input_dict[chrom])
        cur_loop_lost, cur_loop_gain, cur_loop_same = compare_loops(control_loc, input_loc,max_scope)
        count_loop_lost += cur_loop_lost
        count_loop_gain += cur_loop_gain
        count_loop_same += cur_loop_same
    total_amount = count_loop_lost + count_loop_gain + count_loop_same
    print("Loop Lost: %d (%.2f)"%(count_loop_lost, count_loop_lost/total_amount))
    print("Loop Gain: %d (%.2f)"%(count_loop_gain, count_loop_gain/total_amount))
    print("Loop Same: %d (%.2f)"%(count_loop_same, count_loop_same/total_amount))
    #plot the figure make venn diagram
    plt.figure()
    plattee = sns.color_palette("pastel")
    plt.pie([count_loop_lost, count_loop_gain, count_loop_same], 
            labels=["Loop Lost", "Loop Gain", "Loop Unchanged"],
            colors=plattee,autopct='%1.1f%%', textprops={'fontsize': 14})
    #plt.legend(fontsize=16,loc="upper right")
    plt.tight_layout()
    plt.savefig(output_pdf, dpi=600)


"""
This script is used to compare the loop change between two bed files.
```
python3 cmp_loop_report.py [control.bed] [input.bed] [resolution] [output.pdf]
```
[control.bed]: the control bed file. <br>
[input.bed]: the input bed file. <br>
[resolution]: the resolution of the data. <br>
[output.pdf]: the output pdf/png file. It is a pie chart showing the loop change. <br>
"""





if __name__ == '__main__':
    if len(sys.argv)!=5:
        print("Usage: python3 cmp_loop_report.py [control.bed] [input.bed] [resolution] [output.pdf]")
        print("This script is used to compare the loop change between two bed files.")
        print("[control.bed]: the control bed file.")
        print("[input.bed]: the input bed file.")
        print("[resolution]: the resolution of the data.")
        print("[output.pdf]: the output pdf/png file. It is a pie chart showing the loop change.")
        sys.exit(1)
    control_bed = os.path.abspath(sys.argv[1])
    input_bed = os.path.abspath(sys.argv[2])
    resolution = int(sys.argv[3])
    output_pdf = os.path.abspath(sys.argv[4])
    output_dir = os.path.dirname(output_pdf)
    os.makedirs(output_dir, exist_ok=True)
    cmp_loop_report(control_bed, input_bed, output_pdf,max_scope=5*resolution)