
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
def write_bed(output_path, loc_dict):
    with open(output_path, "w") as f:
        for chrom in loc_dict:
            for loc in loc_dict[chrom]:
                f.write(chrom + "\t" + str(loc[0]) + "\t" + str(loc[1]) + "\n")
def generate_peak_overlap(input1_bed, input2_bed, overlap_ratio, output_dir):
    input_dict1 = extract_loc(input1_bed)
    input_dict2 = extract_loc(input2_bed)
    overlap_dict1 = defaultdict(list)
    overlap_dict2 = defaultdict(list)
    independent_dict1 = defaultdict(list)
    independent_dict2 = defaultdict(list)
    for chrom in input_dict1:
        if chrom not in input_dict2:
            independent_dict1[chrom] = input_dict1[chrom]
            continue
        loc1 = np.array(input_dict1[chrom])
        loc2 = np.array(input_dict2[chrom])
        #sort based on the 0 axis
        loc1 = loc1[np.argsort(loc1[:,0])]
        loc2 = loc2[np.argsort(loc2[:,0])]
        i = 0
        j = 0
        while i < len(loc1) and j < len(loc2):
            start1 = loc1[i][0]
            end1 = loc1[i][1]
            start2 = loc2[j][0]
            end2 = loc2[j][1]
            if (start1 <= start2 and end1 >= start2) or (start2 <= start1 and end2 >= start1):
                #calculate current overlap ratio
                overlap_len = min(end1, end2) - max(start1, start2)
                cur_overlap_ratio = max(overlap_len/(end1-start1), overlap_len/(end2-start2))
                if cur_overlap_ratio >= overlap_ratio:
                    overlap_dict1[chrom].append([start1, end1])
                    overlap_dict2[chrom].append([start2, end2])
                    i += 1
                    j += 1
                else:
                    if start1 < start2:
                        independent_dict1[chrom].append([start1, end1])
                        i += 1
                    else:
                        independent_dict2[chrom].append([start2, end2])
                        j += 1
            else:
                if start1 < start2:
                    independent_dict1[chrom].append([start1, end1])
                    i += 1
                else:
                    independent_dict2[chrom].append([start2, end2])
                    j += 1
        print("Chromosome", chrom, "done")
        print("Overlap peaks in input1.bed:", len(overlap_dict1[chrom]))
        print("Overlap peaks in input2.bed:", len(overlap_dict2[chrom]))
        print("Independent peaks in input1.bed:", len(independent_dict1[chrom]))
        print("Independent peaks in input2.bed:", len(independent_dict2[chrom]))
    write_bed(os.path.join(output_dir, "overlap1.bed"), overlap_dict1)
    write_bed(os.path.join(output_dir, "overlap2.bed"), overlap_dict2)
    write_bed(os.path.join(output_dir, "independent1.bed"), independent_dict1)
    write_bed(os.path.join(output_dir, "independent2.bed"), independent_dict2)
"""
This script is to compare two bed files and find the overlapping peaks. <br>
```
python3 peak_overlap.py [input1.bed] [input2.bed] [overlap_ratio] [output_dir]
```
- input1.bed: the first input bed file  <br>
- input2.bed: the second input bed file  <br>
- overlap_ratio: the ratio of overlap to consider as overlap  <br>
- output_dir: the output directory  <br>
- The output files are overlap1.bed, overlap2.bed, independent1.bed, independent2.bed, indicating the overlap peaks in input1.bed, overlap peaks in input2.bed, independent peaks in input1.bed, independent peaks in input2.bed.
"""

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python3 peak_overlap.py [input1.bed] [input2.bed] [overlap_ratio] [output_dir]")
        print("input1.bed: the first input bed file")
        print("input2.bed: the second input bed file")
        print("overlap_ratio: the ratio of overlap to consider as overlap")
        print("The ratio is calculated based on the length of the shorter peak.")
        print("output_dir: the output directory")
        print("The output files are overlap1.bed, overlap2.bed, independent1.bed, independent2.bed, indicating the overlap peaks in input1.bed, overlap peaks in input2.bed, independent peaks in input1.bed, independent peaks in input2.bed")
        sys.exit(1)
    input1_bed = os.path.abspath(sys.argv[1])
    input2_bed = os.path.abspath(sys.argv[2])
    overlap_ratio = float(sys.argv[3])
    output_dir = os.path.abspath(sys.argv[4])
    os.makedirs(output_dir, exist_ok=True)
    generate_peak_overlap(input1_bed, input2_bed, overlap_ratio, output_dir)
