import sys 
import numpy as np
import os 
from collections import defaultdict
from scipy.spatial.distance import cdist
def extract_loc(pred_detect_path):
    overall_dict = defaultdict(list)
    with open(pred_detect_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            if "x1" in line or "start" in line:
                print("skip header",line)
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
def get_cluster_loc(map_dict, map_dict_reverse, key1,key2,visit_set1,visit_set2):
    
    if key1 is not None:
        for idx in map_dict[key1]:
            if idx not in visit_set2:
                visit_set2.add(idx)
                visit_set1,visit_set2 = get_cluster_loc(map_dict, map_dict_reverse, None, idx,visit_set1,visit_set2)
    if key2 is not None:
           
        for idx in map_dict_reverse[key2]:
            if idx not in visit_set1:
                visit_set1.add(idx)
                visit_set1,visit_set2 = get_cluster_loc(map_dict, map_dict_reverse, idx, None,visit_set1,visit_set2)
                
    return visit_set1,visit_set2
    

def compare_loops(control_loc, input_loc, max_dist_allowance):
    dist_matrix = cdist(control_loc, input_loc)
    independent1 = []
    independent2 = []
    map_dict = defaultdict(set)
    map_dict_reverse = defaultdict(set)
    for i in range(len(control_loc)):
        min_dist = np.min(dist_matrix[i])
        min_idx = int(np.argmin(dist_matrix[i]))
        if min_dist > max_dist_allowance:
            independent1.append(control_loc[i])
        else:
            all_min_idx = np.argwhere(dist_matrix[i]<=max_dist_allowance).flatten()
            
            for idx in all_min_idx:
                idx = int(idx)
                map_dict[i].add(idx)
                map_dict_reverse[idx].add(i)
    for i in range(len(input_loc)):
        if i not in map_dict_reverse:
            independent2.append(input_loc[i])
    common_loc = []
    visit_dict={}
    for key in map_dict:
        if key in visit_dict:
            continue
        key_set1,key_set2 = get_cluster_loc(map_dict, map_dict_reverse, key,None,set(),set())
        for tmp in key_set1:
            visit_dict[tmp] = 1
        #get average positions
        location_list=[]
        for idx in key_set1:
            location_list.append(control_loc[idx])
        for idx in key_set2:
            location_list.append(input_loc[idx])
        location_list = np.array(location_list)#change to n*3
        avg_loc = np.mean(location_list, axis=0)
        common_loc.append(avg_loc)
    return independent1, independent2, common_loc
def write_bed(output_path, loc_dict,resolution):
    #write loop bed file
    with open(output_path, "w") as f:
        for chrom in loc_dict:
            cur_list = loc_dict[chrom]
            for loc in cur_list:
                start1 = loc[0]
                start2 = loc[1]
                end1 = loc[0]+resolution
                end2 = loc[1]+resolution
                f.write(chrom + "\t" + str(start1) + "\t" + str(end1) + "\t" + chrom + "\t" + str(start2) + "\t" + str(end2) + "\n")

def count_loop(input_dict):
    count = 0
    for chrom in input_dict:
        count += len(input_dict[chrom])
    return count
def generate_loop_overlap(control_bed, input_bed,resolution, output_dir,max_dist=5):
    contro_dict = extract_loc(control_bed)
    input_dict = extract_loc(input_bed)
    overlap_dict = defaultdict(list)
    independent1_dict = defaultdict(list)
    independent2_dict = defaultdict(list)
    for chrom in contro_dict:
        if chrom not in input_dict:
            independent1_dict[chrom] = contro_dict[chrom]
            continue
        control_loc = np.array(contro_dict[chrom])
        input_loc = np.array(input_dict[chrom])
        #get cdist matrix
        indepedent_loc1, indepedent_loc2,common_loc = compare_loops(control_loc, input_loc, 
                                                                    max_dist_allowance=resolution*max_dist)

        overlap_dict[chrom] = common_loc
        independent1_dict[chrom] = indepedent_loc1
        independent2_dict[chrom] = indepedent_loc2
        print("Chromosome", chrom, "done")
    #write the output
    write_bed(os.path.join(output_dir, "overlap_loop.bed"), overlap_dict,resolution)
    write_bed(os.path.join(output_dir, "control_independent.bed"), independent1_dict,resolution)
    write_bed(os.path.join(output_dir, "input_independent.bed"), independent2_dict,resolution)
    #output each category number
    print("Overlap loops:", count_loop(overlap_dict))
    print("Independent loops in control:", count_loop(independent1_dict))
    print("Independent loops in input:", count_loop(independent2_dict))
"""
This script is used to compare the loop change between two bed files and outputs independent/overlap loops.
```
python3 loop_overlap.py [control.bed] [input.bed] [resolution] [output_dir]
```
[control.bed]: the control bed file recording the control loop location. <br>
[input.bed]: the input bed file recording the input loop location. <br>
[resolution]: the resolution of the Hi-C data. <br>
[output_dir]: the output directory. <br>
The output files are overlap.bed, independent1.bed, independent2.bed, indicating the overlap loops, independent loops in control, independent loops in input. <br>

"""



if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python3 loop_overlap.py [control.bed] [input.bed] [resolution] [output_dir]")
        print("This script is used to compare the loop change between two bed files and outputs independent/overlap loops.")
        print("[control.bed]: the control bed file recording the loop location.")
        print("[input.bed]: the input bed file recording the loop location.")
        print("[resolution]: the resolution of the Hi-C data.")
        print("[output_dir]: the output directory.")
        print("The output files are overlap.bed, independent1.bed, independent2.bed, indicating the overlap loops, independent loops in control, independent loops in input.")
        sys.exit(1)
    control_bed = os.path.abspath(sys.argv[1])
    input_bed = os.path.abspath(sys.argv[2])
    resolution = int(sys.argv[3])
    output_dir = os.path.abspath(sys.argv[4])
    os.makedirs(output_dir, exist_ok=True)
    generate_loop_overlap(control_bed, input_bed, resolution, output_dir)
