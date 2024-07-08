import sys
import os
from collections import defaultdict
import pyBigWig
import numpy as np
def extract_loc(pred_detect_path):
    record_list=[]
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
            
            if "chr" not in chrom2: 
                chrom2 = "chr"+chrom2

            record_list.append([chrom1, start1, end1, chrom2, start2, end2])
    print("Extracted",len(record_list)," loop detection records")
    return record_list

def clean_loops(input_bed, mappability_bw, output_bed, threshold):
    #read the input bed file get all records
    record_info = extract_loc(input_bed)
    #read the mappability bigwig file
    final_record_list = []
    with pyBigWig.open(mappability_bw) as input_bw:
       
        for k in range(len(record_info)):
            #chrom,start,end =peak_region[k]
            chrom1, start1, end1, chrom2, start2, end2 = record_info[k]
            # Get the signal values (reads) within the peak
            signal_values = input_bw.values(chrom1, start1, end1)

            # Sum the signal values (reads) within the peak
            signal_values=np.array(list(signal_values))
            signal_values = np.nan_to_num(signal_values)
            mappable_percentage1 = np.sum(signal_values>0.5)/len(signal_values)
        

            signal_values = input_bw.values(chrom2, start2, end2)
            signal_values=np.array(list(signal_values))
            signal_values = np.nan_to_num(signal_values)
            mapable_percentage2 = np.sum(signal_values>0.5)/len(signal_values)
            if mappable_percentage1 < threshold or mapable_percentage2 < threshold:
                continue

            final_record_list.append([chrom1, start1, end1, chrom2, start2, end2])
    print("After filtering, there are",len(final_record_list),"loops left")
    with open(output_bed, 'w') as f:
        f.write("#chrom1\tstart1\tend1\tchrom2\tstart2\tend2\n")
        for record in final_record_list:
            f.write("\t".join(map(str, record)) + "\n")

"""
This script is for filter out the loops on low mappability regions
```
python3 loop_cleaner.py [input.bed] [mappablility.bw] [output.bed] [threshold]
```
- input.bed: the input bed file <br>
- mappablility.bw: the mappablility bigwig file <br>
- output.bed: the output bed file <br>
- threshold: the mappablility threshold used to clean loops <br>
"""
if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python3 loop_cleaner.py [input.bed] [mappablility.bw] [output.bed] [threshold]")
        print("This script is for filter out the loops on low mappability regions")
        print("input.bed: the input bed file")
        print("mappablility.bw: the mappablility bigwig file")
        print("output.bed: the output bed file")
        print("threshold: the mappablility percentage in the loop region to clean loops")
        sys.exit(1)
    input_bed = os.path.abspath(sys.argv[1])
    mappability_bw = os.path.abspath(sys.argv[2])
    output_bed = os.path.abspath(sys.argv[3])
    output_dir = os.path.dirname(output_bed)
    os.makedirs(output_dir,exist_ok=True)
    threshold = float(sys.argv[4])

    clean_loops(input_bed, mappability_bw, output_bed, threshold)