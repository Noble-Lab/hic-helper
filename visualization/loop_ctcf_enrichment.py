import sys
import numpy as np
import os 
from collections import defaultdict
import pyBigWig


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
def identify_key(chrom,key_list):
    if chrom in key_list:
        return chrom
    if "chr" in chrom:
        chrom = chrom.replace("chr", "")
        if chrom in key_list:
            return chrom
    chrom = "chr"+chrom
    if chrom in key_list:
        return chrom
    return None


def extract_enrichment(loop_dict, ctcf_chip_bw, resolution):
    
    with pyBigWig.open(ctcf_chip_bw) as bw:
        chrom_sizes = bw.chroms()
        chrom_list = list(chrom_sizes.keys())
        
        output_list = []
        for chrom in chrom_list:
            current_chrom = identify_key(chrom, loop_dict.keys())
            if current_chrom is None:
                continue
            cur_size = chrom_sizes[current_chrom]
            for loop in loop_dict[current_chrom]:
                start1 = loop[0]//resolution*resolution
                end1 = (loop[0]+resolution)//resolution*resolution
                signal_start1 = max(0, start1-resolution)
                signal_end1 = min(cur_size, end1+resolution) #check 3pixel peak regions
                if signal_start1 > cur_size:
                    continue

                start2 = loop[1]//resolution*resolution
                end2 = (loop[1]+resolution)//resolution*resolution
                signal_start2 = max(0, start2-resolution)
                signal_end2 = min(cur_size, end2+resolution)
                if signal_start2 > cur_size:
                    continue
                signal1 = bw.stats(chrom, signal_start1, signal_end1, type="max")
                signal2 = bw.stats(chrom, signal_start2, signal_end2, type="max")
                signal1 = float(signal1[0]) if signal1[0] is not None else 0
                signal2 = float(signal2[0]) if signal2[0] is not None else 0
                output_list.append([current_chrom, start1, end1, start2, end2, signal1, signal2])
    return output_list
                
def write_bed(output_list,output_bed):
    with open(output_bed, "w") as f:
        for line in output_list:
            f.write("\t".join([str(x) for x in line])+"\n")
                
"""
This script is used to calculate the CTCF enrichment in the loop regions. 
```
python3 loop_ctcf_enrichment.py [loop.bed] [ctcf_chip.bw] [resolution] [output.bed]
```
[loop.bed]: the chromatin loop coordinate, in bed format <br>
[ctcf_chip.bw]: the CTCF ChIP-seq signal in bigwig format  <br>
[resolution]: the resolution of the Hi-C data <br>
[output.bed]: the output bed file. Format: chr1 start1 end1 start2 end2 enrichment1 enrichment2 <br>


"""


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python3 loop_ctcf_enrichment.py [loop.bed] [ctcf_chip.bw] [resolution] [output.bed]")
        print("loop.bed: the chromatin loop coordinate, in bed format")
        print("ctcf_chip.bw: the CTCF ChIP-seq signal in bigwig format")
        print("resolution: the resolution of the Hi-C data")
        print("output.bed: the output bed file. Format: chr1 start1 end1 start2 end2 enrichment1 enrichment2")
        sys.exit(1)
    loop_bed = os.path.abspath(sys.argv[1])
    ctcf_chip_bw = os.path.abspath(sys.argv[2])
    resolution = int(sys.argv[3])
    output_bed = os.path.abspath(sys.argv[4])
    output_dir = os.path.dirname(output_bed)
    os.makedirs(output_dir, exist_ok=True)
    loop_dict = extract_loop_loc(loop_bed)
    enrichment_list = extract_enrichment(loop_dict, ctcf_chip_bw, resolution)
    write_bed(output_list=enrichment_list, output_bed=output_bed)
    print("The enrichment file is saved in", output_bed)