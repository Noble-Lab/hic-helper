import sys
import os 
import sys
import os
import numpy as np
import pyfaidx
import pickle
from collections import defaultdict
import numpy as np


def parse_bed(input_bed):
    input_dict = defaultdict(list)
    with open(input_bed, 'r') as f:
        for line in f:
            line = line.strip().split()
            try:
                chrom = line[0]
                start = int(line[1])
                end = int(line[2])
                input_dict[chrom].append([start, end])
            except:
                print("skip line:", line)
                print("should be format: chrom start end")
    return input_dict

def count_peak_num(input_dict):
    peak_num = 0
    for chrom in input_dict.keys():
        peak_num += len(input_dict[chrom])
    return peak_num

def revise_input_dict(input_dict, window_region):
    revised_dict = defaultdict(list)
    half_window = window_region // 2
    for chrom in input_dict.keys():
        for peak in input_dict[chrom]:
            start = peak[0]
            end = peak[1]
            if window_region == 0:
                revised_dict[chrom].append([start, end])
            else:
                center = (start + end) // 2
                new_start = max(0, center - half_window)
                new_end = center + half_window
                revised_dict[chrom].append([new_start, new_end])
    return revised_dict
def extract_peak_sequence(input_dict, genome_fasta, output_fasta):
    output_seqlist = []
    genome = pyfaidx.Fasta(genome_fasta)
    all_chrom_keys = list(genome.keys())
    for chrom in input_dict.keys():
        if chrom in all_chrom_keys:
            chrom_key= chrom
        elif "chr"+chrom in all_chrom_keys:
            chrom_key = "chr"+chrom
        elif chrom.replace("chr", "") in all_chrom_keys:
            chrom_key = chrom.replace("chr", "")
        else:
            print("Chromosome not found in genome fasta file:", chrom)
            continue
        chrom_seq = genome[chrom]
        for peak in input_dict[chrom]:
            start = peak[0]
            end = peak[1]
            seq = chrom_seq[start:end].seq
            peak_seq = ">%s:%d-%d(+)\n"%(chrom_key,start,end)+str(seq)
            output_seqlist.append(peak_seq)
    with open(output_fasta, 'w') as f:
        for seq in output_seqlist:
            f.write(seq + '\n')
"""
This script is used to extract peak sequences from genome fasta file.
```
python3 extract_peak_sequence.py [peak_bed] [genome_fasta] [output_fasta] [window_region]
```
[peak_bed]: the bed file containing peak regions. <br>
[genome_fasta]: the genome fasta file. <br>
[output_fasta]: the output fasta file containing peak sequences. <br>
[window_region]: the window region to extract sequence around peak regions. <br>
If set to 0, then use the peak region itself. <br>
"""
if __name__ == '__main__':
    if len(sys.argv)!=5:
        print("Usage: python extract_peak_sequence.py [peak_bed] [genome_fasta] [output_fasta] [window_region]")
        print("peak_bed: the bed file containing peak regions")
        print("genome_fasta: the genome fasta file")
        print("output_fasta: the output fasta file containing peak sequences")
        print("window_region: the window region to extract sequence around peak regions. \n \
              If set to 0, then use the peak region itself")
        sys.exit(1)
    peak_bed = os.path.abspath(sys.argv[1])
    genome_fasta = os.path.abspath(sys.argv[2])
    output_fasta = os.path.abspath(sys.argv[3])
    window_region = int(sys.argv[4])
    output_dir = os.path.dirname(output_fasta)
    os.makedirs(output_dir, exist_ok=True)
    input_dict = parse_bed(peak_bed)
    print("Finished parsing peak bed file! Total number of peaks:", count_peak_num(input_dict))
    #revise the input_dict to include window region
    input_dict = revise_input_dict(input_dict, window_region)
    print("Finished revising peak regions!")
    #extract peak sequences
    extract_peak_sequence(input_dict, genome_fasta, output_fasta)
    print("Finished extracting peak sequences! Output fasta file:", output_fasta)





