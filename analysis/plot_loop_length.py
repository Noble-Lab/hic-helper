
import os
import sys 
import matplotlib.pyplot as plt
from collections import defaultdict
import pandas as pd
import seaborn as sns 
def parse_input_bed(input_bed):
    """
    input_bed: the input bed file.
    """
    loop_dict = defaultdict(list)
    with open(input_bed, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line)==0:
                continue
            if line.startswith("#"):
                continue
            if "start" in line or "chrom" in line:
                continue
            arr = line.split('\t')
            chrom1 = arr[0]
            start1 = int(arr[1])
            end1 = int(arr[2])
            chrom2 = arr[3]
            start2 = int(arr[4])
            end2 = int(arr[5])
            if chrom1!=chrom2:
                print("Warning: the input bed file contains inter-chromosomal loops, only intra-chromosomal loops are considered.") 
                continue
            loop_dict[chrom1].append((start1,start2))
    return loop_dict

def plot_loop_length(input_bed, output_pdf):
    #parse the input bed file
    input_dict = parse_input_bed(input_bed)
    #get all the distance
    final_df = {"distance": [],"Count": []}
    for chrom in input_dict:
        cur_loop_list = input_dict[chrom]
        for loop in cur_loop_list:
            start1 = loop[0]
            start2 = loop[1]
            distance = abs(start1-start2)
            final_df["distance"].append(distance)
            final_df["Count"].append(1)
    final_df = pd.DataFrame(final_df)
    plt.figure(figsize=(5,5))
    #plot the figure
    sns.histplot(
        data=final_df, x="distance", 
        fill=True, 
       alpha=0.5, linewidth=0,
    )
    plt.xlabel("Loop Length",fontsize=24)
    plt.ylabel("Density",fontsize=24)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    #plt.yscale('log')
    plt.tight_layout()
    plt.savefig(output_pdf, dpi=600)
"""
This script is used to plot the loop strength report.
```
python3 plot_loop_strengt.py [input.bed] [output.pdf]
```
[input.bed]: the input bed file. <br>
[output.pdf]: the output pdf/png file. <br>
"""
if __name__ == '__main__':
    if len(sys.argv)!=3:
        print("Usage: python3 plot_loop_strengt.py [input.bed] [output.pdf]")
        print("This script is used to plot the loop strength report.")
        print("[input.bed]: the input bed file.")
        print("[output.pdf]: the output pdf/png file.")
        sys.exit(1)

    input_bed = os.path.abspath(sys.argv[1])
    output_pdf = os.path.abspath(sys.argv[2])
    output_dir = os.path.dirname(output_pdf)
    os.makedirs(output_dir, exist_ok=True)
    plot_loop_length(input_bed, output_pdf)