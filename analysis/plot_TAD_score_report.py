
import os
import sys 
import matplotlib.pyplot as plt
from collections import defaultdict
import pickle
import pandas as pd
import seaborn as sns 
def parse_boundary(input_bed):
    """
    input_bed: the input bed file.
    """
    boundary_dict = defaultdict(list)
    with open(input_bed, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line)==0:
                continue
            arr = line.split('\t')
            chrom = arr[0]
            start = int(arr[1])
            boundary_dict[chrom].append(start)
    return boundary_dict

def plot_TAD_score_report(input_bed, input_pkl, output_pdf, window_size):
    """
    input_bed: the input bed file.
    input_pkl: the input pickle file containing the insulation score.
    output_pdf: the output pdf/png file.
    window_size: the window size for the plot.
    """
    input_dict = parse_boundary(input_bed)
    insulation_score_dict = pickle.load(open(input_pkl, 'rb'))
    average_dict = defaultdict(list) #key: distance relative to boundary, value: corresponding insulation score
    for chrom in input_dict:
        #ignore chrX,chrY,chrM, chrUn and random chromosomes
        if "chrX" in chrom or "chrY" in chrom or "chrM" in chrom or "chrUn" in chrom or "random" in chrom:
            continue
        cur_bound_list = input_dict[chrom]
        cur_insulation_score = insulation_score_dict[chrom]
        for bound in cur_bound_list:
            for k in range(-window_size, window_size+1):
                if bound+k>=0 and bound+k<len(cur_insulation_score):
                    average_dict[k].append(cur_insulation_score[bound+k])
    final_df = {"distance": [], "insulation_score": []}
    for k in average_dict:
        cur_dist= k
        cur_score_list = average_dict[k]
        final_df["distance"].extend([cur_dist]*len(cur_score_list))
        final_df["insulation_score"].extend(cur_score_list)
    final_df = pd.DataFrame(final_df)

    #plot the figure
    fig, ax = plt.subplots(figsize=(5,5))
    sns.lineplot(data=final_df, x="distance", y="insulation_score", orient="x",alpha=0.5)
    plt.xlabel("Distance to boundary",fontsize=24)
    plt.ylabel("Insulation score",fontsize=24)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.tight_layout()
    plt.savefig(output_pdf,dpi=600)







if __name__ == '__main__':
    if len(sys.argv)!=5:
        print("Usage: python3 plot_TAD_score_report.py [input.bed] [insulation_score.pkl] [output.pdf] [window_size]")
        print("This script is used to plot the TAD score report.")
        print("[input.bed]: the input bed file.")
        print("[insulation_score.pkl]: the input pickle file containing the insulation score.")
        print("[output.pdf]: the output pdf/png file.")
        print("[window_size]: the window size for the plot.")
        sys.exit(1)
    input_bed = os.path.abspath(sys.argv[1])
    input_pkl = os.path.abspath(sys.argv[2])
    output_pdf = os.path.abspath(sys.argv[3])
    window_size = int(sys.argv[4])
    plot_TAD_score_report(input_bed, input_pkl, output_pdf, window_size)