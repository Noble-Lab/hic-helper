import pyBigWig
import numpy as np
import sys 
import os
import matplotlib.pyplot as plt
import seaborn as sns
def find_key(chromosome, chromosome_sizes):
    """
    This function is used to find the key of the chromosome in the chromosome sizes dictionary.
    chromosome: the chromosome name [String].
    chromosome_sizes: the chromosome sizes dictionary [Dictionary].
    """
    chromosome_choice1= chromosome
    if "chr" in chromosome:
        chromome_choice2= chromosome.replace('chr','')
    else:
        chromome_choice2= 'chr'+chromosome

    for key in chromosome_sizes.keys():
        if key == chromosome_choice1 or key == chromome_choice2:
            return key
    return None 
def plot_bigwig_signal(input_bigwig, output_png, chromosome, start, end):
     
    with  pyBigWig.open(input_bigwig) as bw:
        chromosome_sizes = bw.chroms()
        use_key = find_key(chromosome, chromosome_sizes)
        current_size = chromosome_sizes[use_key]
        start = max(0, start)
        end = min(current_size, end)
        signal = bw.values(use_key, start, end,numpy=True)
    plt.figure()

    select_range = np.arange(start, end)
    plt.plot(select_range,signal,alpha= 0.5,color="#a1c9f4")
    plt.axis('off')
    plt.fill_between(
        x= select_range, 
        y1= signal,
        color="#a1c9f4",
        alpha= 0.5)
    plt.tight_layout()
    plt.savefig(output_png, dpi=600)
"""
This script is used to plot the region of signal of a bigwig file.
```
python3 plot_bigwig_signal.py [input.bw] [output.png] [chromosome] [start] [end]
```
input.bw: the path to the bigwig file [String]. <br>
output.png: the name of the output png file [String]. <br>
chromosome: the chromosome name [String]. <br>
start: the start position of the region [Integer]. <br>
end: the end position of the region [Integer]. <br>

"""
if __name__ == '__main__':
    if len(sys.argv)!=6:
        print('Usage: python3 plot_bigwig_signal.py [input.bw] [output.png] [chromosome] [start] [end]')
        print('input.bw: the path to the bigwig file [String].')
        print('output.png: the name of the output png file [String].')
        print('chromosome: the chromosome name [String].')
        print('start: the start position of the region [Integer].')
        print('end: the end position of the region [Integer].')
        sys.exit(1)
    input_bigwig = os.path.abspath(sys.argv[1])
    output_png = os.path.abspath(sys.argv[2])
    output_dir = os.path.dirname(output_png)
    os.makedirs(output_dir, exist_ok=True)
    chromosome = sys.argv[3]
    start = int(sys.argv[4])
    end = int(sys.argv[5])
    plot_bigwig_signal(input_bigwig, output_png, chromosome, start, end)