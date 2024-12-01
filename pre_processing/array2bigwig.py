import numpy as np
import os
import shutil
from collections import defaultdict
import json
import pyBigWig
import sys

def parse_genome_size(genome_size_file):
    chrom_size = {}
    with open(genome_size_file) as f:
        for line in f:
            chrom, size = line.strip("\n").split()
            chrom_size[chrom] = int(size)
    return chrom_size

def array2bigwig(genome_path,input_file,output_bigwig,chrom_list=None, resolution=1000000):
    """
    convert the 1D array to bigwig file
    """

    #read genome info
    chrom_size = parse_genome_size(genome_path)
    
    # read the input file
    if input_file.endswith(".npy"):
        data = np.load(input_file)
    elif input_file.endswith(".txt") or input_file.endswith(".Ev1") or input_file.endswith(".Ev2"):
        data = np.loadtxt(input_file)
    else:
        print("The input file should be in npy or txt format")
        sys.exit(1)
    data = np.nan_to_num(data)
    data = data.reshape(-1)

    
    with pyBigWig.open(output_bigwig, "w") as bw:
        if not chrom_list:
            chromosomes = [(key,chrom_size[key]) for key in chrom_size]
        else:
            chromosomes = [(key,chrom_size[key]) for key in chrom_size if key in chrom_list]
        print(chromosomes)
        # Add chromosome information to the BigWig file
        bw.addHeader(chromosomes)
        for info in chromosomes:
            chrom = info[0]
            print("start chrom",chrom,"size",chrom_size[chrom])
            pos_embedding = np.arange(0, len(data), dtype=int)
            start_embedding = np.arange(0, chrom_size[chrom], resolution)
            end_embedding = np.arange(resolution-1, chrom_size[chrom]+resolution, resolution)
            end_embedding = np.clip(end_embedding, 0, chrom_size[chrom]-1)
            # bw.addEntries("chr1", [500, 600, 635], values=[-2.0, 150.0, 25.0], span=20)
            bw.addEntries([chrom]*len(data), start_embedding, ends=end_embedding, values=list(data))
            # bw.addEntries(str(chrom), pos_embedding, values=list(data), span=resolution)
            print("finished",chrom)
    return output_bigwig

"""
This script is used to merge bigwig files into one bigwig file.
```
python3 array2bigwig.py [genome_path] [input_file] [output_bigwig] [chrom_list] [resolution]
```
[genome_path]: the path to the genome size file. <br>
[input_file]: the file to be converted, can be npy or txt or Ev. <br>
[output_bigwig]: the output bigwig file. <br>
[chrom_list]: the chromosome list to be included, if None, using full chromosome list. <br>
[resolution]: the span. <br>

"""


if __name__ == '__main__':
    
    if len(sys.argv) != 4:
        print("Usage: python3 array2bigwig.py [genome_path] [input_file] [output_bigwig]")
        print("genome_path: the path to the genome size file")
        print("input_file: the directory containing all the bigwig files")
        print("output_bigwig: the output bigwig file")
        sys.exit(1)
    genome_path = os.path.abspath(sys.argv[1])
    input_file = os.path.abspath(sys.argv[2])
    output_bigwig = os.path.abspath(sys.argv[3])
    chrom_list = [os.path.split(input_file)[1].split(".")[0]]
    output_dir = os.path.basename(output_bigwig)
    os.makedirs(output_dir, exist_ok=True)
    output_bw = array2bigwig(genome_path, input_file, output_bigwig, chrom_list)

    print("Finished converting saved to %s" % output_bw)