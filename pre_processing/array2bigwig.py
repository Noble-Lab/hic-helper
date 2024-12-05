import numpy as np
import os
import shutil
from collections import defaultdict
import pyBigWig
import sys
import pickle 
# def parse_genome_size(genome_size_file):
#     chrom_size = {}
#     with open(genome_size_file) as f:
#         for line in f:
#             chrom, size = line.strip("\n").split()
#             chrom_size[chrom] = int(size)
#     return chrom_size

def array2bigwig(input_file,output_bigwig,resolution=1000):
    """
    convert the 1D array to bigwig file
    """

    #read genome info
    #chrom_size = parse_genome_size(genome_path)
    
    # read the input file
    # if input_file.endswith(".npy"):
    #     data = np.load(input_file)
    # elif input_file.endswith(".txt") or input_file.endswith(".Ev1") or input_file.endswith(".Ev2"):
    #     data = np.loadtxt(input_file)
    # else:
    #     print("The input file should be in npy or txt format")
    #     sys.exit(1)
    with open(input_file, 'rb') as f:
        data = pickle.load(f)
    #generate chromosomes size dict
    chrom_size={}
    for chrom in data:
        cur_list=data[chrom]
        chrom_size[chrom]=len(cur_list)*resolution
    
    with pyBigWig.open(output_bigwig, "w") as bw:
        chromosomes = [(key,chrom_size[key]) for key in chrom_size]
        print(chromosomes)
        # Add chromosome information to the BigWig file
        bw.addHeader(chromosomes)
        for info in chromosomes:
            chrom = info[0]
            print("start chrom",chrom,"size",chrom_size[chrom])
            #pos_embedding = np.arange(0, len(data), dtype=int)
            start_embedding = np.arange(0, chrom_size[chrom], resolution)
            start_embedding = start_embedding.astype(int)
            end_embedding= start_embedding+resolution
            end_embedding = np.clip(end_embedding, 0, chrom_size[chrom]-1)
            # bw.addEntries("chr1", [500, 600, 635], values=[-2.0, 150.0, 25.0], span=20)
            chrom_data = data[chrom]
            chrom_data= np.nan_to_num(chrom_data)
            start_embedding = [int(x) for x in start_embedding]
            end_embedding = [int(x) for x in end_embedding]
            chrom_data= [float(x) for x in chrom_data]
            print(set([type(x) for x in chrom_data]))
            print(len(chrom_data),len(start_embedding),len(end_embedding))
            assert all(isinstance(c, str) for c in chrom), "Chromosomes must be strings"
            assert all(isinstance(s, int) for s in start_embedding), "Start positions must be integers"
            assert all(isinstance(e, int) for e in end_embedding), "End positions must be integers"
            assert all(isinstance(v, (int, float)) for v in chrom_data), "Values must be int or float"
            assert all(s >= 0 for s in start_embedding), "Start positions must be non-negative"
            assert all(e >= s for s, e in zip(start_embedding, end_embedding)), "End must be >= Start"
            for s, e, v in zip(start_embedding[:10], end_embedding[:10], chrom_data[:10]):
                print(f"Chromosome: {chrom}, Start: {s}, End: {e}, Value: {v}")
            bw.addEntries([chrom]*len(chrom_data), list(start_embedding), ends=list(end_embedding), values=list(chrom_data))
            print("finished",chrom)
    return output_bigwig

"""
This script is used to merge bigwig files into one bigwig file.
```
python3 array2bigwig.py [input_file] [output_bigwig] [resolution]
```
[input_file]: the pkl file to be converted, should be a dict in format of [chr]:[array]. <br>
[output_bigwig]: the output bigwig file. <br>
[resolution]: the resolution stored in the pkl file. <br>

"""


if __name__ == '__main__':
    
    if len(sys.argv) != 4:
        print("Usage: python3 array2bigwig.py [input_file] [output_bigwig] [resolution]")
        print("input_file: the pkl file to be converted, should be a dict in format of [chr]:[array].")
        print("output_bigwig: the output bigwig file.")
        print("resolution: the resolution stored in the pkl file.")
        sys.exit(1)
    input_file = os.path.abspath(sys.argv[1])
    output_bigwig = os.path.abspath(sys.argv[2])
    output_dir = os.path.basename(output_bigwig)
    os.makedirs(output_dir, exist_ok=True)
    resolution = int(sys.argv[3])
    output_bw = array2bigwig(input_file, output_bigwig, resolution)

    print("Finished converting saved to %s" % output_bw)
