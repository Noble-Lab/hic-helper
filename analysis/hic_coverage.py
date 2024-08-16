import sys
import os
import pickle
def calculate_coverage(input_pkl):
    data = pickle.load(open(input_pkl, 'rb'))
    count_reads=0
    count_length=0
    for chrom in data:
        cur_data = data[chrom]
        count_length += cur_data.shape[0]
        count_reads += cur_data.sum()
    count_length = count_length
    count_reads = count_reads
    coverage = count_reads/count_length
    return coverage
"""
This script calculates the coverage of the Hi-C data.
```
python3 hic_coverage.py [input.pkl]
```
[input.pkl]: the input pkl file containing the Hi-C data <br>

"""


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python3 hic_coverage.py [input.pkl]")
        print("[input.pkl]: the input pkl file containing the Hi-C data")
        # print("[fragment_size]: the size of the fragment for building Hi-C")
        sys.exit(1)
    input_pkl = os.path.abspath(sys.argv[1])
    #resolution = int(sys.argv[2])
    #fragment_size = int(sys.argv[3])
    coverage = calculate_coverage(input_pkl)
    print("Hi-C Coverage: ", coverage)