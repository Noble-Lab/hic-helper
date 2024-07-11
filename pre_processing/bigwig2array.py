import sys
import os
import pyBigWig
import pickle 
import numpy as np
def bigwig2array(input_bw, output_pkl, resolution):
    bw = pyBigWig.open(input_bw)
    chroms = bw.chroms()
    signal_dict = {}
    for chrom in chroms:
        print("Processing", chrom)
        chrom_size = chroms[chrom]
        print("Chrom size:", chrom_size)
        # signal = bw.values(chrom, 0, chrom_size, numpy=True)
        #each resolution interval should sum to get the overall signal
        cutoff_length = chrom_size // resolution * resolution
        # signal = signal[:cutoff_length]
        
        # signal = np.nan_to_num(signal)
        # if len(signal)==0:
        #     print("Empty signal for", chrom)
        #     continue
        # if resolution > 1:
        #     signal = signal.reshape(-1, resolution).mean(axis=1)
        # signal_dict[chrom] = signal
        value_list = []
        for i in list(range(0, cutoff_length, resolution)):
            value_list.append(bw.stats(chrom, i, i + resolution)[0])

        value_list = [0 if v is None else v for v in value_list]
        signal = np.array(value_list)

        signal_dict[chrom] = signal
        print("Finished procssing! Signal shape:", signal.shape)
        #output signal stats
        print("Signal stats: mean ",np.mean(signal), "std ", np.std(signal), "max ", np.max(signal), "min ", np.min(signal))
    bw.close()
    with open(output_pkl, 'wb') as f:
        pickle.dump(signal_dict, f)
"""
This script converts bigwig file to array format specified by resolution.
```
python3 bigwig2array.py [input_bw] [output_pkl] [resolution]
```
[input_bw]: the input bigwig file. <br>
[output_pkl]: the output pkl file with [chrom]:[signal] format. <br>
[resolution]: the resolution of the signal. <br>

"""


if __name__ == '__main__':
    if len(sys.argv)!=4:
        print("Usage: python bigwig2array.py [input_bw] [output_pkl] [resolution]")
        print("input_bw: the input bigwig file")
        print("output_pkl: the output pkl file with [chrom]:[signal] format")
        print("resolution: the resolution of the signal")
        sys.exit(1)
    input_bw = os.path.abspath(sys.argv[1])
    output_pkl = os.path.abspath(sys.argv[2])
    resolution = int(sys.argv[3])

    #convert bigwig to array
    bigwig2array(input_bw, output_pkl, resolution)
