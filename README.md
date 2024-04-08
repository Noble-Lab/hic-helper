# HiC Helper

Welcome to the HiC Helper repository! This repository contains a collection of tools and scripts for processing HiC (High-throughput Chromosome Conformation Capture) data.

HiC Helper provides a set of utilities to preprocess, analyze, and visualize HiC data. Whether you are working with raw HiC data or processed matrices, this toolkit aims to simplify the analysis pipeline and provide useful functions for exploring chromatin interactions.

Please refer to the documentation and code examples in this repository to learn more about the available tools and how to use them effectively.

Happy HiC processing!

Copyright (C) 2024 Xiao Wang, Anupama Jha , Tangqi Fang and the University of Washington.

License: GPL v3. 

Contact: William Stafford Noble (wnoble@uw.edu) and Sheng Wang (swang@cs.washington.edu).

For technical problems or questions, please reach to Xiao Wang (wang3702@uw.edu) and Anupama Jha (anupamaj@uw.edu).

## HiC Data Format Conversion
#### 1. cool2array.py
[cool2array.py](pre_processing/cool2array.py) <br> 
Usage
```
python3 cool2array.py [input.cool] [output.pkl] [mode]
```
This is the full cool2array script, converting both intra, inter chromosome regions to array format. <br>
The output array is saved in a pickle file as dict: [chrom1_chrom2]:[array] format. <br>
Two modes are supported: 
```
0: scipy coo_array format output; 
1: numpy array format output;
2: normed scipy coo_array format output; 
3: normed numpy array format output.
```
For different resolution, please simply input cool path as [xx.cool::resolutions/5000], here resolution specified is 5Kb. You can modify it to support different resolutions that you want to focus.

#### 2. hic2array.py
[hic2array.py](pre_processing/hic2array.py) <br>
Usage
```
python3 hic2array.py [input.hic] [output.pkl] [resolution] [normalization_type] [mode]
```

This is the full cool2array script, converting both intra, inter chromosome regions to array format. <br>
The output array is saved in a pickle file as dict: [chrom1_chrom2]:[array] format. <br>
[resolution] is used to specify the resolution that stored in the output array. <br>
[normalization_type] supports the following type: <br>
```
0: NONE normalization applied, save the raw data to array.
1: VC normalization; 
2: VC_SQRT normalization; 
3: KR normalization; 
4: SCALE normalization.
```
Two modes are supported for different format saving: 
```
0: scipy coo_array format output; 
1: numpy array format output.
```

#### 3. array2hic.py
[array2hic.py](pre_processing/array2hic.py) <br>

Usage
```
python3 array2hic.py [input.pkl] [output.hic] [resolution] [refer_genome_name] [mode]
```
The input pickle should be in a pickle file as dict: [chrom1_chrom2]:[array] format for common mode. Here array should be scipy sparce array. <br>
For intra-chromsome only, the dict format can be [chrom]:[array] in pickle files.<br>
[output.hic] is the name of the output hic file. <br>
[resolution] is used to specify the resolution that stored in the output array. <br>
[refer_genome_name] is used to specify the reference genome name. For example, "hg38","hg19","mm10" are valid inputs. <br>
[mode]: 0: all chromosome mode; 1: intra-chromosome mode. <br>


#### 4. array2png.py
[array2png.py](visualization/array2png.py) <br>

Usage
```
python3 array2png.py [input.pkl] [output.png] [chrom1] [start_index1] [end_index1] [chrom2] [start_index2] [end_index2] [resolution] [max_value] [mode]
```
This is the full array2png script. <br>
[input.pkl] is the path to the pickle file containing the array. <br>
[input.pkl] format: [chrom1_chrom2]:[array] format for common mode. Here array should be scipy sparce array. <br>
For intra-chromsome only, the dict format can be [chrom]:[array] in pickle files. <br>
[output.png] is the name of the output png file. <br>
[chrom1] is the name of the first chromosome. <br>
[start_index1] is the start index of the first chromosome. <br>
[end_index1] is the end index of the first chromosome. <br>
[chrom2] is the name of the second chromosome. <br>
[start_index2] is the start index of the second chromosome. <br>
[end_index2] is the end index of the second chromosome. <br>
[resolution] is the resolution of the input array. <br>
[max_value] is the maximum threshold of the input array for figures. <br>
[mode] is 0 for all chromosome mode and 1 for intra-chromosome mode. <br>
All index input should be absolute index counted by base. <br>
