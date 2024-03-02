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
This is the intra-chromosome cool2array scripts, can be easily extended to full cool2array. <br>
The output array is saved in a pickle file as dict: [chrom]:[array] format. <br>
Two modes are supported: 
```
0: scipy coo_array format output; 
1: numpy array format output;
2: normed scipy coo_array format output; 
3: normed numpy array format output.
```


