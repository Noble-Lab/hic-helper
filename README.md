# HiC Helper

Welcome to the HiC Helper repository! This repository contains a collection of tools and scripts for processing HiC (High-throughput Chromosome Conformation Capture) data.

HiC Helper provides a set of utilities to preprocess, analyze, and visualize HiC data. Whether you are working with raw HiC data or processed matrices, this toolkit aims to simplify the analysis pipeline and provide useful functions for exploring chromatin interactions.

Please refer to the documentation and code examples in this repository to learn more about the available tools and how to use them effectively.

Happy HiC processing!

Copyright (C) 2024 Xiao Wang, Anupama Jha , Tangqi Fang and the University of Washington.

License: GPL v3. 

Contact: William Stafford Noble (wnoble@uw.edu) and Sheng Wang (swang@cs.washington.edu).

For technical problems or questions, please reach to Xiao Wang (wang3702@uw.edu) and Anupama Jha (anupamaj@uw.edu).

## Installation
Configure software dependency
```
conda env create -f environment.yml
```

Configure data dependency (optional)

```
bash set_up.sh
```

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
The input pickle should be in a pickle file as dict: [chrom1_chrom2]:[array] format for common mode. Here array should be scipy sparce array/numpy array. <br>
For intra-chromsome only, the dict format can be [chrom]:[array] in pickle files.<br>
[output.hic] is the name of the output hic file. <br>
[resolution] is used to specify the resolution that stored in the output array. <br>
[refer_genome_name] is used to specify the reference genome name. For example, "hg38","hg19","mm10" are valid inputs. <br>
[mode]: 0: all chromosome mode (scipy sparce array); 1: intra-chromosome mode(scipy sparce array); 2: all chromosome mode (numpy array); 3: intra-chromosome mode(numpy array). <br>


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

#### 5.hiccups_loop.py
[hiccups_loop.py](post_processing/hiccups_loop.py) <br>
Use HiCCUPs to detect loop from Hi-C input
```
python3 hiccups_loop.py [hicFile] [output_dir] [resolution]
```
[hicFile]: the path to the input hic file [String]. <br>
[output_dir]: the directory to the output loops [String]. <br>
[resolution]: the resolution of the input hic file [Integer]. <br>
Currently only support 5000,10000,25000 resolutions. <br>
The output loop bedpe file will be saved in [output_dir]/merged_loops.bedpe. <br>

#### 6. loop_f1.py
[loop_f1.py](analysis/loop_f1.py)
Compute F1 metrics of predicted loop and ground truth loop.
```
python3 loop_f1.py [true.bed] [pred.bed] [resolution]
```
[true.bed]: the true peaks, in bed format <br>
[pred.bed]: the predicted peaks, in bed format <br>
[resolution]: the resolution of the Hi-C data <br>

#### 7. fastq2bam.py
[fastq2bam.py](pre_processing/fastq2bam.py)

```
python3 fastq2bam.py [fastq_file1] [fastq_file2] [refer.fa] [output_dir]
```
[fastq_file1]: the first fastq file. <br>
[fastq_file2]: the second fastq file. <br>
[refer.fa]: the reference genome file. You can download the reference genome file from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/ for human.<br>
You can also run set_up.sh to download the reference genome files for human and mouse. <br>
[output_dir]: the output directory. The output file will be named as 4DN.sorted.bam under this direcotry. <br>

#### 8. fastq_4dn.py
This script is used to convert fastq files to cool or hic files following 4DN's pipeline.
```
python3 fastq_4dn.py [fastq_file1] [fastq_file2] [refer.fa] [chrom_size_file] [output_dir] [mode] [number_cpu] [max_memory] [resume_flag]
```
[fastq_file1]: the first fastq file. <br>
[fastq_file2]: the second fastq file. <br>
[refer.fa]: the reference genome file. You can download the reference genome file from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/ for human. <br>
You can also run set_up.sh to download the reference genome files for human and mouse. <br>
[chrom_size_file]: the chromosome size file. Example file: hg38.chrom.sizes for human genome build GRCh38 <br>
[output_dir]: the output directory. The output file will be named as 4DN.cool under this direcotry. <br>
[mode]: the mode of the conversion. 0: convert to cool file; 1：convert to hic file <br>
[number_cpu]: the number of cpu to use <br>
[max_memory]: the max memory to use (GB) <br>
[resume_flag]: 0: do not resume; 1: resume from the previously generated files. default should be 0. <br>
Recommended running with 8 cores and 64GB memory. <br>

##### Example Command:
After you run set_up.sh to download the reference genome files, you can run the following command to convert example fastq files to hic/cool files.<br>
1. Convert fastq files to cool files: <br>
```
python3 fastq_4dn.py ../reference_data/sample_data/GM12878_SRR1658581_1pc_1_R1.h10000.fastq.gz ../reference_data/sample_data/GM12878_SRR1658581_1pc_1_R2.h10000.fastq.gz ../reference_data/hg19.fa ../reference_data/hg19.chrom.sizes output_test 0 32 128 0
```
The output file will be named as 4DN.cool under the "output_test" direcotry. <br>

2. Convert fastq files to hic files: <br>
```
python3 fastq_4dn.py ../reference_data/sample_data/GM12878_SRR1658581_1pc_1_R1.h10000.fastq.gz ../reference_data/sample_data/GM12878_SRR1658581_1pc_1_R2.h10000.fastq.gz ../reference_data/hg19.fa ../reference_data/hg19.chrom.sizes output_test 1 32 128 0
```
The output file will be named as 4DN.hic under the "output_test" direcotry. <br>