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
This script is to convert .cool format to dict of arrays format.
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
This script is to convert .hic format to dict of arrays format.
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
Four modes are supported for different format saving: 
```
0: scipy coo_array format output; 
1: numpy array format output;
2: scipy csr_array format output (only include intra-chromsome region).
3: numpy array format output (only include intra-chromsome region).
```


#### 3. array2hic.py
[array2hic.py](pre_processing/array2hic.py) <br>
This script is used to convert dict of arrays format to .hic format.
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
This script is used to visualize Hi-C images from  format.
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
[loop_f1.py](analysis/loop_f1.py) <br>
Compute F1 metrics of predicted loop and ground truth loop.
```
python3 loop_f1.py [true.bed] [pred.bed] [resolution]
```
[true.bed]: the true peaks, in bed format <br>
[pred.bed]: the predicted peaks, in bed format <br>
[resolution]: the resolution of the Hi-C data <br>

#### 7. fastq2bam.py
[fastq2bam.py](pre_processing/fastq2bam.py) <br>
This script is used to convert fastq file to bam file.
```
python3 fastq2bam.py [fastq_file1] [fastq_file2] [refer.fa] [output_dir]
```
[fastq_file1]: the first fastq file. <br>
[fastq_file2]: the second fastq file. <br>
[refer.fa]: the reference genome file. You can download the reference genome file from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/ for human.<br>
You can also run set_up.sh to download the reference genome files for human and mouse. <br>
[output_dir]: the output directory. The output file will be named as 4DN.sorted.bam under this direcotry. <br>

#### 8. fastq_4dn.py
[fastq_4dn.py](pre_processing/fastq_4dn.py) 
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

#### 9. bam_4dn.py
[bam_4dn.py](pre_processing/bam_4dn.py) <br>
This script is used to convert bam file to cool or hic files following 4DN's pipeline.
```
python3 bam_4dn.py [input.bam] [chrom_size_file] [output_dir] [mode] [number_cpu] [max_memory] [resume_flag]
```
[input.bam]: the input bam file. <br>
[chrom_size_file]: the chromosome size file. Example file: hg38.chrom.sizes for human genome build GRCh38 <br>
[output_dir]: the output directory. The output file will be named as 4DN.cool under this direcotry. <br>
[mode]: the mode of the conversion. 0: convert to cool file; 1：convert to hic file <br>
[number_cpu]: the number of cpu to use <br>
[max_memory]: the max memory to use (GB) <br>
[resume_flag]: 0: do not resume; 1: resume from the previously generated files. default should be 0. <br>
Recommended running with 8 cores and 64GB memory. <br>

#### 10. pairs_4dn.py
[pairs_4dn.py](pre_processing/pairs_4dn.py) <br>
This script is used to convert pairs file to cool or hic files following 4DN's pipeline.
```
python3 pairs_4dn.py [input.pairs.gz] [chrom_size_file] [output_dir] [mode] [number_cpu] [max_memory] [resume_flag]
```
[input.pairs]: the input pairs.gz file <br>.
[chrom_size_file]: the chromosome size file. Example file: hg38.chrom.sizes for human genome build GRCh38 <br>
[output_dir]: the output directory. The output file will be named as 4DN.cool under this direcotry. <br>
[mode]: the mode of the conversion. 0: convert to cool file; 1：convert to hic file <br>
[number_cpu]: the number of cpu to use <br>
[max_memory]: the max memory to use (GB) <br>
[resume_flag]: 0: do not resume; 1: resume from the previously generated files. default should be 0. <br>
Recommended running with 8 cores and 64GB memory. <br>

#### 11. run-fastqc.sh
[run-fastqc.sh](pre_processing/run-fastqc.sh)  <br>
This script is used for fastq file's quality analysis.
```
./run-fastqc.sh [input_fastq] [num_threads] [output_dir]
```
[input_fastq]: an input fastq file, either gzipped or not.<br>
[num_threads]: number of threads to use.<br>
[output_dir] : output directory, will be created automatically if not exists.<br>

#### 12. pairs_qc.py
[pairs_qc.py](pre_processing/pairs_qc.py)  <br>
This script is used to perform quality control on the pairs file.<br>
**You must sussessfully run ./set_up.sh before you can use this script.** <br>
```
python3 pairs_qc.py [input.pairs.gz] [chrom_size_file] [output_dir] [enzyme]
```
- input.pairs.gz: input pairs file. <br>
- chrom_size_file: chrom size file. <br>
- output_dir: output directory. <br>
- enzyme: enzyme used for Hi-C experiment, either 4 or 6. <br>

#### 13. bamqc
[bamqc](bin/bamqc)  <br>
This script is used to perform quality control on the bam files.<br>
**You must sussessfully run ./set_up.sh before you can use this script.** <br>
```
./bin/bamqc --outdir=[output_dir] --noextract -t 8 [bam_file/bam_file_dir]
```
[output_dir] : output directory, you should create it before running this. <br>
[bam_file/bam_file_dir]: the file or directory that includes bam files. <br>


#### 14. array2cool
[array2cool.py](pre_processing/array2cool.py)  <br>
This script is used to convert dict of arrays format to .cool format.<br>
Usage
```
python3 array2cool.py [input.pkl] [output.cool] [resolution] [refer_genome_name] [mode]
```
The input pickle should be in a pickle file as dict: [chrom1_chrom2]:[array] format for common mode. Here array should be scipy sparce array. <br>
For intra-chromsome only, the dict format can be [chrom]:[array] in pickle files.<br>
[output.cool] is the name of the output cool file. <br>
[resolution] is used to specify the resolution that stored in the output array. <br>
[refer_genome_name] is used to specify the reference genome name. For example, "hg38","hg19","mm10" are valid inputs. <br>
[mode]:  0: all chromosome mode (scipy sparce array); 1: intra-chromosome mode(scipy sparce array); 2: all chromosome mode (numpy array); 3: intra-chromosome mode(numpy array). <br>

#### 15. bam_align_quality.py
[bam_align_quality.py](analysis/bam_align_quality.py)  <br>
This script calculates the alignment quality of a given bam file.
```
python3 bam_align_quality.py [input.bam] [output_dir] [number_cpu] [mode]
```
[input.bam]: the input bam file. <br>
[output_dir]: the output directory. <br>
[number_cpu]: the number of cpu used. <br>
[mode]: 0 for unsorted bam file, 1 for sorted bam file. <br>
The output includes stats of Unmapped, Low quality (mapq), Singleton, Multimapped, Duplicate, Other, Unique, and Total.


#### 16. pairs_quality.py
[pairs_quality.py](analysis/pairs_quality.py) <br>
This script is used to analyze the quality of the pairs file.
```
python3 pairs_quality.py [input.pairs.gz] [output_dir] [number_cpu]
```
[input.pairs.gz]: the input pairs gz file from 4DN pipeline (*.marked.sam.pairs.gz).<br>
[output_dir]: the output directory<br>
[number_cpu]: the number of cpu used<br>
The script will generate a report file in the output directory.<br>
The report file will contain the following information:<br>
1. Unmapped sequences: the number of unmapped sequences and the percentage.<br>
2. Singleton sequences: the number of singleton sequences and the percentage.<br>
3. Multimapped sequences: the number of multimapped sequences and the percentage.<br>
4. Duplicate sequences: the number of duplicate sequences and the percentage.<br>
5. Unique sequences: the number of unique sequences and the percentage.<br>
6. Total sequences: the total number of sequences.<br>
7. Detailed Unique Sequences Information:<br>
    RU sequences: the number of RU sequences and the percentage.<br>
    UR sequences: the number of UR sequences and the percentage.<br>
    UU sequences: the number of UU sequences and the percentage.<br>
    For detailed definition of RU/UR/UU, please see https://pairtools.readthedocs.io/en/latest/formats.html#pair-types<br>

#### 17. addnorm2hic.py
[addnorm2hic.py](pre_processing/addnorm2hic.py) <br>
This script is used to add norm to the hic files.
```
 python3 addnorm2hic.py [input.hic] [resolution] [num_cpu] [memory]
```
[input.hic]: input hic path. <br>
[resolution]: minimum resolution that normalization works. <br>
[num_cpu]: number of cpus used to normalize. <br>
[memory]: maximum memory (GB) that this script is allowed to use. <br>
The calculated norm vectors will be automatically saved in the input.hic file. <br>

#### 18. count_hic_read.py
[count_hic_read.py](analysis/count_hic_read.py) <br>
This script is to count the total/total non-diag reads of cis/all.
```
python3 count_hic_read.py [input.hic] [resolution] [normalization_type] 
```
[input.hic]: input hic path. <br>
[resolution] is used to specify the resolution that stored in the output array. <br>
[normalization_type]: should be an integer 0-4, corresponds the following type: <br>
```
0: NONE normalization applied, save the raw data to array.
1: VC normalization; 
2: VC_SQRT normalization; 
3: KR normalization; 
4: SCALE normalization.
```

#### 19. extract_hicnorms.py
[extract_hicnorms.py](pre_processing/extract_hicnorms.py) <br>
This script is to extract the normalization vectors from a .hic file. <br>
```
python3 extract_hicnorms.py [input.hic] [resolution] [normalization_type] [output_pkl]
```
[input.hic]: input hic path. <br>
[resolution]: resolution to extract the normalization vector, [Integer]. <br>
[normalization_type]: should be one of the following: NONE, VC, VC_SQRT, KR, SCALE, [string]. <br>
[output_pkl]: output pickle file path. <br> 
The normalization vector is saved in dict format, where the key is the chromosome name and the value is the normalization vector.

#### 20. loop_cleaner.py
[loop_cleaner.py](post_processing/loop_cleaner.py) <br>
This script is for filter out the loops on low mappability regions
```
python3 loop_cleaner.py [input.bed] [mappablility.bw] [output.bed] [threshold]
```
- input.bed: the input bed file <br>
- mappablility.bw: the mappablility bigwig file <br>
- output.bed: the output bed file <br>
- threshold: the mappablility threshold used to clean loops <br>

#### 21. merge_bigwig.py
[merge_bigwig.py](pre_processing/merge_bigwig.py) <br>
This script is used to merge bigwig files into one bigwig file.
```
python3 merge_bigwig.py [input_dir] [output_bw] [refer_genome.sizes]
```
[input_dir]: the directory containing all the bigwig files. <br>
[output_bw]: the output bigwig file. <br>
[refer_genome.sizes]: the chromosome sizes of the reference genome. <br>

#### 22. bigwig2array.py
[bigwig2array.py](pre_processing/bigwig2array.py) <br>
This script converts bigwig file to array format specified by resolution.
```
python3 bigwig2array.py [input_bw] [output_pkl] [resolution]
```
[input_bw]: the input bigwig file. <br>
[output_pkl]: the output pkl file with [chrom]:[signal] format. <br>
[resolution]: the resolution of the signal. <br>


#### 23. merge_bam.py
[merge_bam.py](pre_processing/merge_bam.py) <br>
This script is used to merge bam files into one bam file.
```
python3 merge_bam.py [input_dir] [output_bam]
```
[input_dir]: the directory containing all the bam files. <br>
[output_bam]: the output merged bam file. <br>

#### 24. bed_cleaner.py
[bed_cleaner.py](pre_processing/bed_cleaner.py) <br>
This script merges overlapping regions from bed file.
```
python3 bed_cleaner.py [input_bed] [output_bed]
```
[input_bed]: the input bed file. <br>
[output_bed]: the output bed file without overlapping regions. <br>


#### 25. pkl_contact_stat.py
[pkl_contact_stat.py](visualization/pkl_contact_stat.py) <br>
This script is used to plot contact frequency vs. genomic distance.
```
python3 pkl_contact_stat.py [input.pkl] [output.png] [genomic_dist]
```
[input.pkl]: the path to the pickle file containing the contact matrix [String]. <br>
[output.png]: the name of the output png file [String]. <br>
[genomic_dist]: the genomic distance for the plot [Integer]. <br>

#### 26. merge_hic.py
[merge_hic.py](pre_processing/merge_hic.py) <br>
This script is to merge two hic files to a new merged hic file with specified resolution.
```
python3 merge_hic.py [hic_file1] [hic_file2] [output_hic] [resolution] [refer_genome]
```
[hic_file1]: the first hic file to be merged. <br>
[hic_file2]: the second hic file to be merged. <br>
[output_hic]: the output merged hic file. <br>
[resolution]: the resolution of the output merged hic file. <br>
[refer_genome]: the name of the reference genome. For example, hg38, hg19, mm10. <br>

#### 26, merge_pkl.py
[merge_pkl.py](pre_processing/merge_pkl.py) <br>
This script is to merge two pkl files to a new merged pkl file.
```
python3 merge_pkl.py <pkl_file1> <pkl_file2> <output_pkl>
```
[pkl_file1]: the first pkl file to be merged. <br>
[pkl_file2]: the second pkl file to be merged. <br>
[output_pkl]: the output merged pkl file. <br>

#### 27. loop_apa.py
[loop_apa.py](visualization/loop_apa.py) <br>
This script is for plot the loop average peak analysis (APA) on the hic matrix.
```
python3 loop_apa.py [hic.pkl] [input.bed] [output.png] [resolution] [window_size]
```
- hic.pkl: the hic matrix file <br>
- input.bed: the input bed file including the loop regions <br>
- output.png: the output loop APA png file <br>
- resolution: the resolution of the hic matrix <br>
- window_size: the window size of the loop region <br>


#### 28. annotate_loop_gene.py
[annotate_loop_gene.py](post_processing/annotate_loop_gene.py) <br>
This script is for annotate the loop with gene information.
```
python3 annotate_loop_gene.py [input.bed] [gene_annotation] [output.bed]
```
input.bed: the input bed file that contains the loop information. <br>
gene_annotation: the gene annotation file format:.gtf, like hg38.ncbiRefSeq.gtf. <br>
output.bed: the output bed file that contains the annotated loop information. <br>
The last two columns in the output.bed file are the closest gene and the distance to the loop (corresponds to x and y).

#### 29. hic_coverage.py
[hic_coverage.py](analysis/hic_coverage.py) <br>
This script calculates the coverage of the Hi-C data.
```
python3 hic_coverage.py [input.pkl]
```
[input.pkl]: the input pkl file containing the Hi-C data <br>

#### 30. count_1d_read.py
[count_1d_read.py](analysis/count_1d_read.py)
This script is to calculate the average read count per base and total read count in the bigwig file.
```
python3 bigwig2count.py [input_bw]
```
[input_bw]: the input bigwig file. <br>


#### 31. peak_f1.py
[peak_f1.py](analysis/peak_f1.py)
```
python3 peak_f1.py [true.bed] [pred.bed] [max_dist]
```
[true.bed]: the true peaks, in bed format <br>
[pred.bed]: the predicted peaks, in bed format <br>
[max_dist]: the maximum distance to match the peaks <br>