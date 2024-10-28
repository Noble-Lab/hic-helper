import numpy as np
from scipy.sparse import coo_matrix
import hicstraw
import os
def read_chrom_count(chr1,chr2, normalization,hic_file, resolution):
    """
    The read_chrom_count function takes in a chromosome, normalization method, hic file and resolution.
    It then returns the total number of reads for that chromosome as well as the number of non-diagonal reads.
    
    :param chr1: Specify the chromosome of interest
    :param chr2: Specify the chromosome to be compared with chr_name
    :param normalization: Normalize the hic matrix
    :param hic_file: Specify the hi-c file
    :param resolution: Determine the resolution of the hi-c matrix
    :return: The total number of reads in the chromosome
    :doc-author: Trelent
    """
    chr1_name = chr1.name
    chr2_name = chr2.name
    infos = []
    infos.append('observed')
    infos.append(normalization)
    infos.append(hic_file)
    infos.append(chr1_name)
    infos.append(chr2_name)
    infos.append('BP')
    infos.append(resolution)
    print(infos)
    val = []
    nondiag_val =[]
    rets = hicstraw.straw(*infos)
    print('\tlen(rets): {:3e}'.format(len(rets)))
    for ret in rets:
        cur_row=(int)(ret.binX // resolution)
        cur_col = (int)(ret.binY // resolution)
        val.append(ret.counts)
        if cur_row!=cur_col:
            nondiag_val.append(ret.counts)
    total_read = np.sum(val)
    nondiag_read = np.sum(nondiag_val)
    return total_read,nondiag_read
def count_cistotal(input_dict):
    """
    The count_cistotal function takes a dictionary as input and returns the total number of interactions in the dictionary,
    and the number of cis interactions. The function does this by splitting each key into two parts (split_key 1 and split_key 2)
    and then checking if they are equal. If they are equal, it adds one to cis count for every interaction that is found between 
    the same chromosome.
    
    :param input_dict: Pass in the dictionary that will be used to count cis and total interactions
    :return: A tuple of the total number of interactions and the cis interactions
    :doc-author: Trelent
    """
    cis_count=0
    total_count=0
    for key in input_dict:
        split_key1,split_key2=key.split("_")
        if split_key1==split_key2:
            cis_count+=input_dict[key]
        total_count+=input_dict[key]
    return total_count,cis_count
def count_hic_read(input_hic_path,resolution,normalization_type):
    """
    The count_hic_read function takes in a HiC file and returns the number of reads for each chromosome.
    
    :param input_hic_path: Specify the path of the hi-c file
    :param normalization_type: Determine how the data is normalized
    :return: A dictionary of chromosomes and their lengths
    :doc-author: Trelent
    """
    hic = hicstraw.HiCFile(input_hic_path)
    chrom_list=[]
    chrom_dict={}
    for chrom in hic.getChromosomes():
        print(chrom.name, chrom.length)
        if "all" in chrom.name.lower():
            continue
        chrom_list.append(chrom)
        chrom_dict[chrom.name]=chrom.length
    resolution_list = hic.getResolutions()
    resolution_list = list(resolution_list)
    #resolution = np.max(resolution_list)
    if resolution not in resolution_list:
        print("Resolution not found in the hic file, please choose from the following list:")
        print(resolution_list)
        exit()
    read_dict={}
    nondiag_read_dict={}
    for i in range(len(chrom_list)):
        for j in range(i,len(chrom_list)):
            chrom1 = chrom_list[i]
            chrom1_name = chrom_list[i].name
            chrom2 = chrom_list[j]
            chrom2_name = chrom_list[j].name
            cur_read,cur_nondiag_read = read_chrom_count(chrom1,chrom2, normalization_type, input_hic_path, resolution)
            if chrom1_name!=chrom2_name:
                cur_nondiag_read=0 #if they are not same chrom, diag is meaningless
            print(f"chrom-{chrom1_name} chrom-{chrom2_name} read {cur_read}, non-diag read {cur_nondiag_read}")
            read_dict["%s_%s"%(chrom1_name,chrom2_name)]=cur_read
            nondiag_read_dict["%s_%s"%(chrom1_name,chrom2_name)]=cur_nondiag_read
    print("chromosome-wise total read:")
    print(read_dict)
    print("chromosome-wise non-diagonal read:")
    print(nondiag_read_dict)

    #sum all the reads
    total_count,cis_count = count_cistotal(read_dict)
    print(f"Total Read {total_count}, Total Cis Read {cis_count}")

    #sum all non-diag reads

    total_count,cis_count = count_cistotal(nondiag_read_dict)
    print(f"Intra-chromo stat: Total Non-Diag Read {total_count}, Total Non-Diag Cis Read {cis_count}")
"""
This script is to count the total/total non-diag reads of cis/all. <br>
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
"""
if __name__ == '__main__':
    import os 
    import sys
    if len(sys.argv) != 4:
        print('Usage: python3 count_hic_read.py [input.hic] [resolution] [normalization_type]')
        print("This is the script to count the total reads, total non-diagonal reads from a input hic. ")
        print("[input.hic]: the input hic path.")
        print("[resolution]: specify the resolution to decide if the reads belong to diagonal.")
        print("[normalization_type]: 0: None normalization; 1: VC normalization; 2: VC_SQRT normalization; 3: KR normalization; 4: SCALE normalization")
        sys.exit(1)

    input_hic_path = os.path.abspath(sys.argv[1])
    resolution = int(sys.argv[2])
    normalization_type = int(sys.argv[3])
    normalization_dict={0:"NONE",1:"VC",2:"VC_SQRT",3:"KR",4:"SCALE"}
    if normalization_type not in normalization_dict:
        print('normalization type should be 0,1,2,3,4')
        print("normalization type: 0: None normalization; 1: VC normalization; 2: VC_SQRT normalization; 3: KR normalization; 4: SCALE normalization")
        sys.exit(1)
    normalization_type = normalization_dict[normalization_type]
    count_hic_read(input_hic_path,resolution,normalization_type)

    