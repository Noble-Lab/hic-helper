import numpy as np
from scipy.sparse import coo_matrix
import hicstraw
import os
def read_chrom_count(chr1,chr2, normalization,hic_file, resolution):
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
        if cur_row==cur_col:
            nondiag_val.append(ret.counts)
    total_read = np.sum(val)
    nondiag_read = np.sum(nondiag_val)
    return total_read,nondiag_read
def count_cistotal(input_dict):
    cis_count=0
    total_count=0
    for key in input_dict:
        split_key1,split_key2=key.split("_")
        if split_key1==split_key2:
            cis_count+=input_dict[key]
        total_count+=input_dict[key]
    return total_count,cis_count
def count_hic_read(input_hic_path,normalization_type):
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
    print(f"Total Non-Diag Read {total_count}, Total Non-Diag Cis Read {cis_count}")
    
if __name__ == '__main__':
    import os 
    import sys
    if len(sys.argv) != 3:
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

    