

import numpy as np
from scipy.sparse import coo_matrix
import hicstraw
import os
import pickle 
def write_pkl(data, path):
    with open(path, 'wb') as f:
        pickle.dump(data, f)
def read_chrom_array(chr1, chr2, normalization, hic_file, resolution):
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
    row, col, val = [], [], []
    rets = hicstraw.straw(*infos)
    print('\tlen(rets): {:3e}'.format(len(rets)))
    for ret in rets:
        row.append((int)(ret.binX // resolution))
        col.append((int)(ret.binY // resolution))
        val.append(ret.counts)
    print('\tsum(val): {:3e}'.format(sum(val)))
    if chr1_name==chr2_name:
        max_shape =max(max(row),max(col))+1
        mat_coo = coo_matrix((val, (row, col)), shape = (max_shape,max_shape),dtype=np.float32)
    else:
        max_row = max(row)+1
        max_column = max(col)+1
        mat_coo = coo_matrix((val, (row, col)), shape = (max_row,max_column),dtype=np.float32)

    mat_coo = mat_coo #+ triu(mat_coo, 1).T #no below diagonaline records

    return mat_coo


def hic2array(input_hic,output_pkl=None,resolution=25000):
    """
    input_hic: str, input hic file path
    output_pkl: str, output pickle file path
    resolution: int, resolution of the hic file
    """

    hic = hicstraw.HiCFile(input_hic)
    chrom_list=[]
    chrom_dict={}
    for chrom in hic.getChromosomes():
        print(chrom.name, chrom.length)
        if "All" in chrom.name or "all" in chrom.name:
            continue
        chrom_list.append(chrom)
        chrom_dict[chrom.name]=chrom.length
    normalization="NONE"
    output_dict={}
    for i in range(len(chrom_list)):
        for j in range(i,len(chrom_list)):
            chrom1 = chrom_list[i]
            chrom1_name = chrom_list[i].name
            chrom2 = chrom_list[j]
            chrom2_name = chrom_list[j].name
            read_array=read_chrom_array(chrom1,chrom2, normalization, input_hic, resolution)
            output_dict[chrom1_name+"_"+chrom2_name]=read_array
    if output_pkl is not None:
        output_dir = os.path.dirname(os.path.realpath(output_pkl))
        os.makedirs(output_dir, exist_ok=True)
        write_pkl(output_dict,output_pkl)

    return output_dict
