import sys
import os 
import pickle 
import numpy as np
import pandas as pd
import cooler 
def array2sparse(array):
    """
    The array2sparse function converts a numpy array to a scipy sparce array.
    
    :param array: Specify the numpy array
    :return: A scipy sparce array
    :doc-author: Trelent
    """
    from scipy.sparse import coo_matrix
    row, col = np.where(array)
    data = array[row, col]
    return coo_matrix((data, (row, col)), shape=array.shape)
def array2cool(input_array_pickle,output_cool,resolution,refer_genome_name,mode):
    """
    The array2cool function converts a dict of numpy array to cool file.
    
    :param juicer_tools: Specify the location of the juicer_tools
    :param input_array_pickle: Specify the path to the pickle file containing the array
    :param output_cool: Specify the name of the output cool file
    :param resolution: Set the resolution of the hic file
    :param refer_genome_name: Specify the reference genome name
    :return: A hic file
    :doc-author: Trelent
    """
    #load array
    with open(input_array_pickle, 'rb') as f:
        data = pickle.load(f)
    output_dir = os.path.dirname(output_hic)
    os.makedirs(output_dir, exist_ok=True)

    #set each chromosome's length
    chromsizes={"name":[],"length":[]}
    chromosize_add_dict ={}#[chr_name:add_size]
    sort_keys = sorted(data.keys())
    accumulate_index = 0
    for chrom_name in sort_keys:
        
        if mode == 0 or mode == 2:
            chrom1, chrom2 = chrom_name.split('_')
        else:
            chrom1 = chrom_name
            chrom2 = chrom_name
        if chrom1!=chrom2:
            continue
        if "chr" not in chrom1:
            chrom1 = "chr"+chrom1
        cur_array = data[chrom_name]
        chrom_size_total=resolution*cur_array.shape[0]
        #chromsizes={"name":chrom_name,"length":[chrom_size_total]}
        chromsizes['name'].append(chrom1)
        chromsizes['length'].append(chrom_size_total)
        chromosize_add_dict[chrom1]=accumulate_index
        accumulate_index += cur_array.shape[0]
    print("collecting bin dict size",chromsizes)
    chrom_dict=pd.DataFrame.from_dict(chromsizes).set_index("name")['length']
    bins = cooler.binnify(chrom_dict, resolution)
    #then convert data array to index raw column, count array
    
    
    data_dict = {"bin1_id":[],"bin2_id":[],"count":[]}
    
    for key in data:
        chrom_name = key
        if mode == 0 or mode == 2:
            chrom1, chrom2 = chrom_name.split('_')
        else:
            chrom1 = chrom_name
            chrom2 = chrom_name
        if "chr" not in chrom1:
            chrom1 = "chr"+chrom1
        if "chr" not in chrom2:
            chrom2 = "chr"+chrom2
        print("processing",chrom1,chrom2,"...")
        matrix = data[key]
        if mode>=2:
            matrix = array2sparse(matrix)
        matrix_row = matrix.row
        matrix_col = matrix.col
        matrix_data = matrix.data

        matrix_row += chromosize_add_dict[chrom1]
        matrix_col += chromosize_add_dict[chrom2]
        data_dict['bin1_id']+=list(matrix_row)
        data_dict["bin2_id"]+=list(matrix_col)
        data_dict['count'] +=list(matrix_data)
        accumulate_index += matrix.shape[0]
    print("creating cool file...")
    #cooler.create_cooler(hic_path, bins,data_dict, dtypes={"count":"int"}, assembly="hg38")
    cooler.create_cooler(output_cool, bins=pd.DataFrame.from_dict(bins), pixels=pd.DataFrame.from_dict(data_dict), dtypes={'count': float},assembly=refer_genome_name)
"""
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
"""
if __name__ == '__main__':
    
    if len(sys.argv)!=6:
        print('Usage: python3 array2cool.py [input.pkl] [output.cool] [resolution] [refer_genome_name] [mode]')
        print("This is the full array2cool script. ")
        print("input.pkl: the path to the pickle file containing the array [String].")
        print("input.pkl format: [chrom1_chrom2]:[array] format for common mode. Here array should be scipy sparce array. For intra-chromsome only, the dict format can be [chrom]:[array] in pickle files.")
        print("output.cool: the name of the output cool file [String].")
        print("resolution: resolution of the input array [Integer].")
        print("refer_genome_name: the name of the reference genome [String]. Example: hg38, hg19, mm10.")
        print("mode: 0: all chromosome mode (scipy sparce array); 1: intra-chromosome mode(scipy sparce array); 2: all chromosome mode (numpy array); 3: intra-chromosome mode(numpy array).")
        sys.exit(1)
    
    script_dir = os.path.dirname(os.path.realpath(__file__))
    input_array_pickle = os.path.abspath(sys.argv[1])
    output_hic = os.path.abspath(sys.argv[2])
    resolution = int(sys.argv[3])
    refer_genome_name = str(sys.argv[4])
    mode = int(sys.argv[5])
    array2cool(input_array_pickle,output_hic,resolution,refer_genome_name,mode)