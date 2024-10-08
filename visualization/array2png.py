import pickle
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
import numpy as np
from PIL import Image
from numpy import triu
def convert_rgb(data,max_value):
    """
    The convert_rgb function takes in a 2D array and converts it to a 3D RGB array.
    :param data: Specify the input 2D array
    :param max_value: Specify the maximum threshold of the input array for figures
    :return: A 3D RGB array of size (data.shape[0], data.shape[1], 3)
    """
    data_red = np.ones(data.shape)
    data = np.minimum(data,max_value)
    data = (max_value-data)/max_value
    data_rgb = np.stack([data_red,data,data],axis=0,dtype=np.float32)#transform only accept channel last case
    #data_rgb = data_rgb - imagenet_mean #put normalize to transform
    #data_rgb = data_rgb / imagenet_std
    data_rgb = data_rgb.transpose(1,2,0)
    data_rgb = (data_rgb*255).astype(np.uint8)
    return data_rgb
def find_key(data,input_chrom1,input_chrom2):
    if "chr" in input_chrom1:
        input_chrom1 = input_chrom1.replace("chr","")
    if "chr" in input_chrom2:
        input_chrom2 = input_chrom2.replace("chr","")
    input_key = input_chrom1 + "_" + input_chrom2
    if input_key in data.keys():
        return input_key
    input_key = "chr" + input_chrom1 + "_" + "chr" + input_chrom2
    if input_key in data.keys():
        return input_key
    input_key = "chr" + input_chrom2 + "_" + "chr" + input_chrom1
    if input_key in data.keys():
        return input_key
    if input_chrom1==input_chrom2:
        input_key = "chr" + input_chrom1
        if input_key in data.keys():
            return input_key
        input_key = input_chrom1
        if input_key in data.keys():
            return input_key
    else:
        print('The input chromosomes are not found in the input array.')
        print('The available chromosomes are:',data.keys())
        raise ValueError('The input chromosomes are not found in the input array.')
def array2png(input_array_pickle,output_png,input_chrom1,start_index1,
              end_index1,input_chrom2,start_index2,end_index2,max_value,run_mode):
    #load array
    """
    The array2png function takes in a pickled array, and outputs a png image.
    
    :param input_array_pickle: Specify the input pickle file
    :param output_png: Specify the output file name
    :param chrom1: Specify the chromosome of the first region
    :param start_index1: Specify the start index of the first chromosome
    :param end_index1: Specify the end of the first chromosome
    :param chrom2: Specify the second chromosome in inter-chromosome mode
    :param start_index2: Specify the start index of the second chromosome
    :param end_index2: Specify the end index of the second chromosome
    :param max_value: Specify the maximum threshold of the input array for figures
    :param run_mode: Specify the mode of the visualization. 0:raw visualization; 1: log visualization.
    :return: A figure of size (end_index2-start_index2) * (end_index2-start_index2)
    :doc-author: Trelent
    """
    with open(input_array_pickle, 'rb') as f:
        data = pickle.load(f)
    
    print("allow input data with chromosomes: ",data.keys())    
    output_dir = os.path.dirname(output_png)
    os.makedirs(output_dir, exist_ok=True)
    input_key = find_key(data,input_chrom1,input_chrom2)
    matrix = data[input_key]
    print('The matrix has been loaded.',input_key,'shape:',matrix.shape)
    matrix_row = np.concatenate([matrix.row,matrix.col])
    matrix_col = np.concatenate([matrix.col,matrix.row])
    matrix_data = np.concatenate([matrix.data,matrix.data])
    select_index1 = (matrix_row >= start_index1) & (matrix_row < end_index1)
    select_index2 = (matrix_col >= start_index2) & (matrix_col < end_index2)
    select_index = select_index1 & select_index2
    matrix_row = matrix_row[select_index]
    matrix_col = matrix_col[select_index]
    matrix_data = matrix_data[select_index]
    matrix_row = matrix_row - start_index1
    matrix_col = matrix_col - start_index2
    #gen coo_matrix with the selected data
    output_data = coo_matrix((matrix_data, (matrix_row, matrix_col)), shape=(end_index1-start_index1, end_index2-start_index2))
    output_data = output_data.toarray()
    #output_data = output_data + triu(output_data,1).T #limit its application to rectangular matrix
    if run_mode == 1:
        output_data = np.log(output_data+1)
        max_value = np.log(max_value+1)
    output_data = convert_rgb(output_data,max_value)
    #image=np.array(output_data,dtype=np.uint8)
    
    img = Image.fromarray(output_data, 'RGB')
    expect_min_size=224
    if img.size[0] < expect_min_size or img.size[1] < expect_min_size:
        ratio=max(expect_min_size/img.size[0],expect_min_size/img.size[1])
        img = img.resize((int(img.size[0]*ratio),int(img.size[1]*ratio)),Image.BICUBIC)
    
    img.save(output_png)
    print('The image has been saved to', output_png + '.')
    return output_png

"""
Usage of array2png.py
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
[mode]: 0:raw visualization; 1: log visualization. <br>
All index input should be absolute index counted by base. <br>

"""

if __name__ == '__main__':
    import os
    import sys
    #take the overall array as input
    if len(sys.argv) != 12:
        print('Usage: python3 array2png.py [input.pkl] [output.png] \
              [chrom1] [start_index1] [end_index1] \
              [chrom2] [start_index2] [end_index2] [resolution] [max_value] [mode]')
        print("This is the full array2png script. ")
        print("input.pkl: the path to the pickle file containing the array [String].")
        print("input.pkl format: [chrom1_chrom2]:[array] format for common mode. Here array should be scipy sparce array."\
              "For intra-chromsome only, the dict format can be [chrom]:[array] in pickle files.")
        print("output.png: the name of the output png file [String].")
        print("chrom1: the name of the first chromosome [String].")
        print("start_index1: the start index of the first chromosome [Integer].")
        print("end_index1: the end index of the first chromosome [Integer].")
        print("chrom2: the name of the second chromosome [String].")
        print("start_index2: the start index of the second chromosome [Integer].")
        print("end_index2: the end index of the second chromosome [Integer].")
        print("resolution: resolution of the input array [Integer].")
        print("max_value: the maximum threshold of the input array for figures[Float].")
        print("mode: 0:raw visualization; 1: log visualization. [Integer].")
        print("All index input should be absolute index counted by base.")
        sys.exit(1)
    input_array_pickle = os.path.abspath(sys.argv[1])
    output_png = os.path.abspath(sys.argv[2])
    chrom1 = sys.argv[3]
    start_index1 = int(sys.argv[4])
    end_index1 = int(sys.argv[5])
    chrom2 = sys.argv[6]
    start_index2 = int(sys.argv[7])
    end_index2 = int(sys.argv[8])
    resolution = int(sys.argv[9])
    max_value = float(sys.argv[10])
    mode = int(sys.argv[11])

    start_index1 = start_index1//resolution
    end_index1 = end_index1//resolution
    start_index2 = start_index2//resolution
    end_index2 = end_index2//resolution
    array2png(input_array_pickle,output_png,chrom1,start_index1,
              end_index1,chrom2,start_index2,end_index2,max_value,mode)
    
    