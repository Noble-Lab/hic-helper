import pickle
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
import numpy as np
from PIL import Image
from numpy import triu
from collections import defaultdict
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
def load_pickle(input_array_pickle):
    with open(input_array_pickle, 'rb') as f:
        data = pickle.load(f)
    return data
def load_input_array(input_path,input_chrom1,input_chrom2):
    data = load_pickle(input_path)
    input_key = find_key(data,input_chrom1,input_chrom2)
    input_array = data[input_key]
    print("loaded the corresponding matrix shape: ",input_array.shape)
    return input_array
def parse_bed_file(input_bed):
    loop_dict=defaultdict(list)
    with open(input_bed, 'r') as f:
        
        for line in f:
            try:
                line = line.strip().split('\t')
                chrom1 = line[0]
                start1 = int(line[1])
                end1 = int(line[2])
                chrom2 = line[3]
                start2 = int(line[4])
                end2 = int(line[5])
                loop_dict[chrom1].append((start1,end1,start2,end2))
            except:
                print("skip line: ", line)
    return loop_dict
def locate_array_region(matrix,start_index1,end_index1,start_index2,end_index2):
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
    return output_data

def annotate_loop(input_array,loop_list,output_png,start_index1,end_index1,
                   start_index2,end_index2,max_value,run_mode):
    #locate the array region
    # if array is sparse version
    if isinstance(input_array,coo_matrix):
        input_array = locate_array_region(input_array,start_index1,end_index1,start_index2,end_index2)
    # if array is dense numpy version
    elif isinstance(input_array,np.ndarray):
        input_array = input_array[start_index1:end_index1,start_index2:end_index2]
    else:
        print("Current input array type: ",type(input_array))
        raise ValueError("The input array is not supported.")
    print("select matrix stat, min:",np.min(output_data),"max:",np.max(output_data),"mean:",np.mean(output_data))
    #output_data = output_data + triu(output_data,1).T #limit its application to rectangular matrix
    if run_mode == 1:
        output_data = np.log(output_data+1)
        max_value = np.log(max_value+1)
    output_data = convert_rgb(output_data,max_value)
    #image=np.array(output_data,dtype=np.uint8)
    #read the output data to directly set the loop region to blue
    # check loop list and filter the remained loops
    count_mark_loop=0
    for item in loop_list:
        start1,end1,start2,end2 = item
        #judge if it is in the current region
        if start1 >= end_index1 or end1 <= start_index1 or start2 >= end_index2 or end2 <= start_index2:
            continue
        revise_start1 = start1 - start_index1
        revise_end1 = end1 - start_index1
        revise_start2 = start2 - start_index2
        revise_end2 = end2 - start_index2
        revise_start1 = max(0,revise_start1)
        revise_end1 = min(end_index1-start_index1,revise_end1)
        revise_start2 = max(0,revise_start2)
        revise_end2 = min(end_index2-start_index2,revise_end2)
        #set the loop region to blue
        output_data[revise_start1:revise_end1,revise_start2:revise_end2,0] = 0
        output_data[revise_start1:revise_end1,revise_start2:revise_end2,1] = 0
        output_data[revise_start1:revise_end1,revise_start2:revise_end2,2] = 255
        count_mark_loop+=1
    print("total marked loops in the region: ",count_mark_loop)
    img = Image.fromarray(output_data, 'RGB')
    expect_min_size=224
    if img.size[0] < expect_min_size or img.size[1] < expect_min_size:
        ratio=max(expect_min_size/img.size[0],expect_min_size/img.size[1])
        img = img.resize((int(img.size[0]*ratio),int(img.size[1]*ratio)),Image.BICUBIC)
    img.save(output_png)
    print('The image has been saved to', output_png + '.')
    return output_png
def update_loop_to_mask(loop_list,resolution):
    update_list=[]
    for item in loop_list:
        start1,end1,start2,end2 = item
        #make sure the start1 is smaller than start2, otherwise swap
        if start1 > start2:
            start1,start2 = start2,start1
            end1,end2 = end2,end1
        difference_margin= end1-start1
        difference_margin = difference_margin//resolution
        difference_margin = max(1,difference_margin)
        difference_margin = difference_margin
        start1 = start1//resolution
        end1 = end1//resolution
        start2 = start2//resolution
        end2 = end2//resolution
        
        start1 = start1 - difference_margin
        end1 = end1 + difference_margin
        start2 = start2 - difference_margin
        end2 = end2 + difference_margin
        update_list.append((start1,end1,start2,end2))
    return update_list
"""
This script is to annotate and visualize the input array with loop information.
```
python3 annotate_array_loop.py [input.pkl] [loop.bed] [output.png] [chrom1] [start_index1] [end_index1] [chrom2] [start_index2] [end_index2] [resolution] [max_value] [mode]
```
[input.pkl] is the path to the pickle file containing the Hi-C array. <br>
[input.pkl] format: [chrom1_chrom2]:[array] format for common mode. Here array should be scipy sparce array. <br>
For intra-chromsome only, the dict format can be [chrom]:[array] in pickle files. <br>
[loop.bed] is the path to the bed file containing the loop information. <br>
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
"""




if __name__ == '__main__':
    import os
    import sys

    #take the overall array as input
    if len(sys.argv) != 13:
        print('Usage: python3 array2png.py [input.pkl] [loop.bed] [output.png] \
              [chrom1] [start_index1] [end_index1] \
              [chrom2] [start_index2] [end_index2] [resolution] [max_value] [mode]')
        print("This is the full array2png script. ")
        print("input.pkl: the path to the pickle file containing the array [String].")
        print("input.pkl format: [chrom1_chrom2]:[array] format for common mode. Here array should be scipy sparce array."\
              "For intra-chromsome only, the dict format can be [chrom]:[array] in pickle files.")
        print("loop.bed: the path to the bed file containing the loop information [String].")
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
    loop_file = os.path.abspath(sys.argv[2])
    output_png = os.path.abspath(sys.argv[3])
    output_dir = os.path.dirname(output_png)
    os.makedirs(output_dir, exist_ok=True)
    chrom1 = sys.argv[4]
    start_index1 = int(sys.argv[5])
    end_index1 = int(sys.argv[6])
    chrom2 = sys.argv[7]
    start_index2 = int(sys.argv[8])
    end_index2 = int(sys.argv[9])
    resolution = int(sys.argv[10])
    max_value = float(sys.argv[11])
    mode = int(sys.argv[12])
    start_index1 = start_index1//resolution
    end_index1 = end_index1//resolution
    start_index2 = start_index2//resolution
    end_index2 = end_index2//resolution
    #first parse loop bed to get the loop regions
    loop_dict = parse_bed_file(loop_file)
    loop_key = find_key(loop_dict,chrom1,chrom2)
    loop_list = loop_dict[loop_key]
    print("current focus chromosome loop: ", len(loop_list))
    loop_list= update_loop_to_mask(loop_list,resolution)
    #load the input array
    input_array = load_input_array(input_array_pickle,chrom1,chrom2)
    annotate_loop(input_array,loop_list,output_png,start_index1,end_index1,
                   start_index2,end_index2,max_value,mode)
