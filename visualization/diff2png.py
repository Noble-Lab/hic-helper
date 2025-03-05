import pickle
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
import numpy as np
from PIL import Image

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

def diff2png(input_array_pickle,output_png,chrom1,start_index1,
              end_index1,chrom2,start_index2,end_index2,vmin,vmax):
    """
      The diff2png function takes in a pickled array, and outputs a png image.
      input_array_pickle: Specify the input pickle file
      output_png: Specify the output file name
      chrom1: Specify the chromosome of the first region
      start_index1: Specify the start index of the first chromosome
      end_index1: Specify the end of the first chromosome
      chrom2: Specify the second chromosome in inter-chromosome mode
      start_index2: Specify the start index of the second chromosome
      end_index2: Specify the end index of the second chromosome
      vmin: Specify the minimum threshold of the input array for figures
      vmax: Specify the maximum threshold of the input array for figures
    
    """
    with open(input_array_pickle, 'rb') as f:
        data = pickle.load(f)
    
    print("allow input data with chromosomes: ",data.keys())    
    output_dir = os.path.dirname(output_png)
    os.makedirs(output_dir, exist_ok=True)
    input_key = find_key(data,chrom1,chrom2)
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

    print("select matrix stat, min:",np.min(output_data),"max:",np.max(output_data),"mean:",np.mean(output_data))
    plt.imshow(output_data, cmap='coolwarm',vmin=vmin,vmax=vmax)
    #do not show the axis
    plt.axis('off')
    #remove all outer white space
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    plt.savefig(output_png, bbox_inches='tight', pad_inches=0,dpi=600)
    plt.close()
    print('The figure has been saved to',output_png)
    

"""
This script is to visualize the difference comparison array in png format.
```
python3 diff2png.py [input.pkl] [output.png] [chrom1] [start_index1] [end_index1] [chrom2] [start_index2] [end_index2] [resolution] [vmin] [vmax]
```
input.pkl: the path to the pickle file containing the difference comparison array [String]. <br>
input.pkl format: [chrom1_chrom2]:[array] format for common mode. Here array should be scipy sparce array."For intra-chromsome only, the dict format can be [chrom]:[array] in pickle files.<br>
output.png: the name of the output png file [String].<br>
chrom1: the name of the first chromosome [String].<br>
start_index1: the start index of the first chromosome [Integer].<br>
end_index1: the end index of the first chromosome [Integer].<br>
chrom2: the name of the second chromosome [String].<br>
start_index2: the start index of the second chromosome [Integer].<br>
end_index2: the end index of the second chromosome [Integer].<br>
resolution: resolution of the input array [Integer].<br>
vmin: the minimum threshold of the input array for figures[Float].<br>
vmax: the maximum threshold of the input array for figures[Float].<br>

"""


if __name__ == '__main__':
    import os
    import sys
    #take the overall array as input
    if len(sys.argv) != 12:
        print('Usage: python3 diff2png.py [input.pkl] [output.png] \
              [chrom1] [start_index1] [end_index1] \
              [chrom2] [start_index2] [end_index2] [resolution] [vmin] [vmax]')
        print("This script is to visualize the difference comparison array in png format. ")
        print("input.pkl: the path to the pickle file containing the difference comparison array [String].")
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
        print("vmin: the minimum threshold of the input array for figures[Float].")
        print("vmax: the maximum threshold of the input array for figures[Float].")
        print("All index input should be absolute index counted by base.")
        sys.exit(1)
    input_array_pickle = os.path.abspath(sys.argv[1])
    output_png = os.path.abspath(sys.argv[2])
    output_dir = os.path.dirname(output_png)
    os.makedirs(output_dir, exist_ok=True)
    chrom1 = sys.argv[3]
    start_index1 = int(sys.argv[4])
    end_index1 = int(sys.argv[5])
    chrom2 = sys.argv[6]
    start_index2 = int(sys.argv[7])
    end_index2 = int(sys.argv[8])
    resolution = int(sys.argv[9])
    vmin = float(sys.argv[10])
    vmax = float(sys.argv[11])

    start_index1 = start_index1//resolution
    end_index1 = end_index1//resolution
    start_index2 = start_index2//resolution
    end_index2 = end_index2//resolution
    diff2png(input_array_pickle,output_png,chrom1,start_index1,
              end_index1,chrom2,start_index2,end_index2,vmin,vmax)