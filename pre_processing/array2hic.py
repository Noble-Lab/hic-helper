import os 
import sys
import pickle
import numpy as np
#assume at least run on a machine with 8 CPUs+64G memory

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

def array2hic(juicer_tools,input_array_pickle,
              output_hic,resolution,refer_genome_name,mode=0):
    """
    The array2hic function converts a numpy array to hic file.
    
    :param juicer_tools: Specify the location of the juicer_tools
    :param input_array_pickle: Specify the path to the pickle file containing the array
    :param output_hic: Specify the name of the output hic file
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
    raw_path = output_hic.replace('.hic','.raw')
    with open(raw_path, 'w') as wfile:
        for key in data:
            if mode == 0 or mode == 2:
                chrom1, chrom2 = key.split('_')
            else:
                chrom1 = key
                chrom2 = key
            matrix = data[key]
            if mode>=2:
                matrix = array2sparse(matrix)
            #matrix merge records in the same loc
            matrix.eliminate_zeros()
            matrix.sum_duplicates()
            matrix_row = matrix.row
            matrix_col = matrix.col
            matrix_data = matrix.data
            if "chr" not in chrom1:
                chrom1 = "chr"+chrom1
            if "chr" not in chrom2:
                chrom2 = "chr"+chrom2
            for i in range(len(matrix_row)):
                wfile.write(f'{0} {chrom1} {int(matrix_row[i]*resolution+1)} {0} {0} {chrom2} {matrix_col[i]*resolution+1} {1} {matrix_data[i]:.2f}\n')
    code_path = os.path.dirname(juicer_tools)
    root_path = os.getcwd()
    os.chdir(code_path)
    os.system(f'java -Xmx64g -Xmx64g -jar juicer_tools.jar pre -j 8 -d -r {resolution} {raw_path} {output_hic}  {refer_genome_name}')
    os.remove(raw_path)

    os.chdir(root_path)

"""
Usage
```
python3 array2hic.py [input.pkl] [output.hic] [resolution] [refer_genome_name] [mode]
```
The input pickle should be in a pickle file as dict: [chrom1_chrom2]:[array] format for common mode. Here array should be scipy sparce array. <br>
For intra-chromsome only, the dict format can be [chrom]:[array] in pickle files.<br>
[output.hic] is the name of the output hic file. <br>
[resolution] is used to specify the resolution that stored in the output array. <br>
[refer_genome_name] is used to specify the reference genome name. For example, "hg38","hg19","mm10" are valid inputs. <br>
[mode]:  0: all chromosome mode (scipy sparce array); 1: intra-chromosome mode(scipy sparce array); 2: all chromosome mode (numpy array); 3: intra-chromosome mode(numpy array). <br>
"""

if __name__ == '__main__':
    
    #get current script directory 

    if len(sys.argv) != 6:
        print('Usage: python3 array2hic.py [input.pkl] [output.hic] [resolution] [refer_genome_name] [mode]')
        print("This is the full array2hic script. ")
        print("input.pkl: the path to the pickle file containing the array [String].")
        print("input.pkl format: [chrom1_chrom2]:[array] format for common mode. Here array should be scipy sparce array. For intra-chromsome only, the dict format can be [chrom]:[array] in pickle files.")
        print("output.hic: the name of the output hic file [String].")
        print("resolution: resolution of the input array [Integer].")
        print("refer_genome_name: the name of the reference genome [String]. Example: hg38, hg19, mm10.")
        print("mode: 0: all chromosome mode (scipy sparce array); 1: intra-chromosome mode(scipy sparce array); 2: all chromosome mode (numpy array); 3: intra-chromosome mode(numpy array).")
        sys.exit(1)
    script_dir = os.path.dirname(os.path.realpath(__file__))
    juicer_tools = os.path.join(script_dir, 'juicer_tools.jar')
    input_array_pickle = os.path.abspath(sys.argv[1])
    output_hic = os.path.abspath(sys.argv[2])
    resolution = int(sys.argv[3])
    refer_genome_name = str(sys.argv[4])
    mode = int(sys.argv[5])
    array2hic(juicer_tools,input_array_pickle,output_hic,resolution,refer_genome_name,mode)


