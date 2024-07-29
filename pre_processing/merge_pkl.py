import os 
import sys
import pickle
import numpy as np
def merge_pkl(pkl_file1, pkl_file2, output_pkl):
    """
    Merge two pkl files to a new merged pkl file.
    Args:
        pkl_file1: the first pkl file to be merged.
        pkl_file2: the second pkl file to be merged.
        output_pkl: the output merged pkl file.
    return:
        None
    """
    
    with open(pkl_file1, 'rb') as rfile:
        data1 = pickle.load(rfile)
    with open(pkl_file2, 'rb') as rfile:
        data2 = pickle.load(rfile)
    final_data = {}
    for key in data1.keys():
        if key not in data2.keys():
            print(f"Warning: {key} not found in the second pkl file.")
            continue   
        cur_data1 = data1[key]
        cur_data2 = data2[key]
        #judge if cur_data1 is scipy.sparse.coo matrix or numpy
        if hasattr(cur_data1, 'toarray') and hasattr(cur_data2, 'toarray'):
            cur_data1.row = np.concatenate((cur_data1.row, cur_data2.row))
            cur_data1.col = np.concatenate((cur_data1.col, cur_data2.col))
            cur_data1.data = np.concatenate((cur_data1.data, cur_data2.data))
        else:
            cur_data1 = cur_data1 + cur_data2
        final_data[key] = cur_data1 
        print(f"Merge {key} successfully.")
    with open(output_pkl, 'wb') as wfile:
        pickle.dump(final_data, wfile)

"""
This script is to merge two pkl files to a new merged pkl file.
```
python3 merge_pkl.py <pkl_file1> <pkl_file2> <output_pkl>
```
[pkl_file1]: the first pkl file to be merged. <br>
[pkl_file2]: the second pkl file to be merged. <br>
[output_pkl]: the output merged pkl file. <br>
"""

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python3 merge_pkl.py <pkl_file1> <pkl_file2> <output_pkl>")
        print("This script is to merge two pkl files to a new merged pkl file.")
        print("[pkl_file1]: the first pkl file to be merged.")
        print("[pkl_file2]: the second pkl file to be merged.")
        print("[output_pkl]: the output merged pkl file.")
        exit(0)
    
    pkl_file1 = os.path.abspath(sys.argv[1])
    pkl_file2 = os.path.abspath(sys.argv[2])
    output_pkl = os.path.abspath(sys.argv[3])
    merge_pkl(pkl_file1, pkl_file2, output_pkl)