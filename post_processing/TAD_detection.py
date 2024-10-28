import argparse
import numpy as np
import os
import logging
import pickle 
import sys
eps=1e-20
from scipy.sparse import triu,coo_matrix
from collections import defaultdict


def argparser():
    parser = argparse.ArgumentParser('TAD detection reports and insulation score reports.', add_help=True)
    
    parser.add_argument('--input', type=str, help='Input pkl file path', required=True)
    parser.add_argument('--output', type=str, required=True, help='Output directory, bound file will be saved in TAD.bed, insulation score will be saved in insulation_score.pkl')
    parser.add_argument('--window-size', type=int, default=50)
    parser.add_argument('--delta-smooth-size', type=int, default=10)
    parser.add_argument('--bound-strength', type=float, default=0.1)
    args = parser.parse_args()
    return args 

def filter_sparse_rectangle(input_row,input_col,input_data,
                            start_row,end_row,start_col,end_col):
    """
    input_row: the row index of the sparse matrix
    input_col: the column index of the sparse matrix
    input_data: the data of the sparse matrix
    start_row: the start row index of the rectangle to filter
    end_row: the end row index of the rectangle to filter
    start_col: the start column index of the rectangle to filter
    end_col: the end column index of the rectangle to filter
    """
    select_index1 = (input_row>=start_row) & (input_row<end_row)
    select_index2 = (input_col>=start_col) & (input_col<end_col)
    select_index = select_index1 & select_index2
    new_row = input_row[select_index]-start_row
    new_col = input_col[select_index]-start_col
    new_data = input_data[select_index]
    new_shape = (end_row-start_row,end_col-start_col)
    final_array = coo_matrix((new_data,(new_row,new_col)),shape=new_shape)
    return final_array
def compute_insulation_score_sparse(matrix, window_size):
   
    L, _ = matrix.shape
    scores = []
    for i in range(L):
        if i<window_size or i+window_size >= L: scores.append(None)
        else:
            current_region = filter_sparse_rectangle(matrix.row, matrix.col, matrix.data,
                                                     i-window_size, i,i+1, i+window_size+1)
            #Returns the average of the array/matrix elements. The average is taken over all elements in the array/matrix by default, otherwise over the specified axis. float64 intermediate and return values are used for integer inputs.

            scores.append(current_region.sum()/window_size/window_size)
    return scores
def compute_insulation_score(matrix, window_size):
    #no need to sym matrix, since it only calcualte upper diagonal region.
    L, _ = matrix.shape
    scores = []
    for i in range(L):
        if i<window_size or i+window_size >= L: scores.append(None)
        else:
            current_region = matrix[i-window_size:i,i:i+window_size]
            scores.append(current_region.sum()/window_size/window_size)
    return scores

def compute_bounds(insulation_scores, delta_smooth_size, bound_strength):
    L = len(insulation_scores)

    mean_score = np.mean([s for s in insulation_scores if s is not None])

    normalized_scores = []
    for score in insulation_scores:
        if score is not None:
            normalized_scores.append((np.log(score+eps) - np.log(mean_score+eps))/np.log(2))
        else:
            normalized_scores.append(None)

    delta = []
    for i in range(L):
        left = []
        right = []
        for j in range(delta_smooth_size):
            if i+1+j<L and normalized_scores[i+1+j] is not None:
                right.append(normalized_scores[i+1+j])
            if i-j>=0 and normalized_scores[i-j] is not None:
                left.append(normalized_scores[i-j])
        if len(left) == 0 or len(right) == 0:
            delta.append(None)
        else:
            delta.append(np.mean(right) - np.mean(left))

    minimas = []
    for i in range(L):
        if i==0 or delta[i-1] is None or delta[i] is None: continue
        if delta[i-1] < 0 and delta[i] >0:
            minimas.append(i)

    # print(np.around(insulation_scores[5420:5470], decimals=4))
    # print(np.around(normalized_scores[5420:5470], decimals=4))
    # print(np.around(delta[5420:5470], decimals=4))

    bounds = []
    for minima in minimas:
        l = minima-1
        while delta[l-1] is not None and delta[l-1] <= delta[l]:
            l -= 1

        r = minima
        while delta[r+1] is not None and delta[r+1] >= delta[r]:
            r += 1

        if delta[r] - delta[l] >= bound_strength:
            bounds.append(minima)

    return bounds
"""
This script is used to detect TADs from the input pickle file. <br>
```
python3 TAD_detection.py --input [input.pkl] --output [output_dir] 
```
[input.pkl]: the input pickle file. <br>
[output_dir]: the output directory. <br>
In the output dir, the bound file will be saved in TAD.bed, insulation score will be saved in insulation_score.pkl.
"""

if __name__ == '__main__':
    args = argparser()
    input_pkl = os.path.abspath(args.input)
    data=pickle.load(open(input_pkl,'rb'))
    bound_dict={}
    insulation_score_dict={}
    for chrom in data:
        current_data = data[chrom]
        try:
            current_data = current_data.toarray()
            cur_insulation_scores = compute_insulation_score(current_data,args.window_size)
        except:
            print("Matrix is too large, use sparse matrix instead!")
            cur_insulation_scores = compute_insulation_score_sparse(current_data,args.window_size)
        cur_bounds= compute_bounds(cur_insulation_scores, args.delta_smooth_size, args.bound_strength)
        bound_dict[chrom]=cur_bounds
        insulation_score_dict[chrom]=cur_insulation_scores
        print(f"Finish processing {chrom}", "detected bounds:", len(cur_bounds))
    output_dir = os.path.abspath(args.output)
    os.makedirs(output_dir, exist_ok=True)
    bound_bed = os.path.join(output_dir, 'TAD.bed')
    insulation_score_pkl = os.path.join(output_dir, 'insulation_score.pkl')
    with open(bound_bed,'w') as f:
        for chrom in bound_dict:
            for bound in bound_dict[chrom]:
                f.write(f"{chrom}\t{bound}\t{bound+1}\n")
    pickle.dump(insulation_score_dict, open(insulation_score_pkl, "wb"))
    

