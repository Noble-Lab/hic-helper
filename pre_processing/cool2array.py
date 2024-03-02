import numpy as np
import cooler
from scipy.sparse import coo_matrix
import pickle 
def write_pkl(data, path):
    with open(path, 'wb') as f:
        pickle.dump(data, f)
def cool2array(cooler_path,normalize=False,tondarray=False):
    """
    # *** This is only for the intra-chromosome processing!!!!  only work for fixed bin size!!!!***
    # *** This is only for the intra-chromosome processing!!!!  only work for fixed bin size!!!!***
    #you can easily modify it for all chromosome processing
    cooler_path: the path to the cooler file
    normalize: if True, the matrix will be normalized by the norm matrix saved in the cooler file
    tondarray: if True, the return a numpy array dict
    return a numpy/scipy.sparse array dict
    [chromosome]:sparce matrix
    """
    c = cooler.Cooler(cooler_path)
    binsize= c.info['bin-size']
    chromosomes= c.chromnames
    chromosome_sizes = c.chromsizes
    bins = c.bins()[:] #including the chromosome staring ending info for each bin
    bins['bin_id'] = bins.index
    #column of bins[['chrom', 'start', 'end','weight']
    pixels = c.pixels()[:] #including the chromosome staring ending info for each bin
    #including bin1_id,bin2_id,count columns
    return_dict={}
    for k,chromsome in enumerate(chromosomes):
        cur_chromosome_size=chromosome_sizes[k]
        cur_array_length = int(np.ceil(cur_chromosome_size/binsize))
        #filter bins first via the chromosome name
        cur_bins = bins[bins['chrom']==chromsome] #most important is the bin_id
        #filter pixels first via the bin_id 
        cur_bin_ids = cur_bins['bin_id'].tolist()
        select_index1 = pixels['bin1_id'].isin(cur_bin_ids)
        select_index2 = pixels['bin2_id'].isin(cur_bin_ids)
        cur_pixels = pixels[(select_index1)&(select_index2)]
        #merge join first based on bin_id1 
        cur_pixels['bin_id'] = cur_pixels['bin1_id']
        cur_pixels = cur_pixels.merge(cur_bins,on='bin_id',how='left',suffixes=('','1')) #chrom1, start1, end1
        cur_pixels['bin_id'] = cur_pixels['bin2_id']
        cur_pixels = cur_pixels.merge(cur_bins,on='bin_id',how='left',suffixes=('','2'))
        #get the matrix
        current_table = cur_pixels
        row = np.array(current_table['start'].tolist())/binsize
        column = np.array(current_table['start2'].tolist())/binsize
        #convert to int
        row = row.astype(np.int32)
        column = column.astype(np.int32)
        count = np.array(current_table['count'].tolist())
        #can be easily extended to support inter-cross-chromosome
        if normalize:
            # apply the balancing weights
            weight_row = np.array(current_table['weight'].tolist())
            weight_col = np.array(current_table['weight2'].tolist())
            count = weight_row*weight_col* count
        #coo_matrix will automatically accumulate for same row/col with different count
        final_mat = coo_matrix((count, (row,column)), shape = (cur_array_length,cur_array_length),dtype=np.float32)
        if tondarray:
            final_mat = final_mat.toarray()
        return_dict[chromsome]=final_mat
    return return_dict


"""
Usage
```
python3 cool2array.py [input.cool] [output.pkl] [mode]
```
This is the intra-chromosome cool2array scripts, can be easily extended to full cool2array. <br>
The output array is saved in a pickle file as dict: [chrom]:[array] format. <br>
Two modes are supported: 
```
0: scipy coo_array format output; 
1: numpy array format output;
2: normed scipy coo_array format output; 
3: normed numpy array format output.
```
"""

if __name__ == '__main__':
    import os 
    import sys
    if len(sys.argv) != 4:
        print('Usage: python3 cool2array.py [input.cool] [output.pkl] [mode]')
        print("This is the intra-chromosome cool2array scripts, can be easily extended to full cool2array. ")
        print("mode: 0 for sparse matrix, 1 for dense matrix, 2 for normed sparse matrix, 3 for normed dense matrix")
        sys.exit(1)

    cooler_path = os.path.abspath(sys.argv[1])
    output_pkl_path = os.path.abspath(sys.argv[2])
    output_dir = os.path.dirname(output_pkl_path)
    os.makedirs(output_dir,exist_ok=True)
    mode = int(sys.argv[3])
    if mode not in [0,1,2,3]:
        print('mode should be 0,1,2,3')
        sys.exit(1)
    if mode == 0:
        normalize = False
        tondarray = False
    elif mode == 1:
        normalize = False
        tondarray = True
    elif mode == 2:
        normalize = True
        tondarray = False
    elif mode == 3: 
        normalize = True
        tondarray = True
    return_dict = cool2array(cooler_path,normalize=normalize,tondarray=tondarray)
    write_pkl(return_dict,output_pkl_path)
    