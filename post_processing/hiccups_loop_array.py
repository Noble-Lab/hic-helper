import os
import sys
import numpy as np
from scipy.ndimage import correlate
from scipy.stats import poisson
import scipy
import logging
import pickle
def print_info(s:str):
    print(s)
    logging.info(s)

def donut_kernel(donut_size, peak_size):
    R1 = donut_size
    R2 = peak_size

    kernel = np.ones((R1*2+1, R1*2+1))
    center = (R1, R1)

    kernel[center[0] - R2 : center[0]+R2+1, center[1] - R2 : center[1]+R2+1] = 0
    kernel[center[0], :] = 0
    kernel[:, center[1]] = 0

    return kernel

def lowerleft_kernel(donut_size, peak_size):
    R1 = donut_size
    R2 = peak_size

    kernel = np.ones((R1*2+1, R1*2+1))
    center = (R1, R1)

    kernel[center[0] - R2 : center[0]+R2+1, center[1] - R2 : center[1]+R2+1] = 0
    kernel[:center[0]+1, :] = 0
    kernel[:, center[1]:] = 0

    return kernel

def horizontal_kernel(donut_size, peak_size):
    R1 = donut_size
    R2 = peak_size

    kernel = np.zeros((3, R1*2+1))
    center = (1, R1)

    kernel[ : , : center[1] - R2 ] = 1
    kernel[ : , center[1] + R2 + 1 : ] = 1

    return kernel

def vertical_kernel(donut_size, peak_size):
    R1 = donut_size
    R2 = peak_size

    kernel = np.zeros((R1*2+1, 3))
    center = (R1, 1)

    kernel[ : center[0] - R2, : ] = 1
    kernel[ center[0] + R2 + 1 : , : ] = 1

    return kernel


def loop_clustering( peak_cands,clustering_boundary, singleton_qvalue):
    num_cands = len(peak_cands)
    peaks_final = []
    while len(peak_cands) > 0:
        top_peak = max(peak_cands)
        peak_cands.remove(top_peak)
        peaks_cluster = [top_peak[1]]
        centroid = top_peak[1]
        r = 0
        find = True
        while find:
            find = False

            def dis(x, y):
                return np.linalg.norm((x[0]-y[0], x[1]-y[1]))

            centroid = np.mean(peaks_cluster, axis = 0)
            r = max([dis(peak, centroid) for peak in peaks_cluster ])

            for peak in peak_cands:
                if dis(peak[1], centroid) - r < clustering_boundary:
                    peaks_cluster.append(peak[1])
                    peak_cands.remove(peak)
                    find = True
                    break
        if r>0 or top_peak[2] <= singleton_qvalue:
            peaks_final.append((top_peak[1], centroid, r))

    print_info(f'Found {len(peaks_final)} peaks from {num_cands} candidates')

    return peaks_final

def get_oe_matrix(matrix, bounding = 100000000, oe=True):
    """
    matrix: numpy array
    bounding: int, the maximum distance to calculate the expected value
    oe: bool, whether to return the oe matrix
    """
    max_offset = min(matrix.shape[0], bounding)

    expected = [np.mean(np.diagonal(matrix, offset)) for offset in range(max_offset)]

    if oe:
        e_matrix = np.zeros_like(matrix, dtype=np.float32)
        oe_matrix = np.zeros_like(matrix, dtype=np.float32)
        for i in range(matrix.shape[0]):
            for j in range(max(i-bounding+1, 0), min(i+bounding, matrix.shape[1])):
                e_matrix[i][j] = expected[abs(i-j)]
                oe_matrix[i][j] = matrix[i][j]/expected[abs(i-j)] if expected[abs(i-j)] != 0 else 0

        return oe_matrix, e_matrix
    else:
        return expected


def identify_ignore_regions(matrix):
    """
    scipy array input
    """
    input_shape = matrix.shape
    matrix = np.nan_to_num(matrix)
    matrix_sum = np.sum(matrix, axis=1)
    compact_index = set()
    for i in range(input_shape[0]):
        if matrix_sum[i] != 0:
            compact_index.add(i)
    return compact_index




def find_peaks(full_matrix, donut_size,peak_size,clustering_boundary, 
               singleton_qvalue,lambda_step,
               FDR=0.1, thresholds=[0.02,1.5,1.75,2], bounding_size=400,gap_filter_range=5):
    compact_ids = identify_ignore_regions(full_matrix)
    print("compact %d/%d"%(len(compact_ids),len(full_matrix))) 
    kernels = [donut_kernel(donut_size,peak_size), lowerleft_kernel(donut_size,peak_size),
               horizontal_kernel(donut_size,peak_size), vertical_kernel(donut_size,peak_size)]
    kernel_names = ['donut', 'lowerleft', 'horizontal', 'vertical']
    matrix_shape = full_matrix.shape
    l = matrix_shape[0]
    B = min(bounding_size, l)
    window_size = min(2*B, l)
    upper_triangle = np.triu(np.ones((window_size, window_size)), 0)

    expect_vector = get_oe_matrix(full_matrix, bounding = window_size, oe=False)

    expect = np.zeros((window_size, window_size))
    for i in range(window_size):
        for j in range(window_size):
            if abs(i-j) < len(expect_vector):
                expect[i][j] = expect_vector[abs(i-j)]
    esums = []
    for kernel in kernels:
        esum = correlate(expect, kernel, mode='constant')
        esums.append(esum)

    enriched_pixels = []

    first_patch_qvalues = np.copy(upper_triangle)
    pbar = range(0, l, B)
    for s0 in pbar:
        pbar.set_description(f'Currently find {len(enriched_pixels)} enriched pixels')

        s = min(s0, l-window_size)
        matrix = full_matrix[s:s+window_size, s:s+window_size]
        
        observed = matrix 

        observed = (np.rint(observed)).astype(int)

        log_lambda_step = np.log(lambda_step)

        pixel_scores = {}

        # print(observed)

        for kid, kernel in enumerate(kernels):
            msum = correlate(matrix, kernel, mode='constant')
            esum = esums[kid]

            Ek = np.nan_to_num(msum/esum*expect)
            Ek = Ek 

            #lambda-chunk FDR

            logEk = np.nan_to_num(np.log(Ek))

            bin_id = np.maximum(0, np.ceil(logEk/log_lambda_step).astype(int))
            pvalues = poisson.sf(observed, np.exp(bin_id*log_lambda_step))

            max_bin = bin_id.max()+1

            for id in range(max_bin):
                bin_pos = np.where((bin_id == id) & (upper_triangle == 1))
                p = pvalues[bin_pos]

                bin = sorted(zip(p.tolist(), bin_pos[0].tolist(), bin_pos[1].tolist()))
                size = len(bin)

                qvalue = 1
                for rank in range(len(bin), 0, -1):
                    pvalue, i, j = bin[rank-1]
                    qvalue = min(qvalue, pvalue /(rank / size))

                    if s==0:
                        first_patch_qvalues[i,j] = min(first_patch_qvalues[i,j], 1-qvalue if qvalue <= FDR else 0)

                    if qvalue <= FDR and observed[i][j]/Ek[i][j] > thresholds[kid]:
                        # pass BHFDR, check ratio
                        flag = (kid in (0,1)) and (observed[i][j]/Ek[i][j] > thresholds[-1])
                        if (i,j) not in pixel_scores: pixel_scores[(i,j)] = [0, 0]
                        pixel_scores[(i,j)][0] += 2 + (1 if flag else 0)
                        pixel_scores[(i,j)][1] += qvalue

        for p, v in pixel_scores.items():
            if v[0]>=9 and abs(p[0]-p[1]) <= bounding_size:
                enriched_pixels.append((observed[p[0], p[1]], (p[0]+s, p[1]+s), v[1]))

    gaps = set(range(l)) - set(compact_ids)
    near_gap = [False for _ in range(l)]
    for gap in gaps:
        for i in range(gap_filter_range):
            if gap-i >= 0:
                near_gap[gap-i] = True
            if gap+i < l:
                near_gap[gap+i] = True

    filtered_enriched_pixels = []
    for pixels in enriched_pixels:
        if not near_gap[pixels[1][0]] and not near_gap[pixels[1][1]]:
            filtered_enriched_pixels.append(pixels)

    peaks_final = loop_clustering(filtered_enriched_pixels,clustering_boundary, singleton_qvalue)

    return peaks_final

def annotate_peaks(annotate_path,peaks,resolution,chr="chr1"):
    with open(annotate_path, 'a+') as wfile:
        
        for top, center, radius in peaks:
            x1 = int(top[0] * resolution)
            x2 = int(top[0] * resolution + resolution)
            y1 = int(top[1] * resolution)
            y2 = int(top[1] * resolution + resolution)
            wfile.write(f'{chr}\t{x1}\t{x2}\t{chr}\t{y1}\t{y2}\t.\t.\t.\t.\t0,0,255\n')
def hiccups_loop(input_pkl,output_loop,resolution):
    output_dir = os.path.dirname(output_loop)
    os.makedirs(output_dir,exist_ok=True)

    if resolution not in [5000,10000,25000]:
        raise ValueError('The resolution is not supported. The supported resolutions are 5000,10000,25000')
    with open(output_loop, 'w') as wfile:
        wfile.write(f'chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\tcolor\n')
    
    if resolution == 5000:
        peak_size = 4
        donut_size = 7
        clustering_boundary = 20000//resolution
        
    elif resolution == 10000:
        peak_size = 2
        donut_size = 5
        clustering_boundary = 20000//resolution
        
    elif resolution == 25000:
        peak_size = 1
        donut_size = 3
        clustering_boundary = 50000//resolution
    singleton_qvalue = 0.02
    lambda_step = 2**(1/3)
    #load array
    with open(input_pkl, 'rb') as f:
        data = pickle.load(f)
    for key in data:
        if "_" in key:
            chrom1, chrom2 = key.split('_')
        else:
            chrom1 = key
            chrom2 = key
        if chrom1 != chrom2:
            print_info(f'Skip {chrom1} {chrom2}, loop only detect in the same chromosome')
            continue
        matrix = data[key]
        #if it is sparcse matrix, convert to dense matrix
        if type(matrix) == scipy.sparse.coo_matrix:
            matrix = matrix.toarray()
        print(f"start identifying peaks from {chrom1}",matrix.shape)
        peaks_final =find_peaks(matrix, donut_size,peak_size,clustering_boundary, 
               singleton_qvalue,lambda_step,
               FDR=0.1, thresholds=[0.02,1.5,1.75,2], bounding_size=400,gap_filter_range=5)
        annotate_peaks(output_loop,peaks_final,resolution,chr=chrom1)
        print(f"for chrom {chrom1}, collected {len(peaks_final)} peaks!")    

# Path: post_processing/hiccups_loo_array.py
"""
```
python3 hiccups_loop_array.py [input_pkl] [output_bed] [resolution]
```
[input_pkl]: the path to the input pkl file [String]. <br>
the pkl file should be a dictionary with keys as the cell line names and values as the path to the hic file [String]. <br>
[output_bed]: the bed files to save the output loops [String]. <br>
[resolution]: the resolution of the input hic file [Integer]. <br>
Currently only support 5000,10000,25000 resolutions. <br>
The output loop bedpe file will be saved in [output_dir]/merged_loops.bedpe. <br>
"""

if __name__ == '__main__':
    if len(sys.argv) !=4:
        print('Usage: python3 hiccups_loop_array.py [input_pkl] [output_bed] [resolution]')
        print("input_pkl: the path to the input pkl file [String].")
        print("output_bed: the bed files to save the output loops [String].")
        print("resolution: the resolution of the input hic file [Integer].")
        sys.exit(1)
    input_pkl = os.path.abspath(sys.argv[1])
    output_bed = os.path.abspath(sys.argv[2])
    resolution = int(sys.argv[3])
    hiccups_loop(input_pkl,output_bed,resolution)
