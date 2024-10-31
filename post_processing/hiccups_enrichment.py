import pickle
from scipy.sparse import coo_matrix
from scipy.stats import poisson
import numpy as np
import pandas as pd
from tqdm import tqdm
import argparse

import numba


# Four Helper routines for creating different kernels
# used to compute expectations
# See Rao 2014 Cell paper for visualizations
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
    
    kernel = np.zeros((R1*2+1, R1*2+1))
    center = (R1, R1)

    kernel[center[0] - 1 : center[0] + 2, :center[1] - R2] = 1
    kernel[center[0] - 1 : center[0] + 2, center[1] + R2 + 1:] = 1
    
    return kernel

def vertical_kernel(donut_size, peak_size):
    R1 = donut_size
    R2 = peak_size
    
    kernel = np.zeros((R1*2+1, R1*2+1))
    center = (R1, R1)

    kernel[:center[0] - R2, center[1] - 1 : center[1] + 2] = 1
    kernel[center[0] + R2 + 1:, center[1] - 1 : center[1] + 2] = 1

    return kernel

def compute_expvec(coo_mat, savepath=None, savememory=False):
    """
    Compute the expected diagonal sum vector from a COO sparse matrix.

    Parameters
    ----------
    coo_mat : scipy.sparse.coo_matrix
        The input sparse matrix in COO format from which the expected diagonal 
        sum vector is to be computed.
    
    savepath : str or None, optional
        File path to save the computed `exp_vec` as a .npy file.
    
    savememory : bool, optional
        If True, performs a memory-efficient computation using the sparse matrix 
        directly, which is slower but uses less memory. If False (default), 
        converts the matrix to a dense format for faster computation.

    Returns
    -------
    exp_vec : numpy.ndarray
        The expected diagonal sum vector.

    """
    l = coo_mat.shape[0]

    if not savememory:
        # Convert the sparse matrix to a dense matrix
        dense_mat = coo_mat.toarray()
        exp_vec = [np.diag(dense_mat, d).sum() / (l - d) for d in tqdm(range(l), desc="Computing exp_vec")]

    else:
        exp_vec = np.zeros(l)
        num_nonzeros = coo_mat.nnz
        for row, col, value in tqdm(zip(coo_mat.row, coo_mat.col, coo_mat.data),
                                    total=num_nonzeros, desc="Computing exp_vec"):
            d = abs(row - col)
            if d < l:
                exp_vec[d] += value

        exp_vec = exp_vec / (l - np.arange(l))
    
    # Save to a file if a savepath is provided
    if savepath is not None:
        np.save(savepath, np.array(exp_vec))

    return exp_vec

def compute_exp_submat(exp_vec, bin1, bin2, donut_size):
    """
    Create a submatrix of expected values centered at (bin1, bin2).

    Parameters
    ----------
    exp_vec : numpy.ndarray
        Vector of expected values for each diagonal.
    bin1 : int
        Row index of the center.
    bin2 : int
        Column index of the center.
    donut_size : int
        Radius of the submatrix.

    Returns
    -------
    expsubmat : numpy.ndarray
        A (2 * donut_size + 1) x (2 * donut_size + 1) matrix of expected values.
    """
    expsubmat = np.zeros((2 * donut_size + 1, 2 * donut_size + 1))

    diff = bin2 - bin1

    np.fill_diagonal(expsubmat, exp_vec[diff])

    for i in range(1, 2 * donut_size + 1):
        np.fill_diagonal(expsubmat[:, i:], exp_vec[diff + i])
        np.fill_diagonal(expsubmat[i:, :], exp_vec[diff - i])
    
    return expsubmat

def center_within_scope(max_scope):
    """
    Generate integer offsets within a circular scope.

    Parameters
    ----------
    max_scope : int
        Maximum radius of the circular area.

    Returns
    -------
    offsets : numpy.ndarray
        Array of (x, y) offsets as int32, representing coordinates 
        within the specified circular scope.
    """
    x = np.arange(-max_scope, max_scope + 1)
    y = np.arange(-max_scope, max_scope + 1)
    xx, yy = np.meshgrid(x, y)
    
    distances = np.sqrt(xx**2 + yy**2)
    mask = distances <= max_scope
    offsets = np.column_stack((xx[mask], yy[mask]))
    
    return offsets.astype(np.int32)

@numba.njit(parallel=True)
def compute_loop_ratio_with_wobble(bin1_arr, bin2_arr, obs_arr, 
                                   exp_d_arr, exp_ll_arr, exp_v_arr, exp_h_arr,
                                   donut_size, kernel_d, kernel_ll, kernel_h, kernel_v, 
                                   expmat_list, obsmat_norm_list, obsmat_raw_list,
                                   wobble_scope, center_offset_list, zeros_thresh):
    """
    Compute loop enrichment ratios with wobble for given loci.

    Parameters
    ----------
    bin1_arr, bin2_arr : numpy.ndarray
        Arrays of initial bin coordinates for loci.
    obs_arr : numpy.ndarray
        Array to store observed contact counts at the final coordinates.
    exp_d_arr, exp_ll_arr, exp_v_arr, exp_h_arr : numpy.ndarray
        Arrays to store expected contact counts using donut, lower-left, 
        vertical, and horizontal kernels.
    donut_size : int
        Radius of the donut kernel.
    kernel_d, kernel_ll, kernel_h, kernel_v : numpy.ndarray
        Kernels used to compute expected contact counts.
    expmat_list : numpy.ndarray
        List of expected submatrices for each locus.
    obsmat_norm_list : numpy.ndarray
        List of normalized observed submatrices for each locus.
    obsmat_raw_list : numpy.ndarray
        List of raw observed submatrices for each locus.
    wobble_scope : int
        Maximum allowed wobble around the initial coordinates.
    center_offset_list : numpy.ndarray
        List of (x, y) offsets within the wobble scope.
    zeros_thresh : int
        Maximum number of zeros allowed in the normalized submatrix.

    Notes
    -----
    Updates input arrays in place with the best loop enrichment ratios 
    and adjusted coordinates based on wobble.
    """
    R1 = donut_size
    center_list_len = center_offset_list.shape[0]
    for idx1 in numba.prange(expmat_list.shape[0]):
        obsmat_norm = obsmat_norm_list[idx1]
        obsmat_raw = obsmat_raw_list[idx1]
        expmat = expmat_list[idx1]
        
        loop_ratios_candidates = np.full(center_list_len, -1.0, dtype=np.float32)
        obs_candidates = np.full(center_list_len, 0, dtype=np.float32)
        exp_d_candidates = np.full(center_list_len, 0, dtype=np.float32)
        exp_ll_candidates = np.full(center_list_len, 0, dtype=np.float32)
        exp_h_candidates = np.full(center_list_len, 0, dtype=np.float32)
        exp_v_candidates = np.full(center_list_len, 0, dtype=np.float32)
        
        binx_offset_candidates = np.full(center_list_len, 0, dtype=np.int32)
        biny_offset_candidates = np.full(center_list_len, 0, dtype=np.int32)
    
        for idx2 in range(center_offset_list.shape[0]):
            center_offset_x, center_offset_y = center_offset_list[idx2]
            center_x = R1 + wobble_scope + center_offset_x
            center_y = R1 + wobble_scope + center_offset_y
            
            obssubmat_norm = obsmat_norm[center_x - R1: center_x + R1 + 1, center_y - R1: center_y + R1 + 1]
            expsubmat = expmat[center_x - R1: center_x + R1 + 1, center_y - R1: center_y + R1 + 1]
            obs_raw_center = obsmat_raw[wobble_scope + center_offset_x, wobble_scope + center_offset_x]
     
            if obs_raw_center == 0.0:
                continue
                
            # Check for correct shape
            if obssubmat_norm.shape != (2*R1+1, 2*R1+1):
                continue

            # Check zeros threshold
            zeros_count = np.sum((obssubmat_norm == 0).astype(np.int32))
            if zeros_count > zeros_thresh:
                continue
            
            exp_center = expsubmat[R1, R1]
            center_obs = obssubmat_norm[R1, R1]
            
            # Compute numerators and denominators for exp values
            numerator_d = np.sum(obssubmat_norm * kernel_d)
            denominator_d = np.sum(expsubmat * kernel_d)
            numerator_ll = np.sum(obssubmat_norm * kernel_ll)
            denominator_ll = np.sum(expsubmat * kernel_ll)
            numerator_h = np.sum(obssubmat_norm * kernel_h)
            denominator_h = np.sum(expsubmat * kernel_h)
            numerator_v = np.sum(obssubmat_norm * kernel_v)
            denominator_v = np.sum(expsubmat * kernel_v)
            
            ratio_d = -1.0
            ratio_ll = -1.0
            ratio_h = -1.0
            ratio_v = -1.0
            
            # Compute ratios if denominators are not zero
            if exp_center != 0:
                if denominator_d != 0.0 and numerator_d != 0.0:
                    exp_d = numerator_d / denominator_d * exp_center
                    ratio_d = center_obs / exp_d

                if denominator_ll != 0.0 and numerator_ll != 0.0:
                    exp_ll = numerator_ll / denominator_ll * exp_center
                    ratio_ll = center_obs / exp_ll

                if denominator_h != 0.0 and numerator_h != 0.0:
                    exp_h = numerator_h / denominator_h * exp_center
                    ratio_h = center_obs / exp_h

                if denominator_v != 0.0 and numerator_v != 0.0:
                    exp_v = numerator_v / denominator_v * exp_center
                    ratio_v = center_obs / exp_v
                
            loop_ratios_candidates[idx2] = min(ratio_d, ratio_ll, ratio_h, ratio_v)
            
            if loop_ratios_candidates[idx2] <= 0.0:
                continue
            else:
                # Convert to raw counts for poisson stats
                obs_candidates[idx2] = obs_raw_center
                exp_d_candidates[idx2] = exp_d * (obs_raw_center / center_obs)
                exp_ll_candidates[idx2] = exp_ll * (obs_raw_center / center_obs)
                exp_h_candidates[idx2] = exp_h * (obs_raw_center / center_obs)
                exp_v_candidates[idx2] = exp_v * (obs_raw_center / center_obs)
                
                binx_offset_candidates[idx2] = center_offset_x
                biny_offset_candidates[idx2] = center_offset_y
            
        if np.max(loop_ratios_candidates) <= 0.0:
            bin1_arr[idx1] = -1
            bin2_arr[idx1] = -1
        else:
            max_index = np.argmax(loop_ratios_candidates)
            bin1_arr[idx1] += binx_offset_candidates[max_index]
            bin2_arr[idx1] += biny_offset_candidates[max_index]
            obs_arr[idx1] = obs_candidates[max_index]
            exp_d_arr[idx1] = exp_d_candidates[max_index]
            exp_ll_arr[idx1] = exp_ll_candidates[max_index]
            exp_h_arr[idx1] = exp_h_candidates[max_index]
            exp_v_arr[idx1] = exp_v_candidates[max_index]           

def read_bedpe(filepath):
    """
    Read a .bedpe file into a DataFrame.

    Parameters
    ----------
    filepath : str
        Path to the .bedpe file.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing bedpe data with inferred column names.
    """
    with open (filepath, 'r') as file:
        for line in file:
            if line.startswith("#chr1"):
                header = line.strip().lstrip("#").split('\t')
                break
   
    df = pd.read_csv(filepath, sep = "\t", 
                     comment = "#", names = header)
    return df

def get_coords_df(bedpe_df, resolution=5000):
    """
    Compute bin coordinates from a bedpe DataFrame.

    Parameters
    ----------
    bedpe_df : pandas.DataFrame
        DataFrame containing bedpe data with 'chr1', 'chr2', 'x1', and 'y1' columns.
    resolution : int, optional
        Bin size for calculating bin coordinates, default is 5000.

    Returns
    -------
    df_subset : pandas.DataFrame
        DataFrame with 'chr1', 'chr2', 'x_bin', and 'y_bin' columns.
    """
    df_subset = bedpe_df[["chr1", "chr2"]].copy()
    df_subset["x_bin"] = bedpe_df["x1"] // resolution
    df_subset["y_bin"] = bedpe_df["y1"] // resolution
    return df_subset
    
def read_input_pkl(pkl_filepath, norm, drop_chroms=["chrY", "chrM"]):
    """
    Load raw and normalized data from a pickle file, filtering chromosomes.

    Parameters
    ----------
    pkl_filepath : str
        Path to the input pickle file.
    norm : str
        Key for the desired normalization type (e.g., 'VC', 'KR').
    drop_chroms : list of str, optional
        List of chromosomes to exclude, default is ['chrY', 'chrM'].

    Returns
    -------
    raw_data_dict : dict
        Dictionary with raw data (COO matrices) for each chromosome pair.
    norm_data_dict : dict
        Dictionary with normalized data (COO matrices) for each chromosome pair.
    """
    with open(pkl_filepath, 'rb') as f:
        data = pickle.load(f)
        
    raw_data_dict = {}
    norm_data_dict = {}
    for key, value in data.items():
        chrom1, chrom2 = key.split("_")
        
        if chrom1 in drop_chroms or chrom2 in drop_chroms:
            continue
        
        assert "NONE" in value, f"{key}: not containing raw data"
        assert isinstance(data[key]["NONE"], coo_matrix), f"{key}: raw data is not a COO sparse matrix."
        assert norm in value, f"{key}: not containing {norm} data"
        assert isinstance(data[key][norm], coo_matrix), f"{key}: norm data is not a COO sparse matrix."
        
        raw_data_dict[key] = value["NONE"]
        norm_data_dict[key] = value[norm]
    
    return raw_data_dict, norm_data_dict

def write_peak_output(output_bedpe_path, resolution,
                       chroms_arr, xbin_arr, ybin_arr,
                       match_xbin_arr, match_ybin_arr,
                       obs_arr, 
                       exp_d_arr, exp_ll_arr, exp_v_arr, exp_h_arr, 
                       pval_d_arr, pval_ll_arr, pval_v_arr, pval_h_arr):
    """
    Write peak data to a .bedpe file.

    Parameters
    ----------
    output_bedpe_path : str
        File path for the output bedpe file.
    resolution : int
        Resolution for converting bin indices to genomic coordinates.
    chroms_arr : numpy.ndarray
        Array of chromosome pair strings (e.g., 'chr1_chr2').
    xbin_arr, ybin_arr : numpy.ndarray
        Arrays of bin coordinates query peak locations in input .bedpe file.
    match_xbin_arr, match_ybin_arr : numpy.ndarray
        Arrays of bin coordinates for matched peaks found in .pkl file.
    obs_arr : numpy.ndarray
        Array of observed contact counts (raw counts).
    exp_d_arr, exp_ll_arr, exp_v_arr, exp_h_arr : numpy.ndarray
        Arrays of expected contact counts using donut, lower-left, 
        vertical, and horizontal kernels, respectively.
    pval_d_arr, pval_ll_arr, pval_v_arr, pval_h_arr : numpy.ndarray
        Arrays of p-values for the corresponding expected contact counts.
    """
    chrom1 = [chrom.split("_")[0] for chrom in chroms_arr]
    chrom2 = [chrom.split("_")[1] for chrom in chroms_arr]
    x1 = match_xbin_arr * resolution
    x2 = x1 + resolution
    y1 = match_ybin_arr * resolution
    y2 = y1 + resolution
    ref_x1 = xbin_arr * resolution
    ref_x2 = ref_x1 + resolution
    ref_y1 = ybin_arr * resolution
    ref_y2 = ref_y1 + resolution
    peak_df = pd.DataFrame({
        "#chr1": chrom1, "x1": x1, "x2": x2,
        "chr2": chrom2, "y1": y1, "y2": y2,
        "observed": obs_arr,
        "expectedBL": exp_ll_arr, 
        "expectedDonut": exp_d_arr,
        "expectedH": exp_h_arr,
        "expectedV": exp_v_arr,
        "pvalBL": pval_ll_arr,
        "pvalDonut": pval_d_arr,
        "pvalH": pval_h_arr,
        "pvalV": pval_v_arr,
        "ref_chr1": chrom1, "ref_x1": ref_x1, "ref_x2": ref_x2, 
        "ref_chr2": chrom2, "ref_y1": ref_y1, "ref_y2": ref_y2
    })
    peak_df.to_csv(output_bedpe_path, index=False, sep="\t")

def run_loop_ratio_wobble(input_bedpe_path, pkl_filepath, output_bedpe_path, norm, 
                          wobble_scope=4, zeros_thresh=225, donut_size=7, peak_size=4, 
                          resolution=5000, drop_chroms=["chrY", "chrM"], savememory=False):
    """
    Calculate loop enrichment with wobble and write results to a .bedpe file.

    Parameters
    ----------
    input_bedpe_path : str
        Path to the input .bedpe file containing loci information.
    pkl_filepath : str
        Path to the pickle file with raw and normalized Hi-C data.
    output_bedpe_path : str
        Path for saving the output BEDPE file with computed loop ratios.
    norm : str
        Key for the normalization type (e.g., 'VC', 'KR').
    wobble_scope : int, optional
        Maximum allowed wobble around initial bin coordinates, default is 4.
    zeros_thresh : int, optional
        Maximum number of zeros allowed in the observed submatrix, default is 225.
    donut_size : int, optional
        Radius of the donut kernel, default is 7.
    peak_size : int, optional
        Radius of the peak region, default is 4.
    resolution : int, optional
        Resolution for binning genomic coordinates, default is 5000.
    drop_chroms : list of str, optional
        Chromosomes to exclude from analysis, default is ['chrY', 'chrM'].
    savememory : bool, optional
        If True, uses memory-efficient calculations, default is False.

    Notes
    -----
    This function reads Hi-C data from a pickle file,
    using coordinates from input .bedpe file as queries, 
    calculates loop enrichment ratios with wobble, 
    filters results based on sparsity, 
    computes poisson p-values, and writes results to an output .bedpe file.
    """  
    raw_data_dict, norm_data_dict = read_input_pkl(pkl_filepath, norm, drop_chroms = drop_chroms)
    
    df_bedpe = read_bedpe(input_bedpe_path)
    df_coords = get_coords_df(df_bedpe, resolution=resolution)

    chrom_full_list = ["chr1", "chr2", "chr3", "chr4", "chr5", 
                  "chr6", "chr7", "chr8", "chr9", "chr10",
                  "chr11", "chr12", "chr13", "chr14", "chr15", 
                  "chr16", "chr17", "chr18", "chr19", "chr20", 
                  "chr21", "chr22", "chrX", "chrY", "chrM"]
    
    chrom_list = [chrom for chrom in chrom_full_list 
                  if chrom not in drop_chroms]
    
    chroms_arr = (df_coords["chr1"] + "_" + df_coords["chr2"]).values.astype("str")
    xbin_arr = df_coords["x_bin"].values
    ybin_arr = df_coords["y_bin"].values
    
    mask = [
        (chrom.split("_")[0] in chrom_list) and (chrom.split("_")[1] in chrom_list) 
        for chrom in chroms_arr
    ]
    
    chroms_arr = chroms_arr[mask]
    xbin_arr = xbin_arr[mask]
    ybin_arr = ybin_arr[mask]
    
    exp_vec_dict = {}
    for chrom in tqdm(chrom_list):
        sparse_norm_mat = norm_data_dict[f"{chrom}_{chrom}"]
        exp_vec = compute_expvec(sparse_norm_mat, savememory=savememory)
        exp_vec_dict[chrom] = exp_vec
    print(f"Finish computing exp_vec_dict")

    csr_norm_dict = {}
    for chrom in tqdm(chrom_list):
        sparse_norm_mat = norm_data_dict[f"{chrom}_{chrom}"]
        csr_norm_dict[chrom] = sparse_norm_mat.tocsr()
    print(f"Finish computing csr_mat_dict")
    
    csr_raw_dict = {}
    for chrom in tqdm(chrom_list):
        sparse_raw_mat = raw_data_dict[f"{chrom}_{chrom}"]
        csr_raw_dict[chrom] = sparse_raw_mat.tocsr()
    print(f"Finish computing csr_raw_dict")     
    
    # Initialize arrays to be passed into numba function
    # loop_ratios = np.zeros(len(chroms_arr))

    R1 = donut_size
    R2 = peak_size
    kernel_d = donut_kernel(R1, R2)
    kernel_ll = lowerleft_kernel(R1, R2)
    kernel_h = horizontal_kernel(R1, R2)
    kernel_v = vertical_kernel(R1, R2)

    # Initialize arrays with correct shape
    looplist_len = len(chroms_arr)
    expmat_list = np.zeros((looplist_len, 2 * (R1 + wobble_scope) + 1, 2 * (R1 + wobble_scope) + 1))
    obsmat_norm_list = np.zeros((looplist_len, 2 * (R1 + wobble_scope) + 1, 2 * (R1 + wobble_scope) + 1))
    obsmat_raw_list = np.zeros((looplist_len, 2 * (wobble_scope) + 1, 2 * (wobble_scope) + 1))

    valid_indices = []
    for index, (chrom, xbin, ybin) in enumerate(tqdm(zip(chroms_arr, xbin_arr, ybin_arr), total=looplist_len)):
        chrom1 = chrom.split("_")[0]
        # Extract submatrices
        obsmat_norm = csr_norm_dict[chrom1][xbin - (R1 + wobble_scope): xbin + (R1 + wobble_scope) + 1, 
                                        ybin - (R1 + wobble_scope): ybin + (R1 + wobble_scope) + 1].todense()
        obsmat_raw = csr_raw_dict[chrom1][xbin - (wobble_scope): xbin + (wobble_scope) + 1, 
                                    ybin - (wobble_scope): ybin + (wobble_scope) + 1].todense()
        if obsmat_norm.shape != (2 * (R1 + wobble_scope) + 1, 2 * (R1 + wobble_scope) + 1):
            continue
        if obsmat_raw.shape != (2 * wobble_scope + 1, 2 * wobble_scope + 1):
            continue
        
        # Assign the submatrices
        obsmat_norm_list[index] = obsmat_norm
        obsmat_raw_list[index] = obsmat_raw 
              
        exp_vec = exp_vec_dict[chrom1]
        expmat_list[index] = compute_exp_submat(exp_vec, xbin, ybin, (R1 + wobble_scope))
        
        # Store valid indices
        valid_indices.append(index)

    expmat_list = expmat_list[valid_indices]
    obsmat_norm_list = obsmat_norm_list[valid_indices]
    obsmat_raw_list = obsmat_raw_list[valid_indices]
    xbin_arr = xbin_arr[valid_indices]
    ybin_arr = ybin_arr[valid_indices]
    chroms_arr = chroms_arr[valid_indices]

    match_xbin_arr = xbin_arr.copy()
    match_ybin_arr = ybin_arr.copy()
    obs_arr = np.full(len(xbin_arr), 0, dtype=np.float32)
    exp_d_arr = np.full(len(xbin_arr), 0, dtype=np.float32)
    exp_ll_arr = np.full(len(xbin_arr), 0, dtype=np.float32)
    exp_v_arr = np.full(len(xbin_arr), 0, dtype=np.float32)
    exp_h_arr = np.full(len(xbin_arr), 0, dtype=np.float32)
    
    expmat_list = np.array(expmat_list, dtype=np.float32)
    obsmat_norm_list = np.array(obsmat_norm_list, dtype=np.float32)
    obsmat_raw_list = np.array(obsmat_raw_list, dtype=np.float32)
    print(f"Finish computing expmat_list, obsmat_norm_list, obsmat_raw_list\n")
    
    center_offset_list = center_within_scope(wobble_scope)
    
    compute_loop_ratio_with_wobble(match_xbin_arr, match_ybin_arr, obs_arr, 
                                   exp_d_arr, exp_ll_arr, exp_v_arr, exp_h_arr,
                                   donut_size, kernel_d, kernel_ll, kernel_h, kernel_v, 
                                   expmat_list, obsmat_norm_list, obsmat_raw_list,
                                   wobble_scope, center_offset_list, zeros_thresh)

    locus_pair_mask = (match_xbin_arr >= 0) & (match_ybin_arr >= 0)
    chroms_arr = chroms_arr[locus_pair_mask]
    xbin_arr = xbin_arr[locus_pair_mask]
    ybin_arr = ybin_arr[locus_pair_mask]
    match_xbin_arr = match_xbin_arr[locus_pair_mask]
    match_ybin_arr = match_ybin_arr[locus_pair_mask]
    obs_arr = obs_arr[locus_pair_mask]
    exp_d_arr = exp_d_arr[locus_pair_mask]
    exp_ll_arr = exp_ll_arr[locus_pair_mask]
    exp_v_arr = exp_v_arr[locus_pair_mask]
    exp_h_arr = exp_h_arr[locus_pair_mask]
    pval_d_arr = poisson.sf(obs_arr, exp_d_arr)
    pval_ll_arr = poisson.sf(obs_arr, exp_ll_arr)
    pval_v_arr = poisson.sf(obs_arr, exp_v_arr)
    pval_h_arr = poisson.sf(obs_arr, exp_h_arr)
    
    write_peak_output(output_bedpe_path, resolution,
                       chroms_arr, xbin_arr, ybin_arr,
                       match_xbin_arr, match_ybin_arr,
                       obs_arr, 
                       exp_d_arr, exp_ll_arr, exp_v_arr, exp_h_arr, 
                       pval_d_arr, pval_ll_arr, pval_v_arr, pval_h_arr)
                                   
    return 0    

def main():
    """
    Main function to run loop enrichment calculation with wobble.
    
    sample usage
    python3 hiccups_enrichment.py \
    --input_bedpe input_data.bedpe \
    --input_pkl hic_data.pkl \
    --output_bedpe output_results.bedpe \
    --norm KR \
    --num_cpus 8 \
    --wobble_scope 4 \
    --zeros_thresh 225 \
    --donut_size 7 \
    --peak_size 4 \
    --resolution 5000 \
    --drop_chroms chrY chrM \
    --savememory
    """
    parser = argparse.ArgumentParser(description="Calculate loop enrichment with wobble and output to BEDPE file.")
    # Required parameters
    parser.add_argument("--input_bedpe", type=str, required=True, 
                        help="Path to the input .bedpe file.")
    parser.add_argument("--input_pkl", type=str, required=True, 
                        help="Path to the input pickle file with Hi-C data.")
    parser.add_argument("--output_bedpe", type=str, required=True, 
                        help="Path to save the output .bedpe file.")
    parser.add_argument("--norm", type=str, required=True, 
                        help="Normalization type for Hi-C data (e.g., 'VC', 'KR').")
    parser.add_argument("--num_cpus", type=int, required=True, 
                        help="Number of CPUs for parallelization.")
    # Optional parameters   
    parser.add_argument("--wobble_scope", type=int, default=4,
                        help="Maximum allowed wobble around initial coordinates (default: 4).")
    parser.add_argument("--zeros_thresh", type=int, default=225,
                        help="Max zeros allowed in observed submatrix (default: 225).")
    parser.add_argument("--donut_size", type=int, default=7,
                        help="Radius of the donut kernel (default: 7).")
    parser.add_argument("--peak_size", type=int, default=4,
                        help="Radius of the peak region (default: 4).")
    parser.add_argument("--resolution", type=int, default=5000,
                        help="Resolution for binning genomic coordinates (default: 5000).")
    parser.add_argument("--drop_chroms", nargs="+", default=["chrY", "chrM"],
                        help="Chromosomes to exclude (default: ['chrY', 'chrM']).")
    parser.add_argument("--savememory", action="store_true",
                        help="Use memory-efficient calculations if set.")
    
    args = parser.parse_args()

    numba.set_num_threads(args.num_cpus)

    run_loop_ratio_wobble(
        input_bedpe_path=args.input_bedpe,
        pkl_filepath=args.input_pkl,
        output_bedpe_path=args.output_bedpe,
        norm=args.norm,
        wobble_scope=args.wobble_scope,
        zeros_thresh=args.zeros_thresh,
        donut_size=args.donut_size,
        peak_size=args.peak_size,
        resolution=args.resolution,
        drop_chroms=args.drop_chroms,
        savememory=args.savememory
    )
"""
This script is to calculate loop enrichment and output to BEDPE file.
```
python3 hiccups_enrichment.py --input_bedpe [input.bed] --input_pkl [hic.pkl] \
--output_bedpe [output.bed] --norm [norm_type] --num_cpus [int] --wobble_scope [int] \
--zeros_thresh [int] --donut_size [int] --peak_size [int] --resolution [int] \
--drop_chroms [chrY chrM] --savememory
```
[input.bed]: Path to the input .bedpe file, which records the loci information. <br>
[hic.pkl]: Path to the input pickle file stored Hi-C data. <br>
[output.bed]: Path to save the output .bedpe file. <br>
[norm_type]: Normalization type for Hi-C data (e.g., 'VC', 'KR'). Required for peak enrichment calculation. <br>
[num_cpus]: Number of CPUs for parallelization. <br>
[wobble_scope]: Maximum allowed wobble around the initial bin coordinates. Default is 4(5Kb),2(10Kb),2(25Kb). <br>
[zeros_thresh]: Maximum number of zeros allowed in the observed submatrix. Should set to (2*donut_size+1)^2. <br>
[donut_size]: Radius of the donut kernel. Default is 7(5Kb),5(10Kb),3(25Kb). <br>
[peak_size]: Radius of the peak region. Default is 4(5Kb),2(10Kb),1(25Kb). <br>
[resolution]: Resolution for binning genomic coordinates. Default is 5000. <br>
[drop_chroms]: Chromosomes to exclude from analysis. Default is ['chrY', 'chrM']. For example, ``--drop_chroms chrY chrM`` <br>
[savememory]: Use memory-efficient calculations if set. <br>
"""
if __name__ == "__main__":
    main()
