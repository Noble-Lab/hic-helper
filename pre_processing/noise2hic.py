import numpy as np
import time
import random
from numba import jit
import os
import sys
import pickle
from scipy.sparse import coo_matrix

def load_pkl(file_path):
    """
    Load a .pkl file containing Hi-C data.

    Parameters
    ----------
    file_path : str
        Path to the pickle file.

    Returns
    -------
    data : object
        Data loaded from the pickle file.
    """
    with open(file_path, "rb") as f:
        data = pickle.load(f)
    return data


def save_pkl(data, file_path):
    """
    Save data to a .pkl file.

    Parameters
    ----------
    data : object
        Data to be saved.
    file_path : str
        Path to the output pickle file.
    """
    with open(file_path, "wb") as f:
        pickle.dump(data, f)


@jit(nogil=True, nopython=True)
def weighted_choice(weights):
    """
    Return a single index sampled from 0..len(weights)-1
    according to the 1D array of probabilities `weights` (which should sum to >0).

    Parameters
    ----------
    weights : np.ndarray
        1D array of probabilities.

    Returns
    -------
    int
        Selected index.
    """
    cum = np.empty(weights.shape[0], dtype=weights.dtype)
    total = 0.0
    for i in range(weights.shape[0]):
        total += weights[i]
        cum[i] = total

    u = np.random.random() * total

    for i in range(cum.shape[0]):
        if u < cum[i]:
            return i
    return weights.shape[0] - 1


@jit(nogil=True, nopython=True)
def choice_nb(a, size, replace, p):
    """
    A simple re-implementation of numpy.random.choice for 1D array `a`,
    sampling `size` elements **with** or **without** replacement,
    using probability vector `p`.

    Parameters
    ----------
    a : np.ndarray
        1D array of values to choose from.
    size : int
        Number of draws.
    replace : bool
        Sample with replacement if True.
    p : np.ndarray
        1D probability array same length as `a`, sum(p) > 0.

    Returns
    -------
    np.ndarray
        Sampled values.
    """
    out = np.empty(size, dtype=a.dtype)
    w = p.copy()

    for j in range(size):
        idx = weighted_choice(w)
        out[j] = a[idx]
        if not replace:
            w[idx] = 0.0
    return out


@jit(nogil=True, nopython=True)
def list_to_array(lst):
    """
    Convert a list to a numpy array of type float64.

    Parameters
    ----------
    lst : list
        List of values.

    Returns
    -------
    np.ndarray
        Array of float64.
    """
    n = len(lst)
    arr = np.empty(n, dtype=np.float64)
    for i in range(n):
        arr[i] = lst[i]
    return arr


@jit(nogil=True, nopython=True)
def stratifiedSample(V, F, strataSize=100):
    """
    Perform stratified sampling on V using F as the stratification variable.

    Parameters
    ----------
    V : list or np.ndarray
        Values to sample from.
    F : list or np.ndarray
        Stratification variable.
    strataSize : int
        Number of bins per stratum.

    Returns
    -------
    list
        Stratified sample.
    """
    N = len(V)
    V = list_to_array(V)
    F = list_to_array(F)
    strataCount = int(np.ceil(float(N) / strataSize))
    sortInd = np.argsort(F)
    strata = []
    strataMax = []

    for i in range(strataCount):
        stratum = V[sortInd[(strataSize * (i)) : (strataSize * (i + 1))]]
        stratumF = F[sortInd[(strataSize * (i)) : (strataSize * (i + 1))]]
        strata.append(stratum)
        strataMax.append(max(stratumF))

    sample = []
    for i in range(len(V)):
        if F[i] == 0:
            sample.append(0)
        else:
            stratumInd = 0
            for k in range(strataCount):
                if F[i] <= strataMax[k]:
                    stratumInd = k
                    break
            if stratumInd == 0:
                stratumInd = k
            sample.append(np.random.choice(strata[k], size=1)[0])

    return sample


@jit(nogil=True, nopython=True)
def shuffleMatrix(CM, stratumSize=50):
    """
    Shuffle a contact matrix using stratified sampling.

    Parameters
    ----------
    CM : np.ndarray
        Contact matrix.
    stratumSize : int
        Number of bins per stratum.

    Returns
    -------
    np.ndarray
        Shuffled matrix.
    """
    contactSum = np.sum(CM, 1)
    N = len(CM)
    countByDist = []
    matrixByDist = []
    for i in range(0, N):
        for k in range(i, N):
            dist = k - i
            if len(countByDist) - 1 < dist:
                countByDist.append([float(contactSum[i]) * contactSum[k]])
                matrixByDist.append([int(CM[i, k])])
            else:
                countByDist[dist].append(float(contactSum[i]) * contactSum[k])
                matrixByDist[dist].append(int(CM[i, k]))

    noiseMatrix = np.zeros((N, N))

    for i in range(len(matrixByDist)):
        thisSample = stratifiedSample(matrixByDist[i], countByDist[i], stratumSize)
        for k in range(len(thisSample)):
            noiseMatrix[k, k + i] = thisSample[k]

    for i in range(0, N):
        for k in range(i, N):
            noiseMatrix[k, i] = noiseMatrix[i, k]

    return noiseMatrix


@jit(nogil=True, nopython=True)
def uniformMatrix(CM, subSampleCount=1000000, bias=False):
    """
    Generate a symmetric uniformly sampled contact matrix.

    Parameters
    ----------
    CM : np.ndarray
        Input contact matrix.
    subSampleCount : int
        Number of samples.
    bias : bool
        Whether to weight sampling by marginal products.

    Returns
    -------
    np.ndarray
        Symmetric matrix of sample counts.
    """
    (R, C) = np.shape(CM)
    marginal = np.sum(CM, 1)
    uniSampleCM = np.zeros((R, C))

    indexMap = []
    indexProb = []
    for i in range(R):
        for k in range(i, R):
            if marginal[i] != 0 and marginal[k] != 0:
                indexMap.append([i, k])
                if bias:
                    indexProb.append(marginal[i] * marginal[k])

    if bias:
        totalProb = float(sum(indexProb))
        indexProb = [iP / totalProb for iP in indexProb]
        triuSample = choice_nb(
            np.arange(len(indexMap)), subSampleCount, True, list_to_array(indexProb)
        )
    else:
        triuSample = np.random.choice(len(indexMap), subSampleCount)

    for s in triuSample:
        i, k = indexMap[s]
        uniSampleCM[i, k] += 1
    uniSampleCM += np.transpose(np.triu(uniSampleCM, 1))

    return uniSampleCM

def uniformMatrix_fast(CM, subSampleCount=1_000_000, bias=False):
    """
    Generate a symmetric “uniformly sampled” contact matrix of the same shape as CM.

    If bias=False, each upper‐triangle pair (i,k) with nonzero marginals is equally likely.
    If bias=True, pairs are drawn with probability proportional to marginal[i] * marginal[k].

    Parameters
    ----------
    CM : np.ndarray, shape (R,R)
        Input contact matrix.
    subSampleCount : int
        Total number of (i,k) draws to make.
    bias : bool
        Whether to weight sampling by marginal products.

    Returns
    -------
    uniSampleCM : np.ndarray, shape (R,R)
        Symmetric matrix of sample counts.
    """
    R = CM.shape[0]
    marginal = CM.sum(axis=1)
    i_idx, k_idx = np.triu_indices(R)
    if bias:
        probs = marginal[i_idx] * marginal[k_idx]
        probs = probs / probs.sum()
        choices = np.random.choice(i_idx.size, size=subSampleCount, replace=True, p=probs)
    else:
        choices = np.random.randint(0, i_idx.size, size=subSampleCount)

    uniSampleCM = np.zeros((R, R), dtype=np.int64)
    selected_i = i_idx[choices]
    selected_k = k_idx[choices]
    np.add.at(uniSampleCM, (selected_i, selected_k), 1)
    uniSampleCM = uniSampleCM + np.triu(uniSampleCM, 1).T

    return uniSampleCM

@jit(nogil=True, nopython=True)
def SubSampleMatrix(CM, subSampleN=1000000, symmetric=True):
    """
    Subsample entries from a contact matrix.

    Parameters
    ----------
    CM : np.ndarray
        Input contact matrix.
    subSampleN : int
        Number of samples to draw.
    symmetric : bool
        Whether to symmetrize the output.

    Returns
    -------
    np.ndarray
        Subsampled matrix.
    """
    if subSampleN >= np.sum(np.triu(CM)):
        print(
            "Asked for "
            + str(subSampleN)
            + "entries, matrix has "
            + str(np.sum(np.triu(CM)))
        )
        print("Sampling more entries than available, returning original matrix")
        return CM

    index1 = []
    index2 = []
    subCM = np.zeros((len(CM), len(CM)))

    for i in range(len(CM)):
        for k in range(i, len(CM)):
            count = int(CM[i, k])
            v1 = np.empty(count)
            v1.fill(i)
            v2 = np.empty(count)
            v2.fill(k)
            index1.extend(v1)
            index2.extend(v2)

    index1 = list_to_array(index1)
    index2 = list_to_array(index2)

    subSampleIndex = np.random.choice(len(index1), size=subSampleN, replace=False)

    for i in range(len(subSampleIndex)):
        cur_index = int(subSampleIndex[i])
        a = int(index1[cur_index])
        b = int(index2[cur_index])
        subCM[a, b] = subCM[a, b] + 1

    subCM = subCM + np.triu(subCM, 1).T
    return subCM.astype(CM.dtype)


def array2sparse(array):
    """
    Convert a numpy array to a scipy sparse COO matrix.

    Parameters
    ----------
    array : np.ndarray
        Numpy array to convert.

    Returns
    -------
    coo_matrix
        Scipy sparse COO matrix.
    """
    row, col = np.where(array)
    data = array[row, col]
    return coo_matrix((data, (row, col)), shape=array.shape)


def inject_noise(
    input_pkl, output_pkl, noise_percent, ligation_noise_percent, stratum_size
):
    """
    Inject noise into Hi-C data.

    Parameters
    ----------
    input_pkl : str
        Path to input pickle file.
    output_pkl : str
        Path to output pickle file.
    noise_percent : float
        Fraction of noise to inject.
    ligation_noise_percent : float
        Fraction of noise that is random ligation noise.
    stratum_size : int
        Stratum size for stratified sampling.
    """
    input_data = load_pkl(input_pkl)
    output_dict = {}
    for key in input_data:
        #skip Un, Alt, and other non-standard chromosomes
        if "Un" in key or "Alt" in key or "chrM" in key or "alt" in key or "random" in key:
            print(f"Skipping {key} due to non-standard chromosome name.")
            continue
        input_mat = input_data[key]
        if hasattr(input_mat, "toarray"):
            input_mat = input_mat.toarray()
        if not np.array_equal(input_mat, input_mat.T):
            input_mat = np.triu(input_mat) + np.triu(input_mat, 1).T
        
        time_start = time.time()
        inputCoverage = int(np.sum(np.triu(input_mat)))
        GDnoiseMatrix = shuffleMatrix(input_mat, stratum_size)
        RLnoiseMatrix = uniformMatrix_fast(input_mat, inputCoverage, bias=True)
        realSampleCount = int(inputCoverage * (1 - float(noise_percent)))
        sinputMatrix = SubSampleMatrix(input_mat, subSampleN=realSampleCount)

        GDSampleCount = int(
            inputCoverage * float(noise_percent) * float(1 - ligation_noise_percent)
        )
        sGDnoiseMatrix = SubSampleMatrix(GDnoiseMatrix, subSampleN=GDSampleCount)

        RLSampleCount = int(
            inputCoverage * float(noise_percent) * float(ligation_noise_percent)
        )
        sRLnoiseMatrix = SubSampleMatrix(RLnoiseMatrix, subSampleN=RLSampleCount)
        noisedMatrix = sinputMatrix + sGDnoiseMatrix + sRLnoiseMatrix
        if hasattr(input_mat, "toarray"):
            noisedMatrix = array2sparse(noisedMatrix)
            noisedMatrix.eliminate_zeros()
        time_end = time.time()
        print(f"Processed {key}: Time taken = {time_end - time_start:.2f} seconds")
        output_dict[key] = noisedMatrix
    save_pkl(output_dict, output_pkl)


if __name__ == "__main__":
    """
    Script for injecting noise into Hi-C data.

    Usage:
        python3 noise2hic.py [input.pkl] [output.pkl] [noise_percent] [ligation_noise_percent] [stratum_size]

    Arguments:
        [input.pkl]: input .pkl that saved Hi-C data
        [output.pkl]: output .pkl that saved Hi-C data with noise injected
        [noise_percent]: percent of noise to inject (float between 0 and 1)
        [ligation_noise_percent]: percent of simulated noise that is random ligation noise (float between 0 and 1)
        [stratum_size]: number of bins per stratum for stratified sampling (e.g., 100 for 25/40kb)
    """
    if len(sys.argv) != 6:
        print(
            "Usage: python3 noise2hic.py [input.pkl] [output.pkl] [noise_percent] [ligation_noise_percent] [stratum_size]"
        )
        print(
            "The function is to inject noise to hic data. Here hic data is saved in .pkl files by hic2array.py. "
        )
        print("[input.pkl]: input .pkl that saved Hi-C data")
        print("[output.pkl]: output .pkl that saved Hi-C data with noise injected")
        print(
            "[noise_percent]: this argument specifies what percent of noise will be injected to the real Hi-C contact matrix. "
        )
        print(
            "[ligation_noise_percent]: this argument specifies what percent of simulated noise should consist of random ligation noise."
        )
        print(
            "[stratum_size]: the number of bins that make up each strata for stratified sampling. For 25/40kb, 100 bin is a reasonable choice."
        )
        sys.exit(1)
    input_pkl = os.path.abspath(sys.argv[1])
    output_pkl = os.path.abspath(sys.argv[2])
    output_dir = os.path.dirname(output_pkl)
    os.makedirs(output_dir, exist_ok=True)
    noise_percent = float(sys.argv[3])
    ligation_noise_percent = float(sys.argv[4])
    stratum_size = int(sys.argv[5])
    if not os.path.exists(input_pkl):
        print(f"Input file {input_pkl} does not exist.")
        sys.exit(1)
    inject_noise(
        input_pkl, output_pkl, noise_percent, ligation_noise_percent, stratum_size
    )
