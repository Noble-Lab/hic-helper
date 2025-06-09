import numpy as np
import time 
import random
from numba import jit
def load_pkl(file_path):
    """
    Load a .pkl file containing Hi-C data.
    """
    import pickle
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
    return data

def save_pkl(data, file_path):
    """
    Save data to a .pkl file.
    """
    import pickle
    with open(file_path, 'wb') as f:
        pickle.dump(data, f)
@jit(nogil=True,nopython=True)
def weighted_choice(weights):
    """
    Return a single index sampled from 0..len(weights)-1
    according to the 1D array of probabilities `weights` (which should sum to >0).
    """
    # 1) build a cumulative distribution
    cum = np.empty(weights.shape[0], dtype=weights.dtype)
    total = 0.0
    for i in range(weights.shape[0]):
        total += weights[i]
        cum[i] = total

    # 2) draw a uniform in [0, total)
    u = np.random.random() * total

    # 3) find the first index where cum[i] > u
    for i in range(cum.shape[0]):
        if u < cum[i]:
            return i
    # due to rounding, fall back to last index
    return weights.shape[0] - 1


@jit(nogil=True,nopython=True)
def choice_nb(a, size, replace, p):
    """
    A simple re-implementation of numpy.random.choice for 1D array `a`,
    sampling `size` elements **with** or **without** replacement,
    using probability vector `p`.

    a       : 1D array of values to choose from
    size    : number of draws
    replace : True/False
    p       : 1D probability array same length as `a`, sum(p) > 0
    """

    out = np.empty(size, dtype=a.dtype)
    # if sampling without replacement, we will zero out used weights
    w = p.copy()

    for j in range(size):
        idx = weighted_choice(w)
        out[j] = a[idx]
        if not replace:
            # zero out and renormalize
            w[idx] = 0.0
    return out
@jit(nogil=True,nopython=True)
def list_to_array(lst):
    n = len(lst)
    arr = np.empty(n, dtype=np.float64)
    for i in range(n):
        arr[i] = lst[i]
    return arr
@jit(nogil=True,nopython=True)
def stratifiedSample ( V, F, strataSize = 100 ):
	N = len(V)
	# V = np.array(V)
	# F = np.array(F)
	V= list_to_array(V)
	F= list_to_array(F)
	strataCount = int(np.ceil(float(N) / strataSize))
	sortInd = np.argsort(F)
	strata = []
	strataMax = []

	#print '%d to stratify, %d strata to be filled' % (N,strataCount)

	
	for i in range(strataCount) :
		stratum = V [ sortInd[ (strataSize*(i) ) : (strataSize*(i+1)) ] ]
		stratumF = F [ sortInd[ (strataSize*(i) ) : (strataSize*(i+1)) ] ]
		strata.append( stratum )
		strataMax.append(max(stratumF))
		#print str(strataSize*(i)) + ' ' + str(strataSize*(i+1)) + ' ' + str(len(stratum))


	sample = []
	for i in range(len(V) ):
		if ( F[i] == 0 ) :
			sample.append (0)
		else :
			stratumInd = 0
			for k in range(strataCount) :
				#if ( F[i] >= strataMax[k] ):
				if F[i] <= strataMax[k]:
					stratumInd = k
					break
			if stratumInd == 0:
				stratumInd = k
			sample.append ( np.random.choice(strata[k],size=1)[0] )

	return sample
@jit(nogil=True,nopython=True)
def shuffleMatrix ( CM, stratumSize = 50 ):
	#Convert to integer
	#CM = CM.astype(int)
	#Get marginals and number of rows
	contactSum = np.sum(CM,1)#np.sum(np.array(CM),1)
	N = len(CM)

	# For matrix entry Mik, store Marginal i * Marginal k in CountByDist
	# and the Mik itself in matrixByDist
	countByDist = []
	matrixByDist = []
	for i in range(0,N):
		for k in range(i,N):
			dist = k-i
			if ( len(countByDist)-1 < dist ):
				countByDist.append( [ float(contactSum[i]) * contactSum[k] ] )
				matrixByDist.append( [ int( CM[i,k] ) ] )
			else:
				countByDist[dist].append( float(contactSum[i]) * contactSum[k] )
				matrixByDist[dist].append( int( CM[i,k] ) )
	
	noiseMatrix = np.zeros((N,N))
	

	for i in range(len(matrixByDist)):
	#for i in range(1):
		#print "dist is %d" % (i)
		thisSample = stratifiedSample(matrixByDist[i],countByDist[i],stratumSize)
		for k in range(len(thisSample)):				
			noiseMatrix[k,k+i] = thisSample[k]
		
	for i in range(0,N):
		for k in range(i,N):
			noiseMatrix[k,i] = noiseMatrix[i,k]
	
	
	return noiseMatrix


@jit(nogil=True,nopython=True)
def uniformMatrix ( CM, subSampleCount = 1000000, bias = False ):
	(R,C) = np.shape(CM)
	marginal = np.sum(CM,1)#np.sum(np.array(CM),1)
	uniSampleCM = np.zeros((R,C)) 
	#triuSum = sum(np.arange(R)+1)
	
	indexMap = []
	indexProb = []
	for i in range(R) :
		for k in range(i,R) :
			if marginal[i] != 0 and marginal[k] != 0 :
				indexMap.append([i,k])
			if bias :
				indexProb.append(marginal[i] * marginal[k])

	if bias :
		totalProb = float(sum(indexProb))
		indexProb = [ iP / totalProb for iP in indexProb ]
		triuSample = choice_nb(np.arange(len(indexMap)),subSampleCount,True,list_to_array(indexProb))
	else :
		triuSample = np.random.choice(len(indexMap),subSampleCount)
        	
	for s in triuSample :
		(i,k) = indexMap[s]
		uniSampleCM[i,k] += 1
	uniSampleCM += np.transpose(np.triu(uniSampleCM,1))

	return uniSampleCM

@jit(nogil=True,nopython=True)
def SubSampleMatrix(CM, subSampleN = 1000000, symmetric = True):

	if subSampleN >= np.sum(np.triu(CM)) : 
		print('Asked for ' + str(subSampleN) + 'entries, matrix has ' + str(np.sum(np.triu(CM))))
		print('Sampling more entries than available, returning original matrix')
		return CM

	index1 = []
	index2 = []
	subCM = np.zeros((len(CM),len(CM)))

	for i in range(len(CM)):
		for k in range(i,len(CM)):
			count = int(CM[i,k])
			v1=np.empty(count); v1.fill(i)				
			v2=np.empty(count); v2.fill(k)				
			index1.extend(v1); index2.extend(v2)
	
	# index1 = np.array(index1)
	# index2 = np.array(index2)
	index1 = list_to_array(index1)
	index2 = list_to_array(index2)
	shufIndex = range(0,len(index1))
	random.shuffle(shufIndex)
	subSampleIndex = np.random.choice(shufIndex,size=subSampleN,replace=False)
	index1 = index1[subSampleIndex]
	index2 = index2[subSampleIndex]

	for i in range(len(index1)):
		a = int(index1[i]); b = int(index2[i])
		subCM[a,b] = subCM[a,b] + 1

	subCM = subCM + np.triu(subCM,1).T		
	return subCM


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

def inject_noise(input_pkl, output_pkl, noise_percent, ligation_noise_percent, stratum_size):
	"""
	Inject noise into Hi-C data.
	"""
	input_data = load_pkl(input_pkl)
	output_dict={}
	for key in input_data:
		input_mat = input_data[key]
		if hasattr(input_mat, 'toarray'):
			#if scipy sparse, convert to dense
			input_mat = input_mat.toarray()
		if not np.array_equal(input_mat, input_mat.T):
			#judge if the array is symmetric, if not, change it to symmetric
			input_mat = np.triu(input_mat) + np.triu(input_mat, 1).T
		time_start = time.time()
		inputCoverage =  int(np.sum(np.triu(input_mat)))
		GDnoiseMatrix = shuffleMatrix(input_mat,stratum_size)
		RLnoiseMatrix = uniformMatrix(input_mat,inputCoverage,bias=True)
		realSampleCount = int(inputCoverage * ( 1 - float(noise_percent) ) )
		sinputMatrix = SubSampleMatrix(input_mat, subSampleN = realSampleCount )

		GDSampleCount = int(inputCoverage * float(noise_percent) * float(1-ligation_noise_percent) )
		sGDnoiseMatrix = SubSampleMatrix(GDnoiseMatrix, subSampleN = GDSampleCount )

		RLSampleCount = int(inputCoverage * float(noise_percent) * float(ligation_noise_percent) )
		sRLnoiseMatrix = SubSampleMatrix(RLnoiseMatrix, subSampleN = RLSampleCount )
		noisedMatrix = sinputMatrix + sGDnoiseMatrix + sRLnoiseMatrix
		#if input matrix is sparse, convert to sparse
		if hasattr(input_mat, 'toarray'):
			noisedMatrix = array2sparse(noisedMatrix)
		time_end = time.time()
		print(f"Processed {key}: Time taken = {time_end - time_start:.2f} seconds")
		#save the noisedMatrix to output_dict
		output_dict[key] = noisedMatrix
	save_pkl(output_dict, output_pkl)


if __name__ == '__main__':
    import os 
    import sys
    """
    from paper: Yardımcı, Galip Gürkan, et al. "Measuring the reproducibility and quality of Hi-C data." Genome biology 20 (2019): 1-19.

    noisePercentage, percent noise to be injected
    This argument (a float between 0 and 1) specifies what percent of noise will be injected to the real Hi-C contact matrix. For example, if this argument is 0.5, 50% noise is injected.

    Random ligation noise percentage
    Simualted noise is a mixture of genomic distance noise and random ligation noise.
    This argument (a float between 0 and 1) specifies what percent of simulated noise should consist of random ligation noise.

    Stratum Size, number of bins to make up each strata
    The number of bins that make up each strata for stratified sampling. For 25/40kb, 100 bin is a reasonable choice.

    """
    if len(sys.argv) != 6:
        print('Usage: python3 noise2hic.py [input.pkl] [output.pkl] [noise_percent] [ligation_noise_percent] [stratum_size]')
        print("The function is to inject noise to hic data. Here hic data is saved in .pkl files by hic2array.py. ")
        print("[input.pkl]: input .pkl that saved Hi-C data")
        print("[output.pkl]: output .pkl that saved Hi-C data with noise injected")
        print("[noise_percent]: this argument specifies what percent of noise will be injected to the real Hi-C contact matrix. ")
        print("[ligation_noise_percent]: this argument specifies what percent of simulated noise should consist of random ligation noise.")
        print("[stratum_size]: the number of bins that make up each strata for stratified sampling. For 25/40kb, 100 bin is a reasonable choice.")
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
    inject_noise(input_pkl, output_pkl, noise_percent, ligation_noise_percent, stratum_size)