

import mmap
import os
import sys
import struct
import numpy as np
import pickle

def print_stderr(message):
    """
    Simply print str message to stderr
    """
    print(message, file=sys.stderr)


def force_exit(message, req=None):
    """
    Exit the program due to some error. Print out message and close the given
    input files.
    """
    if req:
        req.close()
    print_stderr(message)
    sys.exit(1)

# read function
def readcstr(f):
    # buf = bytearray()
    buf = b""
    while True:
        b = f.read(1)
        if b is None or b == b"\0":
            # return buf.encode("utf-8", errors="ignore")
            return buf.decode("utf-8")
        elif b == "":
            raise EOFError("Buffer unexpectedly empty while trying to read null-terminated string")
        else:
            buf += b
def read_header(req):
    """
    Takes in a .hic file and returns a dictionary containing information about
    the chromosome. Keys are chromosome index numbers (0 through # of chroms
    contained in file) and values are [chr idx (int), chr name (str), chrom
    length (str)]. Returns the masterindex used by the file as well as the open
    file object.

    """
    chrs = {}
    resolutions = []
    magic_string = struct.unpack(b'<3s', req.read(3))[0]
    req.read(1)
    if magic_string != b"HIC":
        error_string = '... This does not appear to be a HiC file; magic string is incorrect'
        force_exit(error_string, req)
    global version
    version = struct.unpack(b'<i', req.read(4))[0]
    footerPosition = struct.unpack(b'<q', req.read(8))[0]
    genome = b""
    c = req.read(1)
    while c != b'\0':
        genome += c
        c = req.read(1)
    genome = genome.decode('ascii')
    if version >= 9:
        frag_resolutions = []
        normVectorIndexPosition = struct.unpack('<q', req.read(8))[0]
        normVectorIndexLength = struct.unpack('<q', req.read(8))[0]
    nattributes = struct.unpack(b'<i', req.read(4))[0]
    metadata = {}
    for _ in range(nattributes):
        key = readcstr(req)
        value = readcstr(req)
        metadata[key] = value

    nChrs = struct.unpack(b'<i', req.read(4))[0]
    for i in range(nChrs):
        name = readcstr(req)
        if version >= 9:
            length = struct.unpack(b'<q', req.read(8))[0]
        else:
            length = struct.unpack(b'<i', req.read(4))[0]
        if name and length:
            chrs[i] = [i, name, length]

    nBpRes = struct.unpack(b'<i', req.read(4))[0]
    for _ in range(nBpRes):
        res = struct.unpack(b'<i', req.read(4))[0]
        resolutions.append(res)
        
    nBpResFrag = struct.unpack(b'<i', req.read(4))[0]
    for _ in range(nBpResFrag):
        res = struct.unpack(b'<i', req.read(4))[0]
        frag_resolutions.append(res)

    return chrs, resolutions, footerPosition, genome, metadata
def read_footer(f, buf, footerPosition,NORMS):
    f.seek(footerPosition)

    cpair_info = {}
    if version >= 9:
        nBytes = struct.unpack(b'<q', f.read(8))[0]
    else:
        nBytes = struct.unpack(b'<i', f.read(4))[0]
    nEntries = struct.unpack(b'<i', f.read(4))[0]
    for _ in range(nEntries):
        key = readcstr(f)
        fpos = struct.unpack(b'<q', f.read(8))[0]
        sizeinbytes = struct.unpack(b'<i', f.read(4))[0]
        cpair_info[key] = fpos

    expected = {}
    factors = {}
    norm_info = {}
    # raw (norm == 'NONE')
    nExpectedValues = struct.unpack(b'<i', f.read(4))[0]
    for _ in range(nExpectedValues):
        unit = readcstr(f)
        binsize = struct.unpack(b'<i', f.read(4))[0]
        if version >= 9:
            nValues = struct.unpack(b'<q', f.read(8))[0]
            expected['RAW', unit, binsize] = np.frombuffer(
                buf,
                dtype='<f',
                count=nValues,
                offset=f.tell())
            f.seek(nValues * 4, 1)
            nNormalizationFactors = struct.unpack(b'<i', f.read(4))[0]
            factors['RAW', unit, binsize] = np.frombuffer(
                buf,
                dtype={'names':['chrom','factor'], 'formats':['<i', '<f']},
                count=nNormalizationFactors,
                offset=f.tell())
            f.seek(nNormalizationFactors * 8, 1)
        else:
            nValues = struct.unpack(b'<i', f.read(4))[0]
            expected['RAW', unit, binsize] = np.frombuffer(
                buf,
                dtype=np.dtype('<d'),
                count=nValues,
                offset=f.tell())
            f.seek(nValues * 8, 1)
            nNormalizationFactors = struct.unpack(b'<i', f.read(4))[0]
            factors['RAW', unit, binsize] = np.frombuffer(
                buf,
                dtype={'names':['chrom','factor'], 'formats':['<i', '<d']},
                count=nNormalizationFactors,
                offset=f.tell())
            f.seek(nNormalizationFactors * 12, 1)
            
    # normalized (norm != 'NONE')
    possibleNorms = f.read(4)
    if not possibleNorms:
        print_stderr('!!! WARNING. No normalization vectors found in the hic file.')
        return cpair_info, expected, factors, norm_info
    
    
    nExpectedValues = struct.unpack(b'<i', possibleNorms)[0]
    for _ in range(nExpectedValues):
        normtype = readcstr(f)
        if normtype not in NORMS:
            NORMS.append(normtype)
        unit = readcstr(f)
        binsize = struct.unpack(b'<i', f.read(4))[0]
        if version >= 9:
            nValues = struct.unpack(b'<q', f.read(8))[0]
            expected[normtype, unit, binsize] = np.frombuffer(
                buf,
                dtype='<f',
                count=nValues,
                offset=f.tell())
            f.seek(nValues * 4, 1)
            nNormalizationFactors = struct.unpack(b'<i', f.read(4))[0]
            factors[normtype, unit, binsize] = np.frombuffer(
                buf,
                dtype={'names':['chrom','factor'], 'formats':['<i', '<f']},
                count=nNormalizationFactors,
                offset=f.tell())
            f.seek(nNormalizationFactors * 8, 1)
        else:
            nValues = struct.unpack(b'<i', f.read(4))[0]
            expected['RAW', unit, binsize] = np.frombuffer(
                buf,
                dtype=np.dtype('<d'),
                count=nValues,
                offset=f.tell())
            f.seek(nValues * 8, 1)
            nNormalizationFactors = struct.unpack(b'<i', f.read(4))[0]
            factors['RAW', unit, binsize] = np.frombuffer(
                buf,
                dtype={'names':['chrom','factor'], 'formats':['<i', '<d']},
                count=nNormalizationFactors,
                offset=f.tell())
            f.seek(nNormalizationFactors * 12, 1)

    nEntries = struct.unpack(b'<i', f.read(4))[0]
    for _ in range(nEntries):
        normtype = readcstr(f)
        chrIdx = struct.unpack(b'<i', f.read(4))[0]
        unit = readcstr(f)
        resolution = struct.unpack(b'<i', f.read(4))[0]
        filePosition = struct.unpack(b'<q', f.read(8))[0]
        if version >= 9:
            sizeInBytes = struct.unpack(b'<q', f.read(8))[0]
        else:
            sizeInBytes = struct.unpack(b'<i', f.read(4))[0]
        norm_info[normtype, unit, resolution, chrIdx] = {
            'filepos': filePosition,
            'size': sizeInBytes
        }

    return cpair_info, expected, factors, norm_info,NORMS

def read_normalization_vector(f, buf, entry):
    filepos = entry['filepos']
    f.seek(filepos)
    nValues = struct.unpack(b'<i', f.read(4))[0]
    return np.frombuffer(buf, dtype=np.dtype('<d'), count=nValues, offset=filepos+4)

def extract_hicnorms_binary(input_hic,norm_type,use_resolution):
    """
    input_hic: .hic file
    """
    hic_norms = {}
    unit = 'BP' # only using base pair unit for now
    with open(input_hic, 'r') as req:
        buf = mmap.mmap(req.fileno(), 0, access=mmap.ACCESS_READ)
        used_chrs, resolutions, masteridx,genome,metadata = read_header(buf)
        #chrs, resolutions, frag_resolutions, footerPosition, normVectorIndexPosition, normVectorIndexLength, genome, metadata
       
        #chrs[i] = [i, name, length] dict
        #resolutions = list of resolutions
        #masteridx = master index
        #genome = genome name
        #metadata in a dict
        if use_resolution not in resolutions:
            output_string='... Resolution not found in the hic file. '
            output_string+='Available resolutions are: '+', '.join([str(x) for x in resolutions])
            force_exit(output_string, req)
        pair_footer_info, expected, factors, norm_info,NORMS = read_footer(req, buf, masteridx,[])
        if norm_type not in NORMS:
            output_string='... Normalization type not found in the hic file. '
            output_string+='Available normalization types are: '+', '.join(NORMS)
            force_exit(output_string, req)
        for chr_val in [uc for uc in used_chrs.values() if uc[1].lower() != 'all']:
            chr_idx, chr_name, chr_len = chr_val
            chr_num_bins = int(np.ceil(chr_len / use_resolution))
            try:
                norm_key = norm_info[norm_type,unit,use_resolution,chr_idx]
            except KeyError:
                print_stderr('!!! WARNING. No normalization vectors found for chr %s in the hic file.' % chr_name)
                norm_vector = [np.nan]*chr_num_bins
            else:
                norm_vector = read_normalization_vector(req, buf, norm_key)[:chr_num_bins]
            hic_norms[chr_name] = norm_vector   
    
    return hic_norms
def write_pickle(output_pkl,hic_norms):
    with open(output_pkl, 'wb') as f:
        pickle.dump(hic_norms, f)
    print(f"Normalization vectors saved to {output_pkl}")
import hicstraw
def extract_hicnorms(input_hic,norm_type,resolution):
    """
    input_hic: .hic file
    norm_type: normalization type
    resolution: resolution to extract the normalization vector
    """
    hic = hicstraw.HiCFile(input_hic)
    chrom_list=[]
    chrom_dict={}
    for chrom in hic.getChromosomes():
        print(chrom.name, chrom.length)
        if "all" in chrom.name.lower():
            continue
        chrom_list.append(chrom)
        chrom_dict[chrom.name]=chrom.length
    resolution_list = hic.getResolutions()
    if resolution not in resolution_list:
        print("Resolution not found in the hic file, please choose from the following list:")
        print(resolution_list)
        exit()
    hicnorms = {}
    for i in range(len(chrom_list)):
        for j in range(i,len(chrom_list)):
            if i!=j:
                #skip inter-chromosome region
                #since shared the norm vector
                continue
            chrom1 = chrom_list[i]
            chrom1_name = chrom_list[i].name
            chrom2 = chrom_list[j]
            chrom2_name = chrom_list[j].name
            mzd = hic.getMatrixZoomData(chrom1_name, chrom2_name, 'observed', norm_type, "BP", resolution)
            norm_vector = mzd.getNormVector(chrom1.index)
            hicnorms[chrom1_name] = norm_vector
    return hicnorms


"""
This script is to extract the normalization vectors from a .hic file. <br>
```
python3 extract_hicnorms.py [input.hic] [resolution] [normalization_type] [output_pkl]
```
[input.hic]: input hic path. <br>
[resolution]: resolution to extract the normalization vector, [Integer]. <br>
[normalization_type]: should be one of the following: NONE, VC, VC_SQRT, KR, SCALE, [string]. <br>
[output_pkl]: output pickle file path. <br> 
The normalization vector is saved in dict format, where the key is the chromosome name and the value is the normalization vector. <br>
"""

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python3 extract_hicnorms.py [input.hic] [resolution] [normalization_type] [output_pkl]")
        print("[input.hic]: input hic path. [str]")
        print("[resolution]: resolution to extract the normalization vector, [Integer].")
        print("[normalization_type]: should be one of the following: NONE, VC, VC_SQRT, KR, SCALE. [str]")
        print("[output_pkl]: output pickle file path. [str]")
        print("The normalization vector is saved in dict format, where the key is the chromosome name and the value is the normalization vector.")
        sys.exit(1)
    input_hic = os.path.abspath(sys.argv[1])
    resolution = int(sys.argv[2])
    norm_type = sys.argv[3]
    output_pkl = os.path.abspath(sys.argv[4])
    output_dir = os.path.dirname(output_pkl)
    os.makedirs(output_dir, exist_ok=True)
    if norm_type not in ['NONE', 'VC', 'VC_SQRT', 'KR', 'SCALE']:
        print("Normalization type should be one of the following: NONE, VC, VC_SQRT, KR, SCALE.")
        sys.exit(1)
    #hic_norms = extract_hicnorms_binary(input_hic,norm_type,resolution)
    #switch to hicstraw to get the normalization vectors
    hic_norms = extract_hicnorms(input_hic,norm_type,resolution)
    write_pickle(output_pkl,hic_norms)
        