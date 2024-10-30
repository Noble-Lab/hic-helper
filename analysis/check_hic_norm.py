import argparse
import numpy as np
import hicstraw

def comma_sep_ints(value):
    """Converts a comma-separated string of integers into a list of integers."""
    return [int(item) for item in value.split(',')]

def comma_sep_strs(value):
    """Converts a comma-separated string into a list of strings."""
    return value.split(',')

def get_hic_chrom_list(hic):
    """
    Returns a list of chromosome names from the .hic file.
    """
    # return [chrom.name for chrom in hic.getChromosomes() if chrom.name != "All"]
    return [chrom.name for chrom in hic.getChromosomes()]
    
def check_hic_norm(hic_file, chrom_list, resolutions, normalizations):
    """
    Checks the normalization vectors for intra-chromosomal matrices in a .hic file.
    
    Parameters:
    - hic_file: Path to the .hic file.
    - chrom_list: List of chromosomes to check.
    - resolutions: List of resolutions to check.
    - normalizations: List of normalization types to check.
    
    Raises:
    - Exception if a chromosome is not found in the .hic file's genome assembly.
    """
    hic = hicstraw.HiCFile(hic_file)
    valid_chrom_list = get_hic_chrom_list(hic)
    # Map chrom name to their indices in the list.
    valid_chrom_index_dict = {entry: index for index, entry in enumerate(valid_chrom_list)}
    
    # Validate that all chromosomes in the list are present in the .hic file
    for chrom in chrom_list:
        if chrom not in valid_chrom_list:
            raise Exception(f"{chrom} not in genome assembly of .hic file!")
    
    # Iterate through resolutions, chromosomes, and normalizations
    for resolution in resolutions:
        for idx, chrom in enumerate(chrom_list):
            for norm in normalizations:
                try:
                    # Attempt to get matrix zoom data for the given parameters
                    mzd = hic.getMatrixZoomData(chrom, chrom, "observed", norm, "BP", resolution)
                    
                    # Check if the normalization vector is NaN
                    chrom_idx = valid_chrom_index_dict[chrom]
                    # use chrom.index
                    norm_vector = mzd.getNormVector(chrom_idx)
                    if norm_vector is None or np.isnan(norm_vector).all():
                        print(f"{norm} {chrom} {resolution}: norm vector all NaN")
                
                except MemoryError:
                    # Handle MemoryError gracefully
                    print(f"MemoryError: Could not process {chrom} at {resolution} with {norm} normalization.")
                
                except Exception as e:
                    # General exception handling with more context
                    print(f"Unexpected error for {chrom} at {resolution} with {norm}: {e}")

    print("Checking on normalization vectors complete!")

def main():
    """
    Main function for parsing command-line arguments and checking normalization vectors.
    
    Example usage:
    python check_hic_norm.py --hic_path [hic_filepath] --chrom_list [comma-sep chrom list]
        --resolutions [comma-sep resolution list] --normalizations [comma-sep normalization methods]
    
    Note: This script currently only checks intra-chromosomal matrices, not inter-chromosomal or genome-wide matrices.
    """
    parser = argparse.ArgumentParser(description="Check existing normalizations for a .hic file")
    parser.add_argument("--hic_path", type=str, required=True,
                        help="Path to the .hic file")
    parser.add_argument("--chrom_list", type=comma_sep_strs, 
                        default=['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
                                 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 
                                 'chr20', 'chr21', 'chr22', 'chrX'],
                        help="Comma-separated list of chromosomes to check (e.g., chr1,chr2,chr3)")
    parser.add_argument("--res", type=comma_sep_ints, 
                        default=[2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000],
                        help="Comma-separated list of resolutions to check (e.g., 25000, 10000, 5000)")
    parser.add_argument("--norm", type=comma_sep_strs, 
                        default=['VC', 'VC_SQRT', 'KR', 'SCALE'],
                        help="Comma-separated list of normalization types to check (e.g., VC,VC_SQRT,KR,SCALE)")

    args = parser.parse_args()
    check_hic_norm(args.hic_path, args.chrom_list, args.res, args.norm)

if __name__ == "__main__":
    main()