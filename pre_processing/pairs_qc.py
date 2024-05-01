
"""
This script is used to perform quality control on the pairs file.
```
python3 pairs_qc.py [input.pairs.gz] [chrom_size_file] [output_dir] [enzyme]
```
- input.pairs.gz: input pairs file. <br>
- chrom_size_file: chrom size file. <br>
- output_dir: output directory. <br>
- enzyme: enzyme used for Hi-C experiment, either 4 or 6. <br>
"""

if __name__ == '__main__':
    import sys
    import os
    if len(sys.argv) !=5:
        print('Usage: python3 pairs_qc.py [input.pairs.gz] [chrom_size_file] [output_dir] [enzyme]')
        print("[input.pairs.gz]: input pairs file.")
        print("[chrom_size_file]: chrom size file.")
        print("[output_dir]: output directory.")
        print("[enzyme]: enzyme used for Hi-C experiment, either 4 or 6.")
        sys.exit(1)
    input_pairs = sys.argv[1]
    chromsize_file = sys.argv[2]
    output_dir = sys.argv[3]
    enzyme = sys.argv[4]    
    input_pairs = os.path.abspath(input_pairs)
    chromsize_file = os.path.abspath(chromsize_file)
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    script_dir = os.path.dirname(os.path.realpath(__file__))
    script_dir = script_dir.replace('post_processing','bin')
    script_dir = os.path.join(script_dir, "pairsqc")
    script_path = os.path.join(script_dir, 'pairsqc.py')
    output_prefix = os.path.join(output_dir,"4DN_QC")
    os.system(f'python3 {script_path} -p {input_pairs} -c {chromsize_file} -tP -s {output_prefix} -O {output_prefix}')
    #python3 $scriptdir/pairsqc.py -p $input_pairs -c $chromsize -tP -s $sample_name -O $sample_name

    script_path = os.path.join(script_dir, 'plot.r')
    #Rscript $scriptdir/plot.r $enzyme $sample_name\_report
    os.system(f'Rscript {script_path} {enzyme} {output_prefix}_report ')
    
    os.system(f'zip -r {output_prefix}_report.zip {output_prefix}_report')
    