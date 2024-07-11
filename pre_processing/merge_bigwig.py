
import sys
import os 
import pyBigWig
import numpy as np
def parse_genome_size(genome_size_file):
    chrom_size = {}
    with open(genome_size_file) as f:
        for line in f:
            chrom, size = line.strip("\n").split()
            chrom_size[chrom] = int(size)
    return chrom_size


def merge_bigwig(input_dir, output_bw, refer_genome_size):
    all_bw_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) 
                 if f.endswith('.bigwig') or f.endswith('.bw') or f.endswith('.bigWig') or f.endswith('.bigWig')]
    #read genome info
    chrom_size = parse_genome_size(refer_genome_size)
    
    #read bigwig file
    bw_handler_list=[]
    for bw_file in all_bw_files:
        bw_handler = pyBigWig.open(bw_file)
        bw_handler_list.append(bw_handler)
    
    with pyBigWig.open(output_bw, "w") as bw:
        chromosomes = [(key,chrom_size[key]) for key in chrom_size]
        # Add chromosome information to the BigWig file
        bw.addHeader(chromosomes)
        for chrom in chrom_size:
            print("start chrom",chrom,"size",chrom_size[chrom])
            cur_chrom_size =chrom_size[chrom]
            value_list=[]
            try:
                for bw_handle in bw_handler_list:
                    tmp_value = bw_handle.values(chrom, 0, cur_chrom_size, numpy=True)
                    value_list.append(tmp_value)
            except:
                print("error",chrom)
                continue
            values = value_list[0]
            for i in range(1,len(value_list)):
                values = values + value_list[i]
            pos_embedding = np.arange(0, len(values), dtype=int)
            pos_embedding = pos_embedding.tolist()
            values = values.tolist()
            new_pos_embedding = []
            new_values = []
            for i in range(len(pos_embedding)):
                if values[i] != 0:
                    new_pos_embedding.append(pos_embedding[i])
                    new_values.append(values[i])
            # bw.addEntries("chr1", [500, 600, 635], values=[-2.0, 150.0, 25.0], span=20)
            bw.addEntries(str(chrom), new_pos_embedding, values=new_values, span=1)
            print("finished",chrom)
    for bw_handle in bw_handler_list:
        bw_handle.close()
    return output_bw

"""
This script is used to merge bigwig files into one bigwig file.
```
python3 merge_bigwig.py [input_dir] [output_bw] [refer_genome.sizes]
```
[input_dir]: the directory containing all the bigwig files. <br>
[output_bw]: the output bigwig file. <br>
[refer_genome.sizes]: the chromosome sizes of the reference genome. <br>

"""


if __name__ == '__main__':
    
    if len(sys.argv) != 4:
        print("Usage: python merge_bigwig.py [input_dir] [output_bw] [refer_genome.sizes]")
        print("input_dir: the directory containing all the bigwig files")
        print("output_bw: the output bigwig file")
        print("refer_genome.sizes: the chromosome sizes of the reference genome")
        sys.exit(1)
    input_dir = os.path.abspath(sys.argv[1])
    output_bw = os.path.abspath(sys.argv[2])
    output_dir = os.path.dirname(output_bw)
    os.makedirs(output_dir, exist_ok=True)
    refer_genome_size = os.path.abspath(sys.argv[3])

    output_bw = merge_bigwig(input_dir, output_bw, refer_genome_size)
    print("Finished merging bigwig files saved to %s" % output_bw)





