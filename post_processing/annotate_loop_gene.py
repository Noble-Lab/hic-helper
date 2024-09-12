import sys
import os
from collections import defaultdict

def extract_loc(pred_detect_path):
    record_list=[]
    with open(pred_detect_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip().split()
            chrom1 = line[0]
            try:
                start1 = int(float(line[1]))
                end1 = int(float(line[2]))
            except:
                continue
            chrom2 = line[3]
            try:
                start2 = int(float(line[4]))
                end2 = int(float(line[5]))
            except:
                continue
            if "chr" not in chrom1:
                chrom1 = "chr"+chrom1
            
            if "chr" not in chrom2: 
                chrom2 = "chr"+chrom2

            record_list.append([chrom1, start1, end1, chrom2, start2, end2])
    print("Extracted",len(record_list)," loop detection records")
    return record_list

def extract_gene_dict(gene_annotation):
    """
    This function is for extract gene information from gene annotation file
    """
    gene_dict = defaultdict(list)
    with open(gene_annotation) as f:
        for line in f:
            if line.startswith("#"):
                continue
            split_res = line.split("\t")
            chrom = split_res[0]
            start_index = int(split_res[3])
            end_index = int(split_res[4])
            gene_name = split_res[8].split(";")[0].split(" ")[1].replace("\"","")
            gene_dict[chrom].append([start_index, end_index, gene_name])
    #sort the gene list
    for chrom in gene_dict:
        gene_dict[chrom] = sorted(gene_dict[chrom], key=lambda x: x[0])
    return gene_dict
def reorganize_loop(loop_list):
    loop_dict = {}
    for loop in loop_list:
        chrom1, start1, end1, chrom2, start2, end2 = loop
        if chrom1!=chrom2:
            print("loop chromsome not match")
            sys.exit(1)
        if chrom1 not in loop_dict:
            loop_dict[chrom1] = []
        loop_dict[chrom1].append(loop)
    #sort the loop list
    for chrom in loop_dict:
        loop_dict[chrom] = sorted(loop_dict[chrom], key=lambda x: x[1])
    return loop_dict

def annotate_loop_gene(input_bed, gene_annotation, output_bed):
    """
    This function is for annotate the loop with gene information
    input_bed: the input bed file
    gene_annotation: the gene annotation file format:.gtf, like hg38.ncbiRefSeq.gtf
    output_bed: the output bed file
    """
    loop_list = extract_loc(input_bed)
    total_loop = len(loop_list)
    gene_dict = extract_gene_dict(gene_annotation)
    #reorganize loop_list to chrom-based loop
    loop_dict = reorganize_loop(loop_list)
    #check the loop_list add gene information
    output_list = []
    loop_label_count = 0
    for chrom in loop_dict:
        current_loop_list = loop_dict[chrom]
        current_gene_list = gene_dict[chrom]
        for loop in current_loop_list:
            chrom1, start1, end1, chrom2, start2, end2 = loop
            #make sure start1<start2, otherwise swap
            if start1>start2:
                start1, start2 = start2, start1
                end1, end2 = end2, end1
            close1_distance = 100000000
            close2_distance = 100000000
            close1_gene = ""
            close2_gene = ""
            inside_gene_list=[]
            for gene in current_gene_list:
                gene_start, gene_end, gene_name = gene
                #if the gene is in the loop
                if gene_start>=start1 and gene_end<=end1:
                    close1_distance = 0
                    close1_gene = gene_name
                    inside_gene_list.append(gene_name)
                    continue
                if gene_start>=start2 and gene_end<=end2:
                    close2_distance = 0
                    close2_gene = gene_name
                    inside_gene_list.append(gene_name)
                    break
                #check the distance
                cur_distance1 = min(abs(gene_start-start1), abs(gene_end-end1))
                cur_distance2 = min(abs(gene_start-start2), abs(gene_end-end2))
                if cur_distance1<close1_distance:
                    close1_distance = cur_distance1
                    close1_gene = gene_name
                if cur_distance2<close2_distance:
                    close2_distance = cur_distance2
                    close2_gene = gene_name
            tmp_list=[chrom1, start1, end1, chrom2, start2, end2]
            #[chrom1, start1, end1, chrom2, start2, end2, close1_gene+"-%d"%close1_distance, close2_gene+"-%d"%close2_distance])
            if len(inside_gene_list)>0:
                for gene_name in inside_gene_list:
                    tmp_list.append(gene_name+"-0")
            else:
                tmp_list.append(close1_gene+"-%d"%close1_distance)
                tmp_list.append(close2_gene+"-%d"%close2_distance)
            output_list.append(tmp_list)
            loop_label_count+=1
            if loop_label_count%100==0:
                print("Annotated",loop_label_count,"/%d loops"%total_loop)
                print("Example gene annotation:",output_list[-1])
    with open(output_bed, "w") as f:
        for record in output_list:
            f.write("\t".join([str(x) for x in record])+"\n")
"""
This script is for annotate the loop with gene information.
```
python3 annotate_loop_gene.py [input.bed] [gene_annotation] [output.bed]
```
input.bed: the input bed file that contains the loop information. <br>
gene_annotation: the gene annotation file format:.gtf, like hg38.ncbiRefSeq.gtf. <br>
output.bed: the output bed file that contains the annotated loop information. <br>
The last two columns in the output.bed file are the closest gene and the distance to the loop (corresponds to x and y). <br>

"""



if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python3 annotate_loop_gene.py [input.bed] [gene_annotation] [output.bed] ")
        print("This script is for annotate the loop with gene information")
        print("input.bed: the input bed file that contains the loop information")
        print("gene_annotation: the gene annotation file format:.gtf, like hg38.ncbiRefSeq.gtf")
        print("output.bed: the output bed file that contains the annotated loop information")
        sys.exit(1)
    input_bed = os.path.abspath(sys.argv[1])
    gene_annotation = os.path.abspath(sys.argv[2])
    output_bed = os.path.abspath(sys.argv[3])
    output_dir = os.path.dirname(output_bed)
    os.makedirs(output_dir, exist_ok=True)
    #read the input bed file get all records
    annotate_loop_gene(input_bed, gene_annotation, output_bed)