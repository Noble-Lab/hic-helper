import sys
import os
#requirement to run: samtools, pysam
import pysam
from collections import defaultdict
def add_dict(dict,key):
    if key in dict:
        dict[key] += 1
    else:
        dict[key] = 1
    return dict
def calcuate_alignment_name_stat(sorted_bam_file,output_dir):
    samfile = pysam.AlignmentFile(sorted_bam_file,"rb")

    all_chroms = list(samfile.references)
    all_chroms_length = list(samfile.lengths)
    output_record = os.path.join(output_dir,"bam_alignments_count.tsv")
    final_stats = {}

    for k,chrom in enumerate(all_chroms):
        print("Processing %s" % chrom)
        chrom_length = all_chroms_length[k]
        alignments = samfile.fetch(chrom, 0, chrom_length)
        for aln in alignments:
            search_key = aln.query_name 
            add_dict(final_stats,search_key)
        print("accumulated %d alignments" % len(final_stats))
    total_alignments = len(final_stats)
    mapping_dict = {}
    for k,v in final_stats.items():
        if v in mapping_dict:
            mapping_dict[v] += 1
        else:
            mapping_dict[v] = 1
    with open(output_record,"w") as f:
        f.write("Total_alignments\t%d\n" % total_alignments)
        f.write("#Mapping Count\t#Alignment\n")
        for k,v in mapping_dict.items():
            f.write("%d\t%d\n" % (k,v))


if __name__ == '__main__':
    if len(sys.argv) !=5:
        print("Usage: python raw_quality.py [input.bam] [output_dir] [number_cpu] [mode]")
        print("input.bam: the input bam file")
        print("output_dir: the output directory")
        print("number_cpu: the number of cpu used")
        print("mode: 0: for unsorted bam file, 1: for sorted bam file")
        sys.exit(1)
    run_mode = int(sys.argv[4])
    output_dir = os.path.abspath(sys.argv[2])
    os.makedirs(output_dir, exist_ok=True)
    input_file = os.path.abspath(sys.argv[1])
    if not os.path.exists(input_file):
        print("Input file does not exist")
        sys.exit(1)
    # reference_fasta = os.path.abspath(sys.argv[2])
    # if not os.path.exists(reference_fasta):
    #     print("Reference fasta file does not exist")
    #     sys.exit(1)
    input_bam = os.path.join(output_dir,"input.bam")
    if os.path.exists(input_bam):
        print("Warning: input bam file already exists, will overwrite it.")
        os.remove(input_bam)
    os.symlink(input_file, input_bam)
    number_cpu = int(sys.argv[3])
    #sort first 
    os.chdir(output_dir)
    if run_mode==0:
        sorted_bam_file = os.path.join(output_dir,"input.sorted.bam")
        os.system("samtools sort -o %s -T %s --threads %d %s" % (sorted_bam_file,"input.tmp",number_cpu,input_bam))
        os.system("samtools index --threads %d %s" % (number_cpu,sorted_bam_file))
        if not os.path.exists(sorted_bam_file):
            print("sorted bam generation failed, please check if samtools is properly installed.")
            sys.exit(1)
    else:
        sorted_bam_file = input_bam 
        sorted_bam_index_file = input_file+".bai"
        if not os.path.exists(sorted_bam_index_file):
            print("Index file does not exist, please run with mode 0.")
            sys.exit(1)
        #link file to output_dir
        if os.path.exists(os.path.join(output_dir,"input.bam.bai")):
            print("Warning: input bam index file already exists, will overwrite it.")
            os.remove(os.path.join(output_dir,"input.bam.bai"))
        os.symlink(input_file+".bai",os.path.join(output_dir,"input.bam.bai"))
    
    #read sam file to get the alignment name and to the count of each alignment
    calcuate_alignment_name_stat(sorted_bam_file,output_dir)