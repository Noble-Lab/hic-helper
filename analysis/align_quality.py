import sys
import os
#requirement to run: samtools, pysam
import pysam
from collections import defaultdict
def cigar_pass_checking(current_cigar):
    """
    cigar operation meaning
    operations = {
    0: 'M',  # match or mismatch
    1: 'I',  # insertion
    2: 'D',  # deletion
    3: 'N',  # skipped region
    4: 'S',  # soft clipping
    5: 'H',  # hard clipping
    6: 'P',  # padding
    7: '=',  # sequence match
    8: 'X'   # sequence mismatch
    """
    #soft clipping:  Number of reads where there is some softclipping at some point in the read's alignment
    if current_cigar is None:
        return -1
    elif len(current_cigar) == 0:
        return -1
    fail_flag_list = [1,2,3,4,5,6,8]
    for current_cigar_op in current_cigar:
        if current_cigar_op[0] in fail_flag_list:
            return current_cigar_op[0]
    return 0
def md_flag_checking(aln):
    """
    MD tag is a string that describes the differences between the read sequence and the reference sequence.
    The MD tag is used to identify SNPs and indels.
    """
    try:
        md_tag = aln.get_tag("MD")
    except KeyError:
        md_tag = None
    if md_tag is None:
        return -1
    if md_tag.isnumeric():
        #only no snp change can be numeric
        return 1
    else:
        return 0
def add_dict(dict,key):
    if key in dict:
        dict[key] += 1
    else:
        dict[key] = 1
    return dict
def calculate_chrom_stat(alignments,min_mapq=0):
    alignments = [a for a in alignments]
    count_all = len(alignments)
    #remove unmapped alignment
    alignments = [a for a in alignments if not a.is_unmapped]
    count_mapped = len(alignments)
    count_unmapped_remove = count_all - count_mapped
    print("Chrom total number of sequences: %d" % count_all)
    print("Chrom total number of mapped sequences: %d" % count_mapped)
    #remove low map quality alignment
    alignments = [a for a in alignments if a.mapping_quality> min_mapq]
    count_mapq = len(alignments)
    count_mapq_remove = count_mapped - count_mapq
    print("Chrom total number of sequences with mapq > %d: %d" % (min_mapq,count_mapq))
    #remove duplicate alignment
    alignments = [a for a in alignments if not a.is_duplicate]
    count_nodup = len(alignments)
    count_dup_remove = count_mapq - count_nodup
    print("Chrom total number of sequences without duplicate: %d" % count_nodup)

    #remove mapped multiple times
    alignments = [a for a in alignments if not a.is_secondary]
    count_primary = len(alignments)
    count_secondary_remove = count_nodup - count_primary
    print("Chrom total number of sequences with only primary alignment: %d" % count_primary)

    #count singletons
    count_singletons = 0
    final_alignments =[]
    count_singletons ={}
    for aln in alignments:
        """
        It maps to a unique position in the genome.
        Its mate read does not map to the expected location or does not map at all.
        It does not form a valid Hi-C pair with any other read.

        Criteria for Unique Alignment
        Single Mapping Location: A read must map to one and only one location in the reference genome.
        High Mapping Quality: The alignment must have a high mapping quality score, indicating confidence in the alignment. This is often represented by the MAPQ score in BAM files.
        No Secondary Alignments: Reads with secondary alignments (alternative mappings to other locations) are generally not considered unique alignments.
        """
        
        if aln.mate_is_unmapped:
            add_dict(count_singletons,"mate_unmapped")
        elif aln.reference_id != aln.next_reference_id:
            """
            reference_id: the reference sequence number as defined in the header
            next_reference_id: the reference id of the mate/next read.
            """
            add_dict(count_singletons,"mate2other_chrom")
            #count_singletons += 1
        elif aln.is_reverse == aln.mate_is_reverse:
            add_dict(count_singletons,"mate_same_strand")
        elif (
            # mapped to reverse strand but leftmost
            (aln.is_reverse and aln.template_length > 0)
            # mapped to fwd strand but rightmost
            or (not aln.is_reverse and aln.template_length < 0)
        ):
            #reads_faceaway - Number of reads where the read and its mate are mapped facing away from each other.
            #count_singletons += 1
            add_dict(count_singletons,"reads_faceaway")
        else:
            final_alignments.append(aln)
    print("possible singleton categories: ",count_singletons)
    count_remain = len(final_alignments)
    print("Chrom total number of sequences with mate correctly mapped: %d" % count_remain)
    count_singletons = count_primary - count_remain
    alignments = final_alignments
    count_other = {}
    final_alignments = []
    for aln in alignments:
        cigar_flag = cigar_pass_checking(aln.cigar)
        md_flag = md_flag_checking(aln)
        if cigar_flag!=0:
            if 'cigar_%d'%cigar_flag not in count_other:
                count_other['cigar_%d'%cigar_flag] = 1
            else:
                count_other['cigar_%d'%cigar_flag] += 1
            #1,2,3,4,5,6,8 flags are not allowed
        # """
        # cigar operation meaning
        # operations = {
        # 0: 'M',  # match or mismatch
        # 1: 'I',  # insertion
        # 2: 'D',  # deletion
        # 3: 'N',  # skipped region
        # 4: 'S',  # soft clipping
        # 5: 'H',  # hard clipping
        # 6: 'P',  # padding
        # 7: '=',  # sequence match
        # 8: 'X'   # sequence mismatch
        # """
        elif aln.is_qcfail:
            if "qcfail" not in count_other:
                count_other["qcfail"] = 1
            else:
                count_other["qcfail"] += 1
        elif not aln.is_paired:
            if 'unpaired' not in count_other:
                count_other['unpaired'] = 1
            else:
                count_other['unpaired'] += 1
        elif md_flag!=1:
            #snp checking
            if 'md_%d'%md_flag not in count_other:
                count_other['md_%d'%md_flag] = 1
            else:
                count_other['md_%d'%md_flag] += 1
            """
            Identifying SNPs:

            For each read, the script compares the query sequence (read sequence) to the reference sequence.
            Mismatches between the query base and the reference base are counted as potential SNPs.
            Identifying Indels:

            The CIGAR string of each read is processed to identify insertions (INS) and deletions (DEL).
            Indels are stored in the indels dictionary within the variants dictionary.
            """

        else:
            final_alignments.append(aln)
    count_proper = len(final_alignments)
    print("Detailed removed in others: ",count_other)
    print("Chrom total number of sequences with proper cigar and md tag: %d" % count_proper)
    alignments = final_alignments
    count_other_total = sum([count_other[key] for key in count_other])
    stats = {"proper":count_proper,
            "unmapped":count_unmapped_remove,
            "low quality (mapq)":count_mapq_remove,
            "duplicate":count_dup_remove,
            "Multimapped":count_secondary_remove,
            "singleton":count_singletons,
            "other":count_other_total,
            "all":count_all,}
    return stats



def calculate_stat(sorted_bam_file,output_dir):
    #open bam file
    """
    The calculate_stat function takes a sorted bam file and an output directory as input.
    It then opens the bam file using pysam, and creates a list of all chromosomes in the reference genome.
    The function then iterates through each chromosome, calculating statistics for each one individually. 
    These statistics are written to an output file called &quot;chromosome_stats&quot; in the specified output directory.
    
    :param sorted_bam_file: Open the bam file
    :param output_dir: Specify the output directory for the chromosome_stats
    :return: The stats of each chromosome in the bam file
    :doc-author: Trelent
    """
    samfile = pysam.AlignmentFile(sorted_bam_file,"rb")

    all_chroms = list(samfile.references)
    all_chroms_length = list(samfile.lengths)
    output_record = os.path.join(output_dir,"chromosome_stats.tsv")
    final_stats = defaultdict(list)
    with open(output_record,"w") as f:
        f.write("chromosome\tproper\tunmapped\tlow quality (mapq)\tduplicate\tMultimapped\tsingleton\tother\tall\n")
    for k,chrom in enumerate(all_chroms):
        print("Processing %s" % chrom)
        chrom_length = all_chroms_length[k]
        alignments = samfile.fetch(chrom, 0, chrom_length)
        stats = calculate_chrom_stat(alignments)
        for key in stats:
            final_stats[key].append(stats[key])
        with open(output_record,"a+") as f:
            f.write("%s\t%d\t%d\t%d\t%d \
                    \t%d\t%d\t%d\t%d\n" % (chrom,stats["proper"],
                                            stats["unmapped"],
                                            stats["low quality (mapq)"],
                                            stats["duplicate"],
                                            stats["Multimapped"],
                                            stats["singleton"],
                                            stats["other"],
                                            stats["all"]))
            #calculate the percentage in the chromsome report
            chrom_stats= {}
            for key in stats:
                if key == "all":
                    continue
                chrom_stats[key] = stats[key]/stats["all"]*100
            chrom_stats["all"] = 100
            f.write("%s\t%d\t%d\t%d\t%d \
                    \t%d\t%d\t%d\t%d\n" % (chrom,chrom_stats["proper"],
                                            chrom_stats["unmapped"],
                                            chrom_stats["low quality (mapq)"],
                                            chrom_stats["duplicate"],
                                            chrom_stats["Multimapped"],
                                            chrom_stats["singleton"],
                                            chrom_stats["other"],
                                            chrom_stats["all"]))
    #calculate total stats, with percentage level
    total_stats = defaultdict(int)
    for key in final_stats:
        total_stats[key] = sum(final_stats[key])
    
    percent_stats = {}
    for key in total_stats:
        if key == "all":
            continue
        percent_stats[key] = total_stats[key]/total_stats["all"]*100
    
    final_record = os.path.join(output_dir,"output_report.txt")
    with open(final_record,"w") as wfile:
        wfile.write("Total stats for all chromosomes\n")
        #write the total number of sequences
        wfile.write("Total number of sequences: %d\n" % total_stats["all"])
        for key in total_stats:
            if key == "all":
                continue
            wfile.write("%s: %d (%.6f%%)\n" % (key,total_stats[key],percent_stats[key]))
    return final_record





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

    #start to do quality check
    output_report=calculate_stat(sorted_bam_file,output_dir)
    print("Please check the final report at %s"%output_report)