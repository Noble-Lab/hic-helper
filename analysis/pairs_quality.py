
import sys
import os

def access_count(pair_dict,key):
    if key in pair_dict:
        count = pair_dict[key]
    else:
        count = 0
    return count
def parser_stats(output_stats_tsv, output_report):
    map_dict={}
    with open(output_stats_tsv, "r") as f:
        for line in f:
            line = line.strip()
            fields = line.split("\t")
            pair_type = fields[0]
            pair_count = float(fields[1])
            map_dict[pair_type] = pair_count
    total_counts = map_dict["total"]
    total_unmapped = map_dict["total_unmapped"]#including WW,NN,XX,NM,MM type
    total_multimapped = access_count(map_dict,"pair_types/MM")
    total_unmapped = total_unmapped-total_multimapped

    total_singleton = map_dict["total_single_sided_mapped"]#NU,NR,MU,MR
    mu_mr_count = access_count(map_dict,"pair_types/MU") + access_count(map_dict,"pair_types/MR")
    total_singleton = total_singleton - mu_mr_count
    total_multimapped = total_multimapped + mu_mr_count
    total_duplicate = map_dict["total_dups"]#UU
    final_remain = map_dict['total_nodups']#RR,RU,UR,
    #convert int for all total count
    total_unmapped = int(total_unmapped)
    total_singleton = int(total_singleton)
    total_duplicate = int(total_duplicate)
    total_multimapped = int(total_multimapped)
    final_remain = int(final_remain)
    with open(output_report, "w") as f:
        f.write("Unmapped sequences: \t" + str(total_unmapped) + "(%.5f%%)\n"%(total_unmapped/total_counts*100))
        f.write("Singleton sequences: \t" + str(total_singleton) + "(%.5f%%)\n"%(total_singleton/total_counts*100))
        f.write("Duplicate sequences: \t" + str(total_duplicate) + "(%.5f%%)\n"%(total_duplicate/total_counts*100))
        f.write("Multimapped sequences: \t" + str(total_multimapped) + "(%.5f%%)\n"%(total_multimapped/total_counts*100))
        f.write("Final remain sequences: \t" + str(final_remain) + "(%.5f%%)\n"%(final_remain/total_counts*100))
        f.write("Total sequences: \t" + str(total_counts) + "\t(%.5f%%)\n"%100)



if __name__ == '__main__':
    if len(sys.argv) !=4:
        #4DN.marked.sam.pairs.gz
        print("Usage: python raw_quality.py [input.marked.sam.pairs.gz] [output_dir] [number_cpu]")
        print("input.bam: the input bam file")
        print("output_dir: the output directory")
        print("number_cpu: the number of cpu used")
        sys.exit(1)
    #pairtools stats 4DN.ff.pairs -o 4DN.ff.stats.tsv --no-yaml --nproc-in 32 --nproc-out 32
    input_pairs = os.path.abspath(sys.argv[1])
    output_dir = os.path.abspath(sys.argv[2])
    os.makedirs(output_dir, exist_ok=True)
    number_cpu = int(sys.argv[3])
    prefix_name = "input"
    output_stats_tsv = os.path.join(output_dir, prefix_name + ".stats.tsv")
    if not os.path.exists(output_stats_tsv):
        command="pairtools stats " + input_pairs + " -o " + output_stats_tsv + " --no-yaml --nproc-in " + str(number_cpu) + " --nproc-out " + str(number_cpu)
        print(command)
        os.system(command)
    #analyze the stats
    output_report = os.path.join(output_dir, prefix_name + ".report.txt")
    parser_stats(output_stats_tsv, output_report)