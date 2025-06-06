
import sys
import os

def access_count(pair_dict,key):
    if key in pair_dict:
        count = pair_dict[key]
    else:
        count = 0
    return count

def int2scientific(num):
    if num == 0:
        return "0"
    output_str = str.format("{:.3e}", num)
    return output_str
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

    total_ru_count = access_count(map_dict,"pair_types/RU")
    total_ur_count = access_count(map_dict,"pair_types/UR")
    total_uu_count = access_count(map_dict,"pair_types/UU")
    with open(output_report, "w") as f:
        f.write("Unmapped sequences: \t" + int2scientific(total_unmapped) + "(%.3f%%)\n"%(total_unmapped/total_counts*100))
        f.write("Singleton sequences: \t" + int2scientific(total_singleton) + "(%.3f%%)\n"%(total_singleton/total_counts*100))
        f.write("Multimapped sequences: \t" + int2scientific(total_multimapped) + "(%.3f%%)\n"%(total_multimapped/total_counts*100))
        f.write("Duplicate sequences: \t" + int2scientific(total_duplicate) + "(%.3f%%)\n"%(total_duplicate/total_counts*100))
        f.write("Unique sequences: \t" + int2scientific(final_remain) + "(%.3f%%)\n"%(final_remain/total_counts*100))
        f.write("Total sequences: \t" + int2scientific(total_counts) + "(%.3f%%)\n"%100)
        f.write("-"*50 + "\n")
        f.write("Detailed Unique Sequences Information:\n")
        f.write("RU sequences: \t" + int2scientific(total_ru_count) + "(%.3f%%)\n"%(total_ru_count/final_remain*100))
        f.write("UR sequences: \t" + int2scientific(total_ur_count) + "(%.3f%%)\n"%(total_ur_count/final_remain*100))
        f.write("UU sequences: \t" + int2scientific(total_uu_count) + "(%.3f%%)\n"%(total_uu_count/final_remain*100))
        f.write("For detailed definition of RU/UR/UU, please see https://pairtools.readthedocs.io/en/latest/formats.html#pair-types\n")


    return output_report
"""
This script is used to analyze the quality of the pairs file.
```
python3 pairs_quality.py [input.pairs.gz] [output_dir] [number_cpu]
```
[input.pairs.gz]: the input pairs gz file from 4DN pipeline (*.marked.sam.pairs.gz).<br>
[output_dir]: the output directory<br>
[number_cpu]: the number of cpu used<br>
The script will generate a report file in the output directory.<br>
The report file will contain the following information:<br>
1. Unmapped sequences: the number of unmapped sequences and the percentage.<br>
2. Singleton sequences: the number of singleton sequences and the percentage.<br>
3. Multimapped sequences: the number of multimapped sequences and the percentage.<br>
4. Duplicate sequences: the number of duplicate sequences and the percentage.<br>
5. Unique sequences: the number of unique sequences and the percentage.<br>
6. Total sequences: the total number of sequences.<br>
7. Detailed Unique Sequences Information:<br>
    RU sequences: the number of RU sequences and the percentage.<br>
    UR sequences: the number of UR sequences and the percentage.<br>
    UU sequences: the number of UU sequences and the percentage.<br>
    For detailed definition of RU/UR/UU, please see https://pairtools.readthedocs.io/en/latest/formats.html#pair-types<br>
"""


if __name__ == '__main__':
    if len(sys.argv) !=4:
        #4DN.marked.sam.pairs.gz
        print("Usage: python pairs_quality.py [input.marked.sam.pairs.gz] [output_dir] [number_cpu]")
        print("input.pairs.gz: the input pairs gz file from 4DN pipeline.")
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
    print(f"Output stat is saved in {output_report}")