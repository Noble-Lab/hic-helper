import os 
import sys
import hicstraw

def get_chrom_list(hic):
    """
    Get the chromosome list from a hic file.
    Args:
        hic: the hic handler
    return:
        chrom_list: the chromosome list in the hic file.
    """
    chrom_list = []
    for chrom in hic.getChromosomes():
        print(chrom.name, chrom.length)
        if "all" in chrom.name.lower():
            continue
        chrom_list.append(chrom)
    return chrom_list
def read_chrom_data(chr1_name, chr2_name, normalization, hic_file, resolution):
    #chr1_name = chr1.name
    #chr2_name = chr2.name
    infos = []
    infos.append('observed')
    infos.append(normalization)
    infos.append(hic_file)
    infos.append(chr1_name)
    infos.append(chr2_name)
    infos.append('BP')
    infos.append(resolution)
    print(infos)
    row, col, val = [], [], []
    rets = hicstraw.straw(*infos)
    print('\tlen(rets): {:3e}'.format(len(rets)))
    for ret in rets:
        row.append((int)(ret.binX // resolution))
        col.append((int)(ret.binY // resolution))
        val.append(ret.counts)
    print('\tsum(val): {:3e}'.format(sum(val)))
    return row,col,val

def write_record_txt(chrom1, chrom2,row, col, val, output_txt):
    """
    Write the hic data to a txt file.
    Args:
        row: the row index of the hic data.
        col: the column index of the hic data.
        val: the value of the hic data.
        output_txt: the output txt file.
    return:
        None
    """
    with open(output_txt, 'a+') as wfile:
        for i in range(len(row)):
            wfile.write(f'{0} {chrom1} {int(row[i]*resolution+1)} {0} {0} {chrom2} {col[i]*resolution+1} {1} {val[i]:.2f}\n')
    
def merge_hic(juicer_tools,hic_file1, hic_file2, output_hic, resolution,refer_genome_name):
    """
    Merge two hic files to a new merged hic file with specified resolution.
    Args:
        juicer_tools: the juicer tools jar file path.
        hic_file1: the first hic file to be merged.
        hic_file2: the second hic file to be merged.
        output_hic: the output merged hic file.
        resolution: the resolution of the output merged hic file.
        refer_genome_name: the name of the reference genome
    return:
        None
    """
    hic1 = hicstraw.HiCFile(hic_file1)
    hic2 = hicstraw.HiCFile(hic_file2)
    chrom_list1 = get_chrom_list(hic1)
    chrom_list2 = get_chrom_list(hic2)
    chrom_name_list1 = [chrom.name for chrom in chrom_list1]
    chrom_name_list2 = [chrom.name for chrom in chrom_list2]
    chrom_list = list(set(chrom_name_list1)  | set(chrom_name_list2))
    if len(chrom_list)>len(chrom_name_list1):
        print("Warning: The chromosome list in the two hic files are not the same.")
        more_chrom = list(set(chrom_list) - set(chrom_name_list1))
        print("The following chromosomes are only in the second hic file.")
        print(more_chrom)
    if len(chrom_list)>len(chrom_name_list2):
        print("Warning: The chromosome list in the two hic files are not the same.")
        more_chrom = list(set(chrom_list) - set(chrom_name_list2))
        print("The following chromosomes are only in the first hic file.")
        print(more_chrom)

    resolution_list = hic1.getResolutions()
    if resolution not in resolution_list:
        print("Resolution not found in the 1st hic file, please choose from the following list:")
        print(resolution_list)
        exit()
    resolution_list = hic2.getResolutions()
    if resolution not in resolution_list:
        print("Resolution not found in the 2nd hic file, please choose from the following list:")
        print(resolution_list)
        exit()
    output_txt = output_hic.replace('.hic','.raw')
    with open(output_txt, 'w') as wfile:
        wfile.write("")
    # merge the two hic files
    for i in range(len(chrom_list)):
        for j in range(i,len(chrom_list)):
            row1, col1, val1 = read_chrom_data(chrom_list[i], chrom_list[j], "NONE", hic_file1, resolution)
            row2, col2, val2 = read_chrom_data(chrom_list[i], chrom_list[j], "NONE", hic_file2, resolution)
            row = row1 + row2
            col = col1 + col2
            val = val1 + val2
            write_record_txt(chrom_list[i], chrom_list[j], row, col, val, output_txt)
    # convert the txt file to hic file
    code_path = os.path.dirname(juicer_tools)
    root_path = os.getcwd()
    os.chdir(code_path)
    os.system(f'java -Xmx64g -Xmx64g -jar juicer_tools.jar pre -j 8 -d -r {resolution} {output_txt} {output_hic}  {refer_genome_name}')
    os.remove(output_txt)

    os.chdir(root_path)

"""
This script is to merge two hic files to a new merged hic file with specified resolution.
```
python3 merge_hic.py [hic_file1] [hic_file2] [output_hic] [resolution] [refer_genome]
```
[hic_file1]: the first hic file to be merged. <br>
[hic_file2]: the second hic file to be merged. <br>
[output_hic]: the output merged hic file. <br>
[resolution]: the resolution of the output merged hic file. <br>
[refer_genome]: the name of the reference genome. For example, hg38, hg19, mm10. <br>

"""

if __name__ == '__main__':
    
    if len(sys.argv) != 6:
        print("Usage: python3 merge_hic.py <hic_file1> <hic_file2> <output_hic> <resolution> <refer_genome>")
        print("This script is to merge two hic files to a new merged hic file with specified resolution.")
        print("[hic_file1]: the first hic file to be merged.")
        print("[hic_file2]: the second hic file to be merged.")
        print("[output_hic]: the output merged hic file.")
        print("[resolution]: the resolution of the output merged hic file.")
        print("[refer_genome]: the name of the reference genome. For example, hg38, hg19, mm10.")
        exit(0)
    hic_file1 = os.path.abspath(sys.argv[1])
    hic_file2 = os.path.abspath(sys.argv[2])
    output_hic = os.path.abspath(sys.argv[3])
    output_dir = os.path.dirname(output_hic)
    os.makedirs(output_dir, exist_ok=True)
    resolution = int(sys.argv[4])
    refer_genome_name = sys.argv[5] 
    script_dir = os.path.dirname(os.path.realpath(__file__))
    juicer_tools = os.path.join(script_dir, 'juicer_tools.jar')
    merge_hic(juicer_tools,hic_file1, hic_file2, output_hic, resolution,refer_genome_name)