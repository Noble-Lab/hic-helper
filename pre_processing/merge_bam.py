


import sys
import os 


def merge_bam(input_dir, output_bam):
    """
    input_dir: the directory containing all the bam files
    output_bam: the output merged bam file
    
    """
    output_unsorted_bam = output_bam.replace(".bam", ".unsorted.bam")
    command_line=f"samtools merge -o {output_unsorted_bam} "
    listfiles = [os.path.join(input_dir,x) for x in os.listdir(input_dir) if x.endswith(".bam")]
    for item in listfiles:
        command_line += item + " "
    os.system(command_line)

    #output_sorted_bam = output_bam.replace(".bam", ".sorted.bam")
    os.system(f"samtools sort -o {output_bam} {output_unsorted_bam}")
    os.system(f"samtools index {output_bam}")
    return output_bam

"""
This script is used to merge bam files into one bam file.
```
python3 merge_bam.py [input_dir] [output_bam]
```
[input_dir]: the directory containing all the bam files. <br>
[output_bam]: the output merged bam file. <br>

"""

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python merge_bam.py [input_dir] [output_bam]")
        print("input_dir: the directory containing all the bam files")
        print("output_bam: the output bam file")
        sys.exit(1)

    input_dir = os.path.abspath(sys.argv[1])
    output_bam = os.path.abspath(sys.argv[2])
    output_dir = os.path.dirname(output_bam)
    os.makedirs(output_dir, exist_ok=True)
    output_bam = merge_bam(input_dir, output_bam)
    print("Finished merging bigwig files saved to %s" % output_bam)