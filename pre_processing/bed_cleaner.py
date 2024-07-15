import sys
import os 


def bed_cleaner(input_bed, output_bed):
    with open(input_bed, 'r') as f:
        lines = f.readlines()
    lines = [line.strip().split() for line in lines]
    lines = [[line[0], int(line[1]), int(line[2])] for line in lines]
    lines = sorted(lines, key=lambda x: (x[0], x[1]))
    output_lines = []
    current_chrom = None
    current_start = None
    current_end = None
    for line in lines:
        chrom = line[0]
        start = line[1]
        end = line[2]
        if current_chrom is None:
            current_chrom = chrom
            current_start = start
            current_end = end
        elif current_chrom == chrom:
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                output_lines.append([current_chrom, current_start, current_end])
                current_start = start
                current_end = end
        else:
            output_lines.append([current_chrom, current_start, current_end])
            current_chrom = chrom
            current_start = start
            current_end = end
    output_lines.append([current_chrom, current_start, current_end])
    with open(output_bed, 'w') as f:
        for line in output_lines:
            f.write('\t'.join([str(x) for x in line]) + '\n')
    print("Finished cleaning bed file! Output bed file:", output_bed)

"""
This script merges overlapping regions from bed file.
```
python3 bed_cleaner.py [input_bed] [output_bed]
```
[input_bed]: the input bed file. <br>
[output_bed]: the output bed file without overlapping regions. <br>


"""

#merge overlapping regions from bed file
if __name__ == '__main__':
    if len(sys.argv)!=3:
        print("Usage: python bed_cleaner.py [input_bed] [output_bed]")
        print("input_bed: the input bed file")
        print("output_bed: the output bed file without overlapping regions")
        sys.exit(1)
    input_bed = os.path.abspath(sys.argv[1])
    output_bed = os.path.abspath(sys.argv[2])
    output_dir = os.path.dirname(output_bed)
    os.makedirs(output_dir, exist_ok=True)
    #clean bed file
    bed_cleaner(input_bed, output_bed)