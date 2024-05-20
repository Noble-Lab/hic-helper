# Path: pre_processing/fastq2bam.py
"""
```
python3 fastq2bam.py [fastq_file1] [fastq_file2] [ref.fa] [output_dir]
```
[fastq_file1]: the first fastq file. <br>
[fastq_file2]: the second fastq file. <br>
[ref.fa]: the reference genome file. You can download the reference genome file from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/ for human. <br>
[output_dir]: the output directory. The output file will be named as 4DN.sorted.bam under this direcotry. <br>
You can also automatically download the human and mouse reference genome by run set_up.sh. <br>
"""

"""

Always keep the HiCPipeline from 4DN in mind: https://github.com/4dn-dcic/docker-4dn-hic/blob/a0d5318c07793cb7d841387d13b026dedc7f0496/HiCPipeline.md
"""

if __name__ == '__main__':
    import os 
    import sys
    if len(sys.argv) != 5:
        print('Usage: python3 fastq2bam.py [fastq_file1] [fastq_file2] [refer_genome_file] [output_dir]')
        print('fastq_file1: the first fastq file')
        print('fastq_file2: the second fastq file')
        print('refer_genome_file: the reference genome file. You can download the reference genome file from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/ for human.')
        print('output_dir: the output directory. The output file will be named as 4DN.sorted.bam under this direcotry.')
        sys.exit(1)
    number_cpu = 8
    prefix='4DN'
    fastq_file1 = os.path.abspath(sys.argv[1])
    fastq_file2 = os.path.abspath(sys.argv[2])
    refer_fa_file = os.path.abspath(sys.argv[3])
    output_dir = os.path.abspath(sys.argv[4])
    os.makedirs(output_dir,exist_ok=True)
    script_path = os.path.dirname(os.path.realpath(__file__))
    bwa_mem_script_path = os.path.join(script_path,'run-bwa-mem-seq.sh')
    os.system('bash %s %s %s %s %s %s %d' % (bwa_mem_script_path,fastq_file1,fastq_file2,
                                          refer_fa_file,output_dir,prefix,number_cpu))
    gen_file = os.path.join(output_dir,'%s.bam'%prefix)
    if not os.path.exists(gen_file):
        print('Error: fastq to bam conversion failed')
        sys.exit(1)

    # sort the bam file
    sorted_bam_file = os.path.join(output_dir,'%s.sorted.bam'%prefix)
    tmp_bam_file = os.path.join(output_dir,'%s.tmp.bam'%prefix)
    os.system("samtools sort -o %s -T %s %s" % (sorted_bam_file,tmp_bam_file,gen_file))
    #can add --threads 32 to use 32 threads to accelerate the process
    os.system("samtools index %s" % sorted_bam_file)
    #samtools index --threads 32 [input.bam]
    #accelerate command
    os.remove(tmp_bam_file)
    print('The converted bam file: %s. Enjoy!' % sorted_bam_file)
    
    