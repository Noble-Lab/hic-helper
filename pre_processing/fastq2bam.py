# Path: pre_processing/fastq2bam.py
"""
```
python3 fastq2bam.py [fastq_file1] [fastq_file2] [bwa_index_file] [output_dir]
```
[fastq_file1]: the first fastq file. <br>
[fastq_file2]: the second fastq file. <br>
[bwa_index_file]: the bwa index file. Example file: 4DNFIZQZ39L9.bwaIndex.tgz for human genome build GRCh38 <br>
[output_dir]: the output directory. The output file will be named as 4DN.sorted.bam under this direcotry. <br>
"""

"""
This tar-gzipped file includes the bwa index for the human genome build GRCh38. We use the UCSC hg38 version of the genome, corresponsing to GRCh38/GCA_000001405.15 including the 25 assembled chromosomes (1-22, X, Y, M), the 127 unplaced contigs, the 42 unlocalized contigs, but excluding the 261 alternative haplotypes. The assembly also includes the Epstein-Barr virus (EBV) genome. This is the same reference as used by the ENCODE Consortium for data processing: https://www.encodeproject.org/data-standards/reference-sequences/ .
Download Link: https://data.4dnucleome.org/files-reference/4DNFIZQZ39L9/@@download/4DNFIZQZ39L9.bwaIndex.tgz
You must have an account to download this link, see 4DN Data Portal for more information.


Always keep the HiCPipeline from 4DN in mind: https://github.com/4dn-dcic/docker-4dn-hic/blob/a0d5318c07793cb7d841387d13b026dedc7f0496/HiCPipeline.md
"""

if __name__ == '__main__':
    import os 
    import sys
    if len(sys.argv) != 5:
        print('Usage: python3 fastq2bam.py [fastq_file1] [fastq_file2] [bwa_index_file] [output_dir]')
        print('fastq_file1: the first fastq file')
        print('fastq_file2: the second fastq file')
        print('bwa_index_file: the bwa index file. Example file: 4DNFIZQZ39L9.bwaIndex.tgz for human genome build GRCh38')
        print('output_dir: the output directory. The output file will be named as 4DN.sorted.bam under this direcotry.')
        sys.exit(1)
    number_cpu = 8
    prefix='4DN'
    fastq_file1 = os.path.abspath(sys.argv[1])
    fastq_file2 = os.path.abspath(sys.argv[2])
    bwa_index_file = os.path.abspath(sys.argv[3])
    output_dir = os.path.abspath(sys.argv[4])
    os.makedirs(output_dir,exist_ok=True)
    script_path = os.path.dirname(os.path.realpath(__file__))
    bwa_mem_script_path = os.path.join(script_path,'run_bwa_mem.sh')
    os.system('bash %s %s %s %s %s %s %d' % (bwa_mem_script_path,fastq_file1,fastq_file2,
                                          bwa_index_file,output_dir,prefix,number_cpu))
    gen_file = os.path.join(output_dir,'%s.bam'%prefix)
    if not os.path.exists(gen_file):
        print('Error: fastq to bam conversion failed')
        sys.exit(1)

    # sort the bam file
    sorted_bam_file = os.path.join(output_dir,'%s.sorted.bam'%prefix)
    tmp_bam_file = os.path.join(output_dir,'%s.tmp.bam'%prefix)
    os.system("samtools sort -o %s -T %s %s" % (sorted_bam_file,tmp_bam_file,gen_file))
    os.system("samtools index %s" % sorted_bam_file)
    os.remove(tmp_bam_file)
    print('The converted bam file: %s. Enjoy!' % sorted_bam_file)
    
    