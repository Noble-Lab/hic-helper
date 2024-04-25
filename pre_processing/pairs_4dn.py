"""
This script is used to convert pairs file to cool or hic files following 4DN's pipeline.
```
python3 pairs_4dn.py [input.pairs.gz] [chrom_size_file] [output_dir] [mode] [number_cpu] [max_memory] [resume_flag]
```
[input.pairs]: the input pairs.gz file <br>.
[chrom_size_file]: the chromosome size file. Example file: hg38.chrom.sizes for human genome build GRCh38 <br>
[output_dir]: the output directory. The output file will be named as 4DN.cool under this direcotry. <br>
[mode]: the mode of the conversion. 0: convert to cool file; 1：convert to hic file <br>
[number_cpu]: the number of cpu to use <br>
[max_memory]: the max memory to use (GB) <br>
[resume_flag]: 0: do not resume; 1: resume from the previously generated files. default should be 0. <br>
Recommended running with 8 cores and 64GB memory. <br>
"""

if __name__ == '__main__':
    import os 
    import sys
    if len(sys.argv) != 8:
        print('Usage: python3 pairs_4dn.py [input.pairs.gz] \
              [chrom_size_file] [output_dir] [mode] [number_cpu] [max_memory] [resume_flag]')
        print('input.pairs.gz: input pairs.gz file.')
        print('chrom_size_file: the chromosome size file. Example file: hg38.chrom.sizes for human genome build GRCh38')
        print('output_dir: the output directory. The output file will be named as 4DN.cool under this direcotry.')
        print('mode: the mode of the conversion. 0: convert to cool file; 1：convert to hic file')
        print('number_cpu: the number of cpu to use')
        print('max_memory: the max memory to use (GB)')
        print("resume_flag: 0: do not resume; 1: resume from the previously generated files. default should be 0.")
        print('Recommended running with 8 cores and 64GB memory.')
        sys.exit(1)
    
    
    prefix='4DN'
    
    input_pairs_file = os.path.abspath(sys.argv[1])
    chrom_size_file = os.path.abspath(sys.argv[2])
    output_dir = os.path.abspath(sys.argv[3])
    run_mode = int(sys.argv[4])
    number_cpu = int(sys.argv[5])
    max_memory = int(float(sys.argv[6]))
    resume_flag = int(sys.argv[7])

    os.makedirs(output_dir,exist_ok=True)
    #run on the same directory as the output dir, there are many temp files may be generated
    os.chdir(output_dir)

    script_path = os.path.dirname(os.path.realpath(__file__))
    
    cur_prefix = os.path.join(output_dir,prefix)

    # 5. merge pairs
    """
    run-merge-pairs.sh <output_prefix> <pairs1> <pairs2> [<pairs3> [...]]  
    # output_prefix : prefix of the output pairs file.
    # pairs1, pairs2, ... : input pairs files
    """
    pairs_merge_script_path = os.path.join(script_path,'run-merge-pairs.sh')
    merged_pairs_file_path = os.path.join(output_dir,'%s.pairs.gz'%prefix)
    if resume_flag == 1 and os.path.exists(merged_pairs_file_path) and os.path.getsize(merged_pairs_file_path) > 1000:
        print('Resume from the previously generated merged pairs files')
    else:
        os.system('bash %s %s %s' % (pairs_merge_script_path,cur_prefix,input_pairs_file))

    
    if not os.path.exists(merged_pairs_file_path):
        print('Error: merge pairs failed')
        sys.exit(1)

    # 6. convert pairs to cool or hic file
    
    ## 6.1 convert to cool file
    if run_mode == 0:
        """
        run-cooler.sh <input_pairs> <chromsize> <binsize> <ncores> <output_prefix> <max_split>
        # input_pairs : a pairs file
        # chromsize : a chromsize file
        # binsize : binsize in bp
        # ncores : number of cores to use
        # output_prefix : prefix of the output cool file
        # max_split : max_split argument for cooler (e.g. 2 which is default for cooler) 
        """
        cooler_script_path = os.path.join(script_path,'run-cooler.sh')
        cool_file_path = os.path.join(output_dir,'%s.cool'%prefix)
        if resume_flag == 1 and os.path.exists(cool_file_path) and os.path.getsize(cool_file_path) > 1000:
            print('Resume from the previously generated cool files')
        else:
            os.system('bash %s %s %s %d %d %s %d' % (cooler_script_path,merged_pairs_file_path,chrom_size_file,
                                                1000,number_cpu,cur_prefix,2))
        
        if not os.path.exists(cool_file_path):
            print('Error: convert pairs to cool file failed')
            sys.exit(1)

        # 7 convert to multi-resolution cool file
        """
        run-cool2multirescool.sh -i <input_cool> [-p <ncores>] [-o <output_prefix>] [-c <chunksize>] [-j] [-u custom_res] [-B]
        # input_cool : a (singe-res) cool file with the highest resolution you want in the multi-res cool file
        # -p ncores: number of cores to use (default: 1)
        # -o output_prefix: prefix of the output multires.cool file (default: out)
        # -c chunksize : chunksize argument of cooler (e.g. default: 10000000)
        # -j : juicer resolutions (default: use HiGlass resolutions)
        # -u custom_res : custom resolutions separated by commas (e.g. 100000,200000,500000). The minimun of this set must match min_res (-r).
        # -B : no balancing/normalization
        """
        cool2multirescool_script_path = os.path.join(script_path,'run-cool2multirescool.sh')
        command_line="bash %s \
        -i %s -p %d  -o %s -c 10000000 \
        -u 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000"%(cool2multirescool_script_path,cool_file_path,number_cpu,cur_prefix)
        multi_res_cool_file_path = os.path.join(output_dir,'%s.multires.cool'%prefix)
        if resume_flag == 1 and os.path.exists(multi_res_cool_file_path) and os.path.getsize(multi_res_cool_file_path) > 1000:
            print('Resume from the previously generated multi-res cool files')
        else:
            os.system(command_line)
        
        if not os.path.exists(multi_res_cool_file_path):
            print('Error: convert to multi-resolution cool file failed')
            sys.exit(1)
        # 8. add norm files from hic to cooler
        """
        run-add-hicnormvector-to-mcool.sh <input_hic> <input_mcool> <outdir>
        # input_hic : a hic file
        # input_mcool : an mcool file
        # outdir : output directory
        """
        print("Please get the cool file in %s"%multi_res_cool_file_path)
        print("If you want to add norm files from hic to cooler, please run the following command:")
        print("run-add-hicnormvector-to-mcool.sh <input_hic> %s %s"%(multi_res_cool_file_path,output_dir))
        print("You can run mode 1 to get the .hic file.")
    elif run_mode==1:
        print("convert to hic file")
        ## 6.2 convert to hic file
        ### adds juicer fragment file to the pairs file
        """
        run-addfrag2pairs.sh <input_pairs> <restriction_site_file> <output_prefix>
        # input_pairs : a gzipped pairs file (.pairs.gz) with its pairix index (.px2)
        # restriction_site_file : a text file containing positions of restriction enzyme sites, separated by space, one chromosome per line (Juicer style).
        # output prefix: prefix of the output pairs file
        """
        #if you know restriction file, you can input here, if not do not need to change this command
        addfrag2pairs_script_path = os.path.join(script_path,'run-addfrag2pairs.sh')
        command_line="bash %s  -i %s -o %s "%(addfrag2pairs_script_path,merged_pairs_file_path,cur_prefix)
        refined_pairs_file_path = os.path.join(output_dir,'%s.ff.pairs.gz'%prefix)
        if resume_flag == 1 and os.path.exists(refined_pairs_file_path) and os.path.getsize(refined_pairs_file_path) > 1000:
            print('Resume from the previously generated refined pairs files')
        else:
            os.system(command_line)
        
        if not os.path.exists(refined_pairs_file_path):
            print('Error: add juicer fragment file to the pairs file failed')
            sys.exit(1)
        ### convert pairs to hic file
        """
        run-juicebox-pre.sh -i <input_pairs> -c <chromsize_file> [-o <output_prefix>] [-r <min_res>] [-g] [-u custom_res] [-m <maxmem>] [-q mapqfilter] [-B]
        # -i input_pairs : a gzipped pairs file (.pairs.gz) with its pairix index (.px2), preferably containing frag information.
        # -c chromsize_file : a chromsize file
        # -o output prefix: prefix of the output hic file
        # -r min_res : minimum resolution for whole-genome normalization (e.g. 5000)
        # -g : higlass-compatible : if this flag is used, zoom levels are set in a Hi-Glass compatible way, if not, default juicebox zoom levels.
        # -u custom_res : custom resolutions separated by commas (e.g. 100000,200000,500000). The minimun of this set must match min_res (-r).
        # -m maxmem : java max mem (e.g. 14g)
        # -q mapqfilter : mapq filter (e.g. 30, default 0)
        # -n : normalization only : if this flag is used, binning is skipped.
        # -B : no balancing/normalization
        """
        juicer_tool_jar_path = os.path.join(script_path,'juicer_tools.jar')
        code_dir = os.path.dirname(juicer_tool_jar_path)
        root_path = os.getcwd()
        os.chdir(code_dir)
        juicebox_pre_script_path = os.path.join(script_path,'run-juicebox-pre.sh')
        #-q threshold can be 30 for some settings, please update this if you need
        command_line ="bash %s -s juicer_tools.jar -i %s -c %s -o %s -r 1000 -m %dg \
            -q 0 -u 1000,2000,5000,10000,25000,50000,100000,250000,\
        500000,1000000,2500000,5000000,10000000"%(juicebox_pre_script_path,
                                                refined_pairs_file_path,chrom_size_file,cur_prefix,
                                              max_memory)
        hic_file_path = os.path.join(output_dir,'%s.hic'%prefix)
        if resume_flag == 1 and os.path.exists(hic_file_path) and os.path.getsize(hic_file_path) > 1000:
            print('Resume from the previously generated hic files')
        else:
            os.system(command_line)

        os.chdir(root_path)
        if not os.path.exists(hic_file_path):
            print('Error: convert pairs to hic file failed')
            sys.exit(1)
        print("Please get the hic file in %s"%hic_file_path)
   
    else:
        print('Error: invalid mode, please input 0 or 1 for [mode]')
        sys.exit(1)
    print("bam 4DN pipeline finished! Enjoy!")
