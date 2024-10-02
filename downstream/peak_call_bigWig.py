import os
import sys
import argparse
import pyBigWig
import numpy as np
from multiprocessing import Pool
import scipy.stats as stats
def _call_peak(signal_window,config_lambda,broad,current_index=0):
    """
    signal_window: numpy array
    config_lambda: float lambda from control sample
    broad: bool, specify if it is broad peak calling
    """
    max_value = np.max(signal_window)
    
    max_value =int(np.floor(max_value))
    p_value =1 - stats.poisson.cdf(max_value,config_lambda)
    print("p-value for the peak at", current_index, "is", p_value)
    if not broad:
        return p_value
    mean_value = np.mean(signal_window)
    mean_value = int(np.floor(mean_value))
    broad_pval = 1-stats.poisson.cdf(mean_value,config_lambda)
    return p_value,broad_pval
def filter_peak(peak_list, pvalue, qvalue,broad,broad_cutoff):
    tmp_peak_list = peak_list
    print("before filtering has", len(tmp_peak_list), "peaks")
    #sort the peak list according to the p-value, small to large
    tmp_peak_list = sorted(tmp_peak_list,key=lambda x:x[1])
    if pvalue is not None:
        tmp_peak_list = [x for x in tmp_peak_list if x[1] <= pvalue]
        if broad:
            tmp_peak_list = [x for x in tmp_peak_list if x[2] <= broad_cutoff]
    elif qvalue is not None:
        refine_peak_list = []
        for i in range(len(tmp_peak_list)):
            if not broad:
                cur_index, cur_pval = tmp_peak_list[i]
                cur_qval = cur_pval * len(tmp_peak_list) / (i+1)
                refine_peak_list.append([cur_index,cur_qval])
            else:
                cur_index, cur_pval,cur_broad_pval = tmp_peak_list[i]
                cur_qval = cur_pval * len(tmp_peak_list) / (i+1)
                refine_peak_list.append([cur_index,cur_qval,cur_broad_pval])
        tmp_peak_list = [x for x in refine_peak_list if x[1] <= qvalue]
        #if broad peak calling, we also need to check the broad q-value
        if broad:
            tmp_peak_list = sorted(tmp_peak_list,key=lambda x:x[2])
            refine_peak_list = []
            for i in range(len(tmp_peak_list)):
                cur_index, cur_qval,cur_broad_pval = tmp_peak_list[i]
                cur_broad_qval = cur_broad_pval * len(tmp_peak_list) / (i+1)
                refine_peak_list.append([cur_index,cur_qval,cur_broad_qval])
            tmp_peak_list = [x for x in refine_peak_list if x[2] <= broad_cutoff]
    else:
        print("Please provide either p-value or q-value for peak calling")
        print("once pvalue is provided, qvalue will be ignored")
        sys.exit(1)
    print("after filtering has", len(tmp_peak_list), "peaks")

    return tmp_peak_list

def merge_peak(peak_list):  
    refined_peak_list=peak_list
    merged_peak_list = []
    current_peak = refined_peak_list[0]
    for kk in range(1,len(refined_peak_list)):
        if current_peak[1]>=refined_peak_list[kk][0]:
            current_peak[1] = refined_peak_list[kk][1]
        else:
            merged_peak_list.append(current_peak)
            current_peak = refined_peak_list[kk]
    merged_peak_list.append(current_peak)
    print("after merging has", len(merged_peak_list), "peaks")
    return merged_peak_list
def call_peak(input_bigwig, control_bigwig, output_dir, qvalue, pvalue, broad, broad_cutoff, window_size, thread):
    bw = pyBigWig.open(input_bigwig)
    control_bw = pyBigWig.open(control_bigwig)
    chroms = bw.chroms()
    output_bed = os.path.join(output_dir,"output_peaks.bed")
    with open(output_bed,'w') as wfile:
        wfile.write("#Peak calling from the bigWig file by hichelper.\n")
    check_window_size=[1000,5000,10000]
    max_window = max(check_window_size)
    for chrom in chroms:
        print("Processing chromosome: ", chrom)
        chrom_size = chroms[chrom]
        print("Chromosome size: ", chrom_size)
        signal = bw.values(chrom, 0, chrom_size, numpy=True)
        signal = np.nan_to_num(signal)
        control_signal = control_bw.values(chrom, 0, chrom_size, numpy=True)
        control_signal = np.nan_to_num(control_signal)
        """
        The lambda parameter is estimated from the control sample and 
        is deduced by taking the maximum value across various window sizes:
        λlocal = max(λBG, λ1k, λ5k, λ10k).
        """
        overall_lambda = np.mean(control_signal)    
        p =Pool(thread)
        res_list = []
        index_list = []
        for i in range(0, chrom_size-window_size,window_size):
            current_window = signal[i:i+window_size]
            if np.sum(current_window) <=1:
                continue
            index_list.append(i)
            cur_lambda_max = overall_lambda
            for tmp_window_size in check_window_size:
                start_index = max(0,i)
                end_index = min(i+tmp_window_size,chrom_size)
                tmp_mean = np.mean(control_signal[start_index:end_index])
                cur_lambda_max = max(cur_lambda_max,tmp_mean)
            res = p.apply_async(_call_peak, args=(current_window,cur_lambda_max,broad,i))
            res_list.append(res)
        p.close()
        p.join()
        tmp_peak_list=[]
        for i in range(len(res_list)):
            current_index = index_list[i]
            if not broad:
                pval=res_list[i].get()
                tmp_peak_list.append([current_index,pval])
            else:
                pval,broad_pval=res_list[i].get()
                tmp_peak_list.append([current_index,pval,broad_pval])
        
        tmp_peak_list = filter_peak(tmp_peak_list, pvalue, qvalue,broad,broad_cutoff)
        #merge nearby peaks and then calculate if it still pass the threshold
        refined_peak_list=[[x[0],x[0]+window_size] for x in tmp_peak_list]
        merged_peak_list = merge_peak(refined_peak_list)
        #then do the pvalue, qvalue calculation again on those merged peaks
        p =Pool(thread)
        res_list = []
        candidate_list=[]
        for kk in range(len(merged_peak_list)):
            start_index = merged_peak_list[kk][0]
            end_index = merged_peak_list[kk][1]
            current_window = signal[start_index:end_index]
            if np.sum(current_window) <=1:
                continue
            cur_lambda_max = overall_lambda
            mid_index = (start_index+end_index)//2
            for tmp_window_size in check_window_size:
            
                start_index2 = max(0,mid_index-tmp_window_size//2)
                end_index2 = min(mid_index+tmp_window_size//2,chrom_size)
                tmp_mean = np.mean(control_signal[start_index2:end_index2])
                cur_lambda_max = max(cur_lambda_max,tmp_mean)
            res = p.apply_async(_call_peak, args=(current_window,cur_lambda_max,broad,start_index))
            res_list.append(res)
            candidate_list.append([start_index,end_index])
        p.close()
        p.join()
        #get new peak list
        tmp_peak_list=[]
        end_index_map_dict={}
        for i in range(len(res_list)):
            current_index,cur_end_index = candidate_list[i]
            if not broad:
                pval=res_list[i].get()
                tmp_peak_list.append([current_index,cur_end_index,pval])
            else:
                pval,broad_pval=res_list[i].get()
                tmp_peak_list.append([current_index,cur_end_index,pval,broad_pval])
            end_index_map_dict[current_index]=cur_end_index   
        filter_peak_input=[[x[0],x[2],x[3]] for x in tmp_peak_list]
        #filter the peak list again
        tmp_peak_list = filter_peak(filter_peak_input, pvalue, qvalue,broad,broad_cutoff)
        with open(output_bed,'a') as peak_file:
            for i in range(len(tmp_peak_list)):
                if not broad:
                    start_index, pval = tmp_peak_list[i]
                else:
                    start_index, pval,broad_pval = tmp_peak_list[i]
                end_index = end_index_map_dict[start_index]
                
                log_pval = -np.log10(pval)
                fold_enrichment = np.mean(signal[start_index:end_index])/np.mean(control_signal[start_index:end_index])
                
                if broad:
                    log_broad_pval = -np.log10(broad_pval)
                    peak_file.write("%s\t%d\t%d\t%.4f\t%.4f\t%.4f\n"%(chrom,start_index,end_index,log_pval,log_broad_pval))
                else:
                    peak_file.write("%s\t%d\t%d\t%.4f\t%.4f\n"%(chrom,start_index,end_index,log_pval,fold_enrichment))
    bw.close()
    control_bw.close()
    return output_bed
def argparser():
    parser = argparse.ArgumentParser('Peak calling from the bigWig', add_help=False)
    parser.add_argument('-t', type=str, required=True, help='Path to the treatment bigWig file')
    parser.add_argument('-c', type=str, required=True, help='Path to the control bigWig file')
    parser.add_argument('-o', type=str, required=True, help='Output directory')
    parser.add_argument("-q", type=float,default=None, help="q-value cutoff (minimum FDR) for peak calling")
    parser.add_argument("-p",type=float, default=None, help="p-value cutoff for peak calling, if provided, q-value will be ignored")
    parser.add_argument("--broad",action="store_true", help="Use broad peak calling")
    #parser.add_argument("--broad-cutoff",type=float, default=0.1, help="Cutoff for broad peak calling")
    parser.add_argument("--min_length",type=int, default=300, help="Minimum length of the peak, can be set with fragment size")
    parser.add_argument("--thread",type=int, help="Number of threads to use",default=8)
    parser.add_argument("--broad_cutoff",type=float, default=0.1, help="Cutoff for broad peak calling (it will be q-value or qvalue based on either you choose -p or -q)\n Please note we still require that the the narrow summit should also pass its own q-value or p-value threshold")
    parser.add_argument('--help', action='help', help='Show this help message and exit')
    return parser
# macs3 callpeak --broad \
# -t ENCFF395DAJ.chr12.MAPQ30.blcklst.rh.sorted.bam \
# -c ENCFF956GLJ.chr12.MAPQ30.blcklst.rh.sorted.bam \
# -f BAMPE  -g 04.9e8 --broad-cutoff 0.1 -n neuroGM23338_macs3_rep1
"""
This script is to call peaks from the bigWig file.
The bigWig file is generated from the bam file by using the deepTools bamCoverage.
```
python3 peak_call_bigWig.py -t input.bw -c control.bw -o output_dir -q [qval] --min_length [min_length] --thread [thread] [--broad] [--broad_cutoff=[broad_cutoff]]
```
- t: path to the treatment bigWig file. <br>
- c: path to the control bigWig file. <br>
- o: output directory. <br>
- q: q-value cutoff (minimum FDR) for peak calling. <br>
- p: p-value cutoff for peak calling, if provided, q-value will be ignored. <br>
- min_length: minimum length of the peak, can be set with fragment size. <br>
- thread: number of threads to use. <br>
- broad: use broad peak calling. <br>
- broad_cutoff: cutoff for broad peak calling (it will be q-value or qvalue based on either you choose -p or -q). <br>
The output file is a bed file with the peak information. <br>
"""

if __name__ == '__main__':
    parser = argparser()
    args = parser.parse_args()
    input_bigwig = os.path.abspath(args.t)
    control_bigwig = os.path.abspath(args.c)
    output_dir = os.path.abspath(args.o)
    os.makedirs(output_dir, exist_ok=True)
    qvalue = float(args.q) if args.q else None
    pvalue = float(args.p) if args.p else None
    broad = args.broad
    broad_cutoff = args.broad_cutoff
    #https://github.com/crazyhottommy/ChIP-seq-analysis/blob/master/part1.3_MACS2_peak_calling_details.md explain broad param
    #broad_cutoff = args.broad_cutoff
    window_size = args.min_length
    thread = int(args.thread)
    output_path=call_peak(input_bigwig, control_bigwig, output_dir, qvalue, pvalue, broad, broad_cutoff, window_size, thread)
    print("Final peak file is saved at: ", output_path)
