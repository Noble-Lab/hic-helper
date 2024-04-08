import sys
import os


# java -jar juicer_tools.jar hiccups
# Usage:   juicer_tools hiccups [-m matrixSize] [-k normalization (NONE/VC/VC_SQRT/KR)] [-c chromosome(s)] [-r resolution(s)] [-f fdr] [-p peak width] [-i window] [-t thresholds] [-d centroid distances] [--ignore-sparsity]<hicFile> <outputDirectory> [specified_loop_list]

"""

 * -r resolution
    * Resolution of the Hi-C data. Multiple resolutions can be used by separating them with commas (e.g. "-r 5000,10000").
 
 * -m <int> Maximum size of the submatrix within the chromosome passed on to GPU (Must be an even number greater than 40
 * to prevent issues from running the CUDA kernel). The upper limit will depend on your GPU. Dedicated GPUs
 * should be able to use values such as 500, 1000, or 2048 without trouble. Integrated GPUs are unlikely to run
 * sizes larger than 90 or 100. Matrix size will not effect the result, merely the time it takes for hiccups.
 * Larger values (with a dedicated GPU) will run fastest.

 * -f <int(s)> FDR values actually corresponding to max_q_val (i.e. for 1% FDR use 0.01, for 10%FDR use 0.1). Different
 * FDR values can be used for each resolution using commas. (e.g "-r 5000,10000 -f 0.1,0.15" would run HiCCUPS at
 * 10% FDR for resolution 5000 and 15% FDR for resolution 10000)

 * -p <int(s)> Peak width used for finding enriched pixels in HiCCUPS. Different peak widths can be used for each
 * resolution using commas. (e.g "-r 5000,10000 -p 4,2" would run at peak width 4 for resolution 5000 and
 * peak width 2 for resolution 10000)

  * -w <int(s)> Window width used for finding enriched pixels in HiCCUPS. Different window widths can be used for each
 * resolution using commas. (e.g "-r 5000,10000 -p 10,6" would run at window width 10 for resolution 5000 and
 * window width 6 for resolution 10000)

* -t <floats> Thresholds for merging loop lists of different resolutions. Four values must be given, separated by
 * commas (e.g. 0.02,1.5,1.75,2). These thresholds (in order) represent:
 * > threshold allowed for sum of FDR values of horizontal, vertical donut mask, and bottom left regions
 * (an accepted loop must stay below this threshold)
 * > threshold ratio of observed value to expected horizontal/vertical value
 * (an accepted loop must exceed this threshold)
 * > threshold ratio of observed value to expected donut mask value
 * (an accepted loop must exceed this threshold)
 * > threshold ratio of observed value to expected bottom left value
 * (an accepted loop must exceed this threshold)

 
 * -d <ints> Distances used for merging centroids across different resolutions. Three values must be given, separated by
 * commas (e.g. 20000,20000,50000). These thresholds (in order) represent:
 * > distance (radius) around centroid used for merging at 5kB resolution (if present)
 * > distance (radius) around centroid used for merging at 10kB resolution (if present)
 * > distance (radius) around centroid used for merging at 25kB resolution (if present)

     * Reasonable Commands
     *
     * -f fdr = 0.10 for all resolutions
     * -p peak width = 1 for 25kb, 2 for 10kb, 4 for 5kb
     * -w window = 3 for 25kb, 5 for 10kb, 7 for 5kb
     *
     * -d cluster radius is 20kb for 5kb and 10kb res and 50kb for 25kb res
     * fdrsumthreshold is 0.02 for all resolutions
     * -t
     * oeThreshold1 = 1.5 for all res
     * oeThreshold2 = 1.75 for all res
     * oeThreshold3 = 2 for all res
     * -k KR normalization for all resolutions, VC norm can be also an option if KR failed (suggested by mail communication of Erez's group)
     * published GM12878 looplist was only generated with 5kb and 10kb resolutions
     * same with published IMR90 looplist
     * published CH12 looplist only generated with 10kb
     */
"""
def hiccups_loop(juicer_tools,hicFile,output_loop,resolution):
    output_dir = os.path.dirname(output_loop)
    os.makedirs(output_dir,exist_ok=True)
    if resolution not in [5000,10000,25000]:
        raise ValueError('The resolution is not supported. The supported resolutions are 5000,10000,25000')
    if resolution == 5000:
        specific_param = "-p 4 -w 7 -d 20000"
    elif resolution == 10000:
        specific_param = "-p 2 -w 5 -d 20000"
    elif resolution == 25000:
        specific_param = "-p 1 -w 3 -d 50000"
    os.system(f'java -jar {juicer_tools} hiccups -r {resolution} -f 0.1 {specific_param} -k KR {hicFile} {output_loop}')
    return
# Path: post_processing/hiccups_loop.py

if __name__ == '__main__':
    if len(sys.argv) !=4:
        print('Usage: python3 hiccups_loop.py [hicFile] [outputloop_file] [resolution]')
        print("hicFile: the path to the input hic file [String].")
        print("outputloop_file: the path to the output loops [String].")
        print("resolution: the resolution of the input hic file [Integer].")
        sys.exit(1)
    script_dir = os.path.dirname(os.path.realpath(__file__))
    script_dir = script_dir.replace('post_processing','pre_processing')
    juicer_tools = os.path.join(script_dir, 'juicer_tools.jar')
    hicFile = os.path.abspath(sys.argv[1])
    output_loop = os.path.abspath(sys.argv[2])
    resolution = int(sys.argv[3])
    
    hiccups_loop(juicer_tools,hicFile,output_loop,resolution)





