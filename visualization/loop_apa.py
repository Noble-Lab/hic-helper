
import sys 
import os
import pickle
import numpy as np
def extract_loc(pred_detect_path):
    record_list=[]
    with open(pred_detect_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip().split()
            chrom1 = line[0]
            try:
                start1 = int(float(line[1]))
                end1 = int(float(line[2]))
            except:
                continue
            chrom2 = line[3]
            try:
                start2 = int(float(line[4]))
                end2 = int(float(line[5]))
            except:
                continue
            if "chr" not in chrom1:
                chrom1 = "chr"+chrom1
            
            if "chr" not in chrom2: 
                chrom2 = "chr"+chrom2

            record_list.append([chrom1, start1, end1, chrom2, start2, end2])
    print("Extracted",len(record_list)," loop detection records")
    return record_list
def find_chrom_key(hic_data,chrom1):
    if chrom1 in hic_data:
        return chrom1
    chrom1_nochr = chrom1.replace("chr","")
    if chrom1_nochr in hic_data:
        return chrom1_nochr
    chrom_merge=chrom1+"_"+chrom1 
    if chrom_merge in hic_data:
        return chrom_merge
    chrom_merge=chrom1_nochr+"_"+chrom1_nochr
    if chrom_merge in hic_data:
        return chrom_merge
def reorganize_loop(loop_list):
    loop_dict = {}
    for loop in loop_list:
        chrom1, start1, end1, chrom2, start2, end2 = loop
        if chrom1!=chrom2:
            print("loop chromsome not match")
            sys.exit(1)
        if chrom1 not in loop_dict:
            loop_dict[chrom1] = []
        loop_dict[chrom1].append(loop)
    return loop_dict
def loop_apa(hic_file, input_bed, output_png, resolution,region_size=21):
    center_size=5
    #read the hic matrix
    hic_data=pickle.load(open(hic_file,'rb'))
    #read the input bed file get all records
    loop_list = extract_loc(input_bed)
    output_array = np.zeros([region_size,region_size])
    #reorganize the loop list to a chrom key based dict
    loop_dict = reorganize_loop(loop_list)
    for chrom in loop_dict:
        loop_list = loop_dict[chrom]
        chrom = find_chrom_key(hic_data,chrom)
        #judge if it is coo matrix
        if hasattr(hic_data[chrom],'toarray'):
            chrom_hic = hic_data[chrom].toarray()
        chrom_hic = chrom_hic + np.triu(chrom_hic,1).T 
        #symmetrize the matrix
        for kk,loop in enumerate(loop_list):
            chrom1, start1, end1, chrom2, start2, end2 = loop
            #chrom2 = find_chrom_key(hic_data,chrom2)
            half_size = region_size//2
            start1 = start1//resolution
            end1 = end1//resolution
            start2 = start2//resolution
            end2 = end2//resolution
            mid1 = (start1+end1)//2
            mid2 = (start2+end2)//2
            start1 = max(0,mid1-half_size)
            end1 = min(chrom_hic.shape[0],mid1+half_size)
            start2 = max(0,mid2-half_size)
            end2 = min(chrom_hic.shape[0],mid2+half_size)
            output_start1 = half_size - (mid1-start1)
            output_end1 = output_start1 + end1-start1
            output_start2 = half_size - (mid2-start2)
            output_end2 = output_start2 + end2-start2
            output_array[output_start1:output_end1,output_start2:output_end2] += chrom_hic[start1:end1,start2:end2]
        print("Processed",chrom)
    output_array = output_array/len(loop_list)
    #calculate the average peak value
    left_bottom_region = output_array[-center_size:,:center_size]
    center_start = (region_size-center_size)//2
    center_region = output_array[center_start:center_start+center_size,center_start:center_start+center_size]
    peak_strength = np.mean(center_region)
    p2ll = np.sum(center_region)/np.sum(left_bottom_region)
    print("Peak strength:",peak_strength)
    print("Peak to lower-left ratios:",p2ll)
    import matplotlib.pyplot as plt
    plt.figure(figsize=(5,5))
    plt.imshow(output_array,cmap='hot',interpolation='nearest')
    plt.colorbar()
    plt.title("Loop APA(Strength:%.2f, P2LL:%.2f)"%(peak_strength,p2ll))
    plt.savefig(output_png,dpi=600)

"""
This script is for plot the loop average peak analysis (APA) on the hic matrix.
```
python3 loop_apa.py [hic.pkl] [input.bed] [output.png] [resolution]
```
- hic.pkl: the hic matrix file <br>
- input.bed: the input bed file including the loop regions <br>
- output.png: the output loop APA png file <br>
- resolution: the resolution of the hic matrix <br>
"""
    




if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python3 loop_apa.py [hic.pkl] [input.bed] [output.png] [resolution]")
        print("This script is for plot the loop average peak analysis (APA) on the hic matrix")
        print("hic.pkl: the hic matrix file")  
        print("input.bed: the input bed file including the loop regions")
        print("output.png: the output loop APA png file")
        print("resolution: the resolution of the hic matrix")
        sys.exit(1)
    hic_file = os.path.abspath(sys.argv[1])
    input_bed = os.path.abspath(sys.argv[2])
    output_png = os.path.abspath(sys.argv[3])
    resolution = int(sys.argv[4])
    output_dir = os.path.dirname(output_png)
    os.makedirs(output_dir,exist_ok=True)
    loop_apa(hic_file, input_bed, output_png, resolution) 