import sys
import os 
import pickle
import numpy as np
import matplotlib.pyplot as plt
def load_pickle(pickle_path):
    """
    This script is used to load the pickle file.
    pickle_path: the path to the pickle file [String].
    """
    with open(pickle_path, 'rb') as f:
        data = pickle.load(f)
    return data

def plot_contact_stat(input_array_pickle, output_png, genomic_dist):
    """
    This script is used to plot contact frequency vs. genomic distance.
    input_array_pickle: the path to the pickle file containing the contact matrix [String].
    output_png: the name of the output png file [String].
    genomic_dist: the genomic distance for the plot [Integer].
    """
    data=load_pickle(input_array_pickle)
    all_column=[]
    all_row = []
    all_data = []
    for key in data.keys():
        matrix = data[key]
        matrix_row = matrix.row
        matrix_col = matrix.col
        matrix_data = matrix.data
        #we only care relative distance between row and column
        #so we can simply concat
        all_column.append(matrix_col)
        all_row.append(matrix_row)
        all_data.append(matrix_data)
    all_column = np.concatenate(all_column)
    all_row = np.concatenate(all_row)
    all_data = np.concatenate(all_data)
    all_relative_dist = np.abs(all_column-all_row)
    contact_freq_list=[]
    genome_dist_list=[]
    for k in range(1,genomic_dist):
        select_index = all_relative_dist==k
        contact_freq = np.sum(all_data[select_index])
        contact_freq_list.append(contact_freq)
        genome_dist_list.append(k)

    plt.plot(genome_dist_list,contact_freq_list,'r-o')
    plt.xlabel('Genomic distance')
    plt.ylabel('Contact frequency')
    plt.yscale('log')
    plt.title('Contact frequency vs. genomic distance')
    plt.savefig(output_png,dpi=600)
"""
This script is used to plot contact frequency vs. genomic distance.
```
python3 pkl_contact_stat.py [input.pkl] [output.png] [genomic_dist]
```
[input.pkl]: the path to the pickle file containing the contact matrix [String]. <br>
[output.png]: the name of the output png file [String]. <br>
[genomic_dist]: the genomic distance for the plot [Integer]. <br>

"""






if __name__ == '__main__':
    if len(sys.argv)!=4:
        print('Usage: python3 pkl_contact_stat.py [input.pkl] [output.png] [genomic_dist]')
        print("This script is used to plot contact frequency vs. genomic distance.")
        print("[input.pkl]: the path to the pickle file containing the contact matrix [String].")
        print("[output.png]: the name of the output png file [String].")
        print("[genomic_dist]: the genomic distance for the plot [Integer].")
        sys.exit(1)
    input_array_pickle = os.path.abspath(sys.argv[1])
    output_png = os.path.abspath(sys.argv[2])
    output_dir = os.path.dirname(output_png)
    os.makedirs(output_dir, exist_ok=True)
    genomic_dist = int(sys.argv[3])
    plot_contact_stat(input_array_pickle, output_png, genomic_dist)


