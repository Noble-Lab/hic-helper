
import sys
import os

def addnorm2hic(juicer_tools,input_hic,resolution,number_cpu,max_memory):
    command_line=f"java -Xmx{max_memory}g -Xms{max_memory}g -jar {juicer_tools} addNorm -j {number_cpu} -w {resolution} -d -F {input_hic}"
    print("exec command:")
    print(command_line)
    os.system(command_line)
"""
This script is used to add norm to the hic files.
```
 python3 addnorm2hic.py [input.hic] [resolution] [num_cpu] [memory]
```
[input.hic]: input hic path. <br>
[resolution]: minimum resolution that normalization works. <br>
[num_cpu]: number of cpus used to normalize. <br>
[memory]: maximum memory (GB) that this script is allowed to use. <br>
The calculated norm vectors will be automatically saved in the input.hic file. <br>

"""
if __name__ == '__main__':
    if len(sys.argv)!=5:
        print('Usage: python3 addnorm2hic.py [input.hic] [resolution] [num_cpu] [memory]')
        print("[input.hic]: input hic path")
        print("[resolution]: minimum resolution that normalization works.")
        print("[num_cpu]: number of cpus used to normalize.")
        print("[memory]: maximum memory (GB) that this script is allowed to use.")
        print("The calculated norm vectors will be automatically saved in the input.hic file.")
        sys.exit(1)
    
    script_dir = os.path.dirname(os.path.realpath(__file__))
    juicer_tools = os.path.join(script_dir, 'juicer_tools.jar')
    input_hic = os.path.abspath(sys.argv[1])
    resolution = int(sys.argv[2])
    number_cpu = int(sys.argv[3])
    max_memory = int(sys.argv[4])

    addnorm2hic(juicer_tools,input_hic,resolution,number_cpu,max_memory)




    

