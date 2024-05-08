#!/bin/bash
shopt -s extglob
mapqfilter=0
juicer_tool_path=''
min_res=5000
maxmem='64g'
higlass=0  # if 1, higlass-compatible aggregation
custom_res=''  # custom resolutions, separated by commas
normonly=0  # if 1, normalization only
balance=1  # if 0, no normalization
number_cpu=8 # number of CPUs to use
output_prefix='out'


printHelpAndExit() {
    echo "Usage: ${0##*/} [-s juicer_tool_path] [-q mapqfilter] [-m maxmem] [-r min_res] [-g] [-o out_prefix] -i input_pairs -c chromsize_file"
    echo "-s juicer_tool_path : default /usr/local/bin/juicer_tools.jar"
    echo "-i input_pairs : input file in pairs.gz format"
    echo "-c chromsize_file : chromsizes file"
    echo "-o out_prefix : default out"
    echo "-q mapqfilter : default 0"
    echo "-m maxmem : default 64g"
    echo "-r min_res : default 5000"
    echo "-g : use HiGlass resolutions (default juicer resolutions)"
    echo "-u : custom resolutions (separated by comman)"
    echo "-n : normalization only"
    echo "-B : no balancing/normalization"
    echo "-t : number of CPUs to use"
    echo "-h : print this help message"
    exit "$1"
}

while getopts "i:s:c:q:r:m:go:nu:t:B" opt; do
    case $opt in
        i) input_pairs=$OPTARG;;
        s) juicer_tool_path=$OPTARG;;
        c) chromsizefile=$OPTARG;;
        t) number_cpu=$OPTARG;;
        q) mapqfilter=$OPTARG;;
        r) min_res=$OPTARG;;
        m) maxmem=$OPTARG;;
        g) higlass=1 ;;
        u) custom_res=$OPTARG;;
        n) normonly=1 ;;
        B) balance=0 ;;
        o) output_prefix=$OPTARG;;
        h) printHelpAndExit 0;;
        [?]) printHelpAndExit 1;;
        esac
done

echo "input_pairs: $input_pairs"
echo "juicer_tool_path: $juicer_tool_path"
echo "chromsizefile: $chromsizefile"
echo "mapqfilter: $mapqfilter"
echo "min_res: $min_res"
echo "maxmem: $maxmem"
echo "higlass: $higlass"
echo "custom_res: $custom_res"
echo "normonly: $normonly"
echo "balance: $balance"
echo "output_prefix: $output_prefix"
echo "number_cpu: $number_cpu"

# error when both higlass and custom_res are set
if [[ $higlass == '1' && ! -z $custom_res ]]
then
    echo "Do you want higlass resolution (-g) or custom-resolution (-u)? Make up your mind! :)"
    exit 1
fi


# error when custom resolution clashes with min_res
if [[ ! -z $custom_res ]]
then
    lowest_custom_res=${custom_res//,*/} 
    if [[ $lowest_custom_res != $min_res ]]
    then
        echo "Lowest custom res (-u) does not match min_res (-r)."
        exit 1
    fi
fi


# error when default juicer resolution clashes with min_res
if [[ $higlass == '0' && -z $custom_res ]]
then
    if [[ $min_res != '5000' ]]
    then
        echo "min_res (-r) should be 5000 for default juicer resolutions."
        exit 1
    fi
fi


# aggregation
if [[ $normonly == '0' ]]
then

    # creating a hic file
    if [[ $higlass == '1' ]]
    then
        reslist=$(python3 -c "from cooler.contrib import higlass; higlass.print_zoom_resolutions('$chromsizefile', $min_res)")
	java -Xmx$maxmem -Xms$maxmem -jar $juicer_tool_path pre -n -j $number_cpu -r $reslist -q $mapqfilter $input_pairs $output_prefix.hic $chromsizefile 
    elif [[ ! -z $custom_res ]]
    then
        java -Xmx$maxmem -Xms$maxmem -jar $juicer_tool_path pre -n -j $number_cpu -r $custom_res -q $mapqfilter $input_pairs $output_prefix.hic $chromsizefile 
    else
        java -Xmx$maxmem -Xms$maxmem -jar $juicer_tool_path pre -n -j $number_cpu -q $mapqfilter $input_pairs $output_prefix.hic $chromsizefile 
    fi
fi


# normalization
if [[ $balance == '1' ]]
then
  java -Xmx$maxmem -Xms$maxmem -jar $juicer_tool_path addNorm -j $number_cpu -w $min_res -d -F $output_prefix.hic
fi
