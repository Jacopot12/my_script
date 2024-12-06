#!/bin/bash

#### PEAKS SELECTION SCRIPT ####

# This script allows the selection of peaks from .narrowPeak file (or .bed file) on the basis of peak length
# Procedure: 1) calculation of the upper and lower outlier
#            2) outlier removal
#            3) use bedtools intersect and merge to make a common peakset
#            4) remove peaks shorter than average fragment length

# Method for calculating outlier: Boxplot method based on the Interquartile Range (IQR)
# The script accepts at least 2 .narrowPeak files (or .bed files) or at most 3 .narrowPeak files (or .bed files)
# The script accepts a value (integer number) of average fragment length (AFL) or the AFL for each input peaks files


### FUNCTIONS SECTION ###

# Function to show how to use the script
usage() {
  echo "Usage: ./peaksSelection.sh -p <file1.narrowPeak> <file2.narroPeak> [file3.narroPeak] -f <number1> [number2] [number3]"
  exit 1
}

# Function to check file existence 
file_existence_check(){
  local file="$1"

  if [ -n "$file" ]; then
    if [ ! -f "$file" ]; then
      echo "Error: The file path '$file' does not exist or is not a file."
      usage
      exit 1
    fi
  fi
}

# Initialize arrays
peak_files=()
fragments=()
output_names=()


while getopts "p:f:" opt; do
  case $opt in

    p)
      peak_files+=("$OPTARG")
      while [[ "$OPTIND" -le "$#" && "${!OPTIND:0:1}" != "-" ]]; do
        peak_files+=("${!OPTIND}")
        OPTIND=$((OPTIND + 1))
      done
      ;;
    
    f)
      if ! [[ "$OPTARG" =~ ^[0-9]+$ ]]; then
        echo "Error: -f argument must be a positive integer."
        usage
      fi
      fragments+=("$OPTARG")
      while [[ "$OPTIND" -le "$#" && "${!OPTIND:0:1}" != "-" ]]; do
        if ! [[ "${!OPTIND}" =~ ^[0-9]+$ ]]; then
          echo "Error: -f argument must be a positive integer."
          usage
        fi
        fragments+=("${!OPTIND}")
        OPTIND=$((OPTIND + 1))
      done
      ;;
    
    *) usage ;;
  esac
done



### CONTROL SECTION ###

# Checking that the required flags have been passed
if [[ ${#peak_files[@]} -eq 0 ]]; then
  echo "Error: -p argument is required."
  usage
fi

if [[ ${#fragments[@]} -eq 0 ]]; then
  echo "Error: -f argument is required."
  usage
fi

# Checking the number of peaks files
if [[ ${#peak_files[@]} -gt 3 || ${#peak_files[@]} -lt 2 ]]; then
  echo "Error: -f argument accepts at least 2 files or at most 3 files."
  usage
fi

# Check the existence of the file
numPeakFile=${#peak_files[@]}
for (( i=0; i<$numPeakFile; i++ )); do
  file_existence_check ${peak_files[i]}
done

# I check the file extension of the peaks. It can be .narrowPeak or .bed
numPeakFile=${#peak_files[@]}
for (( i=0; i<$numPeakFile; i++ )); do
  extension=$(echo ${peak_files[i]} | cut -d . -f 2)
  if [[ ${extension} != "narrowPeak" && ${extension} != "bed" ]]; then
    echo "Error: Invalid extension! Peak files must be in .narrowPeak or .bed extension"
    exit 1
  fi
done

if [[ ${#fragments[@]} != 1 && ${#fragments[@]} != ${#peak_files[@]} ]]; then
  echo "Error: Invalid average fragment length numbers"
  usage
  exit 1

fi

### MAIN SECTION ###

# calculation of average fragment length AFL
lengthF=${#fragments[@]}
sumF=0
for (( i=0; i<$lengthF; i++ )); do
  sumF=$((sumF + ${fragments[i]}))
done

AFL=$(( $sumF / $lengthF))



# calculating and removing outliers

for (( i=0; i<$numPeakFile; i++ )); do

  risultati=()
  while IFS= read -r risultato; do
      risultati+=("$risultato")
  done < <(awk '{print ($3 - $2)}' "${peak_files[i]}")

  # Convert array of results to a comma-separated string
  risultati_stringa=$(IFS=,; echo "${risultati[*]}")


  output=$(R --vanilla --slave << EOF
  
  # Transforming the Bash variable into an R vector
  peaks = c($risultati_stringa)  # use string with commas
  
  Q1 = quantile(peaks, 0.25)
  Q3 = quantile(peaks, 0.75)
  IQR_value = IQR(peaks)
  lower_limit = round(Q1 - 1.5 * IQR_value)
  upper_limit = round(Q3 + 1.5 * IQR_value)
  
  lower_limit = if (lower_limit < 0) 0 else lower_limit
  upper_limit = if (upper_limit < 0) 0 else upper_limit

  # interquartile range limits to pass in bash
  cat(lower_limit, upper_limit, sep=",")

  
EOF
)

lower_limit=$(echo "$output" | cut -d',' -f1)
upper_limit=$(echo "$output" | cut -d',' -f2)

# Visualizzare i risultati
#echo "Il valore di lower limit è: $lower_limit"
#echo "Il valore di upper limit è: $upper_limit"

name=$(echo ${peak_files[i]} | cut -d "." -f 1)

awk -v FS="\t" -v OFS="\t" -v lower_limit=${lower_limit} -v upper_limit=${upper_limit} '($3 - $2 > lower_limit && $3 - $2 < upper_limit) {print $1, $2, $3}' "${peak_files[i]}" | sort -k1,1V -k2,2n > "${name}"_noOutlier.bed

output_names+=("${name}"_noOutlier.bed)

done


if [[ ${#peak_files[@]} -eq 2 ]]; then

  echo "Sono 2 file"

  #intersect the .bed file together to keep consensus peaks
  bedtools intersect -a ${output_names[0]} -b ${output_names[1]} > consensusPeaksA.bed

  #sort the consensus peaks file 
  sort -k1,1V -k2,2n consensusPeaksA.bed > consensusPeaksA_sorted.bed

  #merge peaks that are overlapped 
  bedtools merge -i consensusPeaksA_sorted.bed > consensusPeaksA_sorted_merged.bed

  #excluding peaks shorter than AFL
  awk -v AFL=${AFL} '(($3-$2) >= AFL) {print $0} ' consensusPeaksA_sorted_merged.bed > peaksSet.bed

  #peaks set with length
  awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,($3-$2)}' peaksSet.bed > peaksSet_length.bed

  rm consensusPeaksA.bed consensusPeaksA_sorted.bed consensusPeaksA_sorted_merged.bed 


elif [[ ${#peak_files[@]} -eq 3 ]]; then

  echo "Sono 3 file"

  #intersect the .bed file together to keep consensus peaks
  bedtools intersect -a ${output_names[0]} -b ${output_names[1]} ${output_names[2]} > consensusPeaksA.bed

  #sort the consensus peaks file 
  sort -k1,1V -k2,2n consensusPeaksA.bed > consensusPeaksA_sorted.bed

  #merge peaks that are overlapped 
  bedtools merge -i consensusPeaksA_sorted.bed > consensusPeaksA_sorted_merged.bed


  #find peaks between file2 and file3
  #some are already present in consensusPeaksA_sorted_merged.bed
  bedtools intersect -a ${output_names[1]} -b ${output_names[2]} > consensusPeaksB.bed 

  #sort peaks 
  sort -k1,1V -k2,2n consensusPeaksB.bed > consensusPeaksB_sorted.bed

  #find peaks that are only in file2 and file3
  bedtools intersect -v -a consensusPeaksB_sorted.bed -b consensusPeaksA_sorted_merged.bed > consensusPeaks_onlyB.bed 

  cat consensusPeaksA_sorted_merged.bed consensusPeaks_onlyB.bed > consensusPeaksAB.bed

  sort -k1,1V -k2,2n consensusPeaksAB.bed > consensusPeaksAB_sorted.bed


  #excluding peaks shorter than AFL
  awk -v AFL=${AFL} '(($3-$2) >= AFL) {print $0} ' consensusPeaksAB_sorted.bed > peaksSet.bed

  #peaks set with length
  awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,($3-$2)}' peaksSet.bed > peaksSet_length.bed

  rm consensusPeaksA.bed consensusPeaksA_sorted.bed consensusPeaksA_sorted_merged.bed 
  rm consensusPeaksB.bed consensusPeaksB_sorted.bed consensusPeaks_onlyB.bed consensusPeaksAB.bed consensusPeaksAB_sorted.bed

fi
