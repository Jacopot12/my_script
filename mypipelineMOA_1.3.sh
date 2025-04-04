#!/bin/bash

#Author: Jacopo Tartaglia
#Verson: 1.3
#What the script does:
#In the first part, the script takes the paired-ends reads derived from MOA-seq (after the QC and trimming) and merge them to single-ends reads (using NGmerge)
#In the second part, the script maps the single-ends reads on the genome of interest. It should be specified whether the genome is to be indexed de novo or whether the folder with the indexed genome is already present (and specified)
#In the third part, the script filterd the .bam output in order to keep the sequences <80 bp and uniqued mapped on the genome (MAPQ 255). After that it calculate the average fragment length and the effective genome size
#If you have more than one reads file this script merged the .bam file togather after the mapping step. You can choose the name of the merged file by the flag -n
#If you have only one sample you can put its name in the flag -n
#In the fourth part, the script uses these information to perform the peakcalling using MACS3


#Usage:
# $ bash ./mypipleineMOA.sh -f read_files1_1.fastq.gz read_files2_1.fastq.gz ... -r read_files1_2.fastq.gz read_files2_2.fastq.gz ... -g genome.fasta [-G /index_genome_dir] [-a annot.gtf] [-n sample1_control]
# The only required files are the paired-ends reads files, the genome file.
# The files must be in the same directory of the script (no simblic link).
# Since the operation messages of the programs are shown change the parameters according to what is printed 

#Programs required (docker container):
# - NGmerge: quay.io/biocontainers/ngmerge:0.3--0
# - STAR: quay.io/biocontainers/star:2.7.9a--h9ee0642_0
# - samtools: staphb/samtools:latest / docker.io/staphb/samtools:1.21
# - unique-kmers.py: quay.io/biocontainers/khmer:3.0.0a3--py311heabec7a_7
# - MACS3: quay.io/biocontainers/macs3:3.0.2--py310h397c9d8_2




### FUNCTIONS ###

# Usage function to show how to use the script
usage() {
  echo "Usage: $0 -f read_files1_1.fastq.gz read_files2_1.fastq.gz ... -r read_files1_2.fastq.gz read_files2_2.fastq.gz ... -g /absolute/path/to/genome.fasta [-G /absolute/path/to/genome_dir] [-a /absolute/path/to/annot.gtf] [-n sample1_control]"
  echo "Options:"
  echo "  -f    Specify one or more input file for forward reads."
  echo "  -r    Specify one or more input file for reverse reads."
  echo "  -g    Specify the genome file in fasta format (one file, no compression)."
  echo "  -G    Specify the genome index directory, if any (optional)."
  echo "  -a    Specify the GTF file for STAR, if any (optional)."
  echo "  -n    Specify a name (no file extension) for merged biological replicate of the sample (optional, default: all_sample_merged)."
  echo "  -h    Display this help message."
  
  exit 1
}

# FUNCTIONS

file_existence_check(){
  local file="$1"

  if [ -n "$file" ]; then
    if [ ! -f "$file" ]; then
      echo "Error: The file path '$file' does not exist or is not a file."
      usage
      exit 1
    fi
      #echo "File path: $file"
  fi
}

dir_existence_check(){
  local dir=$1

  if [ -n "$dir" ]; then
    if [ ! -d "$dir" ]; then
      echo "Error: The directory path '$dir' does not exist or is not a directory."
      usage
      exit 1
    fi
      #echo "Directory path: $dir"
  fi
}


check_reads_file_extension() {
  local file="$1"
  if [[ "$file" != *.fastq.gz && "$file" != *.fq.gz ]]; then
    echo "Error: The file '$file' is not in the correct format (.fastq.gz or .fq.gz)."
    usage
    exit 1
  fi
}

check_genome_file_extension() {
  local file="$1"

  if [[ -n "$file" ]]; then
    if [[ "$file" != *.fasta && "$file" != *.fa  && "$file" != *.fna ]]; then
      echo "Error: The file '$file' is not in the correct format (.fasta, .fa or .fna)."
      exit 1
    fi
  fi
}

check_annotation_file_extension() {
  local file="$1"

  if [[ -n "$file" ]]; then
    if [[ "$file" != *.gtf && "$file" != *.gff ]]; then
      echo "Error: The file '$file' is not in the correct format (.gtf or .gff)."
      exit 1
    fi
  fi
}




# initialization of variables for input files/directories
readsF=()
readsR=()
file_genome_path=""
dir_genome_path=""
file_annotation=""
sample_name="all_sample_merged"
main_folder=$(pwd)

echo ". 
.      __  _                           
 /|/| /  )/_| __   _ _ _    '  _ /'  _ 
/   |(__/(  |    _) (-(/ /)//)(-(//)(- 
                      / /  /           "


echo " -----------------------------"
echo " ------STARTING PIPELINE------"

# checks to verify that all flags have been provided
flag_f_provided=false
flag_r_provided=false
flag_g_provided=false



# input file management
while getopts "f:r:g:G:a:n:h" opt; do
  case $opt in
    
    f) readsF+=("$OPTARG")
       flag_f_provided=true
       while [[ "$OPTIND" -le "$#" && "${!OPTIND:0:1}" != "-" ]]; do
         readsF+=("${!OPTIND}")
         OPTIND=$((OPTIND + 1))
       done
       ;;
    
    r) readsR+=("$OPTARG")
       flag_r_provided=true
       while [[ "$OPTIND" -le "$#" && "${!OPTIND:0:1}" != "-" ]]; do
         readsR+=("${!OPTIND}")
         OPTIND=$((OPTIND + 1))
       done
       ;;
    
    g) file_genome_path="$OPTARG"
       flag_g_provided=true
       ;;
    
    G) dir_genome_path="$OPTARG"
       ;;  
    
    a) file_annotation="$OPTARG"
       ;;
       
    n) sample_name="$OPTARG"
       ;;

    h) usage
       ;;

    *) usage
       ;;
  esac
done

# verifies that all required flags have been provided
if ! $flag_f_provided; then
  echo "Error: -f argument is required."
  usage
  exit 1
fi

if ! $flag_r_provided; then
  echo "Error: -r argument is required."
  usage
  exit 1
fi

if ! $flag_g_provided; then
  echo "Error: Either -g argument is required."
  usage
  exit 1
fi

#if [ -n "$file_genome_path" ] && [ -n "$dir_genome_path" ]; then
 # echo "Error: You can only specify one of -g or -G."
 # usage
 # exit 1
#fi

if [[ ${#readsF[@]} -ne ${#readsR[@]} ]]; then
  echo "Error: -f and -r must have the same number of elements."
  usage
  exit 1
fi

# check existence of files or dir and their extensions
file_existence_check "$file_genome_path"
check_genome_file_extension "$file_genome_path"
check_annotation_file_extension "$file_annotation"
dir_existence_check "$dir_genome_path"

lengthF=${#readsF[@]}
for (( i=0; i<$lengthF; i++ )); do
  file_existence_check "${readsF[$i]}"
  check_reads_file_extension "${readsF[$i]}"
done

lengthR=${#readsR[@]}
for (( i=0; i<$lengthR; i++ )); do
  file_existence_check "${readsR[$i]}"
  check_reads_file_extension "${readsR[$i]}"
done



#### FIRST PART ####
# from paired-ends reads to single-ends reads with NGmerge


if [ ! -d "./NGmerge_dir" ]; then
  mkdir NGmerge_dir
fi

cd NGmerge_dir

for (( i=0; i<$lengthF; i++ )); do
  echo "Start running NGmerge, sample $(( $i + 1 )), merging ${readsF[$i]} and ${readsR[$i]}"

  sample=$(basename "${readsF[$i]}" | cut -d '_' -f 1)

  docker run --rm \
    -v "$(pwd):/data" \
    -v "$(dirname "$(realpath ../"${readsF[$i]}")"):/input" \
    -u 1000:1000 \
    quay.io/biocontainers/ngmerge:0.3--0 \
    NGmerge \
      -1 /input/"${readsF[$i]}" \
      -2 /input/"${readsR[$i]}" \
      -o /data/${sample}_merged.fq.gz \
      -p 0.2 \
      -m 15 \
      -d \
      -e 30 \
      -z \
      -n 48 \
      -v

  echo " "
  echo " "
  echo "Termination of NGmerge run for sample $(( $i + 1 ))"
  echo " "
  echo "############################"
  echo " "
done

cd ../




#### SECOND PART ####
# reads mapping onto a reference genome

if [ ! -d "STAR" ]; then
  mkdir STAR
fi

# Controlla se la cartella contenente il genoma indicizzato esiste
if [[ ! -d "$dir_genome_path" ]]; then

  mkdir genomeDir/
  chmod 755 genomeDir/
  dir_genome_path="genomeDir/"

  if [[ -e "$file_annotation" ]]; then
    docker run --rm -u 1000:1000 -v "$main_folder:/data" quay.io/biocontainers/star:2.7.9a--h9ee0642_0 \
      STAR --runThreadN 2 \
           --runMode genomeGenerate \
           --genomeFastaFiles "/data/${file_genome_path}" \
           --genomeDir "/data/${dir_genome_path}" \
           --genomeSAindexNbases 13 \
           --limitGenomeGenerateRAM 6442450944 \
           --sjdbGTFfile "/data/${file_annotation}" \
           --outFileNamePrefix "/data/${file_genome_path}"

  else
    docker run --rm -u 1000:1000 -v "$main_folder:/data" quay.io/biocontainers/star:2.7.9a--h9ee0642_0 \
      STAR --runThreadN 2 \
           --runMode genomeGenerate \
           --genomeFastaFiles "/data/${file_genome_path}" \
           --genomeDir "/data/${dir_genome_path}" \
           --genomeSAindexNbases 13 \
           --limitGenomeGenerateRAM 6442450944 \
           --outFileNamePrefix "/data/${file_genome_path}"
  fi
fi

echo " "
echo " "
echo "Termination of STAR genome generate step"
echo "############################"
echo " "
echo " "

cd STAR

for (( i=0; i<$lengthF; i++ )); do
  sample=$(basename "${readsF[$i]}" | cut -d '_' -f 1)
  mkdir "$sample" 
  cd "$sample"

  echo "Start running STAR for sample ${sample}"

  # IMPORTANT! use "zcat" on Linux or "gzcat" on macOS (zsh)
  docker run --rm -u 1000:1000 -v "$(pwd):/output" -v "$main_folder:/data"  \
    quay.io/biocontainers/star:2.7.9a--h9ee0642_0 \
    STAR --readFilesCommand zcat \
         --genomeDir "/data/${dir_genome_path}" \
         --runThreadN 8 \
         --readFilesIn "/data/NGmerge_dir/${sample}_merged.fq.gz" \
         --outSAMmultNmax 1 \
         --outFilterMultimapNmax 2 \
         --winAnchorMultimapNmax 100 \
         --outMultimapperOrder Random \
         --runRNGseed 100 \
         --outFileNamePrefix "/output/out${sample}_merged" \
         --outBAMsortingBinsN 15 \
         --alignIntronMax 1 \
         --outSAMtype BAM SortedByCoordinate

  echo " "
  echo "Termination of STAR run for sample ${sample}"
  echo "############################"
  echo " "
  echo " "
  
  cd ..

done

cd ..



#### THIRD PART ####
# parsing of the STAR output

echo "Starting analysis of alignment files"
echo " "
mkdir -p "parsing_bam_files"
cd parsing_bam_files

command="docker run --rm -u 1000:1000 -v $(pwd):/data docker.io/staphb/samtools:1.21 samtools merge -f -@ 6 -o /data/${sample_name}.bam "

for (( i=0; i<$lengthF; i++ )); do
  sample=$(basename "${readsF[$i]}" | cut -d '_' -f 1)
  
  cp ../STAR/"${sample}"/out"${sample}_mergedAligned.sortedByCoord.out.bam" ./
  command+="/data/out${sample}_mergedAligned.sortedByCoord.out.bam "
done

eval "$command"
rm *.out.bam

# Keep only mapped reads max 80 bp long
docker run --rm -u 1000:1000 -v $(pwd):/data docker.io/staphb/samtools:1.21 \
  samtools view -@ 6 -h /data/${sample_name}.bam | \
  awk 'length($10) < 81 || $1 ~ /^@/' | \
  docker run -u 1000:1000 -v $(pwd):/data --rm -i docker.io/staphb/samtools:1.21 \
  samtools sort -@ 6 - | \
  docker run --rm -u 1000:1000 -v $(pwd):/data -i docker.io/staphb/samtools:1.21 \
  samtools view -@ 6 -bS - -o /data/${sample_name}_max80.bam

docker run --rm -u 1000:1000 -v $(pwd):/data docker.io/staphb/samtools:1.21 \
  samtools index -@ 6 /data/${sample_name}_max80.bam

# Keep only reads with MAPQ equal to 255 (uniq mapping reads)
docker run --rm -u 1000:1000 -v $(pwd):/data docker.io/staphb/samtools:1.21 \
  samtools view -@ 6 -q 255 /data/${sample_name}_max80.bam -o /data/${sample_name}_max80_255.bam

docker run --rm -u 1000:1000 -v $(pwd):/data docker.io/staphb/samtools:1.21 \
  samtools index -@ 6 /data/${sample_name}_max80_255.bam

rm -f *_max80.bam
rm -f *_max80.bam.bai

docker run --rm -u 1000:1000 -v $(pwd):/data docker.io/staphb/samtools:1.21 \
  samtools stats /data/${sample_name}_max80_255.bam > stats_${sample_name}_max80_255.txt

AFL=$(gawk -v OFS='\t' '{if($2=="average" && $3=="first" && $4=="fragment"){print $6}}' stats_${sample_name}_max80_255.txt)
echo " "
echo "The average fragment length: "${AFL}
echo " "

echo "Termination of alignment files analysis"
echo " "
echo "############################"
echo " "



# calculation of the effective genome size


echo " "
echo "Starting calculation of effective genome size"
echo " "

docker run --rm -u 1000:1000 -v "$main_folder:/data" quay.io/biocontainers/khmer:3.0.0a3--py311heabec7a_7 \
   unique-kmers.py -q -k ${AFL} -R /data/parsing_bam_files/${sample_name}_max80_255_frg${AFL}_egSize.txt /data/$file_genome_path

echo " "

#EGS=$(cut -b 1-10 -z ${sample_name}_max80_255_frg${AFL}_egSize.txt)
EGS=$(grep "number of unique k-mers" ${sample_name}_max80_255_frg${AFL}_egSize.txt | awk '{print $5}')

#echo "Total estimated number of unique ${AFL}-mers: $EGS"

echo " "
echo "Finish calculation of effective genome size"
echo " "
echo "############################"
echo " "




#### FOURTH PART ####

# Peak calling using Docker

echo " "
echo "Starting peak calling process"
echo " "

max_gap=$(printf "scale=10; %d*%d \n" 2 ${AFL} | bc)

# Crea la cartella di output se non esiste
mkdir -p outMACS3

# Esegui MACS3 nel container Docker
docker run --rm -u 1000:1000 -v "$(pwd):/data" quay.io/biocontainers/macs3:3.0.2--py310h397c9d8_2 \
    macs3 callpeak \
    -t /data/${sample_name}_max80_255.bam \
    -n ${sample_name} \
    --outdir /data/outMACS3 \
    -q 0.05 \
    -f BAM \
    -g ${EGS} \
    -s ${AFL} \
    --min-length ${AFL} \
    --max-gap ${max_gap} \
    --nomodel \
    --extsize ${AFL} \
    --keep-dup all \
    --buffer-size 10000000

echo " "
echo "Finish peak calling process"
echo " "
echo " "
