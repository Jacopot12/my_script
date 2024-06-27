#!/bin/bash

#What the script does:
#The script insert a variable number of SNPs in a file fasta with several sequences
#The aim is to generate a mock dataset of sequencing reads from an original one
#The script does not allow an SNP to be entered twice or more at the same location, but it allow a base substitution with the same base

#Usage: ./randomSNPs.sh [file.fasta] [% SNPs mutation] [outprefix] [compression]
#Example: ./randomSNPs.sh mysequneces.fasta 0.2 out c
#The input and output files are .fasta, input file could be compressed
#The output file could be compressed if [compression] parameter is equal to "c". If leave empty no compression occurred
#The % mutation can range from 0 to 1
#For exaple: 0.2 (or 20%) of SNPs mutation is equal to 20 SNPs on a sequence of 100 bp. The output file will be out.fasta

# Input file .fasta
input_file=$1

# % SNPs mutation can range from 0 (0%) to 1 (100%)
SNPrate=$2

#output file name (estension excluded)
outprefix=$3

#output file compression
compr=$4

file_type=$(file --mime-type -b "$input_file")

decompressed_file="$input_file"

# Make a temporary directory
temp_dir=$(mktemp -d)

# Check if the input file is compressed
case "$file_type" in
    application/gzip)
        decompressed_file="$temp_dir/$(basename "$input_file" .gz)"
        gunzip -c "$input_file" > "$decompressed_file"
        input="$decompressed_file"
        ;;
    application/x-bzip2)
        decompressed_file="$temp_dir/$(basename "$input_file" .bz2)"
        bunzip2 -c "$input_file" > "$decompressed_file"
        input="$decompressed_file"
        ;;
    application/x-xz)
        decompressed_file="$temp_dir/$(basename "$input_file" .xz)"
        xzcat "$input_file" > "$decompressed_file"
        input="$decompressed_file"
        ;;
    *)
        input=$1
        ;;
esac


bases=("A" "T" "C" "G")

if (( $(echo "$SNPrate >= 0" | bc -l) )) && (( $(echo "$SNPrate <= 1" | bc -l) )); then #check if the SNPs mutation is valid
    
    if [[ -f "${outprefix}.fasta" ]] || [[ -f "${outprefix}.fasta.gz" ]]; then #check if the output file already exists

        echo "Error: the file '${outprefix}.fasta' or '${outprefix}.fasta.gz' already exists"
        exit 1

    else

        touch $outprefix.fasta
    
    fi


    while read -r line || [ -n "$line" ]; do #if a line is the header of the fasta sequence append to the output file
    
        if [[ $line == \>* ]] || [[ -z $line ]]; then

            echo $line >> $outprefix.fasta

        else

            length=${#line}                                             #if the line contains the sequence add the mutation based on its length
            result=$(echo "$length * $SNPrate" | bc)
            rounded_result=$(echo "scale=0; ($result+0.999999)/1" | bc)

        
            selected_numbers=()

            while [[ ${#selected_numbers[@]} -lt $rounded_result ]]; do
                random_number=$(( RANDOM % length ))
                unique=true

                for i in "${selected_numbers[@]}"; do      #add a mutation twice or more in the same position is not permitted
                    if [[ $i -eq $random_number ]]; then
                        unique=false
                        break
                    fi
                done

                if $unique; then
                    selected_numbers+=($random_number)
                    index=$(( RANDOM % ${#bases[@]} ))
                    new_char=${bases[$index]}
                    prefix=${line:0:random_number}
                    suffix=${line:random_number+1}
                    line="${prefix}${new_char}${suffix}"
                fi
            done

            echo $line >> $outprefix.fasta
        fi
    done < "$input"

else
    echo "SNPs rate is not permitted! Try again."
fi


#compress the output file if needed
if [[ "$compr" == "c" ]]; then

    gzip $outprefix.fasta

fi

