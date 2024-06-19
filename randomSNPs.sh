#!/bin/bash

#What the script does:
#The script insert a variable number of SNPs in a file fasta with several sequences
#The aim is to generate a mock dataset of sequencing reads from an original one

#Usage: ./randomSNPs.sh [file.fasta] [% SNPs mutation] [outprefix]
#Example: ./randomSNPs.sh mysequneces.fasta 0.2 out
#The input and output files are .fasta
#The % mutation can range from 0 to 1
#For exaple: 0.2 (or 20%) of SNPs mutation is equal to 20 SNPs on a sequence of 100 bp. The output file will be out.fasta

# Input file .fasta
input_file=$1

# % SNPs mutation can range from 0 (0%) to 1 (100%)
SNPrate=$2

#outout file name (estension excluded)
outprefix=$3

bases=("A" "T" "C" "G")

if (( $(echo "$SNPrate >= 0" | bc -l) )) && (( $(echo "$SNPrate <= 1" | bc -l) )); then #check if the SNPs mutation is valid

    touch $outprefix.fasta

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
    done < "$input_file"

else
    echo "SNPs rate is not permitted! Try again."
fi
