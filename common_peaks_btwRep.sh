#!/bin/bash

# This script selects common peaks between 2 replicates.
# Pay attention to the reciprocal overlap threshold in bedtools intersect (flags -f 0.9 -r). 
# Change these parameter if you want more/less stringent reciprocal overlap threshold.

# usage: bash common_peaks_btwRep.sh <file1.bed> <file2.bed> <chr.faidx>
# file1.bed is the file with the peaks of the first replicate (it must have only 3 colums [chr \t start \t end])
# file2.bed is the file with the peaks of the second replicate (it must have only 3 colums [chr \t start \t end])
# chr.faidx is the file with the name of the chromosomes in the correct order


peaks1=$1
peaks2=$2
chrs=$3

# Sort files by chromosomes 

bedtools sort -i "$peaks1" -faidx "$chrs" > file1_sorted.bed
bedtools sort -i "$peaks2" -faidx "$chrs" > file2_sorted.bed


# I take the common peaks that have a minimum 90% mutual overlap
# Change the value of the -f flag to choose a higher or lower overlap threshold
# Print the coordinates of the peaks that overlap next to each other  

bedtools intersect -a file1_sorted.bed -b file2_sorted.bed -wa -wb -f 0.9 -r > commonPeaks_temp1.bed

# Separate the common peaks into two different files and then put them into one file to do the merging the pairs of commin peaks
# I do that in order to have the widest common peak 

awk '{print $1"\t"$2"\t"$3}' commonPeaks_temp1.bed | bedtools sort -faidx "$chrs" > commonPeaks_temp2_1.bed
awk '{print $4"\t"$5"\t"$6}' commonPeaks_temp1.bed | bedtools sort -faidx "$chrs" > commonPeaks_temp2_2.bed
cat commonPeaks_temp2_1.bed commonPeaks_temp2_2.bed | bedtools sort -faidx "$chrs" | bedtools merge > commonPeaks_merged.bed

rm file1_sorted.bed file2_sorted.bed commonPeaks_temp1.bed commonPeaks_temp2_1.bed commonPeaks_temp2_2.bed


