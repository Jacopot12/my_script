#!/bin/bash

#usage: bash make_variants_file.sh file.vcf cultivar_name /absolute/path/to/reference_genome.fa /absolute/path/to/annotation.fa 

EXPECTED_ARGS=4
E_BADARGS=4
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: bash ./make_variants_file.sh file.vcf cultivar_name /absolute/path/to/reference_genome.fa /absolute/path/to/annotation.fa"
  exit $E_BADARGS
fi


### DEFINE VARIABLES ###
input_vcf=$1
name=$2
reference=$3
annotation=$4


### PARSING THE VCF FILE ###

#take the header
grep '^#' ${input_vcf} > header.txt

#take the data
grep -v '^#' ${input_vcf} > data.txt

#take veriants that pass the filter
awk '($7 == "PASS") {print $0}' data.txt > data_PASS.txt

#take only 1/1 sites
grep -w "1/1" data_PASS.txt > data_PASS_homo_NoHeader.txt

#make the filtered vcf file
cat header.txt data_PASS_homo_NoHeader.txt > ${name}_PASS_homo.vcf

#number of final variants
numVar=$(wc -l data_PASS_homo_NoHeader.txt)
echo "The number of variants is: ${numVar}"

#remove temporary files
rm *.txt

#sort the file with bcftools
source $HOME/anaconda3/bin/activate bcftools
bcftools sort -o ${name}_PASS_homo_sorted.vcf ${name}_PASS_homo.vcf
conda deactivate



### CREATE THE PSEUDOREFERENCE GENOME ###

source $HOME/g2gtools/bin/activate

bgzip ${name}_PASS_homo_sorted.vcf

echo " "
echo "MAKING THE VCI FILE"
echo " "

g2gtools vcf2vci --vcf ${name}_PASS_homo_sorted.vcf.gz --fasta ${reference} --strain SAMPLE --vci out.vci --gtf ${annotation}

echo " "
echo "PATCHING THE FASTA FILE"
echo " "

g2gtools patch --fasta ${reference} --vci out.vci.gz --out reference.patched.fa

echo " "
echo "TRANSFORM THE FASTA FILE"
echo " "

g2gtools transform --fasta reference.patched.fa --vci out.vci.gz --out ${name}_cultivar.fa

echo " "
echo "CONVERT THE ANNOTATION FILE"
echo " "

g2gtools convert --in ${annotation} --vci out.vci.gz --out ${name}_cultivar.gtf -f GTF

rm reference.patched.fa reference.patched.fa.fai
gunzip ${name}_PASS_homo_sorted.vcf.gz
