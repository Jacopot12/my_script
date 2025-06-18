#!/bin/bash

#Author: Jacopo Tartaglia
#Verson: 1.0
#What the script does:
#The script prepares the output table of STAR countGene to use for the subsequent analysis
#It takes only the reverse strand count (the 3th column in the output table of STAR). Change the code to take other columns

#Usage:
# $ bash ./rna_count_table.sh rep1_control.tab rep2_control.tab rep3_control.tab rep1_stress.tab rep2_stress.tab rep3_stress.tab cross_name parental_genotype
#The name of the input tables MUST be with the name of the sample separated by a "_". Example: M1_countTable.tab
#The order of the samples MUST be rep1_control rep2_control rep3_control rep1_stress rep2_stress rep3_stress
#The name of the parental genotype MUST be the same of the one in the gene_id of the count tables of STAR

EXPECTED_ARGS=8
E_BADARGS=8
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Number of arugment is not correct"	
  echo "Usage: ./rna_count_table.sh rep1_control.tab rep2_control.tab rep3_control.tab rep1_stress.tab rep2_stress.tab rep3_stress.tab cross_name genotype"
  echo "The name of the input tables must be with the name of the sample separated by a "_". Example: M1_countTable.tab"
  echo "The order of the samples MUST be rep1_control rep2_control rep3_control rep1_stress rep2_stress rep3_stress"
  echo "The name of the parental genotype MUST be the same of the one in the gene_id of the count tables of STAR"
  exit $E_BADARGS
fi
echo " "
echo "REMEMBER!"
echo "The name of the input tables MUST be with the name of the sample separated by a "_". Example: M1_countTable.tab"
echo "The order of the samples MUST be rep1_control rep2_control rep3_control rep1_stress rep2_stress rep3_stress"
echo "The name of the parental genotype MUST be the same of the one in the gene_id of the count tables of STAR"
echo " "

#input files
file1=$1
file2=$2
file3=$3
file4=$4
file5=$5
file6=$6

#input info
cross=$7
geno=$8

name1=$(echo ${file1} | cut -d "_" -f 1)
name2=$(echo ${file2} | cut -d "_" -f 1)
name3=$(echo ${file3} | cut -d "_" -f 1)
name4=$(echo ${file4} | cut -d "_" -f 1)
name5=$(echo ${file5} | cut -d "_" -f 1)
name6=$(echo ${file6} | cut -d "_" -f 1)


mkdir output_countTable
cd ./output_countTable
touch check_line.txt


awk 'NR > 4 ' ../${file1} | awk '{print $1}' > gene_id.txt
wc -l gene_id.txt >> check_line.txt

for (( i=0; i<6; i++ )); do

    j=$(( $i + 1))
    eval "var_file=\$file$j"
    eval "var_name=\$name$j"
    awk 'NR > 4 ' ../"${var_file}" | awk '{print $3}' > ${var_name}_countTable_reverse.tsv
    wc -l ${var_name}_countTable_reverse.tsv >> check_line.txt

done

#control samples analysis
paste -d '\t' gene_id.txt ${name1}_countTable_reverse.tsv ${name2}_countTable_reverse.tsv ${name3}_countTable_reverse.tsv > countTable${cross}_control.tsv

grep "${geno}" countTable${cross}_control.tsv > tmp1.tsv

grep "Quncho" countTable${cross}_control.tsv | awk '{print $2"\t"$3"\t"$4}' > tmp2.tsv

paste -d '\t'  tmp1.tsv tmp2.tsv > tmp3.tsv

#awk -v FS="-" '{print $2}' tmp3.tsv > tmp4.tsv
sed "s/${geno}_//g" tmp3.tsv > tmp4.tsv 

echo -e "gene_id\t${name1}_${geno}\t${name2}_${geno}\t${name3}_${geno}\t${name1}_Quncho\t${name2}_Quncho\t${name3}_Quncho" > header.txt

cat header.txt tmp4.tsv > countTable${cross}_control_genotype.tsv

rm tmp* header.txt


#stress samples analysis
paste -d '\t' gene_id.txt ${name4}_countTable_reverse.tsv ${name5}_countTable_reverse.tsv ${name6}_countTable_reverse.tsv > countTable${cross}_stress.tsv

grep "${geno}" countTable${cross}_stress.tsv > tmp1.tsv

grep "Quncho" countTable${cross}_stress.tsv | awk '{print $2"\t"$3"\t"$4}' > tmp2.tsv

paste -d '\t'  tmp1.tsv tmp2.tsv > tmp3.tsv

#awk -v FS="-" '{print $2}' tmp3.tsv > tmp4.tsv
sed "s/${geno}_//g" tmp3.tsv > tmp4.tsv

echo -e "gene_id\t${name4}_${geno}\t${name5}_${geno}\t${name6}_${geno}\t${name4}_Quncho\t${name5}_Quncho\t${name6}_Quncho" > header.txt

cat header.txt tmp4.tsv > countTable${cross}_stress_genotype.tsv

rm tmp* header.txt *_reverse.tsv


#count matrix for WW vs DS conditions

awk -v OFS="\t" '{print $2, $3, $4}' countTable${cross}_stress.tsv > temp1.txt

paste -d '\t' countTable${cross}_control.tsv temp1.txt > temp2.txt

echo -e "gene_id\t${name1}\t${name2}\t${name3}\t${name4}\t${name5}\t${name6}" > header.txt

cat header.txt temp2.txt > countMatrixWW_DS.tsv

rm temp1.txt temp2.txt header.txt

echo "DONE IT"
