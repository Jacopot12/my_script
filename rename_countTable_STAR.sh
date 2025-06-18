#!/bin/bash

#Usage: bash ./rename_countTable_STAR.sh CONTROL_REP1.ReadsPerGene.out.tab CONTROL_REP2.ReadsPerGene.out.tab CONTROL_REP3.ReadsPerGene.out.tab STRESS_REP1.ReadsPerGene.out.tab STRESS_REP2.ReadsPerGene.out.tab STRESS_REP3.ReadsPerGene.out.tab cross_name
#Simple bash script to re-name the output file of STAR countGene table

EXPECTED_ARGS=7
E_BADARGS=7
if [ $# -ne $EXPECTED_ARGS ]
then
  echo " "
  echo "Number of arugment is not correct! Required 7 arguments."
  echo " "
  echo "Usage: ./rename_countTable_STAR.sh CONTROL_REP1.ReadsPerGene.out.tab CONTROL_REP2.ReadsPerGene.out.tab CONTROL_REP3.ReadsPerGene.out.tab STRESS_REP1.ReadsPerGene.out.tab STRESS_REP2.ReadsPerGene.out.tab STRESS_REP3.ReadsPerGene.out.tab cross_name"
  echo "The order of the samples MUST be rep1_control rep2_control rep3_control rep1_stress rep2_stress rep3_stress"
  echo " "
  exit $E_BADARGS
fi

file1=$1
file2=$2
file3=$3
file4=$4
file5=$5
file6=$6

cross=$7

#name the control files
mv ${file1} ${cross}Control1_countTable.tab
mv ${file2} ${cross}Control2_countTable.tab
mv ${file3} ${cross}Control3_countTable.tab

#name the stress files
mv ${file4} ${cross}Stress1_countTable.tab
mv ${file5} ${cross}Stress2_countTable.tab
mv ${file6} ${cross}Stress3_countTable.tab
