#!/bin/bash

# Author: Jacopo Tartaglia
# Verson: 1.0
# What the script does:
# This script carry out the Quality Control on raw paired ends reads, the trimming (with SeqPurge) of them and carry out another QC on the trimmed reads
# The outputs of fastQC before and after the trimming step are merged by MultiQC
# The parameters for the trimming must be modified directly in this script


# Usage: ./readsQC_trimming.sh <reads_foreward.fq.gz> <reads_reverse.fq.gz> <sample_name>
# The files order MUST be respected

#Programs required (each program is located in its own environment conda):
# - fastQC
# - ngs-bits
# - multiQC 

EXPECTED_ARGS=3
E_BADARGS=3
if [ $# -ne $EXPECTED_ARGS ]
	then
  	echo "Usage: ./readsQC_trimming.sh <reads_foreward.fq.gz> <reads_reverse.fq.gz> <sample_name>"
  	echo "The files order MUST be respected!"
  	exit $E_BADARGS
fi


# input files
readsF=$1
readsR=$2
sample=$3


#function to activate / deactivate conda environment
activate_conda(){
  local ENV_NAME=$1
  source /opt/anaconda3/bin/activate $ENV_NAME #CHANGE THIS ACCORDING TO YOUR CONDA PATH

  if [[ $(conda info --envs | grep "*" | awk '{print $1}') == $ENV_NAME ]]; then
    echo " "
    echo "The $ENV_NAME conda environment has been successfully activated"
    echo " "
  else
    echo " "
    echo "Error in activating the environment $ENV_NAME."
    exit 1
  fi
}


deactivate_conda() {
    local ENV_NAME=$1
    conda deactivate

    if [[ $(conda info --envs | grep "*" | awk '{print $1}') == "base" ]]; then
        echo "The $ENV_NAME conda environment has been successfully deactivated."
    else
        echo "Error in deactivating the environment $ENV_NAME."
        exit 1
    fi
}



#### FIRST PART ####

# quality control before trimming

activate_conda "fastqc"

if [ ! -d "./QC_beforeTrim" ]; then
  mkdir QC_beforeTrim
fi

cd QC_beforeTrim

echo "Start running FastQC, quality checking for ${readsF} and ${readsR}"

fastqc ../${readsF} ../${readsR} -o ./

deactivate_conda "fastqc"
echo "Finish of quality check step"
cd ../


#### SECOND PART ####

# trimming of raw reads

activate_conda "ngs-bits"

if [ ! -d "./SeqPurge" ]; then
  mkdir SeqPurge
fi

cd ./SeqPurge

echo "Start running SeqPurge, trimming step for ${readsF} and ${readsR}"

nameF=$(echo ${readsF} | cut -d "." -f 1)
nameR=$(echo ${readsR} | cut -d "." -f 1)

SeqPurge \
	-min_len 20 \
	-threads 4 \
	-qcut 25 \
  -qwin 10 \
	-in1 ../${readsF} \
	-in2 ../${readsR} \
	-out1 ${nameF}_trim.fq.gz \
	-out2 ${nameR}_trim.fq.gz \
	-summary ${sample}.trim.stats \
	-compression_level 5

deactivate_conda "ngs-bits"
echo "Finish of trimming step"
cd ../


#### THIRD PART ####

# quality control after trimming

activate_conda "fastqc"

if [ ! -d "./QC_afterTrim" ]; then
  mkdir QC_afterTrim
fi

cd QC_afterTrim

echo "Start running FastQC, quality checking for ${readsF} and ${readsR} after trimming"

fastqc ../SeqPurge/${nameF}_trim.fq.gz ../SeqPurge/${nameR}_trim.fq.gz -o ./

deactivate_conda "fastqc"
echo "Finish of quality check step after trimming"
cd ../



#### FOURTH PART ####

# quality control summary

activate_conda "multiqc"

if [ ! -d "./multiQC" ]; then
  mkdir multiQC
fi

cd multiQC

echo "Start running MultiQC"

multiqc ../QC_beforeTrim/ ../QC_afterTrim/

deactivate_conda "multiqc"
echo "Finish of QC and trimming pipeline"
cd ../





















