#!/bin/bash
# -------------------------------------------------------------------------------
# Split the A. thaliana data into subsets
# 
# Example: 
# > split_dataset_athaliana.sh ../data/high_order_epistasis/avrB 40
# will split X.dat into 5 x 40 = 200 subdatasets
# First, each of the five chromosomes
# and then, each chromosome downsampled by a factor of 40, with different 40 offsets
#
# F. Llinares-Lopez
# -------------------------------------------------------------------------------

DATASET=$1
DOWNSAMPLING_FACTOR=$2
N_CHR=5

cut -d ' ' -f1 ${DATASET}/plink.map > ${DATASET}/chr.txt
paste ${DATASET}/chr.txt ${DATASET}/X.dat -d " " > ${DATASET}/X_with_chr.dat
awk -v dataset=${DATASET} -F" " '{ for (i=2; i<NF; i++) printf $i " " >  dataset"/X_"$1".dat"; print $NF >  dataset"/X_"$1".dat"}' ${DATASET}/X_with_chr.dat
rm ${DATASET}/chr.txt ${DATASET}/X_with_chr.dat


for chr in `seq 1 ${N_CHR}`;
do
	for offset in `seq 0 $((DOWNSAMPLING_FACTOR - 1))`;
	do
		awk -v downsampling_factor=${DOWNSAMPLING_FACTOR} -v offset=${offset} 'NR % downsampling_factor == offset' ${DATASET}/X_Chr${chr}.dat > ${DATASET}/X_Chr${chr}_Offset${offset}.dat
	done
done
