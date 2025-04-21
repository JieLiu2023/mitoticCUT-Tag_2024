#!/bin/bash

##################################################
## Project: Mitotic CutTag
## Script purpose: IDR test for Macs2 peakcalling
## Date: 2024-11
## Author: Jie Liu
##################################################

source ~/.bashrc
conda activate ${your_env}

project_dir=${your_project_dir}

process_dir="${project_dir}/data/"
peakCalling_dir="${project_dir}/data/peakCalling/"

mkdir ${peakCalling_dir}/tmp
mkdir ${peakCalling_dir}/IDR

treat_group=${your_treat_groupName}

## Section: Sort narrowPeak based on qValue
##########################################

# If you have produced 3 down-sampled pseudo-replicates
order_seed=(0 1 2)
SEED=(100 150 180)

for order in ${order_seed[@]};
do
	sort -k9,9nr ${peakCalling_dir}/Macs2.${treat_group}.picard_DS_seed${SEED[${order}]}.q005_peaks.narrowPeak >> ${peakCalling_dir}/tmp/Macs2.${treat_group}.picard_DS_seed${SEED[${order}]}.q005_sorted.narrowPeak	

done

## Section: IDR test in pairs
##########################################

idr --samples \
	${peakCalling_dir}/tmp/Macs2.${treat_group}.picard_DS_seed${SEED[0]}.q005_sorted.narrowPeak \
        ${peakCalling_dir}/tmp/Macs2.${treat_group}.picard_DS_seed${SEED[1]}.q005_sorted.narrowPeak \
        --input-file-type narrowPeak --rank q.value --plot \
        --output-file ${peakCalling_dir}/IDR/${treat_group}.rep01.IDR.bed \
        --log-output-file ${peakCalling_dir}/IDR/${treat_group}.rep01.IDR.log

idr --samples \
        ${peakCalling_dir}/tmp/Macs2.${treat_group}.picard_DS_seed${SEED[0]}.q005_sorted.narrowPeak \
        ${peakCalling_dir}/tmp/Macs2.${treat_group}.picard_DS_seed${SEED[2]}.q005_sorted.narrowPeak \
        --input-file-type narrowPeak --rank q.value --plot \
        --output-file ${peakCalling_dir}/IDR/${treat_group}.rep02.IDR.bed \
        --log-output-file ${peakCalling_dir}/IDR/${treat_group}.rep02.IDR.log

idr --samples \
        ${peakCalling_dir}/tmp/Macs2.${treat_group}.picard_DS_seed${SEED[1]}.q005_sorted.narrowPeak \
        ${peakCalling_dir}/tmp/Macs2.${treat_group}.picard_DS_seed${SEED[2]}.q005_sorted.narrowPeak \
        --input-file-type narrowPeak --rank q.value --plot \
        --output-file ${peakCalling_dir}/IDR/${treat_group}.rep12.IDR.bed \
        --log-output-file ${peakCalling_dir}/IDR/${treat_group}.rep12.IDR.log

## Section: Filter based on IDR value
##########################################
# The output file contains the scaled IDR value (min(int(-125*log2(IDR), 1000)) in the 5th field.
# If choose 0.05 as the IDR threshold to identify "reproducible" peaks, then this metric must be at least 540.
# With an IDR threshold defined as 0.5, which means scaled IDR value in 5th field > 125

pairs="01 02 12"

for pair in ${pairs};
do
	awk '$5 >= 125 {print $1"\t"$2"\t"$3}' ${peakCalling_dir}/IDR/${treat_group}.rep${pair}.IDR.bed > ${peakCalling_dir}/IDR/${treat_group}.rep${pair}.IDR_filter0.5.bed
done

## Section: Merge
##########################################

cut -f1-3 ${peakCalling_dir}/IDR/${treat_group}.rep*.IDR_filter0.5.bed | sort -k1,1 -k2,2n | bedtools merge >> ${peakCalling_dir}/IDR/bedtools_merge.${treat_group}.IDR_filter0.5.bed


