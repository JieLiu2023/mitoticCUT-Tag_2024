#!/bin/bash

##################################################
## Project: Cut&Tag 2024
## Script purpose: Down-sample to produce pseudo-replicates
## Date: 2024-11-24
## Author: Jie Liu
##################################################

## Section: input args
#################################################

source ~/.bashrc
conda activate ${your_env}

project_dir=${your_project_dir}
process_dir="${project_dir}/data/"

ifDuplicate="rmDup"
minQualityScore=10

mkdir -p $process_dir/peakCalling/iceberg/tmp/

## Strategy 1: Down-sample from pooled file [optional, if read depths of different replicates vary greatly]
#################################################

# If you want to produce 3 down-sampled pseudo-replicates
group=${this_group_name}
order=(0 1 2)
SEED=(100 150 180)   # set random seeds

# for a group of biological replicates, like proteinA, IgG_for_A, proteinB, IgG_for_B, ...

samtools merge \
        ${process_dir}/iceberg/merge.{group}.process.bam
	${process_dir}/aligned/${sample_1}.${ifDuplicate}.QS${minQualityScore}.map.sorted.bam \
	${process_dir}/aligned/${sample_2}.${ifDuplicate}.QS${minQualityScore}.map.sorted.bam \
	...

## Calculate percent in advance: You can use samtools flagstat
percent=0.2			

for order in ${order[@]};
do

	picard DownsampleSam --ACCURACY 1.0E-5 --RANDOM_SEED ${SEED[${order}]} --STRATEGY ConstantMemory --TMP_DIR ${process_dir}/iceberg/tmp/ --PROBABILITY ${percent} -I ${process_dir}/iceberg/merge.${group}.process.bam -O ${process_dir}/iceberg/downSample.${group}.picard_DS_seed${SEED[${order}]}.bam

done

## Strategy 2: Down-sample from processed replicates
#################################################

# If you have 3 biological replicates in a group
group=${this_group_name}
percent=(0.20 0.25 0.3)  				# Just example, you should calculate in advance
sample=("sample_1" "sample_2" "sample_3")
order_sample=(0 1 2)

# If you want to produce 3 down-sampled pseudo-replicates
order_seed=(0 1 2)
SEED=(100 150 180)   					# set random seeds


for order1 in ${order_seed[@]};
do
	for order2 in ${order_sample[@]};
	do
		picard DownsampleSam --ACCURACY 1.0E-5 --RANDOM_SEED ${SEED[${order1}]} --STRATEGY ConstantMemory --TMP_DIR ${process_dir}/iceberg/tmp/ --PROBABILITY ${percent[${order2}]} -I ${process_dir}/iceberg/${sample[${order2}]}.${ifDuplicate}.QS${minQualityScore}.map.sorted.bam -O ${process_dir}/iceberg/tmp.${sample[${order2}]}.picard_DS_seed${SEED[${order1}]}.bam
	done
	
	samtools merge ${process_dir}/iceberg/downSample.${group}.picard_DS_seed${SEED[${order1}]}.bam ${process_dir}/iceberg/tmp.*.picard_DS_seed${SEED[${order1}]}.bam
done


