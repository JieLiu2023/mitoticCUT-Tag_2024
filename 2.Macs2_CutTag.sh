#!/bin/bash

##################################################
## Project: Cut&Tag 2024
## Script purpose: Macs2 peakCalling
## Date: 2024-11
## Author: Jie Liu
##################################################

source ~/.bashrc
conda activate ${your_env}

project_dir=${your_project_dir}
process_dir="${project_dir}/data/"
mkdir -p $process_dir/peakCalling/

ifDuplicate="rmDup"
minQualityScore=10

treat_group="proteinA"		# example
IgG_group="IgG"			# example

## Section: Macs2 peakcalling based on down-sampled bam files
##################################################

# If you have produced 3 down-sampled pseudo-replicates
order_seed=(0 1 2)
SEED=(100 150 180)

for order in ${order[@]};
do
	macs2 callpeak -t ${process_dir}/iceberg/downSample.${treat_group}.picard_DS_seed${SEED[${order}]}.bam -c ${process_dir}/iceberg/downSample.${IgG_group}.picard_DS_seed${SEED[${order}]}.bam -f BAMPE -g dm -q 0.05 --outdir ${process_dir}/peakCalling/ -n Macs2.${treat_group}.picard_DS_seed${SEED[${order}]}.q005

done

