#!/bin/bash

##################################################
## Project: Cut&Tag for 
## Script purpose: data QC and process
## Date: 2024-11-24
## Author: Jie Liu
##################################################

## Section: input args
#################################################

source ~/.bashrc
conda activate ${your_env}

project_dir=${your_project_dir}

QC_dir=${your_QC_dir}  

ifDuplicate="rmDup"
minQualityScore=10

sample_list=${your sample list}

raw_dir="${project_dir}/raw/"
process_dir="${project_dir}/data/"
log_dir="${project_dir}/logs/"
result_dir="${project_dir}/results/"

dm6_bowtie2_index=${your_dm6_bowtie2_index}
dm6_chromSize=${your_dm6_chromSize}

hg38_bowtie2_index=${your_hg38_bowtie2_index}
hg38_chromSize=${your_hg38_chromSize}

# If spike-In is in need
Ecoli_index=${your_Ecoli_index}

mkdir -p $process_dir
mkdir -p $log_dir
mkdir -p $result_dir
mkdir -p $result_dir/${QC_dir}_afterTrim
mkdir -p $process_dir/cleanData
mkdir -p $process_dir}aligned/bowtie2_summary
mkdir -p $process_dir/aligned/checkDuplicates
mkdir -p $process_dir/aligned/fragmentLen
mkdir -p $process_dir/aligned/bed
mkdir -p $process_dir/aligned/bigwig
mkdir -p $process_dir/peakCalling/

## Section: Data QC 
##################################################

for sample in ${sample_list};
do

	trim_galore --paired --output_dir ${process_dir}/cleanData/ --quality 20 --stringency 3 -j 4 ${raw_dir}/${sample}/${sample}_1.fq.gz ${raw_dir}/${sample}/${sample}_2.fq.gz

done

echo -e "================ Trim_galore Finished  ==============\n"


fastqc -t 8 -o ${result_dir}/${QC_dir}_afterTrim \
    ${process_dir}/cleanData/*.fq.gz

multiqc -n ${QC_dir}_afterTrim.multiQC.html \
        -o ${result_dir}/${QC_dir}_afterTrim \
        ${result_dir}/${QC_dir}_afterTrim

rm ${result_dir}/${QC_dir}_afterTrim/*fastqc.zip
rm ${result_dir}/${QC_dir}_afterTrim/*fastqc.html


echo -e "================ QC & Trim finished ==============\n"


## Section: Bowtie2 Alignment
##################################################

for sample in ${sample_list};
do

	bowtie2 --local --very-sensitive-local --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 8 -x ${dm6_bowtie2_index} -1 ${process_dir}/cleanData/${sample}_1_val_1.fq.gz -2 ${process_dir}/cleanData/${sample}_2_val_2.fq.gz -S ${process_dir}/aligned/${sample}.sam 2>${process_dir}/aligned/bowtie2_summary/${sample}_bowtie2.summary

done

echo -e "================ Bowtie2 Finished! ==============\n"

## Section: Check Duplicates Rate 
##################################################

for sample in ${sample_list};
do
	# Sort by coordinate
	picard SortSam -I ${process_dir}/aligned/${sample}.sam -O ${process_dir}/aligned/${sample}.picard_sorted.sam --SORT_ORDER coordinate
	# Mark Duplicates
	# picard MarkDuplicates -I ${process_dir}/aligned/${sample}.picard_sorted.sam -O ${process_dir}/aligned/${sample}.picard_sorted.MarkDup.sam -M ${process_dir}/aligned/checkDuplicates/${sample}.picard_sorted.MarkDup.summary

	# remove Duplicates, but it is optional
	picard MarkDuplicates -I ${process_dir}/aligned/${sample}.picard_sorted.sam -O ${process_dir}/aligned/${sample}.picard_sorted.rmDup.sam --REMOVE_DUPLICATES true -M ${process_dir}/aligned/checkDuplicates/${sample}.picard_sorted.rmDup.summary

	rm ${process_dir}/aligned/${sample}.sam
	rm ${process_dir}/aligned/${sample}.picard_sorted.sam

done

echo -e "================ Picard check Duplicates End! ==============\n"

## Section: Filtration
##################################################

for sample in ${sample_list};
do
	### Filter and keep the mapped read pairs
	samtools view -bS -q $minQualityScore -F 0x04 ${process_dir}/aligned/${sample}.picard_sorted.${ifDuplicate}.sam > ${process_dir}/aligned/${sample}.${ifDuplicate}.QS${minQualityScore}.map.bam

	samtools sort ${process_dir}/aligned/${sample}.${ifDuplicate}.QS${minQualityScore}.map.bam -o ${process_dir}/aligned/${sample}.${ifDuplicate}.QS${minQualityScore}.map.sorted.bam
	rm ${process_dir}/aligned/${sample}.${ifDuplicate}.QS${minQualityScore}.map.bam

	samtools view -bS ${process_dir}/aligned/${sample}.picard_sorted.${ifDuplicate}.sam > ${process_dir}/aligned/${sample}.picard_sorted.${ifDuplicate}.bam
	rm ${process_dir}/aligned/${sample}.picard_sorted.${ifDuplicate}.sam

done

echo -e "================ Samtools filtration Finish! ===============\n"

## Section: Visualization
##################################################

for sample in ${sample_list};
do

        samtools index ${process_dir}/aligned/${sample}.${ifDuplicate}.QS${minQualityScore}.map.sorted.bam
	bamCoverage --normalizeUsing RPKM --extendReads -b ${process_dir}/aligned/${sample}.${ifDuplicate}.QS${minQualityScore}.map.sorted.bam -o ${process_dir}/aligned/bigwig/${sample}.${ifDuplicate}.QS${minQualityScore}.RPKM.extendReads.bw
        
done


