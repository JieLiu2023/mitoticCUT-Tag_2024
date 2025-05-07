# main scripts for *Protocol for identifying genomic binding sites of mitotic bookmarkers in Drosophila neural stem cells and cultured mammalian cells*

## Purposes
- This repository includes main scripts used for:
  - *Protocol for identifying genomic binding sites of mitotic bookmarkers in Drosophila neural stem cells and cultured mammalian cells*

## Steps
- 0.rawData_process.sh
  - This script is used to process CUT&Tag data (from fastq.gz to bam/bigWig)
- 1.DownSample.sh
  - This script is used for down-sampling to generate pseudo-replicates
- 2.callpeak_Macs2.sh
  - This script is used for peak-calling with Macs2.
- 3.IDRtest.sh
  -  This script is used to perform IDR-test on CUT&Tag peaks.
