# main scripts for *TBP bookmarks and preserves neural stem cell fate memory by orchestrating local chromatin architecture*

## Purposes
- This repo includes main scripts used for:
  - Shen Y, Liu K, Liu J, et al. *TBP bookmarks and preserves neural stem cell fate memory by orchestrating local chromatin architecture*. Mol Cell. 2025
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
