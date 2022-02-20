

##2022-PLOS1-paper

This repository is for a paper under review. Raw data FASTQ files from all cell line samples used in this study are available via BioProject PRJNA803819 in SRA. SRA records will be accessible with the following link upon publication:
https://www.ncbi.nlm.nih.gov/sra/PRJNA803819 Special reviewer link to those data was provided to reviewers via the Journal for immediate access before public release. Two separate Snakemake workflow files were provided for raw data QC and variant calling respectively in the workflows directory. Both could be tested with "snakemake -n" dry run command in the corresponding subdirectories. Human genome reference hg19 should be obtained through GATK bundle. All designated software tools should be obtained via their original providers. A script in R was provided for the creation of all figures, except for Figure 1B, which was composed of IGV screenshots.
