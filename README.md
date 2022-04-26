

2022 PLoS One paper

This repository is a part of the data and method sharing for the research article by Shi et al., "Development and evaluation of ActSeq: a targeted next-generation sequencing panel for clinical oncology use", PLoS ONE 2022, (https://doi.org/10.1371/journal.pone.0266914). Please refer to the full paper for proper use of scripts and data contained here.

Two separate Snakemake (v6.14.0) workflow files were provided for raw data QC and variant calling respectively in the workflows directory. Both workflows take FASTQ files as input. Raw data FASTQ files from all cell line samples used in this study are available via BioProject PRJNA803819 in SRA (https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA803819). Workflow files (Snakefile) can be tested with the "snakemake -n" dry run command in the corresponding subdirectories. The expected final output is a multiqc_report.html file summarizing all samples used for QC workflow. The variant calling workflow produces one VCF file containing the variant calls and two text files containing the NGS panel metrics extracted from BAM file for the individual input sample.

ActSeq panel-specific data, such as bait sequences and BED files are located in the data directory. The human genome reference hg19 can be obtained through GATK bundle (https://gatk.broadinstitute.org/). All designated software tools can be obtained via their original providers.

A script in R was provided for the creation of all figures, except for Figure 1B, which was composed of IGV screenshots. Source data are located in the data directory.
