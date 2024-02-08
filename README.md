# txtools use cases code

This repository contains the Rmd files needed to recreate the use cases presented in
the manuscript "txtools: an R package facilitating analysis of RNA modifications,
structures, and interactions".

It contains two directories with two versions of the code:

- **abridged_version**: Which downloads and loads already processed data by txtools 
and runs the plotting functions. This is ideal for quickly recreating the 
results and inspecting the resulting **txtools' summarized data tables (txDT)**.
- **full_version**: This version of the code fully recreates the processing and 
analysis from the FASTQ files' original data. It will attempt to download the FASTQ
files from the SRA, align the RNA-seq reads to the reference sequence with the
 *Rsubread* package, and process the resulting mapped reads (BAM files) using
txtools. Running this code will require a high amount of computational resources
and time, but is intended for completion and reproducibility of our results.
    - Additionally, a supplementary analysis has been added to show the 
    difference and importance of properly processing paired-end reads to reflect
    accurate coverage values across the transcriptome. This is achieved running 
    the `supp_PairedvsSingle.R` script.
