# Utilities for burden-analysis of short-read sequencing data

This repository contains sequencing scripts that have been used for the paper:

"Exome sequencing identifies rare damaging variants in ATP8B4 and ABCA1 as novel risk factors for Alzheimer’s Disease", Nature Genetics 2022.


## QC scripts 

The QC procedure used can be found in the supplement of this paper. 

- bcftools_split_patches:  These patches adapt the splitting operation implemented by BCFTools, and introduce handling for the PGT (phased-GT), PL (phred-likelihood) and AD (read by allele depth) tags.
- dechimerizer: This script will remove soft-clip overhangs of BAM alignments that are (likely) to have been caused by chimeric reads. Includes also a script to generate detailed statistics for BAM files on supplementary reads.
- posterior_prob: Calculation of posterior likelihoods/dosages and adding these to the VCF
- variant_qc_annot: This script annotates the VCF file with a number of QC annotations, for different variant coverage/quality thresholds, and for genotype 0/1 or posterior probabilistic calls.


## Ordinal Logistic Regression burden analysis

- burden_test:  a stand-alone version of the burden test as used in the publication. Variant prioritization based on VEP/Loftee/transcripts, and with genotype sampling based on posterior probabilistic calls.


### Latest release

[![DOI](https://zenodo.org/badge/440135941.svg)](https://zenodo.org/badge/latestdoi/440135941)

