# QC scripts for Illumina sequencing data

This repository contains sequencing QC utilities that have been used for the paper:

"Exome sequencing identifies rare damaging variants in ATP8B4 and ABCA1 as novel risk factors for Alzheimer’s Disease".

The QC procedure used can be found in the supplement of this paper. 

- bcftools_split_patches:  These patches adapt the splitting operation implemented by BCFTools, and introduce handling for the PGT (phased-GT), PL (phred-likelihood) and AD (read by allele depth) tags.
- dechimerizer: This script will remove soft-clip overhangs of BAM alignments that are (likely) to have been caused by chimeric reads. Includes also a script to generate detailed statistics for BAM files on supplementary reads.
- posterior_prob: Calculation of posterior likelihoods/dosages and adding these to the VCF
- variant_qc_annot: This script annotates the VCF file with a number of QC annotations, for different variant coverage/quality thresholds, and for genotype 0/1 or posterior probabilistic calls.


# Burden analysis script based on Ordinal Logistic Regression

- burden_test:  this script contains a stand-alone version of the burden test as used in our publication, with variant prioritization based on VEP/Loftee/transcripts, and with genotype sampling based on posterior probabilistic calls.


