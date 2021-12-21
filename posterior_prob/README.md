# Calculation of posterior likelihoods/dosages


This script will add three fields to the VCF:

- DS: for each call the posterior dosage
- GP: for each call/genotype pair the posterior probability
- PGQ: for each call the posterior genotype quality

###Problem

Genotype probability vectors cannot be used 'as-is'.  This can be best be explained from the case where no 
reads are available: all genotypes are then equally likely, and therefore a variant caller gives a probability 
of 0.33/0.33/0.33 (diploid setting, reference, heterozygous, homozygous genotype). 

One can easily observe that these probabilities do not reflect our knowledge that this is a rare variant. 


###Solution

This script calculates genotype priors based on estimated population allele frequency, instead of assuming each genotype equally likely.
For example,  a low population frequency means also a lower prior probability for homozygous alternate genotypes. 
The approach is to estimate the allele frequency and to take this through Hardy-Weinberg to calculate a genotype prior. 




###Usage

python get_posterior_dosages.py {input} - | bcftools view -Ob -o {output} --threads 8

Where input and output are a BCF/VCF file. 

Dependencies: cyvcf2, numpy
Limitations: only for haploid and diploid genotype calls. Other calls are skipped.
