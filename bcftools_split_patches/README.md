# BCFTools patches to multi-allelic variant split operation

Multi-allelic variants are usually split into bi-allelic variants before analysis. 

These patches adapt the splitting operation implemented by BCFTools, and introduce
handling for the PGT (phased-GT), PL (phred-likelihood) and AD (read by allele depth) tags. 

The main idea behind these patches is to perform variant splitting such that the resulting bi-allelic variants
represent NON-ALT vs. ALT variants, instead of the regular REF vs. ALT variants model used by BCFTools. 

The disadvantage of the latter model is that in case of multi-allelic variants, a carrier might be neither REF nor ALT.
This results in a situation in which all genotypes are unlikely, but still one of the genotypes is considered the 
'more likely one'. This will affect downstream analyses. A NON-ALT vs. ALT split does not have this problem. 

A more detailed description of the rules followed by these patches is given in the supplement to the paper:
"Exome sequencing identifies rare damaging variants in ATP8B4 and ABCA1 as novel risk factors for Alzheimer’s Disease".

These patches have been developed with respect to BCFTools version 1.8. 

