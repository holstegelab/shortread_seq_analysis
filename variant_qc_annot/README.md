# Variant annotator 

This script annotates the VCF file with a number of annotations, for different
variant coverage/quality thresholds, and for 0/1 or posterior probabilistic calls:

###Currently considered thresholds:

Digital 0/1 calls (GT field):
 - All calls
 - Depth >= 6
 - Depth >= 10
 - Genotype Quality (GQ) >= 2

Posterior probabilistic (see posterior_prob folder):
 - All calls
 - Depth >= 10

This set can be easily adapted as needed in the script.  

###Variant annotations (per setting described above):

- kept: fraction of calls kept by filter
- ac:   count of minor alleles
- gt: string with count of reference, heterozygous and homozygous genotype 
- dpref: string with sum of (number of reference reads times probablity of respectively reference, heterozygous and homozygous genotype). 
- dpalt: string with sum of (number of alternate reads times probablity of respectively reference, heterozygous and homozygous genotype). 
- abhet: allele balance heterozygous calls (0.5 = normal, 1.0 = reference bias)
- abhom: allele balance homozygous calls (1.0 = normal, 0.0 = reference bias)
- fdhet:  reads on het calls vs. all calls.  (ratio)
- abdev: not used
- incoef: inbreeding coefficient
- phred-scaled hardy-weinberg p-value
- control_incoef: inbreeding coefficient of controls
- phred-scaled hardy-weingberg p-value of controls
- althet: average alt read depth per alt call
- total_althet: total alt read depth for alt calls

Usage:
    var_qc_annot.py {output} {METAFILE} {input}

    Or when using a pipe:

    samtools view {input} | PYTHONPATH={PYTHONPATH}  {MYPATH}/var_qc_annot.py {output} {METAFILE} -

Dependencies:
- cyvcf2
- numpy
- rescale_af.py from posterior_prob folder
