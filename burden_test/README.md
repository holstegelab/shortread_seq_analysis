# Logistic ordinal burden test implementation 

Stand-alone version of the burden testing script used in the manuscript "Exome sequencing identifies rare damaging variants in ATP8B4 and ABCA1 as novel risk factors for Alzheimer’s Disease", adapted for stand-alone use.

## Installation:


1) If conda is not yet available: install miniconda and start a new shell

2) If mamba is not yet available:  `conda install mamba`     #mamba increases speed of installation

3) Download repository, move to an appropriate folder and execute:  `git clone https://github.com/holstegelab/shortread_seq_analysis.git`       

4) Enter folder: `cd shortread_seq_analysis/burden_test`

5) Create environment:  `mamba env create -f burden.yaml`

6) Activate environment: `conda activate burden`

7) Download dependency not available in conda:  `git clone https://github.com/mhulsman/ibidas3.git`

8) Enter folder: `cd ibidas3`

9) Install: `python setup.py install`

10) Leave folder: cd ..

11) Excecute script: `./burden_apply.py -h`   #to print usage info.




## Usage
```
usage: burden_apply.py [-h] [--vcf VCF] --meta META [--pca PCA] --genes GENES --variants VARIANTS [--max_maf MAX_MAF] [--max_missingness MAX_MISSINGNESS] [--diff_miss DIFF_MISS] [--exclude_qc] [--sample_removal] [--pheno PHENO] [--causative_mutation]

Process some integers.

optional arguments:
  -h, --help            show this help message and exit
  --vcf VCF             The VCF file or BCF file. If there are multiple files for the different chromosomes, one can use a wildcard (%s) which will be used to automatically fill in the correct chromosome (e.g. result.chr%s.bcf). The VCF/BCF files
                        should have an index.
  --meta META           TSV file with at least a column iid and status, containing respectively the sample id and the phenotype. Phenotype can be coded as case<=65, case>65, control.
  --pca PCA             TSV file with at least a column iid and columns pc1,pc2,pc3,pc4,pc5 and pc6.
  --genes GENES         Comma separated list of genes to analyze.
  --variants VARIANTS   Variant category to analyze. Either "LOF", or REVEL threshold between 0 and 100 (includes also LOF variants).
  --max_maf MAX_MAF     Maximum minor allele frequency of included variants (default: 0.01).
  --max_missingness MAX_MISSINGNESS
                        Maximum missingness (depth < 6) for included variants (default: 0.2).
  --diff_miss DIFF_MISS
                        Differential missingness between cases and controls, as -log10 p-value threshold (default: 20). Note that the default may be overly permissive for smaller studies (and vice-versa).
  --exclude_qc          Remove all variants denoted with the EXCLUDE_QC flag in the VCF (default: True).
  --sample_removal      Remove samples with >80% missingness in the selected variants of the analyzed gene (default: False).
  --pheno PHENO         Name of the column in the meta file that carries the phenotype information
  --causative_mutation  Indicates that there is a column 'causative_mutation_scan' in the pheno file, with samples that carry a causative mutation and need to be removed from the anaysis (default: False).
```


## VCF annotation

VCFs can be annotated with various flags/fields. Most important:

1. EXCLUDE_QC (flag), EXCLUDE_REASON (field):   variants that should be filtered, with the reason (text string)
2. CSQ: variant annotation (VEP), see below
3. GAN_NONNEURO_POPMAX, GAC_NONNEURO_POPMAX, G_NONNEURO_POPMAX, G_FILTER:  annotation copied from Gnomad VCF, describing population with maximum frequency for this variant in non-neuro samples. Used for variant filtering if available.
4. COMB_DEL (flag), COMB_NOTE, COMB_EXCLUDE_QC: variants that should be not considered because they were merged with other variants. COMB_DEL flags variants that were merge dwith other variants. Note describes the merging action, and comb_exclude_qc describes if the variant that was merged into the current variant had qc issues.
5. abhet_p: allele balance annotation (see shortread_seq_analysis/variant_qc_annot scripts). Used to identify burdens affected by somatic mutations.


## CSQ annotation: 

0. Install VEP (e.g. 94.5) + Loftee (e.g. 1.0.2) + dbNSFP (e.g. 4.1)

1. Annotate a VCF without genotype fields


    bcftools view -G {input} --threads 8 |\
    vep --offline --cache --force --stats_text --fork 8 --dir_cache {VEPCACHE} -o stdout --dir_plugins {VEPPLUGINS} --format vcf --fasta {REF_FASTA}
    	--humdiv --polyphen b --sift b --ccds  --uniprot  --hgvs  --symbol  --tsl  --vcf
    	--plugin dbNSFP,{VEPDBNSFP},genename,Ensembl_transcriptid,Ensembl_proteinid,Uniprot_acc,TSL,SIFT_score,SIFT_converted_rankscore,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,FATHMM_score,FATHMM_converted_rankscore,FATHMM_pred,VEST4_score,VEST4_rankscore,REVEL_score,REVEL_rankscore,PrimateAI_score,PrimateAI_rankscore,PrimateAI_pred,Ancestral_allele,CADD_phred,CADD_phred_hg19,clinvar_id,clinvar_clnsig  --plugin LoF,loftee_path:{VEPLOFTEE},human_ancestor_fa:{VEPANCESTOR},conservation_file:{VEPPHYLO} |\
    bcftools view -O b -o {output} --threads 8

Replace indicated variables with following paths (adapt to your situation):

    VEPPLUGINS='/path/to/a/folder/vep_software/share/ensembl-vep-94.5-0/'
    VEPDBNSFP='/path/to/a/folder/vep_resource/dbnsfp4.1/dbNSFP4.1.hg19.gz'
    VEPLOFTEE="/path/to/a/folder/loftee/loftee-1.0.2/"
    VEPANCESTOR="/path/to/a/folder/vep_resource/ancestral_seq/human_ancestor.fa.gz"
    VEPPHYLO="/path/to/a/folder/vep_resource/conservation_loftee/phylocsf_gerp.sql"

2. Transfer annotation to full VCF:

    bcftools annotate -a {input[6]} -h {MYPATH}/vep.header -c CSQ  --threads 8 {input[0]} -Ou




