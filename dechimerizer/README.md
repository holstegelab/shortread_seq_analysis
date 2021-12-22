# DeChimerizer

This script will remove soft-clip overhangs of BAM alignments that are (likely) to have been caused
by chimeric reads. 

It has been developed based on BWA alignments, created with the -Y flag.


For a detailed description, see the supplement to the paper:

"Exome sequencing identifies rare damaging variants in ATP8B4 and ABCA1 as novel risk factors for Alzheimer’s Disease".


Usage:

samtools view -h -@ 8 {input} | python {MYPATH}/bam_clean.py | samtools view -@ 8 -b -o {output}

Tip: for faster performance, use pypy. 


Dependencies:
- swalign

# Bam Stats

Generates BAM statistics for the following fields:

unmapped_ratio:      fraction of alignments unmapped
mqual20_ratio:       fraction of alignments with mapping quality >= 20
secondary_ratio:     fraction of alignments that are secondary
supplementary_ratio: fraction of alignments that are supplementary
sup_diffchrom_ratio: fraction of supplementary alignments that have a primary alignment on a different chromosome
duplicate_ratio:     fraction of duplicate alignments
duplicate_supplement_ratio: fraction of duplicate supplementary alignemnts
soft_clipped_bp_ratio:      fraction of bases that are soft-clipped
aligned_bp_ratio:           fraction of bases that are aligned
inserted_bp_ratio:          fraction of bases that are inserted
poly_a:                     fraction of alignments with a poly-a signature
illumina_adapter:           fraction of alignnment with a illimuna adapter
pcr_adapter_1:              fraction of alignments with a Illumina PCR adapter 1
pcr_adapter_2:              fraction of alignments with a Illumina PCR adapter 2
nextera:                    fraction of alignments with a Nextera adapter


Usage: 

samtools view -s 0.05 -h {input[0]} --threads 8 | {MYPATH}/bam_stats.py stats > {output}

This generates statistics with a sampling ratio = 0.05, this can be adapted as needed. 


Dependencies:
- numpy
