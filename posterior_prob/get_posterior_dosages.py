import sys
import numpy
import cyvcf2
import rescale_af

# usage:  python get_posterior_dosage.py  input.bcf/vcf  chrom:from-start    output.bcf/vcf  [threshold]

# This script will add three fields to the VCF:
# DS: for each call the posterior dosage
# GP: for each call/genotype pair the posterior probability
# PGQ: for each call the posterior genotype quality
#
#
# threshold parameter: only consider calls with depth >= [threshold], consider all other calls missing.
# if threshold parameter is given, only sets GP values currently.
#
#
# For a description of the calculation of the posterior probabilities, see rescale_af.py
# Dependencies:  cyvcf2, numpy
# Limitations: only for haploid and diploid genotype calls. Other calls are skipped.

vcf_reader = cyvcf2.VCF(sys.argv[1])
vcf_reader.add_format_to_header(
    {"ID": "DS", "Number": "1", "Type": "Float", "Description": "Dosage (posterior)"}
)
vcf_reader.add_format_to_header(
    {
        "ID": "GP",
        "Number": "G",
        "Type": "Float",
        "Description": "Probability (posterior)",
    }
)
vcf_reader.add_format_to_header(
    {"ID": "PGQ", "Number": "1", "Type": "Float", "Description": "GQ (posterior)"}
)
if len(sys.argv) > 2:
    out = cyvcf2.Writer(sys.argv[2], vcf_reader)
else:
    out = cyvcf2.Writer(sys.stdout, vcf_reader)

counter = 0
# walk through vcf
if len(sys.argv) > 3:
    threshold = int(sys.argv[3])
else:
    threshold = -1

maxlen = 1
for v in vcf_reader:
    counter += 1

    pls = v.format("PL")
    rad = v.format("AD")
    dp = v.format("DP")

    # identify missing calls
    pl_missing = pls[:, 0] <= -2147483647

    genotypes = v.genotypes
    for pos in numpy.where(pl_missing)[0]:
        genotypes[pos] = [-1, -1, False]
    v.genotypes = genotypes

    pls[pls <= -2147483647] = 0
    rad[rad <= -2147483647] = 0
    depth = rad.sum(axis=1)
    depth[pl_missing] = 0

    ngeno = pls.shape[1]

    if pls.T.shape[0] == 3 or pls.T.shape[0] == 2:
        if threshold == -1:
            afs, posterior_pls = rescale_af.process_pl(pls.T)
        else:
            fil = depth >= threshold
            afs, tposterior_pls = rescale_af.process_pl(pls[fil, :].T)
            posterior_pls = numpy.zeros(
                (tposterior_pls.shape[0], len(depth)), dtype=float
            )
            posterior_pls[:, fil] = tposterior_pls

        if posterior_pls.shape[0] == 3:
            dosage = posterior_pls[1, :] + 2.0 * posterior_pls[2, :]
        else:
            dosage = posterior_pls[1, :]

    else:
        continue

    # ignore runtime warning, it does the right thing

    # determine phred-scaled posterior probabilities
    phred_posterior_pls = -10.0 * numpy.log10(posterior_pls)
    phred_posterior_pls_min = phred_posterior_pls - phred_posterior_pls.min(axis=0)
    phred_posterior_pls_min[phred_posterior_pls_min > 10000.0] = 10000.0

    # determine genotype call that is most likely
    calls = numpy.argmin(phred_posterior_pls_min, axis=0)
    phred_posterior_pls_min = numpy.cast[int](phred_posterior_pls_min)

    # calculate genotype quality
    gq = numpy.sort(phred_posterior_pls_min, axis=0)[1, :]
    gq[gq > 99] = 99

    if (
        threshold >= 0
    ):  # for genotypes with depth < threshold, set genotype quality to 0 (missing genotype).
        gq[~fil] = 0

    v.set_format("PGQ", gq)

    # restore missingness arrays
    rad[pl_missing, 0] = -2147483648
    rad[pl_missing, 1] = -2147483647
    dp[pl_missing] = -2147483648
    if threshold == -1:
        v.set_format("AD", rad)
        v.set_format("DP", dp)
        v.set_format("DS", dosage)
        v.set_format("GP", posterior_pls.T.copy(order="C"))
    else:
        v.set_format("GP", posterior_pls.T.copy(order="C"))

    out.write_record(v)

out.close()
vcf_reader.close()
