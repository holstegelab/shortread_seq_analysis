## Calculates population prior for variant based on genotype likelihoods (PL field)
##
## PROBLEM
## Genotype probability vectors cannot be used 'as-is'.
## This can be best be explained from the case where no reads are available: all genotypes are then equally likely, and therefore GATK gives a probability of 0.33/0.33/0.33.
## One can observe two problems::
## - If this is a rare variant, then using a genotype probability of 0.33 here for an alternate homozygous allele is a large overestimation
## - Giving an equal probability to each genotype does not take into account Hardy-Weinberg.
##
##
## SOLUTION
## The functions in this file calculates genotype priors based on estimated population allele frequency, instead of assuming each genotype equally likely.
## For example,  a low population frequency means also a lower prior probability for homozygous alternate genotypes.
## The approach is to estimate the allele frequency and to take this through Hardy-Weinberg to calculate a genotype prior.
##
## INPUT:
## matrix of (n_samples * n_genotypes) with genotype probabilities (based on PL field)
## - supported n_genotypes: 2 or 3
##
## OUTPUT:
## 1) allele frequency estimates
## 2) genotype priors
##    - vector of n_genotypes with prior per genotype
##
## IMPLEMENTATION NOTES
##  Note that the allele frequency is based on the genotype probabilities. Here we cannot use the direct genotype probabilities from GATK due to the abovementioned problem, as this would severely overstimate/underestimate
##  allele-frequency in many cases. Allele frequency has to be estimated based on the posterior genotype probabilities. As these are not yet available (we need to calculate the population prior for that), this is somewhat
##  of a problem. To solve this chicken/egg problem, we start with an allele frequency estimate of 1/n_allele, and calculate a posterior probability based on that. Then we refine our allele frequency estimate based
##  on these posterior probabilities, and recalculate new psoterior genotype probabilities using this better estimate of the allele frequency. We iterate this until convergence.
##
## Dependencies: numpy

import numpy


def estimate_afs_simple(vectors, start=None):
    nsamples, ngeno = vectors.shape
    assert ngeno == 2 or ngeno == 3
    nallele = 2

    has_info = (vectors != (1.0 / ngeno)).all(axis=1)
    vectors = vectors[has_info, :]

    if start is None:
        afs_estimate = numpy.ones(nallele, dtype=float) / nallele
    else:
        afs_estimate = start
    genotype_priors = numpy.zeros(ngeno, dtype=float)
    max_diff = 1.0

    if ngeno == 2:
        while max_diff > 0.0000001:
            genotype_priors[0] = afs_estimate[0]
            genotype_priors[1] = afs_estimate[1]

            gvectors = genotype_priors * vectors
            gvectors = gvectors / numpy.sum(gvectors, axis=1)[:, numpy.newaxis]
            total_prob = numpy.sum(gvectors, axis=0)

            allele_counts = numpy.zeros(nallele, dtype=float)
            for j in range(nallele):
                idx = j
                allele_counts[j] += total_prob[idx]

            new_afs_estimate = allele_counts / float(len(vectors))

            max_diff = numpy.max(new_afs_estimate - afs_estimate)
            afs_estimate = new_afs_estimate

        genotype_priors[0] = afs_estimate[0]
        genotype_priors[1] = afs_estimate[1]
    else:
        while max_diff > 0.0000001:
            for j in range(nallele):
                for k in range(j, nallele):
                    idx = int(0.5 * k * (k + 1) + j)
                    genotype_priors[idx] = (
                        ((j != k) + 1.0) * afs_estimate[j] * afs_estimate[k]
                    )

            gvectors = genotype_priors * vectors
            gvectors = gvectors / numpy.sum(gvectors, axis=1)[:, numpy.newaxis]
            total_prob = numpy.sum(gvectors, axis=0)

            allele_counts = numpy.zeros(nallele, dtype=float)
            for j in range(nallele):
                for k in range(j, nallele):
                    idx = int(0.5 * k * (k + 1) + j)
                    allele_counts[j] += total_prob[idx]
                    allele_counts[k] += total_prob[idx]

            new_afs_estimate = allele_counts / (len(vectors) * 2.0)

            max_diff = numpy.max(new_afs_estimate - afs_estimate)
            afs_estimate = new_afs_estimate

        for j in range(nallele):
            for k in range(j, nallele):
                idx = int(0.5 * k * (k + 1) + j)
                genotype_priors[idx] = (
                    ((j != k) + 1.0) * afs_estimate[j] * afs_estimate[k]
                )

    return (afs_estimate, genotype_priors[:, numpy.newaxis])


def process_pl(plrow):
    # Returns allele frequency estimate and posterior probabilities.
    #
    # Parameter:
    # - plrow: array of g * n  with phred-scaled likelihood values for n samples and g genotypes.
    #
    # Returns:
    # - afs:
    # - posterior_cor: array of g * n with genotype probabilities.
    plrow = plrow - plrow.min(axis=0)
    plrow = 10.0 ** (plrow / -10.0)
    plrow = plrow / plrow.sum(axis=0)
    afs, gafs = estimate_afs_simple(plrow.T)

    pdata = gafs * plrow
    posterior_cor = pdata / pdata.sum(axis=0)

    return (afs, posterior_cor)
