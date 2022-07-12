#!/usr/bin/env python2
import gzip
import csv
import numpy
import sys

from scipy.stats import beta
from scipy.stats import chi2
from multiprocessing import Pool
import rescale_af
import cyvcf2


def get_csv(filename, fieldnames=None):
    # Read tab delimited csv file.
    # if no fieldnames are given, assumes that first line of file contains fieldnames

    with open(filename) as meta_file:
        meta_reader = csv.reader(meta_file, delimiter="\t")
        if fieldnames is None:
            fieldnames = next(meta_reader)
        meta_data = [row for row in meta_reader]

    return {fieldname: col for fieldname, col in zip(fieldnames, list(zip(*meta_data)))}


def hardy(obs_hom1, obs_hets, obs_hom2):
    # Calculate significance of deviation from Hardy-Weinberg and inbreeding coefficient for n samples with diploid genome.
    # Input:
    # - obs_hom1, obs_hets, obs_hom2: array of size n with for each sample probability of respectively reference, heterozygous and homozygous call
    xx = obs_hom1
    xy = obs_hets
    yy = obs_hom2

    samples = xx + xy + yy
    if samples == 0:
        return (0.0, 0.5)

    observed_heterozygosity = float(xy) / samples

    alleleX_count = (2 * xx) + xy
    alleleX_ratio = float(alleleX_count) / (2 * samples)
    alleleY_count = xy + (2 * yy)
    alleleY_ratio = float(alleleY_count) / (2 * samples)

    expected_xx_n = alleleX_ratio * alleleX_ratio * samples
    expected_xy_n = 2.0 * alleleX_ratio * alleleY_ratio * samples
    expected_yy_n = alleleY_ratio * alleleY_ratio * samples
    expected_total = expected_xx_n + expected_xy_n + expected_yy_n

    expected_xx_ratio = float(expected_xx_n) / expected_total
    expected_xy_ratio = float(expected_xy_n) / expected_total
    expected_yy_ratio = float(expected_yy_n) / expected_total

    if xx == 0.0 and expected_xx_n == 0.0:
        chi_square_xx = 0.0
    else:
        chi_square_xx = ((xx - expected_xx_n) ** 2.0) / expected_xx_n

    if xy == 0.0 and expected_xy_n == 0.0:
        chi_square_xy = 0.0
    else:
        chi_square_xy = ((xy - expected_xy_n) ** 2.0) / expected_xy_n

    if yy == 0.0 and expected_yy_n == 0.0:
        chi_square_yy = 0.0
    else:
        chi_square_yy = ((yy - expected_yy_n) ** 2.0) / expected_yy_n

    chi_square_total = chi_square_xx + chi_square_xy + chi_square_yy

    if expected_xy_n == 0.0:
        icoef = -numpy.inf if xy > 0.0 else numpy.nan
    else:
        icoef = 1.0 - (xy / expected_xy_n)
    hardypval = chi2.sf(chi_square_total, 1)

    if hardypval == 0.0:
        hardypval = 10e-300

    return (icoef, hardypval)


def calc_stats(alleles, probs, fil, controls, hardy_fil, ad):
    # Calculate variant statistics
    # Inputs:
    # - alleles: number of alleles (ploidy)
    # - probs: array with genotype probabilities of size g * n  (g = number of genotypes, n = number of samples)
    # - fil: boolean array of size n, indicating which samples to use. If None, use all samples.
    # - control_fil: boolean array of size n, indicating which samples are control samples.
    # - hardy_fil: boolean array of size n, indicating which samples to use to estimate Hardy-Weinberg statistics (e.g. only males for chrY).
    # - ad: boolean array of size a * n (a = number of alleles), containing read depth for each allele.
    #
    # Outputs:
    #     - filfrac: fraction of calls that pass the filter
    #     - ac: allele count of the minor allele
    #     - gt: string with count of reference, heterozygous and homozygous genotype
    #     - dpref: string with sum of (number of reference reads times probablity of respectively reference, heterozygous and homozygous genotype).
    #     - dpalt: string with sum of (number of alternate reads times probablity of respectively reference, heterozygous and homozygous genotype).
    #     - abhet: allele balance heterozygous calls (0.5 = normal, 1.0 = reference bias)
    #     - abhom: allele balance homozygous calls (1.0 = normal, 0.0 = reference bias)
    #     - dhet_frac:  reads on het calls vs. all calls.  (ratio)
    #     - abdev_rate: not used
    #     - incoef: inbreeding coefficient
    #     - phred-scaled hardy-weinberg p-value
    #     - control_incoef: inbreeding coefficient of controls
    #     - phred-scaled hardy-weingberg p-value of controls
    #     - althet:
    #     - total_althet:

    if fil is None:
        control_fil = controls
        if hardy_fil is None:
            hardy_fil = slice(None, None, None)

        ad_fil = ad
        probs_fil = probs
        filfrac = 1.0
    else:
        if hardy_fil is None:
            hardy_fil = slice(None, None, None)
        else:
            hardy_fil = hardy_fil[fil]
        control_fil = controls[fil]
        ad_fil = ad[:, fil]
        probs_fil = probs[:, fil]
        filfrac = fil.mean()

    # optim hardy: depth 6?
    ac_ref = probs_fil[0, :].sum()
    # post ref
    dref_ref = (probs_fil[0, :] * ad_fil[0, :]).sum()
    dref_alt = (probs_fil[0, :] * ad_fil[1, :]).sum()

    if alleles == 1:
        ac_het = 0
        ac_hom = probs_fil[1, :].sum()

        hardypval = 0.5
        incoef = 0
        control_hardypval = 0.5
        control_incoef = 0
        abhet = 0.5

        vafdev = 0.5
        dhet_ref = 0
        dhet_alt = 0

        dhom_ref = (probs_fil[1, :] * ad_fil[0, :]).sum()
        dhom_alt = (probs_fil[1, :] * ad_fil[1, :]).sum()
        dxhom_ref = dhom_ref + dref_alt
        dxhom_alt = dhom_alt + dref_ref
        if (dxhom_ref + dxhom_alt) == 0.0:
            abhom = 1.0
        else:
            abhom = dxhom_alt / (dxhom_ref + dxhom_alt)

        homprob = probs_fil[1, :].sum()
        if homprob > 0:
            fullref = ad_fil[1, :] == 0.0
            abdev_rate = (numpy.cast[int](fullref) * probs_fil[1, :]).sum() / homprob
        else:
            abdev_rate = 0.0
        ac = min(ac_ref, ac_hom)

    else:
        ac_het = probs_fil[1, :].sum()
        ac_hom = probs_fil[2, :].sum()

        incoef, hardypval = hardy(
            probs_fil[0, hardy_fil].sum(),
            probs_fil[1, hardy_fil].sum(),
            probs_fil[2, hardy_fil].sum(),
        )
        control_incoef, control_hardypval = hardy(
            probs_fil[0, control_fil].sum(),
            probs_fil[1, control_fil].sum(),
            probs_fil[2, control_fil].sum(),
        )

        # post het/hom
        dhet_ref = (probs_fil[1, :] * ad_fil[0, :]).sum()
        dhet_alt = (probs_fil[1, :] * ad_fil[1, :]).sum()
        if (dhet_ref + dhet_alt) == 0.0:
            abhet = 0.5
        else:
            abhet = dhet_ref / (dhet_ref + dhet_alt)

        dhom_ref = (probs_fil[2, :] * ad_fil[0, :]).sum()
        dhom_alt = (probs_fil[2, :] * ad_fil[1, :]).sum()
        dxhom_ref = dhom_ref + dref_alt
        dxhom_alt = dhom_alt + dref_ref
        if (dxhom_ref + dxhom_alt) == 0.0:
            abhom = 1.0
        else:
            abhom = dxhom_alt / (dxhom_ref + dxhom_alt)

        hetprob = probs_fil[1, :].sum()
        if hetprob > 0:
            fullref = ad_fil[1, :] == 0.0
            abdev_rate = (numpy.cast[int](fullref) * probs_fil[1, :]).sum() / hetprob
        else:
            abdev_rate = 0.0
        ac = min(ac_het + 2.0 * ac_ref, ac_het + 2.0 * ac_hom)

    if ac_het == 0:
        dhet_frac = 1.0
    else:
        d_het = (dhet_ref + dhet_alt) / ac_het
        d_all = (dref_ref + dref_alt + dhet_ref + dhet_alt + dhom_ref + dhom_alt) / (
            ac_ref + ac_het + ac_hom
        )

        dhet_frac = d_het / d_all

    gt = "%.2f,%.2f,%.2f" % (ac_ref, ac_het, ac_hom)
    dpref = "%.2f,%.2f,%.2f" % (dref_ref, dhet_ref, dhom_ref)
    dpalt = "%.2f,%.2f,%.2f" % (dref_alt, dhet_alt, dhom_alt)

    if alleles == 1:
        if ac_hom < ac_ref:
            althet = dhom_alt / float(ac_hom)
            total_althet = dhom_alt
        else:
            althet = dref_ref / float(ac_ref)
            total_althet = dref_alt
    else:
        if ac_hom <= ac_ref:
            althet = (dhet_alt + dhom_alt) / (float(ac_het) + 2.0 * float(ac_hom))
            total_althet = dhet_alt + dhom_alt

        else:
            althet = (dhet_ref + dref_ref) / (float(ac_het) + 2.0 * float(ac_ref))
            total_althet = dhet_ref + dref_ref

    return (
        filfrac,
        ac,
        gt,
        dpref,
        dpalt,
        abhet,
        abhom,
        dhet_frac,
        abdev_rate,
        incoef,
        -numpy.log10(hardypval) * 10.0,
        control_incoef,
        -numpy.log10(control_hardypval) * 10.0,
        althet,
        total_althet,
    )


def process_row(args):
    # utility function to process a variant row
    chrom, pos, ref, alt, pl01, plrow, adrow, control_fil, hardy_fil = args

    if chrom == "Y" or chrom == "chrY":
        alleles = 1
        pltypes = 2

    else:
        alleles = 2
        pltypes = 3

    # convert to prob
    plrow = plrow - plrow.min(axis=0)
    gq = numpy.sort(plrow, axis=0)[1, :]
    plrow = 10.0 ** (plrow / -10.0)
    plrow = plrow / plrow.sum(axis=0)

    # to posterior prob
    afs, gafs = rescale_af.estimate_afs_simple(plrow.T)  # determine population prior
    pdata = gafs * plrow
    posterior_cor = pdata / pdata.sum(axis=0)

    ### VAF filter ###
    depth = adrow.sum(axis=0)
    fil_d0 = depth >= 0
    fil_d6 = depth >= 6
    fil_d10 = depth >= 10
    fil_gq20 = gq >= 20

    det_stats_raw_d0 = calc_stats(alleles, pl01, None, control_fil, hardy_fil, adrow)
    det_stats_raw_d6 = calc_stats(alleles, pl01, fil_d6, control_fil, hardy_fil, adrow)
    det_stats_raw_d10 = calc_stats(
        alleles, pl01, fil_d10, control_fil, hardy_fil, adrow
    )
    det_stats_raw_gq20 = calc_stats(
        alleles, pl01, fil_gq20, control_fil, hardy_fil, adrow
    )
    det_stats_p = calc_stats(
        alleles, posterior_cor, None, control_fil, hardy_fil, adrow
    )
    det_stats_p_d10 = calc_stats(
        alleles, posterior_cor, fil_d10, control_fil, hardy_fil, adrow
    )

    result = (
        (chrom, pos, ref, alt)
        + det_stats_raw_d0
        + det_stats_raw_d6
        + det_stats_raw_d10
        + det_stats_raw_gq20
        + det_stats_p
        + det_stats_p_d10
    )
    return result


def parse_matrix(outfile, meta_file, filename, range_description=None, nproc=1):
    # Process variant data
    # Parameters:
    # - outfile: file to output vairant info to
    # - meta_file: file with sample annotation information.Tab-separated. Needed columns: 'iid' (sample iid), 'status' (0 = control), 'gender' (1:male, 2:female).
    # - filename: input VCF/BCF file
    # - range_description: string with chrom:start-stop format.
    # - nproc: number of simultaneous processes

    import var_qc_annot as vqa

    pool = Pool(processes=nproc)

    # get control sample ids
    meta = get_csv(meta_file)
    control_samples = set(
        [iid for iid, status in zip(meta["iid"], meta["status"]) if status == "0"]
    )
    female_samples = set(
        [iid for iid, gender in zip(meta["iid"], meta["gender"]) if gender == "2"]
    )
    male_samples = set(
        [iid for iid, gender in zip(meta["iid"], meta["gender"]) if gender == "1"]
    )

    # open out file and write header
    outf = open(outfile, "w")
    fieldnames = ["chrom", "pos", "ref", "alt"]
    res = [
        "kept",
        "ac",
        "gt",
        "dpref",
        "dpalt",
        "abhet",
        "abhom",
        "fdhet",
        "abdev",
        "incoef",
        "hw",
        "control_incoef",
        "control_hw",
        "althet",
        "total_althet",
    ]
    filters = ["r", "r_d6", "r_d10", "r_gq20", "p", "p_d10"]
    for e in filters:
        nres = [x + "_" + e for x in res]
        fieldnames.extend(nres)

    outf.write("#" + "\t".join(fieldnames) + "\n")

    vcf_reader = cyvcf2.VCF(filename)

    samplenames = vcf_reader.samples
    control_fil = numpy.array([a in control_samples for a in samplenames], dtype=bool)
    female_control_fil = numpy.array(
        [a in control_samples and a in female_samples for a in samplenames], dtype=bool
    )
    female_fil = numpy.array([a in female_samples for a in samplenames], dtype=bool)
    male_control_fil = numpy.array(
        [a in control_samples and a in male_samples for a in samplenames], dtype=bool
    )
    male_fil = numpy.array([a in male_samples for a in samplenames], dtype=bool)

    rownr = 0
    rows = []
    # walk through vcf
    vcf_iter = iter(vcf_reader(range_description))
    try:
        while True:
            v = next(vcf_iter)
            if len(rows) == 0:
                if v.CHROM == "X":
                    control_fil = female_control_fil
                    hardy_fil = female_fil
                elif v.CHROM == "Y":
                    control_fil = male_control_fil
                    hardy_fil = male_fil

                else:
                    hardy_fil = None

            pls = v.format("PL")
            ad = v.format("AD")
            gt = v.gt_types

            if v.CHROM == "Y":
                colnr = 2
            else:
                colnr = 3

            if pls.shape[1] < colnr or ad.shape[1] < 2:
                print((pls.shape))
                xpls = numpy.ones((pls.shape[0], colnr)) * 100
                xpls[:, : pls.shape[1]] = pls
                pls = xpls

                xad = numpy.zeros((ad.shape[0], 2))
                xad[:, 0] = ad[:, 0]
                ad = xad

            pls[pls <= -2147483647] = 0
            ad[ad <= -2147483647] = 0
            gt_new = numpy.zeros((pls.shape[0], colnr))
            gt_new[:, 0] = gt == 0

            if pls.shape[1] > 2:
                gt_new[:, 1] = gt == 1
                gt_new[:, 2] = gt == 3
            else:
                gt_new[:, 1] = gt == 3

            rows.append((v.CHROM, v.POS, v.REF, v.ALT[0], gt_new.T, pls.T, ad.T))

            if len(rows) >= (nproc * 150):
                break

    except StopIteration as e:
        pass
    rownr += len(rows)
    sys.stdout.write("Processing variants:\n")
    fo = (
        "%s\t%d\t%s\t%s"
        + "\t%f\t%f\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f" * 6
        + "\n"
    )
    while rows:
        if nproc == 1:
            for row in rows:
                e = process_row(row + (control_fil, hardy_fil))
                if not e is None:
                    outf.write(fo % e)
        else:
            res = pool.map_async(
                vqa.process_row, [e + (control_fil, hardy_fil) for e in rows], 50
            )

        rows = []

        try:
            while True:
                v = next(vcf_iter)
                pls = v.format("PL")
                ad = v.format("AD")
                gt = v.gt_types

                pls[pls <= -2147483647] = 0
                ad[ad <= -2147483647] = 0

                if v.CHROM == "Y":
                    colnr = 2
                else:
                    colnr = 3

                if pls.shape[1] < colnr or ad.shape[1] < 2:
                    print((pls.shape))
                    xpls = numpy.ones((pls.shape[0], colnr)) * 100
                    xpls[:, : pls.shape[1]] = pls
                    pls = xpls

                    xad = numpy.zeros((ad.shape[0], 2))
                    xad[:, 0] = ad[:, 0]
                    ad = xad

                gt_new = numpy.zeros((pls.shape[0], max(pls.shape[1], 2)))
                gt_new[:, 0] = gt == 0

                if pls.shape[1] > 2:
                    gt_new[:, 1] = gt == 1
                    gt_new[:, 2] = gt == 3
                else:
                    gt_new[:, 1] = gt == 3

                rows.append((v.CHROM, v.POS, v.REF, v.ALT[0], gt_new.T, pls.T, ad.T))

                if len(rows) >= (nproc * 150):
                    break
        except StopIteration as e:
            pass

        rownr += len(rows)
        sys.stdout.write("\rCurrent variant: %d" % rownr)
        sys.stdout.flush()

        if nproc > 1:
            for e in res.get():
                if not e is None:
                    outf.write(fo % e)

    pool.close()


if __name__ == "__main__":
    print((sys.argv))
    if len(sys.argv) == 5:
        parse_matrix(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        parse_matrix(sys.argv[1], sys.argv[2], sys.argv[3])
