#!/usr/bin/env python
import errno    
import os.path
import os
import sys
import argparse

import cyvcf2
import numpy
import fisher

import rescale_af
from ibidas import *
from ibidas.utils import util
import score_models
HOME=os.path.expanduser('~')

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def read_good_transcripts(filename=os.path.join(HOME, 'ades_v29_v19_whitelist_notsl_all_mapped_transcripts.tsv')):
    sys.stderr.write("Loading transcripts...\n"); sys.stderr.flush()
    e = Read(filename,verbose=False).Detect().Copy()
    e = e.To(_.transcript_id, _.gene_id, Do=_.SplitOnPattern('\.')[0]).Get(_.transcript_id, _.gene_id).Copy()
    return e


def get_gene_data(genes_file):
    return Read(genes_file,delimiter='\t', verbose=False).Detect().Copy()

def get_meta_data(filename=os.path.join(HOME,'ades/phenotypes.tsv'), pca_filename = os.path.join(HOME,'ades/pca.tsv')):
    sys.stderr.write("Loading meta-data...\n"); sys.stderr.flush()
    d = Read(filename,delimiter='\t', verbose=False).Detect().Copy()
    
    if not 'wes_plate' in d.Names: 
        d = d.Get('*', Rep('')/'wes_plate').Level(1)

    if pca_filename is None:
        d = d/{'iid':'sample'}
        d = d.Get('*', Rep(None)/'pca_europe').Copy()

    else:        
        sys.stderr.write("Loading PCA components...\n"); sys.stderr.flush()
        pca = Read(pca_filename, delimiter='\t', verbose=False).Detect().Copy()
        d = (d |Match(_.iid, _.individ)| pca).Copy()
        d = d.Get(_.iid/'sample', _.Without(_.iid, _.pc1, _.pc2, _.pc3, _.pc4, _.pc5, _.pc6), \
                HArray(_.pc1, _.pc2, _.pc3, _.pc4, _.pc5, _.pc6).values/'pca_europe').Copy()

        d = prepare_pop_feats(d)

    return d



def prepare_pop_feats(sample_stats):#{{{
    pca_feats  = sample_stats.Get(_.pca_europe)()

    pca_feats = pca_feats - numpy.min(pca_feats, axis=0)[numpy.newaxis,:]
    pca_feats = pca_feats / numpy.max(pca_feats, axis=0)
    pca_feats = pca_feats - numpy.mean(pca_feats, axis=0)[numpy.newaxis,:]

    #bit of a workaround to add to dataset, needs better support
    all_vectors = sample_stats.Without(_.pca_europe)()
    res = tuple(list(all_vectors) + [pca_feats])
    names = tuple(list(sample_stats.Without(_.pca_europe).Names) + ['pca_europe'])
    res = Rep(res)/names

    return res.Copy()
#}}}



def read_ranges(vcf_filename, meta_info, good_transcripts, selected_samples,ranges, **kwargs):
    geneid_ann = kwargs.get('geneid_ann',None)
    selkwargs = kwargs

    keep_variants = kwargs.get('keep_variants',{})

    
    #DETERMINE SELECTION THRESHOLDS
    lof = selkwargs.get('lof',None)
    revel = selkwargs.get('revel',None)
    max_maf = selkwargs.get('max_maf',None)
    max_missingness = selkwargs.get('max_missingness',None)
    
    exclude_qc = selkwargs.get('exclude_qc',None)
    comb_del = selkwargs.get('comb_del', None)
    vbd_tech = selkwargs.get('vbd_tech', None)

    done = False


    processed_variants = set()
    var_identifiers = []


    types = {'chrom':"[variants:*]:string", 
             'pos':"[variants:*]:int64",
             "ref":"[variants:*]:string",
             "alt":"[variants:*]:string",
             "snp_id":"[variants:*]:string",
             "g_filter":"[variants:*]:string?",
             "gac_nonneuro_popmax":"[variants:*]:real64?",
             "gan_nonneuro_popmax":"[variants:*]:real64?",
             "g_nonneuro_popmax":"[variants:*]:string",
             "comb_del":"[variants:*]:int64?",
             "comb_exclude_qc":"[variants:*]:int64?",
             "ma_reflow_count":"[variants:*]:int64?",
             "vbd_tech":"[variants:*]:real64?",
             "abhet_p":"[variants:*]:real64?",
             "pls":"[variants:*]:[genotypes:*]:[samples:*]:real64",
             "posterior_pls":"[variants:*]:[genotypes:*]:[samples:*]:real64",
             "depth":"[variants:*]:[samples:*]:int64",
             "maf":"[variants:*]:real64",
             "ac_prob":"[variants:*]:real64",
             "sample":"[samples:*]:string",
             "csq":"[variants:*]:string",
             "exclude_qc":"[variants:*]:int64?",
             'ad':"[variants:*]:[samples:*]:[alleles:*]:int64",
             "exclude_reason":"[variants:*]:string?",
             "comb_note":"[variants:*]:string?"
             }


    warnings = {
        'csq': 'No variant prioritization information from VEP. ALL VARIANTS WILL BE USED!.',
        'exclude_reason': 'Cannot provide exclusion reasons per variant.',
        'exclude_qc': 'NO VARIANT QC. Will include all variants.',
        'gan_nonneuro_popmax':'Cannot filter on GNOMAD max. population frequency.',
        'gac_nonneuro_popmax':'Cannot filter on GNOMAD max. population frequency.',
        'g_nonneuro_popmax':'Cannot filter on GNOMAD max. population frequency.',
        'g_filter': 'Cannot filter on GNOMAD max. population frequency.',
        'ma_reflow_count': 'Filtering of multi-allelic variants with >= 2 alternate alleses with AF > AF-REF not available. ',
        'comb_del': 'No exclusion of merged variants (variants were probably not merged).',
        'comb_note': 'No explanation of which variants were merged.',
        'comb_exclude_qc': 'No sharing of exclusion information of merged variants.',
        'abhet_p': 'No allele balance information to identify somatic variants.',
        'vbd_tech': 'Variant Batch Detecter filtering not available.'}

    var_infos = {'gac_nonneuro_popmax':[], 'gan_nonneuro_popmax':[], 'g_nonneuro_popmax':[], 'ma_reflow_count':[], 'comb_exclude_qc':[], \
        'exclude_qc':[], 'abhet_p':[],'vbd_tech':[], 'g_filter':[], 'comb_del':[], 'csq':[],'exclude_reason':[], 'comb_note':[]}

    var_posterior_pls = []
    var_afs = []
    var_dps = []
    var_acprob = []


    if '%' in vcf_filename:
        vcf = cyvcf2.VCF(vcf_filename % ranges[0][0])
    else:
        vcf = cyvcf2.VCF(vcf_filename)



    available_fields = set()
    for field in vcf.header_iter():
        if field['HeaderType'] == 'INFO':
            available_fields.add(util.valid_name(field['ID']))

    for f in var_infos:
        if not f in available_fields:
            print(('WARNING: Field %s not available. %s' % (f, warnings.get(f,''))))


    samples = vcf.samples
    set_sel_samples = set(selected_samples)
    sample_fil = numpy.array([i in set_sel_samples for i in samples],dtype=bool)

    sys.stderr.write("Retrieving variants"); sys.stderr.flush()
    for var_range in ranges:
        chrom, start, end = var_range
        if '%' in vcf_filename:
            vcf = cyvcf2.VCF(vcf_filename % chrom)
        else:
            vcf = cyvcf2.VCF(vcf_filename)
        assert tuple(samples) == tuple(vcf.samples), 'VCFs with different samples or sample order for different chroms are not allowed'
        

        i = vcf('%s:%d-%d' % (chrom, start, end)).__iter__()
        while not done:

            try:
                var = next(i)
            except StopIteration:
                done = True
                continue

            vid = (var.CHROM, var.POS, var.REF, var.ALT[0],var.ID if not var.ID is None else '')

            #sys.stderr.write('%s:%d:%s>%s (%s)\n' % vid)
            sys.stderr.write("."); sys.stderr.flush()

            vname = '%s:%d:%s>%s' % (str(var.CHROM), var.POS, var.REF, var.ALT[0]) 

            keep =  var.ID in keep_variants or vname in keep_variants

            if vid in processed_variants:#a check for overlapping variants
                continue

            vinfo = {}
            z = iter(var.INFO)
            while True:
                try:
                    k,v = next(z)
                    name = util.valid_name(k)
                    vinfo[name] = v
                except UnicodeDecodeError:
                    pass
                except StopIteration:
                    break

            if exclude_qc and 'exclude_qc' in vinfo and not keep:
                continue

            if comb_del and 'comb_del' in vinfo and not keep:
                continue

            if not vbd_tech is None and 'vbd_tech' in vinfo:
                if float(vinfo['vbd_tech']) > vbd_tech and not keep:
                    continue

            if (not lof is None or not revel is None) and 'csq' in vinfo:
                vrevel,vis_lof = csq_revel(vinfo['csq'])
                if not lof is None:
                    if not vis_lof and not keep:
                        continue
                
                if not vis_lof and not revel is None:
                    if vrevel < revel and not keep:
                        continue
            
            #get gatk variant probability

            href = var.gt_phred_ll_homref
            hhet = var.gt_phred_ll_het
            halt = var.gt_phred_ll_homalt
               
            pls = numpy.vstack((href, hhet, halt))
            pls = pls - pls.min(axis=0)
            pls = 10**(pls / -10.0)
            pls = pls / pls.sum(axis=0)
            
            afs, gafs = rescale_af.estimate_afs_simple(pls.T)
            maf = afs[1]
            pdata = gafs * pls
            posterior_cor = pdata / pdata.sum(axis=0)
            
            posterior_cor = posterior_cor[:,sample_fil]
            pls = pls[:,sample_fil]

          
            dp = var.format('DP').ravel()[sample_fil]
            geno_missing = var.gt_types[sample_fil] == 2
            dp[dp < 0] = 0
            dp[geno_missing] = -1
            
            if not max_missingness is None and (dp < 6).mean() > max_missingness and not keep:
                continue
           
            if not max_maf is None and min(maf, 1.0 - maf) > max_maf and not keep:
                continue
            
            ac_prob = posterior_cor[1,:].sum() + 2.0 * posterior_cor[2,:].sum()

            processed_variants.add(vid)
            var_identifiers.append(vid)

            for key in list(var_infos.keys()):
                if key in vinfo:
                    if 'real64' in types[key]:
                        var_infos[key].append(float(vinfo[key]))
                    elif 'int64' in types[key]:
                        var_infos[key].append(int(vinfo[key]))
                    else:
                        var_infos[key].append(vinfo[key])
                else:
                    if 'real64' in types[key] or 'int64' in types[key]:
                        var_infos[key].append(Missing)
                    else:
                        var_infos[key].append('')

            var_posterior_pls.append(posterior_cor)
            var_acprob.append(ac_prob)
            var_dps.append(dp)

            var_afs.append(maf)
            
    fields = ['chrom','pos','ref','alt','snp_id']
    cols = [[e[i] for e in var_identifiers] for i in range(len(fields))]
    info_fields = list(var_infos.keys())
    fields.extend(info_fields)
    for ifield in info_fields:
        cols.append(var_infos[ifield])

    fields.extend(['posterior_pls','depth', 'maf', 'ac_prob'])
    cols.extend([var_posterior_pls, var_dps, var_afs, var_acprob])

    fields.append('sample')
    cols.append(numpy.array(samples)[sample_fil])

    if len(var_afs) == 0:
        return None

    dtype = "(" + ", ".join(['%s=%s' % (field, types[field]) for field in fields]) + ")"
    data = Rep(tuple(cols),dtype=dtype)
   
    data = data |Match(_.sample, jointype='left', merge_same=True)| meta_info
    
    data = data.Sort(_.pos, _.alt)
    data = data % ('variants','genotypes', 'samples')
   

    for name in ['comb_del','comb_exclude_qc','exclude_qc']:
        if name in data.Names:
            data = data.To(name, Do=(_ != '').ReplaceMissing(False))
    for name in ['vbd_tech','ma_reflow_count']:
        if name in data.Names:
            data = data.To(name, Do=_.Cast(str).Cast('real64?').ReplaceMissing(-1.0))

    if not geneid_ann is None:
        xgood_transcripts = set(good_transcripts[_.gene_id == geneid_ann].transcript_id())
    else:
        xgood_transcripts = set(good_transcripts.transcript_id())

    
    data = data.To(_.csq, Do=_.Each(lambda x: csq_parse(x, xgood_transcripts),dtype='(f0=(revel_score=real64,mutationtaster_pred=string,cadd38_16=real64,cadd37_16=real64,clinvar_clnsig=string,clinvar_id=string),f1=[transcripts:~]<(exon=string,cds_relpos=string,aa_relpos=string,amino_acids=string,feature_id=string,gene_name=string,gene_id=string,lof=string,lof_filter=string,lof_info=string,sift=string,polyphen=string,annotation=string,annotation_impact=string,bad_transcript=bool,vest4_score=real64,fathmm_score=real64))').Fields().Get(_.f0.Fields(), _.f1.Fields()))


    if not geneid_ann is None:
        data = data[_.gene_id == geneid_ann]
    data = data.To(_.revel_score,  Do=_.ReplaceMissing(-1.0))

    
    data = data.Get("*", _[_.bad_transcript == False][~IsMissing(_.feature_id)].Get(
                                      Any((_.lof == 'HC') & (_.annotation_impact == 'HIGH')).Cast(int)/'is_lof', 
                                      _.Get(_.revel_score, (Any(_.annotation_impact == 'MODERATE'))).Each(lambda a,b: a if b else 0.0,dtype=float)/'revel_modified_score'))


    data = data.Get(_.Without(_.revel_modified_score, _.g_filter), _.Get(_.revel_modified_score, _.is_lof).Each(lambda x,y: max(x,float(y)),dtype='real64')/'lof_revel',
                         _.Get(_.g_filter, _.gac_nonneuro_popmax.Cast('real64?'), _.gan_nonneuro_popmax.Cast('real64?')).Each(lambda g_filter,gac_pop,gan_pop: Missing if (g_filter != '' or gac_pop == -1 or gac_pop is Missing or gan_pop is Missing or float(gan_pop) == 0.0) else float(gac_pop) / float(gan_pop),dtype='real64?') /'gaf_nonneuro_popmax')

    if not 'csq' in available_fields:
        data = data.To(_.lof_revel, Do=_.Each(lambda x: 1.0,dtype='real64'))
        data = data.To(_.is_lof, Do=_.Each(lambda x: 1,dtype='int64'))
        data = data.Get(_.Without(_.bad_transcript), _.bad_transcript.Array().Each(lambda x: [False], dtype='[alt_transcript:*]<bool').Elems()).Copy()

    sys.stderr.write('\n')
    sys.stderr.flush()

    return data.Copy()



def csq_revel(cqs):
    transcript_annotation = [e for e in cqs.split(',') if not e == '']
    is_lof = False
    qrevel_score = 0.0
    for ta in transcript_annotation:
        allele, annotation, annotation_impact, gene_name, gene_id, feature_type, feature_id, transcript_biotype, exon, intron, hgvsc, hgvsp, cdna_relpos, cds_relpos, aa_relpos, amino_acids, codos, existing, distance, strand, flags, symbol_source, hgnc_id, tsl, ccds, swissprot, trembl, uniparc, sift, polyphen, hgvs_offset, ancestral_allele, cadd16, cadd37_16, ensembl_protein_id,\
            ensembl_transcriptid, fathmm_converted_rankscore, fathmm_pred, fathmm_score, mutationtaster_rankscore, mutationtaster_pred, mutationtaster_score, polyphen2_hdiv_pred, polyphen2_hdiv_rank, polyphen2_hdiv_score, primateai_pred, primateai_rankscore, primateai_score, revel_rank, revel_score, sift_rank, sift_pred, sift_score, tsl2, \
            uniprot_acc, vest4_rank, vest4_score, clinvar_clnsig, clinvar_id, genename_dbnsfp, lof, lof_filter, lof_flags, lof_info = ta.split('|')

        #REVEL
        if revel_score != '':
            qrevel_score = float(revel_score)
        
        is_lof = is_lof or ((annotation_impact == 'HIGH') and (lof == 'HC'))


    return (qrevel_score, is_lof)

def csq_parse(cqs, good_transcripts):

    transcript_annotation = [e for e in cqs.split(',') if not e == '']

    result = []

    result_by_transcript = {}
    result_generic = [-1.0, '', -1.0, -1.0,'','']


    for ta in transcript_annotation:
        #(allele=string,annotation=string,annotation_impact=string,gene_name=string,gene_id=string,feature_type=string,feature_id=string,transcript_biotype=string,exon=string,intron=string,hgvsc=string,hgvsp=string,cdna_relpos=string,cds_relpos=string,aa_relpos=string,amino_acids=string,codos=string,existing=string,distance=string,strand=string,
        #    flags=string,symbol_source=string,hgnc_id=string,tsl=string,ccds=string,swissprot=string,trembl=string,uniparc=string,sift=string,polyphen=string,hgvs_offset=string,ancestral_allele=string,ensembl_proteinid=string,
        #    ensembl_transcriptid=string,fathmm_converted_rankscore=string,fathmm_pred=string,fathmm_score=string,mutationtaster_rankscore=string,mutationtaster_pred=string,mutationtaster_score=string,polyphen2_hdiv_pred=string,polyphen2_hdiv_rank=string,polyphen2_hdiv_score=string,primateai_pred=string,primateai_rankscore=string,primateai_score=string,
        #    revel_rank=string,revel_score=string,sift_rank=string,sift_pred=string,sift_score=string,tsl2=string,uniprot_acc=string,vest4_rank=string,vest4_score=string,genename_dbnsfp=string,lof=string,lof_filter=string,lof_flags=string,lof_info=string)'))



        allele, annotation, annotation_impact, gene_name, gene_id, feature_type, feature_id, transcript_biotype, exon, intron, hgvsc, hgvsp, cdna_relpos, cds_relpos, aa_relpos, amino_acids, codos, existing, distance, strand, flags, symbol_source, hgnc_id, tsl, ccds, swissprot, trembl, uniparc, sift, polyphen, hgvs_offset, ancestral_allele, cadd16, cadd37_16, ensembl_protein_id,\
            ensembl_transcriptid, fathmm_converted_rankscore, fathmm_pred, fathmm_score, mutationtaster_rankscore, mutationtaster_pred, mutationtaster_score, polyphen2_hdiv_pred, polyphen2_hdiv_rank, polyphen2_hdiv_score, primateai_pred, primateai_rankscore, primateai_score, revel_rank, revel_score, sift_rank, sift_pred, sift_score, tsl2, \
            uniprot_acc, vest4_rank, vest4_score, clinvar_clnsig, clinvar_id, genename_dbnsfp, lof, lof_filter, lof_flags, lof_info = ta.split('|')


        #VEP
        if not feature_id in result_by_transcript:
            result_by_transcript[feature_id] = {}

        result_by_transcript[feature_id].update({'exon':exon, 'cds_relpos':cds_relpos, 'amino_acids':amino_acids,'aa_relpos':aa_relpos, 'feature_id':feature_id, 'gene_name':gene_name, 'gene_id': gene_id, 'lof':lof, 'lof_filter':lof_filter, 'lof_info':lof_info, 'sift':sift, 'polyphen':polyphen, 'annotation':annotation, 'annotation_impact':annotation_impact, 'bad_transcript':not feature_id in good_transcripts, 'vest4_score':0.0, 'fathmm_score':0.0})

        #REVEL
        if revel_score != '':
            result_generic[0] = float(revel_score)
        

        if mutationtaster_pred != '':
            result_generic[1] = mutationtaster_pred
           

        if cadd16 != '':
            result_generic[2] = float(cadd16)

        if cadd37_16 != '':
            result_generic[3] = float(cadd37_16)

        if clinvar_clnsig != '':
            result_generic[4] = clinvar_clnsig
        
        if clinvar_id != '':
            result_generic[5] = clinvar_id

        #DBNSFP transcript:
        if ensembl_transcriptid != '':
            for tet, tvs, tfs in zip(ensembl_transcriptid.split('&'), vest4_score.split('&'), fathmm_score.split('&')):
                if not tet in result_by_transcript:
                    result_by_transcript[tet] = {'exon':'', 'cds_relpos':'', 'amino_acids':'', 'aa_relpos':'', 'feature_id':tet, 'gene_name':'', 'gene_id':'', 'lof':'', 'lof_filter':'', 'lof_info':lof_info, 'sift':'','polyphen':'', 'annotation':'', 'annotation_impact':'', 'bad_transcript':not tet in good_transcripts, 'fathmm_score':0.0, 'vest4_score':0.0}

                xres = {}
                if tvs != '.' and tvs != '':
                    xres['vest4_score'] = float(tvs)
                if tfs != '.' and tfs != '':
                    xres['fathmm_score'] = float(tfs)
                result_by_transcript[tet].update(xres)

        #result.append((gene_name,annotation, annotation_impact, mutationtaster_pred, sift, polyphen, vest4_score, fathmm_score, revel_score, lof, lof_filter, feature_id, ensembl_transcriptid, tsl, tsl2, feature_id in bad_transcripts))
   
    rt = []
    for row in list(result_by_transcript.values()):
        rt.append(tuple([row[x] for x in ['exon','cds_relpos','aa_relpos','amino_acids','feature_id', 'gene_name', 'gene_id', 'lof','lof_filter','lof_info', 'sift','polyphen', 'annotation','annotation_impact','bad_transcript', 'vest4_score','fathmm_score']]))


    return (tuple(result_generic), rt)
    #r a,b in result_by_transcript.items()])
    #return result



def add_exclude_reason(gene, **kwargs):
    gene = gene.Copy()    
    if not 'max_missingness' in kwargs:
        kwargs['max_missingness'] =0.5
    if not 'max_maf' in kwargs:
        kwargs['max_maf'] = 0.01
    if not 'ac_prob' in kwargs:
        kwargs['ac_prob'] = 0.5
    if not 'vbd_tech' in kwargs:
        kwargs['vbd_tech'] = 25

    exclude_status = kwargs.get('exclude_status','status')

    res = gene.Get(_.Without(_.exclude_reason, _.comb_note), _.Get(_.chrom, _.pos, _.ref, _.alt, _.exclude_qc, _.exclude_reason, _.comb_del, _.comb_exclude_qc, _.ma_reflow_count, _.comb_note, (_.bad_transcript == False).Sum()/'ntranscript', (_.depth < 6).Mean(), 
                                _[_.status_casecontrol!= "NA"].Get(_.posterior_pls[:,1,:] * 0.5 + _.posterior_pls[:,2,:]).Get(_.Sum() * 2.0, _.Count(), _.Mean()),
                                _.gaf_nonneuro_popmax, _.g_nonneuro_popmax,  
                                _.Get((_.depth < 6)/'missing_flag', exclude_status).GroupBy(exclude_status).Sort(exclude_status).Get(_.missing_flag.Count().Array(), _.missing_flag.Sum().Array()),
                                _.vbd_tech.ReplaceMissing(0.0), _.is_lof, _.lof_revel).Each(lambda *x: _add_exclude(x, kwargs), dtype=str)/'exclude_reason').Copy()
    return res

def _add_exclude(params, kwargs):

    chrom, pos, ref, alt, exclude_qc, reason, comb_del, comb_exclude_qc, ma_reflow_count, comb_note, ntranscript,  missingness, ac_prob, ncount, maf, gaf_nn_popmax, g_nn_popmax, total_count, missing_count,  vbd_tech, is_lof, lof_revel = params
    reasons = []
    if kwargs.get('exclude_qc', True) and exclude_qc:
        reasons.append('QC-excluded: %s' %reason)
    if kwargs.get('comb_del', True) and comb_del:
        reasons.append('Combined: %s' % comb_note)
    if kwargs.get('exclude_qc',True) and comb_exclude_qc:
        reasons.append('Combined with QC-excluded')
    if ma_reflow_count >= 2:
        reasons.append('Multi-allelic AF-REF low count >= 2')
        
    if vbd_tech >= kwargs.get('vbd_tech',25):
        reasons.append('VBD-tech: %.2f' % vbd_tech)

    if vbd_tech >= kwargs.get('vbd_tech_lowmaf', 15) and maf < 0.0005:
        reasons.append('VBD-tech-low: %.2f' % vbd_tech)

    if ntranscript == 0:
        reasons.append('Not in high-quality transcript')

    diff_miss_found = False
    if 'diff_miss' in kwargs and (kwargs['diff_miss'] >= 1 or kwargs['diff_miss'] == -1):
        compare_total = total_count[-1]
        compare_missing = missing_count[-1]
        for mc,tc in zip(missing_count[:-1], total_count[:-1]):
            if compare_missing == 0 and mc == 0:
                continue
            
            p = fisher.test1t(mc,tc-mc, compare_missing, compare_total-compare_missing)

            if (kwargs['diff_miss'] < 2 and p <= 1e-20) or (kwargs['diff_miss']>=2 and p <= float('10e-%d' % kwargs['diff_miss'])):
                diff_miss_found = True
                if kwargs['diff_miss'] >= 1:
                    reasons.append('Differential missingness: %g' % p)
                break
    if 'diff_miss' in kwargs and diff_miss_found == False and kwargs['diff_miss'] == -1:
        reasons.append('No differential missingness')
    
    if 'max_missingness' in kwargs and missingness >= kwargs['max_missingness']:
        reasons.append('Missingness: %.2f' % missingness)

    if 'pop_max_maf' in kwargs and (not gaf_nn_popmax is Missing and (min(gaf_nn_popmax, 1.0 - gaf_nn_popmax) >= kwargs['pop_max_maf'])):
        reasons.append('POPMAF is too high: %.3f in ADES, %.3f in GNOMAD-NONNEURO in %s' % (maf, numpy.nan if gaf_nn_popmax is Missing else gaf_nn_popmax, g_nn_popmax))
    
    if 'max_maf' in kwargs and (min(maf, 1.0 - maf) >= kwargs['max_maf']):
        reasons.append('MAF is too high: %.3f in ADES, %.3f in GNOMAD' % (maf, numpy.nan if gaf_nn_popmax is Missing else gaf_nn_popmax))

    if 'ac_prob' in kwargs and min(ac_prob, ncount * 2.0 - ac_prob) < kwargs['ac_prob']:
        reasons.append('Dosage sum: %f' % ac_prob)

    keep = False
    if 'revel' in kwargs and lof_revel < float(kwargs['revel']) and not keep:
        reasons.append("REVEL/LOF score too low: %.2f" % lof_revel)

    if 'revel_max' in kwargs and lof_revel >= float(kwargs['revel_max']) and not keep:
        reasons.append("REVEL/LOF score too high: %.2f" % lof_revel)

    if 'lof' in kwargs and is_lof == 0.0 and not keep:
        reasons.append('Not a LOF variant')

    if 'nolof' in kwargs and is_lof == 1.0:
        reasons.append('A LOF variant')

    return ';'.join(reasons)
    

def vname(name):
    name = name.replace('<=','_min')
    
    name = name.replace('+', 'plus')
    name = name.replace('-_', 'minus_')
    if name.endswith('-'):
        name = name[:-1] + 'minus'
    return util.valid_name(name)


def range_score(gene, score, **kwargs):#{{{
    if gene is None:
        return (0.5,{}, {}, {}, "nvar==0")

    stats = {}
    if kwargs.get('sample_removal',False): #remove samples with missingness > 0.8
        remove_samples = set()
        added = True
        xkwargs =kwargs.copy()
        xkwargs.pop('diff_miss',None)
        xkwargs.pop('max_missingness',None)
        
        while added:
            xgene = gene[~IsMissing(_.Get(kwargs['exclude_status']))][_.Get(kwargs['exclude_status']) != 'NA'][~(_.sample |In| list(remove_samples))]
            xgene = add_exclude_reason(xgene.Get('*', _.status.Each(lambda x: 'case_all' if (x == 1 or (isinstance(x,str) and 'case' in x)) else ('control_all' if (x == 0 or (isinstance(x,str) and 'control' in x)) else 'NA'),dtype=str)/'status_casecontrol', 
                                                   _.Get(kwargs['exclude_status']).Each(lambda x: 'case_all' if (x == 1 or (isinstance(x,str) and 'case' in x)) else ('control_all' if (x == 0 or (isinstance(x,str) and 'control' in x)) else 'NA'),dtype=str)/'exclude_status_casecontrol'), **xkwargs)[_.exclude_reason == ''].Copy()
       

            n = xgene.Get(_.sample, (_.depth < 6).Mean(dim=0)/'sample_missingness')[_.sample_missingness > 0.8].sample()
            added =  len(n) > 0
            remove_samples.update(list(n))
        gene = gene[~(_.sample |In| list(remove_samples))]
        stats['nremoved_samples_lowdepth'] = len(remove_samples)
        stats['removed_samples_lowdepth'] = ';'.join(list(remove_samples))

    gene = gene[~IsMissing(_.status)][_.status.ReplaceMissing('NA') != 'NA']
    egene = add_exclude_reason(gene.Get('*', _.status.Each(lambda x: 'case_all' if (x == 1 or (isinstance(x,str) and 'case' in x)) else ('control_all' if (x == 0 or (isinstance(x,str) and 'control' in x)) else 'NA'),dtype=str)/'status_casecontrol'), **kwargs).Copy()

    if 'keep_variants' in kwargs:
        egene = egene.Get(_.Without(_.exclude_reason),  _.Get(_.exclude_reason, _.Get((_.snp_id |In| kwargs['keep_variants']) | (_.Get(_.chrom.Cast(str) + ":" + _.pos.Cast(str) + ":" + _.ref + ">" + _.alt) |In| kwargs['keep_variants']))).\
                                                            Each(lambda x,y: ';'.join(list(x.split(';')) + ['Manually included']) if y else x, dtype=str)/'exclude_reason').Copy()

    if 'exclude_variants' in kwargs:
        egene = egene.Get(_.Without(_.exclude_reason),  _.Get(_.exclude_reason, _.Get(((_.snp_id |In| kwargs['exclude_variants'])| (_.Get(_.chrom.Cast(str) + ":" + _.pos.Cast(str) + ":" + _.ref + ">" + _.alt) |In| kwargs['exclude_variants'])))).\
                                                            Each(lambda x,y: ';'.join(list(x.split(';')) + ['Manually excluded']) if y else x, dtype=str)/'exclude_reason').Copy()

    gene = egene[(_.exclude_reason == '') | (_.exclude_reason.HasPattern('Manually included'))].Copy()
    if len(gene.chrom()) == 0.0:
        return (0.5,{}, {}, {}, "nvar=0")
    gene = gene.Without(_.maf).Get('*', ((_.posterior_pls[:,1,:] + _.posterior_pls[:,2,:] * 2.0).Sum() / (_.sample.Count() * 2.0)).Each(lambda x: x if x < 0.5 else 1 - x,dtype=float)/'maf').Copy()

    nvar = len(gene.pos())
    
   
    stats['abhet_mean'] = gene.abhet_p.Mean()()
    stats['abhet_std'] = gene.abhet_p.Std()()
    stats['nvar'] = nvar
       
    posterior, status  = gene.Get(_.posterior_pls, _.status)()
    extra_vars = gene.Get(*score._vars)()

    if nvar == 0:
        posterior = posterior.reshape(0, 3, len(status))

    mmaf = (posterior[:,1,:] * 0.5 + posterior[:,2,:] * 1.0).mean(axis=1)

    posterior = posterior / posterior.sum(axis=1)[:,numpy.newaxis,:]
   
    #reverse direction of high freq vairants
    for pos, vmaf  in enumerate(mmaf):
        if vmaf > 0.5:
            sys.stderr.write("Reversing variant %d\n" % pos)
            reflik = posterior[pos,0,:].copy()
            homlik = posterior[pos,2,:].copy()
            posterior[pos,0,:] = homlik
            posterior[pos,2,:] = reflik

    variants = egene[_.bad_transcript == False].Get((_.chrom.Cast(str) + ':' + _.pos.Cast(str) + ':' + _.ref + ':' + _.alt)/'variant_id', (_.depth < 6).Mean()/'missingness', _.lof_revel, _.is_lof, _.revel_score, \
                            _.feature_id.Array().Each(lambda x: ';'.join(x),dtype=str), _.lof.Array().Each(lambda x: ';'.join(x),dtype=str), _.lof_info.Array().Each(lambda x: ';'.join(x),dtype=str), \
                            _.annotation.Array().Each(lambda x: ';'.join(x),dtype=str), _.annotation_impact.Array().Each(lambda x: ';'.join(x),dtype=str),_.exclude_reason,\
                            _.Get('*', (_.posterior_pls[:,1,:] + _.posterior_pls[:,2,:] * 2.0)/'dosage').Get(_.Without(_.maf), (_[_.status_casecontrol != 'NA'].dosage.Mean() / 2.0)/'maf').Get(_.maf,\
                                _[_.status_casecontrol == 'case_all'].Get((_.depth < 6).Mean()/'case_missingness', _.dosage.Sum()/'case_dosage_sum', _[((_.maf < 0.5) & (_.dosage > 0.5)) | ((_.maf > 0.5) & (_.dosage < 1.5))].sample/'cases'),\
                                _[_.status_casecontrol == 'control_all'].Get((_.depth < 6).Mean()/'control_missingness', _.dosage.Sum()/'control_dosage_sum', _[((_.maf < 0.5) & (_.dosage > 0.5)) | ((_.maf > 0.5) & (_.dosage < 1.5))].sample/'controls'))).Array().Dict()()

    if nvar == 0:
        return (0.5,stats, {}, variants, "nvar==0")

    var_weights = numpy.ones(mmaf.shape,dtype=float)

    #prepare sampling
    kwargs = kwargs.copy()
    if 'status' in kwargs:
        del kwargs['status']
    xscore = score(status, extra_vars, **kwargs)
    xscore.correctCovariates()
    
    #prepare sampling
    posterior_ref = posterior[:,0,:] # reference probability
    posterior_het = posterior[:, :2,:].sum(axis=1)

    stats['cmac_all'] = (posterior[:,1,:] + 2.0 * posterior[:,2,:]).sum()
    stats['n'] = posterior.shape[2]

    gt = posterior.argmax(axis=1)
    
    maxdiff = 1.0 / 1000000
    fil_ref = (posterior_ref <= (1.0 - maxdiff)) & (posterior_ref >= maxdiff) #only sample genotypes which have at least some uncertainty
    fil_het = (posterior_het <= (1.0 - maxdiff)) & (posterior_het >= maxdiff) #only sample genotypes which have at least some uncertainty
    fil = fil_ref | fil_het


    filposterior_ref = posterior_ref[fil] #get ref sample probabilities
    filposterior_het = posterior_het[fil] #get het sample probabilities

    xgt = gt * var_weights[:,numpy.newaxis] #weight genotypes by pathogenicity
    filvar_weights = (numpy.ones(xgt.shape,dtype=float) * var_weights[:,numpy.newaxis])[fil] #var_weights for sampled positions


    sys.stderr.write('-> missingness sampling: %d/%d: %.2g\n' % (len(fil), len(xgt.ravel()), float(len(fil)) / float(len(xgt.ravel())) ))
    sys.stderr.write('-> sampling and testing')
    sys.stderr.flush()

    #prepare sum and index for non-sampled positions (speedup)
    xgt[fil] = 0.0
    
    init_sum = xgt.sum(axis=0)
    fil_index = numpy.where(fil)[1]
    pvalue, score_stats = direct_aggregate_score(init_sum, fil_index, filposterior_ref, filposterior_het, filvar_weights, xscore, high_precision=kwargs.get('high_precision',False), perform_permutations=kwargs.get('perform_permutations',0))
       
    if not isinstance(score_stats,dict):
        score_stats = {'beta': score_stats}
    
    score_stats_new = {}
    for key,value in list(score_stats.items()):
        score_stats_new[vname(key)] = value
        if key.endswith('beta'):
            score_stats_new[vname(key[:-4] + 'or')] = numpy.exp(value)

    sys.stderr.write('done\n')
    sys.stderr.flush()
    return (pvalue, stats, score_stats_new,variants, "succes")#}}}

def direct_aggregate_score(init_sum, fil_index, filposterior_ref, filposterior_het, filphred, xscore, high_precision=False, perform_permutations=0):#{{{
    r = numpy.random.uniform(size=filposterior_ref.shape)

    nr = len(r.ravel())#number of samples
    maxroll = int(nr / 25.0) #number of random rolls before resampling
    if maxroll < 1: 
        maxroll = 1
   
    effects = 0

    pvalues = []
    effects = []
    i = 0

    permutes = [numpy.random.permutation(len(xscore.labels)) for p in range(perform_permutations)]
    per_pvalues = numpy.zeros(len(permutes),dtype=float)
    while True:
        sys.stderr.write('.')
        sys.stderr.flush()

        i += 1
        if (i % maxroll) == 0:
            r = numpy.random.uniform(size=filposterior_ref.shape)
        else:    
            r = numpy.roll(r, numpy.random.randint(0, nr))
        
        samples =  numpy.cast[float]((r > filposterior_ref)) * filphred + numpy.cast[float]((r > filposterior_het)) * filphred

        score = init_sum + numpy.bincount(fil_index, weights=samples, minlength=len(init_sum))
        score = score.reshape((len(score),1))
        #score = score / float(max(score.sum(),10e-6))


        
        pvalue, effect = xscore.stat(score)
        pvalues.append(pvalue)
        effects.append(effect)
        sample_mean_std = numpy.std(numpy.log10(pvalues)) / numpy.sqrt(len(pvalues))

        if perform_permutations:
            for pos, pn in enumerate(permutes):
                xpvalue,xstat = xscore.stat(score, perm_idx=pn)
                per_pvalues[pos] += xpvalue


        
        if not high_precision:
            if numpy.mean(pvalues) > 0.1 and i >= 5:
                break

            if numpy.mean(pvalues) > 0.01 and i >= 25:
                break

            if  sample_mean_std < 0.01 and i >= 15:
                break

            if i >= 100:
                break
        else:
            if i >= 250:
                break


    per_pvalues = per_pvalues / len(pvalues)
    if isinstance(effects[0],dict):
        res = {}
        for key in list(effects[0].keys()):
            if key == 'se' or key.endswith('_se'):
                res[key] = numpy.sqrt(numpy.mean([effect[key] ** 2.0 for effect in effects]))
            else:
                res[key] = numpy.mean([effect[key] for effect in effects])
            res[key + '_samples'] = numpy.array([effect[key] for effect in effects if key in effect])
        res['pvalue_samples'] = numpy.array(pvalues,dtype=float)
        res['nsamples'] = len(pvalues)
        res['pvalue_mean_std'] = sample_mean_std
        if perform_permutations:
            res['permutation_pvalues'] = per_pvalues
        effects = res
    else:
        effects = numpy.mean(effects)
    return (numpy.mean(pvalues), effects)#}}}


#}}}






if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--vcf', help='The VCF file or BCF file. If there are multiple files for the different chromosomes, one can use a wildcard (%%s) which will be used to automatically fill in the correct chromosome (e.g. result.chr%%s.bcf). The VCF/BCF files should have an index.')
    parser.add_argument('--meta', required=True, help='TSV file with at least a column iid and status, containing respectively the sample id an the phenotype.  Phenotype can be coded as case<=65, case>65, control.')
    parser.add_argument('--pca', default='-', help='TSV file with at least a column iid and columns pc1,pc2,pc3,pc4,pc5 and pc6.')
    parser.add_argument('--genes', required=True, help='Comma separated list of genes to analyze.')
    parser.add_argument('--variants',required=True, help='Variant category to analyze. Either "LOF", or REVEL threshold between 0 and 100 (includes also LOF variants).')
    parser.add_argument('--max_maf',type=float, default=0.01, help='Maximum minor allele frequency of included variants (default: 0.01).')
    parser.add_argument('--max_missingness',type=float, default=0.2, help='Maximum missingness (depth < 6) for included variants (default: 0.2).')
    parser.add_argument('--diff_miss',type=float, default=20.0, help='Differential missingness between cases and controls, as -log10 p-value threshold (default: 20). Note that the default may be overly permissive for smaller studies (and vice-versa).')
    parser.add_argument('--exclude_qc',action='store_true',default=True, help='Remove all variants denoted with the EXCLUDE_QC flag in the VCF (default: True).')
    parser.add_argument('--sample_removal',action='store_true',default=False, help='Remove samples with >80%% missingness in the selected variants of the analyzed gene (default: False).')
    parser.add_argument('--pheno',default='status',help='Name of the column in the meta file that carries the phenotype information')
    parser.add_argument('--causative_mutation',action='store_true',default=False, help="Indicates that there is a column 'causative_mutation_scan' in the pheno file, with samples that carry a causative mutation and need to be removed from the anaysis (default: False).")
    parser.add_argument('--transcripts_file',default='ades_v29_v19_whitelist_notsl_all_mapped_transcripts.tsv', help='File with transcripts that should be considered. Default: ades_v29_v19_whitelist_notsl_all_mapped_transcripts.tsv (see publication)')
    parser.add_argument('--genes_file',default='gene_data_b37.tsv', help='File with list of genes and their locations. Default: gene_data_b37.tsv (see publication). Adapt when using this cript with a different reference genome (see also --transcripts_file).')
    
    args = parser.parse_args()
    vcf_filename = os.path.expanduser(args.vcf)
    meta_filename = os.path.expanduser(args.meta)
    pca_filename = os.path.expanduser(args.pca)
    genes = args.genes.split(',')
    transcripts_file = args.transcripts_file
    genes_file = args.genes_file

    category = args.variants

    
    kwargs = {'score': score_models.DirectROrderLogitScore, 
              'max_missingness': args.max_missingness, 
              'ac_prob':0.5, 
              'max_maf':args.max_maf,
              'pop_max_maf':args.max_maf,
              'exclude_qc':args.exclude_qc,
              'exclude_status':'status', 
              'diff_miss':args.diff_miss, 
              'causative_mutation':args.causative_mutation,
              'sample_removal':args.sample_removal,
              'status':args.pheno}

    if pca_filename == '-':
        pca_filename = None
        kwargs['use_pca'] = 0

    if category == "LOF":
        kwargs['lof'] = True
    else:
        kwargs['revel'] = float(category) / 100.0


    gene_info = get_gene_data(genes_file)
    good_transcripts = read_good_transcripts(filename=transcripts_file)
    meta_info = get_meta_data(filename=meta_filename, pca_filename = pca_filename)
    
    
    if 'causative_mutation' in kwargs:
        if kwargs['causative_mutation'] == True:
            meta_info = meta_info[_.causative_mutation_scan.Cast("real64?").ReplaceMissing(0) == 0]
    
    selected_samples = meta_info.sample()

    for gene_name in genes:
        gene_stable_id, chrom, gene_start, gene_end = gene_info[_.gene_name == gene_name][0].Get(_.gene_stable_id, _.chrom, _.start, _.end)()
        sys.stderr.write('\n\nLoading %s:\n' % gene_name)
        range_variants = read_ranges(vcf_filename, meta_info, good_transcripts, selected_samples, [(chrom, gene_start, gene_end)], geneid_ann =gene_stable_id, **kwargs)

        if 'status' in kwargs and kwargs['status'] != 'status' and range_variants is not None: #rename pheno field
            range_variants = range_variants.Get(_.Get(kwargs['status'])/'status', _.Without(_.status, kwargs['status'])).Copy()
            #workaround 
            names = list(range_variants.Names)
            names[0] = 'status'
            range_variants = range_variants/tuple(names)


        sys.stderr.write('\nPerforming test on %s:\n' % gene_name)
        pvalue, stats, score_stats, variants, status = range_score(range_variants, **kwargs)

        if variants:
            variant_filename = '%s_variants.tsv' % gene_name
            sys.stderr.write('\nSaving variants to %s\n' % variant_filename)
            Save(Rep(variants), variant_filename)
        
        print(('\n\ngene=%s' % gene_name))
        print(('category=%s' % category))
        print(('status=%s' % status))
        if status == 'succes':
            print(('p-value=%.6g' % pvalue))
            print(('n_variant=%d' % stats['nvar']))
            print(('n_sample=%d' % stats['n']))
            print(('n_selected_alleles=%.1f' % stats['cmac_all']))
            if stats['abhet_mean'] is Missing:
                print('allele_balance=NOT AVAILABLE (provide abhet_p field in VCF)')
            else:
                print(('allele_balance=%.2f' % stats['abhet_mean']))
            print(('beta=%.3f' % score_stats['beta']))
            print(('se=%.3f' % score_stats['se']))



