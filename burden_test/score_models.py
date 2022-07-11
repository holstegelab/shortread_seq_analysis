import numpy
import rpy2
import rpy2.robjects.packages as rpackages
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
lmtest = rpackages.importr('lmtest')
mass = rpackages.importr('MASS')

class Score(object):
    _vars = ['pca_europe']
    _type = 'aggregate'

    def __init__(self, status, args, one_sided=False, use_pca=1, **kwargs):
        if isinstance(args, tuple):
            pca_europe = args[0]
        else:
            pca_europe = args
        labels = status
        self.labels = labels
        
        self.use_pca = use_pca

        if use_pca == 0:
            self.feats = numpy.zeros((len(labels),0),dtype=float)
        else:
            self.feats = numpy.cast[float](pca_europe)
   
        self.intercept = numpy.zeros(len(labels),dtype=float)
        self.one_sided = one_sided
        self.ncovars = self.feats.shape[1]
        
    
    
    def correctCovariates(self):
        #correct for gender, age, crate new age label
        raise NotImplementedError


    def _get(self, perm_idx, sub_index):
        labels = self.labels
        intercept = self.feats

        if not sub_index is None:
            labels = labels[sub_index]
            intercept = intercept[sub_index]

        if not perm_idx is None:
            if perm_idx.dtype == bool:
                assert self.ncovars == 0
                labels = numpy.cast[labels.dtype](perm_idx)
            else:
                labels = labels[perm_idx]
                intercept = intercept[perm_idx]

        return (labels, intercept)




class DirectRLogitScore(Score):#{{{
    _type = 'direct_aggregate'
    
    def __init__(self, *args, **kwargs):
        Score.__init__(self, *args, **kwargs)

    def correctCovariates(self):
        return
    
    def stat(self, score, perm_idx=None, sub_index=None):
        labels, correction = self._get(perm_idx, sub_index)


        s = score.std()
        if s == 0:
            return (0.0, {'beta':0.0, 'se':1.0, 'lik':0.0})

        score = score - score.mean()
        normalization_factor = max(s, 1e-3)
        score = score / normalization_factor
        

        df = {}
        for col in range(correction.shape[1]):
            df['feat_%d' % col] = rpy2.robjects.FloatVector(correction[:,col])
        df['score'] = rpy2.robjects.FloatVector(score)
        df['labels'] = rpy2.robjects.BoolVector(labels)
        df = rpy2.robjects.DataFrame(df)


        fields = 'score + ' +  ' + '.join(['feat_%d' % col for col in range(correction.shape[1])])

        rpy2.robjects.r.assign('data',df)
        rpy2.robjects.r("fit <- glm(labels ~ %s, data=data, family=binomial(link='logit'))" % fields)
        rpy2.robjects.r('p <- coef(summary(fit))[2,4]')
        rpy2.robjects.r('c <- coef(fit)')

        coefs = rpy2.robjects.r['c']
        pvalues = rpy2.robjects.r['p']
        param = numpy.array(coefs)[1]/ normalization_factor


        pvalue = numpy.array(pvalues)[0]

        
        #rpy2.robjects.r('s <- summary(fit)')
        #print rpy2.robjects.r['s']
        rpy2.robjects.r('st <- summary(fit)$coefficients') 
        rpy2.robjects.r('lik <- logLik(fit)[1]')
        #print rpy2.robjects.r('summary(fit)')
        
        st = numpy.array(rpy2.robjects.r['st'])
        lik = float(rpy2.robjects.r['lik'])
        st = st[0,1] / normalization_factor
       
        stats = {'beta': param, 'se': st, 'lik': lik}


        return (pvalue, stats)  
#}}}

class DirectROrderLogitScore(Score):#{{{
    _type = 'direct_aggregate'
    
    def __init__(self,status,  *args, **kwargs):
        if isinstance(status[0],int):
            if numpy.max(status) == 1:
                res = []
                for s in status:
                    if s== 0:
                        res.append('control')
                    else:
                        res.append('case')
                status =  numpy.array(res)
            else:
                res = []
                for s in status:
                    if s== 0:
                        res.append('control')
                    else:
                        res.append(str(s))
                status =  numpy.array(res)
               
        Score.__init__(self, status, *args, **kwargs)

    def correctCovariates(self):
        return
    
    def stat(self, score, perm_idx=None, sub_index=None):
        labels, correction = self._get(perm_idx, sub_index)
        s = score.std()
        if s == 0:
            return (1.0, {'beta':0, 'se':0, 'lik':0})

        #print score.sum()
        score = score - score.mean()
        normalization_factor = max(s, 1e-3)
        score = score / normalization_factor

        df = {}
        for col in range(correction.shape[1]):
            df['feat_%d' % col] = rpy2.robjects.FloatVector(correction[:,col])
        df['score'] = rpy2.robjects.FloatVector(score)
        df['labels'] = rpy2.robjects.StrVector(labels)
        df = rpy2.robjects.DataFrame(df)


        rpy2.robjects.r.assign('data',df)
        rpy2.robjects.r("data$labels <- relevel(factor(data$labels), ref='control')")
        rpy2.robjects.r('l <- levels(data$labels)')
        l = list(rpy2.robjects.r['l'])
        starter = [e for e in l if '<=' in e]
        l = [e for e in l if not '<=' in e]
        l.sort()
        if starter:
            l = [starter[0]] + l
        
        rpy2.robjects.r.assign('l', l)


        #print (score, normalization_factor)
        popfields = ' + '.join(['feat_%d' % col for col in range(correction.shape[1])])
        if correction.shape[1] == 0:
            allfields = 'score'
            popfields = '1'
        else:
            allfields = 'score + ' + popfields
            

        rpy2.robjects.r('data$labels <- factor(data$labels, ordered=TRUE, levels=l)')
        rpy2.robjects.r('l <- levels(data$labels)')
        
        rpy2.robjects.r("nullfit <- polr(labels ~ %s, data=data, Hess=T)" % popfields)
        rpy2.robjects.r("fit <- polr(labels ~ %s, data=data, Hess=T)" % allfields)
        rpy2.robjects.r("res <- lrtest(fit, nullfit)")
        rpy2.robjects.r('c <- coef(fit)')
        rpy2.robjects.r('st <- summary(fit)$coefficients') 
        rpy2.robjects.r('lik <- logLik(fit)[1]')
        pvalue = numpy.array(rpy2.robjects.r['res'])[-1,-1]


        #print rpy2.robjects.r('summary(fit)')
        #print(rpy2.robjects.r('summary(fit)'))
        
        coefs = numpy.array(rpy2.robjects.r['c'])
        st = numpy.array(rpy2.robjects.r['st'])
        coefs = st[0,0] / normalization_factor
        st = st[0,1] / normalization_factor
        lik = float(rpy2.robjects.r['lik'])
        
        stats = {'beta': -coefs, 'se': st, 'lik': lik}

        return (pvalue, stats)  
#}}}


