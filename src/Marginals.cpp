#include "Marginals.h"


void marginalGenotypeLikelihoods(long double posteriorNormalizer, vector<GenotypeComboResult>& genotypeComboProbs, Results& results) {

    // push the marginal likelihoods into the rawMarginals vectors in the results
    for (vector<GenotypeComboResult>::iterator gc = genotypeComboProbs.begin(); gc != genotypeComboProbs.end(); ++gc) {
        for (GenotypeCombo::const_iterator i = gc->combo->begin(); i != gc->combo->end(); ++i) {
            const SampleDataLikelihood& sdl = **i;
            map<Genotype*, vector<long double> >& marginals = results[sdl.name].rawMarginals;
            marginals[sdl.genotype].push_back(gc->priorComboProb);
        }
    }

    // safely add the raw marginal vectors using logsumexp
    for (Results::iterator r = results.begin(); r != results.end(); ++r) {
        ResultData& d = r->second;
        for (map<Genotype*, vector<long double> >::iterator m = d.rawMarginals.begin(); m != d.rawMarginals.end(); ++m) {
            d.marginals[m->first] = logsumexp_probs(m->second) - posteriorNormalizer;
        }
    }

}
