#include "Marginals.h"


void marginalGenotypeLikelihoods(long double posteriorNormalizer, vector<GenotypeCombo>& genotypeCombos, Results& results) {

    // push the marginal likelihoods into the rawMarginals vectors in the results
    for (vector<GenotypeCombo>::iterator gc = genotypeCombos.begin(); gc != genotypeCombos.end(); ++gc) {
        for (GenotypeCombo::const_iterator i = gc->begin(); i != gc->end(); ++i) {
            const SampleDataLikelihood& sdl = **i;
            map<Genotype*, vector<long double> >& marginals = results[sdl.name].rawMarginals;
            marginals[sdl.genotype].push_back(gc->posteriorProb);
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
