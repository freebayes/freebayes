#include "Marginals.h"


void marginalGenotypeLikelihoods(long double posteriorNormalizer, list<GenotypeCombo>& genotypeCombos, Results& results) {


    map<string, map<Genotype*, vector<long double> > > rawMarginals;

    // push the marginal likelihoods into the rawMarginals vectors in the results
    for (list<GenotypeCombo>::iterator gc = genotypeCombos.begin(); gc != genotypeCombos.end(); ++gc) {
        // XXX for testing!!!
        if (gc->isHomozygous()) {
            continue;
        }
        for (GenotypeCombo::const_iterator i = gc->begin(); i != gc->end(); ++i) {
            const SampleDataLikelihood& sdl = **i;
            rawMarginals[sdl.name][sdl.genotype].push_back(gc->posteriorProb);
        }
    }

    // safely add the raw marginal vectors using logsumexp
    for (Results::iterator r = results.begin(); r != results.end(); ++r) {
        ResultData& sample = r->second;
        map<Genotype*, vector<long double> >& rawmgs = rawMarginals[r->first];
        for (map<Genotype*, vector<long double> >::iterator m = rawmgs.begin(); m != rawmgs.end(); ++m) {
            // TODO XXX if we are using only non-polymorphic combos, we probably should adjust the normalizer
            sample.marginals[m->first] = logsumexp_probs(m->second) - posteriorNormalizer;
        }
    }

}

void bestMarginalGenotypeCombo(GenotypeCombo& combo,
        Results& results,
        SampleDataLikelihoods& samples,
        long double theta,
        bool pooled,
        bool binomialObsPriors,
        bool alleleBalancePriors,
        long double diffusionPriorScalar) {

    for (SampleDataLikelihoods::iterator s = samples.begin(); s != samples.end(); ++s) {
        vector<SampleDataLikelihood>& sdls = *s;
        const string& name = sdls.front().name;
        const map<Genotype*, long double>& marginals = results[name].marginals;;
        map<Genotype*, long double>::const_iterator m = marginals.begin();
        long double bestMarginalProb = m->second;
        Genotype* bestMarginalGenotype = m->first;
        ++m;
        for (; m != marginals.end(); ++m) {
            if (abs(m->second) < abs(bestMarginalProb)) {
                bestMarginalProb = m->second;
                bestMarginalGenotype = m->first;
            }
        }
        for (vector<SampleDataLikelihood>::iterator d = sdls.begin(); d != sdls.end(); ++d) {
            if (bestMarginalGenotype == d->genotype) {
                combo.push_back(&*d);
                break;
            }
        }
    }

    combo.init(binomialObsPriors);
    combo.calculatePosteriorProbability(theta, pooled, binomialObsPriors,
            alleleBalancePriors, diffusionPriorScalar);

}
