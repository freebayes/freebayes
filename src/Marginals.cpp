#include "Marginals.h"


/*
// calculate marginals from the combos, store in results
void marginalGenotypeLikelihoods(list<GenotypeCombo>& genotypeCombos, Results& results) {


    map<string, map<Genotype*, vector<long double> > > rawMarginals;

    // push the marginal likelihoods into the rawMarginals vectors in the results
    for (list<GenotypeCombo>::iterator gc = genotypeCombos.begin(); gc != genotypeCombos.end(); ++gc) {
        for (GenotypeCombo::const_iterator i = gc->begin(); i != gc->end(); ++i) {
            const SampleDataLikelihood& sdl = **i;
            rawMarginals[sdl.name][sdl.genotype].push_back(gc->posteriorProb);
        }
    }

    // safely add the raw marginal vectors using logsumexp
    for (Results::iterator r = results.begin(); r != results.end(); ++r) {
        ResultData& sample = r->second;
        map<Genotype*, vector<long double> >& rawmgs = rawMarginals[r->first];
        vector<long double> probs;
        for (map<Genotype*, vector<long double> >::iterator m = rawmgs.begin(); m != rawmgs.end(); ++m) {
            probs.push_back(logsumexp_probs(m->second));
        }
        long double normalizer = logsumexp_probs(probs);
        vector<long double>::iterator p = probs.begin();
        for (map<Genotype*, vector<long double> >::iterator m = rawmgs.begin(); m != rawmgs.end(); ++m, ++p) {
            sample.marginals[m->first] = *p - normalizer;
        }
    }

}
*/

// recompute data likelihoods using marginals from the combos
// assumes that the genotype combos are in the same order as the likelihoods
// assumes that the genotype combos are the same size as the number of samples in the likelihoods
// returns the delta from the previous marginals, informative in the case of EM
long double marginalGenotypeLikelihoods(list<GenotypeCombo>& genotypeCombos, SampleDataLikelihoods& likelihoods) {

    long double delta = 0;

    vector< map<Genotype*, vector<long double> > > rawMarginals;
    rawMarginals.resize(likelihoods.size());
    vector< map<Genotype*, vector<long double> > >::iterator rawMarginalsItr;

    // push the marginal likelihoods into the rawMarginals maps
    for (list<GenotypeCombo>::iterator gc = genotypeCombos.begin(); gc != genotypeCombos.end(); ++gc) {
        rawMarginalsItr = rawMarginals.begin();
        for (GenotypeCombo::const_iterator i = gc->begin(); i != gc->end(); ++i) {
            const SampleDataLikelihood& sdl = **i;
            map<Genotype*, vector<long double> >& rmgs = *rawMarginalsItr++;
            rmgs[sdl.genotype].push_back(gc->posteriorProb);
        }
    }

    // safely add the raw marginal vectors using logsumexp
    // and use to update the sample data likelihoods
    rawMarginalsItr = rawMarginals.begin();
    for (SampleDataLikelihoods::iterator s = likelihoods.begin(); s != likelihoods.end(); ++s) {
        vector<SampleDataLikelihood>& sdls = *s;
        const map<Genotype*, vector<long double> >& rawmgs = *rawMarginalsItr++;
        map<Genotype*, long double> marginals;
        vector<long double> rawprobs;
        for (map<Genotype*, vector<long double> >::const_iterator m = rawmgs.begin(); m != rawmgs.end(); ++m) {
            long double p = logsumexp_probs(m->second);
            marginals[m->first] = p;
            rawprobs.push_back(p);
        }
        long double normalizer = logsumexp_probs(rawprobs);
        for (vector<SampleDataLikelihood>::iterator sdl = sdls.begin(); sdl != sdls.end(); ++sdl) {
            long double newmarginal = marginals[sdl->genotype] - normalizer;
            delta += newmarginal - sdl->marginal;
            sdl->marginal = newmarginal;
        }
    }

    return delta;

}

/*
void bestMarginalGenotypeCombo(GenotypeCombo& combo,
        Results& results,
        SampleDataLikelihoods& samples,
        long double theta,
        bool pooled,
        bool permute,
        bool hwePriors,
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
    combo.calculatePosteriorProbability(theta, pooled, permute, hwePriors, binomialObsPriors,
            alleleBalancePriors, diffusionPriorScalar);

}

// XXX probably not "balanced" :)
void balancedMarginalGenotypeLikelihoods(list<GenotypeCombo>& genotypeCombos, Results& results) {

    map<string, map<Genotype*, vector<long double> > > rawMarginals;

    for (list<GenotypeCombo>::iterator gc = genotypeCombos.begin(); gc != genotypeCombos.end(); ++gc) {
        bool isComboKing = true;
        // if it's not the combo king, we're going to record the probability for just the altered genotype(s)
        for (GenotypeCombo::const_iterator i = gc->begin(); i != gc->end(); ++i) {
            const SampleDataLikelihood& sdl = **i;
            if (sdl.rank != 0) {
                rawMarginals[sdl.name][sdl.genotype].push_back(gc->posteriorProb);
                isComboKing = false;
            }
        }
        // if it's the combo king, then record the marginal probability for each genotype
        if (gc->isHomozygous() || isComboKing) {
            for (GenotypeCombo::const_iterator i = gc->begin(); i != gc->end(); ++i) {
                const SampleDataLikelihood& sdl = **i;
                rawMarginals[sdl.name][sdl.genotype].push_back(gc->posteriorProb);
            }
        }
    }

    // safely add the raw marginal vectors using logsumexp
    for (Results::iterator r = results.begin(); r != results.end(); ++r) {
        ResultData& sample = r->second;
        map<Genotype*, vector<long double> >& rawmgs = rawMarginals[r->first];
        vector<long double> probs;
        for (map<Genotype*, vector<long double> >::iterator m = rawmgs.begin(); m != rawmgs.end(); ++m) {
            probs.push_back(logsumexp_probs(m->second));
        }
        long double normalizer = logsumexp_probs(probs);
        vector<long double>::iterator p = probs.begin();
        for (map<Genotype*, vector<long double> >::iterator m = rawmgs.begin(); m != rawmgs.end(); ++m, ++p) {
            sample.marginals[m->first] = *p - normalizer;
        }
    }

}
*/
