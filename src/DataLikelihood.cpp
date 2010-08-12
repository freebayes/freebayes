#include "DataLikelihood.h"
#include "multichoose.h"
#include "multipermute.h"



// log probability of a given matching of true and observed alleles
// true alleles are the 'real' alleles
// observed alleles are whet we observe
// a mismatch between these two sets implies an error, which is scored accordingly here
long double likelihoodGivenTrueAlleles(vector<Allele*>& observedAlleles, vector<Allele*>& trueAlleles) {
    long double prob = 0;
    vector<Allele*>::iterator o = observedAlleles.begin();
    vector<Allele*>::iterator t = trueAlleles.begin();
    for ( ; o != observedAlleles.end() && t != trueAlleles.end(); ++o, ++t)
    {
        Allele& observedAllele = **o;
        Allele& trueAllele = **t;
        if (observedAllele == trueAllele) {
            prob += log(1 - exp(observedAllele.lncurrentQuality()));
        } else {
            prob += observedAllele.lncurrentQuality();
        }
    }
    return prob;
}

/*
  'Exact' data likelihood, sum of sampling probability * joint Q score for the
  observed alleles over all possible underlying 'true allele' combinations."""
*/
long double
probObservedAllelesGivenGenotype(
        vector<Allele*>& observedAlleles,
        Genotype& genotype) {

    int observationCount = observedAlleles.size();
    vector<long double> alleleProbs = genotype.alleleProbabilities();
    vector<long double> probs;
    vector<Allele> uniqueAlleles = genotype.uniqueAlleles();
    vector<vector<Allele*> > combos = multichoose_ptr(observationCount, uniqueAlleles);
    for (vector<vector<Allele*> >::iterator combo = combos.begin(); combo != combos.end(); ++combo) {
        vector<vector<Allele*> > trueAllelePermutations = multipermute(*combo);
        for (vector<vector<Allele*> >::iterator perm = trueAllelePermutations.begin(); perm != trueAllelePermutations.end(); ++perm) {
            vector<Allele*>& trueAlleles = *perm;
            map<string, int> trueAlleleCounts = countAllelesString(trueAlleles);
            vector<int> observationCounts; // counts of 'observations' of true alleles, ordered according to the genotype's internal ordering
            for (Genotype::const_iterator g = genotype.begin(); g != genotype.end(); ++g) {
                map<string, int>::const_iterator count = trueAlleleCounts.find(g->first.base());
                if (count != trueAlleleCounts.end()) {
                    observationCounts.push_back(count->second);
                } else {
                    observationCounts.push_back(0);
                }
            }
            probs.push_back(multinomialln(alleleProbs, observationCounts) 
                    + likelihoodGivenTrueAlleles(observedAlleles, trueAlleles));
        }
    }
    //cout << "l = " << logsumexp(probs) << endl;
    return logsumexp(probs);
}

/*
 Approximate data likelihoods.
approach: drop lots of terms, normalize
 */
long double
approxProbObservedAllelesGivenGenotype(
        vector<Allele*>& observedAlleles,
        Genotype& genotype
        ) {
        //long double fractionCorrect) {

    // approximation approach
    // condition the probability on a configurable fraction of errors
    // ... you can establish the probability that a given set of observations
    // ... has no errors
    // ... has <N errors
    // ... has enough errors to support the genotype
    //
    // a single term can then be scaled by the dirichlet maximum likelihood
    // ratio for the given observation draws from the specified genotype
    //
    // the gist of this is that we should ask "what is the probability that we
    // have enough errors to support the specific genotype"?
    // an example would be the case in which we have no overlap between genotype and observations
    // then, the only way for the specific genotype to be supported by the observations
    // is that the observations are all wrong
    //
    
    // 1) how many observations support the genotype?
    // 2) how many don't
    // 3) what is the probability that the ones which don't support the genotype 
    //    are all wrong, and the ones that do have no wrongs?
    //    ... alternatively the probability that at least one of the ones which
    //    support is right, which may scale better in accordance with the exact
    //    data likelihood
    //
    // 4) multiply the above by the dirichlet maximum likelihood ratio
    
    int observationCount = observedAlleles.size();
    vector<long double> alleleProbs = genotype.alleleProbabilities();
    //vector<Allele> uniqueAlleles = genotype.uniqueAlleles();
    vector<int> observationCounts = genotype.alleleCountsInObservations(observedAlleles);
    //map<Allele, int> inCounts;
    //map<Allele, int> outCounts;
    //genotype.inOutObservationCounts(observedAlleles, inCounts, outCounts);


    /*
    cerr << "genotype " << genotype << endl;
    for (vector<Allele*>::iterator i = observedAlleles.begin(); i != observedAlleles.end(); ++i) {
        cerr << " " << *i;
    }
    cerr << endl;
    */

    // the approximate observation probability is our rough estimate of the
    // probability that we have some support for the genotype within our
    // observations
    long double probInAllWrong = 0; // 
    long double probOutAreAllErrors = 0; // 
    for (vector<Allele*>::iterator obs = observedAlleles.begin(); obs != observedAlleles.end(); ++obs) {
        Allele& observation = **obs;
        if (genotype.containsAllele(observation)) {
            probInAllWrong += observation.lncurrentQuality();
        } else {
            probOutAreAllErrors += observation.lncurrentQuality();
        }
    }
    
    if (sum(observationCounts) == 0) {
        //cerr << "prob = " << exp(probOutAreAllErrors) << endl;
        return probOutAreAllErrors;
    } else {
        long double probInAllCorrect = log(1 - exp(probInAllWrong));
        long double approxObsProbability = log(1 - exp(probInAllWrong)) + probOutAreAllErrors;
        // this is the probability that we have at least some support for the
        // genotype among the observations, at least one observation
        long double dirichlet = log(dirichletMaximumLikelihoodRatio(alleleProbs, observationCounts));
        /*
        vector<long double>::iterator g = alleleProbs.begin();
        vector<int>::iterator q = observationCounts.begin();
        for (; q != observationCounts.end() && g != alleleProbs.end(); ++q, ++g) {
            cerr << *g << " " << *q << ", ";
        }
        cerr << endl;
        cerr << "dirichlet maximum likelihood ratio = " << exp(dirichlet) << endl;
        cerr << "probInAllWrong = " << exp(probInAllWrong) << endl;
        cerr << "probInAllCorrect = " << exp(probInAllCorrect) << endl;
        cerr << "probOutAreAllErrors = " << exp(probOutAreAllErrors) << endl;
        cerr << "data likelihood of " << genotype << " = "  << exp(approxObsProbability + dirichlet) << endl;
        */
        return approxObsProbability + dirichlet;
    }

}


/*
vector<pair<Genotype, long double> >
probObservedAllelesGivenGenotypes(
        vector<Allele*>& observedAlleles,
        vector<Genotype>& genotypes) {
    vector<pair<Genotype, long double> > results;
    for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        results.push_back(make_pair(*g, approxProbObservedAllelesGivenGenotype(observedAlleles, *g)));
    }
    return results;
}
*/

// uses caching to reduce computation while generating the exact correct result

vector<pair<Genotype*, long double> >
probObservedAllelesGivenGenotypes(
        vector<Allele*>& observedAlleles,
        vector<Genotype>& genotypes) {
    vector<pair<Genotype*, long double> > results;
    // cache the results on the basis of number of alleles in the genotype as a fraction of the genotype
    map<pair<int, int>, long double> cachedProbsGivenAllelesInGenotype;
    bool allSame = (countAlleles(observedAlleles).size() == 1);
    for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        long double prob;
        Genotype& genotype = *g;
        if (allSame) {
            pair<int, int> alleleRatio = make_pair(0, genotype.uniqueAlleles().size());
            //alleleRatio(observedAlleles, genotype);
            for (Genotype::iterator a = genotype.begin(); a != genotype.end(); ++a) {
                if (a->first == *(observedAlleles.front()))
                    alleleRatio.first += 1;
            }
            map<pair<int, int>, long double>::iterator c = cachedProbsGivenAllelesInGenotype.find(alleleRatio);
            if (c != cachedProbsGivenAllelesInGenotype.end()) {
                prob = c->second;
            } else {
                prob = approxProbObservedAllelesGivenGenotype(observedAlleles, genotype);
                cachedProbsGivenAllelesInGenotype[alleleRatio] = prob;
            }
            results.push_back(make_pair(&genotype, prob));
        } else {
            results.push_back(make_pair(&genotype, approxProbObservedAllelesGivenGenotype(observedAlleles, genotype)));
        }
    }
    return results;
}
