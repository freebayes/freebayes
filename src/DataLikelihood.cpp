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
            prob += log(1 - exp(observedAllele.lnquality));
        } else {
            prob += observedAllele.lnquality;
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

    //cout << endl << endl << "genotype: "  << genotype << endl;
    int observationCount = observedAlleles.size();
    vector<long double> alleleProbs = genotype.alleleProbabilities();
    vector<long double> probs;
    vector<Allele> uniqueAlleles = genotype.uniqueAlleles();
    vector<vector<Allele*> > combos = multichoose_ptr(observationCount, uniqueAlleles);
    for (vector<vector<Allele*> >::iterator combo = combos.begin(); combo != combos.end(); ++combo) {
        //cout << endl << "combo: " << *combo << endl;
        vector<vector<Allele*> > trueAllelePermutations = multipermute(*combo);
        for (vector<vector<Allele*> >::iterator perm = trueAllelePermutations.begin(); perm != trueAllelePermutations.end(); ++perm) {
            vector<Allele*>& trueAlleles = *perm;
            //cout << "permu: " << *perm << endl;
            map<string, int> trueAlleleCounts = countAllelesString(trueAlleles);
            long double lnTrueAllelePermutationsCount = log(trueAllelePermutations.size());
            vector<int> observationCounts; // counts of 'observations' of true alleles, ordered according to the genotype's internal ordering
            for (Genotype::const_iterator g = genotype.begin(); g != genotype.end(); ++g) {
                map<string, int>::const_iterator count = trueAlleleCounts.find(g->allele.currentBase);
                if (count != trueAlleleCounts.end()) {
                    observationCounts.push_back(count->second);
                } else {
                    observationCounts.push_back(0);
                }
            }
            //cout << "multinomial: " << exp(multinomialln(alleleProbs, observationCounts)) << endl;
            //cout << "likelihood: " << exp(likelihoodGivenTrueAlleles(observedAlleles, trueAlleles)) << endl;
            probs.push_back(multinomialln(alleleProbs, observationCounts) 
                    + likelihoodGivenTrueAlleles(observedAlleles, trueAlleles)
                    - lnTrueAllelePermutationsCount);
        }
    }
    //cout << "l = " << logsumexp(probs) << endl;
    return logsumexp(probs);
}

// p(out are wrong) * p(in are correct) * (multinomial(in | genotype))
// sum(Q)           * sum(1 - Q)        *          ...
/*
long double
approxProbObservedAllelesGivenGenotype(
        vector<Allele*>& observedAlleles,
        Genotype& genotype
        ) {

    int observationCount = observedAlleles.size();
    vector<long double> alleleProbs = genotype.alleleProbabilities();
    vector<int> observationCounts = genotype.alleleCountsInObservations(observedAlleles);

    long double probInAllCorrect = 0; // 
    long double probOutAllWrong  = 0; // 
    for (vector<Allele*>::iterator obs = observedAlleles.begin(); obs != observedAlleles.end(); ++obs) {
        Allele& observation = **obs;
        if (genotype.containsAllele(observation)) {
            probInAllCorrect += log(1 - exp(observation.lnquality));
        } else {
            probOutAllWrong += observation.lnquality;
        }
    }
    
    if (sum(observationCounts) == 0) {
        return probOutAllWrong;
    } else {
        return probInAllCorrect + probOutAllWrong + multinomialln(alleleProbs, observationCounts);
    }

}
*/

long double
approxProbObservedAllelesGivenGenotype(
        vector<Allele*>& observedAlleles,
        Genotype& genotype
        ) {

    int observationCount = observedAlleles.size();
    vector<long double> alleleProbs = genotype.alleleProbabilities();
    vector<int> observationCounts = genotype.alleleCountsInObservations(observedAlleles);

    long double probInAllCorrect = 0; // 
    long double probOutAllWrong  = 0; // 
    for (vector<Allele*>::iterator obs = observedAlleles.begin(); obs != observedAlleles.end(); ++obs) {
        Allele& observation = **obs;
        if (!genotype.containsAllele(observation)) {
            probOutAllWrong += observation.lnquality;
        }
    }
    
    if (sum(observationCounts) == 0) {
        return probOutAllWrong;
    } else {
        return probOutAllWrong + multinomialln(alleleProbs, observationCounts);
    }

}

long double
bamBayesApproxProbObservedAllelesGivenGenotype(
        vector<Allele*>& observedAlleles,
        Genotype& genotype,
        long double dependenceFactor
        ) {

    int observationCount = observedAlleles.size();

    long double sumQout = 0;
    int countOut = 0;

    for (vector<Allele*>::iterator obs = observedAlleles.begin(); obs != observedAlleles.end(); ++obs) {
        Allele& observation = **obs;
        if (!genotype.containsAllele(observation)) {
            sumQout += observation.lnquality;
            ++countOut;
        }
    }

    if (countOut > 1) {
        sumQout *= (1 + (countOut - 1) * dependenceFactor) / countOut;
        //cerr << "dependencefactor," << observedAlleles.front()->sampleID << "," << IUPAC2GenotypeStr(IUPAC(genotype)) << "," << (1 + (countOut - 1) * dependenceFactor) / countOut << endl;
    }

    long double lnprob = sumQout - countOut * log(3);

    if (!genotype.homozygous()) {
        long double lnBinomialfactor = observationCount * log(0.5);
        //cerr << "lnprob," << observedAlleles.front()->sampleID << "," << IUPAC2GenotypeStr(IUPAC(genotype)) << "," << lnprob + lnBinomialfactor << endl;
        return lnprob + lnBinomialfactor;
    } else {
        //cerr << "lnprob," << observedAlleles.front()->sampleID << "," << IUPAC2GenotypeStr(IUPAC(genotype)) << "," << lnprob << endl;
        return lnprob;
    }

}

vector<pair<Genotype*, long double> >
exactProbObservedAllelesGivenGenotypes(
        vector<Allele*>& observedAlleles,
        vector<Genotype>& genotypes) {
    vector<pair<Genotype*, long double> > results;
    for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        results.push_back(make_pair(&*g, probObservedAllelesGivenGenotype(observedAlleles, *g)));
    }
    return results;
}

vector<pair<Genotype*, long double> >
approxProbObservedAllelesGivenGenotypes(
        vector<Allele*>& observedAlleles,
        vector<Genotype>& genotypes) {
    vector<pair<Genotype*, long double> > results;
    for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        results.push_back(make_pair(&*g, approxProbObservedAllelesGivenGenotype(observedAlleles, *g)));
    }
    return results;
}

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
                if (a->allele == *(observedAlleles.front()))
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


// bambayes style data likelihoods
vector<pair<Genotype*, long double> >
bamBayesProbObservedAllelesGivenGenotypes(
        vector<Allele*>& observedAlleles,
        vector<Genotype>& genotypes,
        long double dependenceFactor) {
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
                if (a->allele == *(observedAlleles.front()))
                    alleleRatio.first += 1;
            }
            map<pair<int, int>, long double>::iterator c = cachedProbsGivenAllelesInGenotype.find(alleleRatio);
            if (c != cachedProbsGivenAllelesInGenotype.end()) {
                prob = c->second;
            } else {
                prob = bamBayesApproxProbObservedAllelesGivenGenotype(observedAlleles, genotype, dependenceFactor);
                cachedProbsGivenAllelesInGenotype[alleleRatio] = prob;
            }
            results.push_back(make_pair(&genotype, prob));
        } else {
            results.push_back(make_pair(&genotype, bamBayesApproxProbObservedAllelesGivenGenotype(observedAlleles, genotype, dependenceFactor)));
        }
    }
    return results;
}
