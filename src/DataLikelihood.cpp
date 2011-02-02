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
            //cout << "multinomial: " << exp(multinomialSamplingProbLn(alleleProbs, observationCounts)) << endl;
            //cout << "likelihood: " << exp(likelihoodGivenTrueAlleles(observedAlleles, trueAlleles)) << endl;
            probs.push_back(multinomialSamplingProbLn(alleleProbs, observationCounts) 
                    + likelihoodGivenTrueAlleles(observedAlleles, trueAlleles)
                    - lnTrueAllelePermutationsCount);
        }
    }
    //cout << "l = " << logsumexp(probs) << endl;
    return logsumexp(probs);
}

long double
approxProbObservedAllelesGivenGenotype(
        Sample& sample,
        Genotype& genotype,
        long double dependenceFactor
        ) {

    int observationCount = sample.observationCount();
    vector<long double> alleleProbs = genotype.alleleProbabilities();
    vector<int> observationCounts = sample.alleleObservationCounts(genotype);
    int countOut = 0;
    long double prodQout = 0;  // the probability that the reads not in the genotype are all wrong

    for (Sample::iterator s = sample.begin(); s != sample.end(); ++s) {
        const string& base = s->first;
        if (!genotype.containsAllele(base)) {
            vector<Allele*>& alleles = s->second;
            for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a)
                prodQout += (*a)->lnquality;
            countOut += alleles.size();
        }
    }

    // read dependence factor, asymptotically downgrade quality values of
    // successive reads to dependenceFactor * quality
    if (countOut > 1) {
        prodQout *= (1 + (countOut - 1) * dependenceFactor) / countOut;
    }
    
    if (sum(observationCounts) == 0) {
        return prodQout;
    } else {
        return prodQout + multinomialSamplingProbLn(alleleProbs, observationCounts);
    }

}

long double
bamBayesApproxProbObservedAllelesGivenGenotype(
        Sample& sample,
        Genotype& genotype,
        long double dependenceFactor
        ) {

    int observationCount = sample.observationCount();

    long double prodQout = 0;
    int countOut = 0;

    for (Sample::iterator s = sample.begin(); s != sample.end(); ++s) {
        const string& base = s->first;
        if (!genotype.containsAllele(base)) {
            vector<Allele*>& alleles = s->second;
            for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a)
                prodQout += (*a)->lnquality;
            countOut += alleles.size();
        }
    }

    if (countOut > 1) {
        prodQout *= (1 + (countOut - 1) * dependenceFactor) / countOut;
        //cerr << "dependencefactor," << observedAlleles.front()->sampleID << "," << IUPAC2GenotypeStr(IUPAC(genotype)) << "," << (1 + (countOut - 1) * dependenceFactor) / countOut << endl;
    }

    long double lnprob = prodQout - countOut * log(3);

    if (!genotype.homozygous) {
        long double lnBinomialfactor = observationCount * log(0.5);
        //cerr << "lnprob," << observedAlleles.front()->sampleID << "," << IUPAC2GenotypeStr(IUPAC(genotype)) << "," << lnprob + lnBinomialfactor << endl;
        return lnprob + lnBinomialfactor;
    } else {
        //cerr << "lnprob," << observedAlleles.front()->sampleID << "," << IUPAC2GenotypeStr(IUPAC(genotype)) << "," << lnprob << endl;
        return lnprob;
    }

}

vector<pair<Genotype*, long double> >
approxProbObservedAllelesGivenGenotypes(
        Sample& sample,
        vector<Genotype>& genotypes,
        long double dependenceFactor) {
    vector<pair<Genotype*, long double> > results;
    for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        results.push_back(make_pair(&*g, approxProbObservedAllelesGivenGenotype(sample, *g, dependenceFactor)));
    }
    return results;
}

vector<pair<Genotype*, long double> >
probObservedAllelesGivenGenotypes(
        Sample& sample,
        vector<Genotype>& genotypes,
        long double dependenceFactor) {
    vector<pair<Genotype*, long double> > results;
    for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        long double prob;
        Genotype& genotype = *g;
        results.push_back(make_pair(&genotype, approxProbObservedAllelesGivenGenotype(sample, genotype, dependenceFactor)));
    }
    return results;
}

vector<pair<Genotype*, long double> >
bamBayesProbObservedAllelesGivenGenotypes(
        Sample& sample,
        vector<Genotype>& genotypes,
        long double dependenceFactor) {
    vector<pair<Genotype*, long double> > results;
    for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        long double prob;
        Genotype& genotype = *g;
        results.push_back(make_pair(&genotype, bamBayesApproxProbObservedAllelesGivenGenotype(sample, genotype, dependenceFactor)));
    }
    return results;
}

// NB: in practice, caching seems to not provide any perf boost
// .... but this needs to be tested now with the revised codebase
/*

// uses caching to reduce computation while generating the exact correct result

vector<pair<Genotype*, long double> >
probObservedAllelesGivenGenotypes(
        Sample& sample,
        vector<Genotype>& genotypes) {
    vector<pair<Genotype*, long double> > results;
    // cache the results on the basis of number of alleles in the genotype as a fraction of the genotype
    map<pair<int, int>, long double> cachedProbsGivenAllelesInGenotype;
    bool allSame = (sample.size() == 1);
    for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        long double prob;
        Genotype& genotype = *g;
        if (allSame) {
            pair<int, int> alleleRatio = make_pair(0, genotype.size());
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
        map<string, vector<Allele*>& observedAlleles,
        vector<Genotype>& genotypes,
        long double dependenceFactor) {
    vector<pair<Genotype*, long double> > results;
    // cache the results on the basis of number of alleles in the genotype as a fraction of the genotype
    map<pair<int, int>, long double> cachedProbsGivenAllelesInGenotype;
    bool allSame = areHomozygous(observedAlleles);
    for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        long double prob;
        Genotype& genotype = *g;
        if (allSame) {
            pair<int, int> alleleRatio = make_pair(0, genotype.size());
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
*/
