#include "GenotypePriors.h"

/*
long double alleleFrequencyProbability(const map<int, int>& alleleFrequencyCounts, long double theta) {

    int M = 0;
    long double p = 1;

    for (map<int, int>::const_iterator f = alleleFrequencyCounts.begin(); f != alleleFrequencyCounts.end(); ++f) {
        int frequency = f->first;
        int count = f->second;
        M += frequency * count;
        p *= (double) pow((double) theta, (double) count) / (double) pow((double) frequency, (double) count) * factorial(count);
    }

    long double thetaH = 1;
    for (int h = 1; h < M; ++h)
        thetaH *= theta + h;

    return factorial(M) / (theta * thetaH) * p;

}

AlleleFrequencyProbabilityCache alleleFrequencyProbabilityCache;

long double alleleFrequencyProbabilityln(const map<int, int>& alleleFrequencyCounts, long double theta) {
    return alleleFrequencyProbabilityCache.alleleFrequencyProbabilityln(alleleFrequencyCounts, theta);
}

// Implements Ewens' Sampling Formula, which provides probability of a given
// partition of alleles in a sample from a population
long double __alleleFrequencyProbabilityln(const map<int, int>& alleleFrequencyCounts, long double theta) {

    int M = 0; // multiplicity of site
    long double p = 0;
    long double thetaln = log(theta);

    for (map<int, int>::const_iterator f = alleleFrequencyCounts.begin(); f != alleleFrequencyCounts.end(); ++f) {
        int frequency = f->first;
        int count = f->second;
        M += frequency * count;
        p += powln(thetaln, count) - powln(log(frequency), count) + factorialln(count);
    }

    long double thetaH = 0;
    for (int h = 1; h < M; ++h)
        thetaH += log(theta + h);

    return factorialln(M) - (thetaln + thetaH) + p;

}
*/


long double probabilityGenotypeComboGivenAlleleFrequencyln(GenotypeCombo& genotypeCombo, Allele& allele) {

    int n = genotypeCombo.numberOfAlleles();
    long double lnhetscalar = 0;

    for (GenotypeCombo::iterator gc = genotypeCombo.begin(); gc != genotypeCombo.end(); ++gc) {
        SampleDataLikelihood& sgp = **gc;
        if (!sgp.genotype->homozygous) {
            lnhetscalar += multinomialCoefficientLn(sgp.genotype->ploidy, sgp.genotype->counts());
        }
    }

    return lnhetscalar - multinomialCoefficientLn(n, genotypeCombo.counts());

}


// core calculation of genotype combination likelihoods
//
GenotypeComboResult
genotypeCombinationPriorProbability(
        GenotypeCombo* combo,
        Allele& refAllele,
        long double theta,
        bool pooled,
        bool binomialObsPriors,
        bool alleleBalancePriors,
        long double diffusionPriorScalar) {

        // when we are operating on pooled samples, we will not be able to
        // ascertain the number of heterozygotes in the pool,
        // rendering P(Genotype combo | Allele frequency) meaningless
        long double priorProbabilityOfGenotypeComboG_Af = 0;
        if (!pooled) {
            priorProbabilityOfGenotypeComboG_Af = probabilityGenotypeComboGivenAlleleFrequencyln(*combo, refAllele);
        }

        long double priorObservationExpectationProb = 0;

        if (binomialObsPriors) {
            // for each alternate and the reference allele
            // calculate the binomial probability that we see the given strand balance and read placement prob
            vector<string> alleles = combo->alleles();
            // cerr << *combo << endl;
            for (vector<string>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                const string& allele = *a;
                map<string, AlleleCounter>::iterator ac = combo->alleleCounters.find(allele);
                if (ac != combo->alleleCounters.end()) {
                    const AlleleCounter& alleleCounter = ac->second;
                    int obs = alleleCounter.observations;
                    /*
                    cerr << allele <<  " counts: " << alleleCounter.frequency
                        << " observations " << alleleCounter.observations
                        << " " << alleleCounter.forwardStrand
                        << "," << alleleCounter.reverseStrand
                        << " " << alleleCounter.placedLeft
                        << "," << alleleCounter.placedRight
                        << " " << alleleCounter.placedStart
                        << "," << alleleCounter.placedEnd
                        << endl;
                    cerr << "priorObservationExpectationProb = " << priorObservationExpectationProb << endl;
                    cerr << "binprobln strand = " << binomialProbln(alleleCounter.forwardStrand, obs, 0.5) << endl;
                    cerr << "binprobln position = " << binomialProbln(alleleCounter.placedLeft, obs, 0.5) << endl;
                    cerr << "binprobln start = " << binomialProbln(alleleCounter.placedStart, obs, 0.5) << endl;
                    cerr << "priorObservationExpectationProb = " << priorObservationExpectationProb << endl;
                    */

                    priorObservationExpectationProb
                        += binomialProbln(alleleCounter.forwardStrand, obs, 0.5)
                        +  binomialProbln(alleleCounter.placedLeft, obs, 0.5)
                        +  binomialProbln(alleleCounter.placedStart, obs, 0.5);

                }
            }
            // ok... now do the same move for the observation counts
            // --- this should capture "Allele Balance"
        }

        if (alleleBalancePriors) {
            priorObservationExpectationProb += multinomialSamplingProbLn(combo->alleleProbs(), combo->observationCounts());
        }

        // with larger population samples, the effect of 
        // P(Genotype combo | Allele frequency) may bias us against reporting
        // true variants which are under selection despite overwhelming evidence
        // for variation.  this allows us to scale the effect of this prior
        if (diffusionPriorScalar != 1) {
            priorProbabilityOfGenotypeComboG_Af /= diffusionPriorScalar;
        }

        // Ewens' Sampling Formula
        long double priorProbabilityOfGenotypeComboAf = 
            alleleFrequencyProbabilityln(combo->countFrequencies(), theta);
        long double priorProbabilityOfGenotypeCombo = 
            priorProbabilityOfGenotypeComboG_Af + priorProbabilityOfGenotypeComboAf;
        long double priorComboProb = priorProbabilityOfGenotypeCombo + combo->prob + priorObservationExpectationProb;

        return GenotypeComboResult(combo,
                    priorComboProb,
                    combo->prob,
                    priorProbabilityOfGenotypeCombo,
                    priorProbabilityOfGenotypeComboG_Af,
                    priorProbabilityOfGenotypeComboAf,
                    priorObservationExpectationProb);

}

void
genotypeCombinationsPriorProbability(
        vector<GenotypeComboResult>& genotypeComboProbs,
        vector<GenotypeCombo>& bandedCombos,
        Allele& refAllele,
        long double theta,
        bool pooled,
        bool binomialObsPriors,
        bool alleleBalancePriors,
        long double diffusionPriorScalar) {

    for (vector<GenotypeCombo>::iterator c = bandedCombos.begin(); c != bandedCombos.end(); ++c) {

        GenotypeCombo* combo = &*c;

        genotypeComboProbs.push_back(
                genotypeCombinationPriorProbability(
                    combo,
                    refAllele,
                    theta,
                    pooled,
                    binomialObsPriors,
                    alleleBalancePriors,
                    diffusionPriorScalar));

    }
}
