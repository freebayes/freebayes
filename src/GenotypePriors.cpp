#include "GenotypePriors.h"


long double alleleFrequencyProbability(map<int, int> alleleFrequencyCounts, long double theta) {

    int M = 0;
    long double p = 1;

    for (map<int, int>::iterator f = alleleFrequencyCounts.begin(); f != alleleFrequencyCounts.end(); ++f) {
        int frequency = f->first;
        int count = f->second;
        M += frequency * count;
        p *= pow(theta, count) / pow(frequency, count) * factorial(count);
    }

    long double thetaH = 1;
    for (int h = 1; h < M; ++h)
        thetaH *= theta + h;

    return factorial(M) / (theta * thetaH) * p;

}

// Implements Ewens' Sampling Formula, which provides probability of a given
// partition of alleles in a sample from a population
long double alleleFrequencyProbabilityln(map<int, int> alleleFrequencyCounts, long double theta) {

    int M = 0; // multiplicity of site
    long double p = 0;
    long double thetaln = log(theta);

    for (map<int, int>::iterator f = alleleFrequencyCounts.begin(); f != alleleFrequencyCounts.end(); ++f) {
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


long double probabilityGenotypeComboGivenAlleleFrequencyln(GenotypeCombo& genotypeCombo, Allele& allele) {

    int n = genotypeCombo.numberOfAlleles();
    long double lnhetscalar = 0;

    for (GenotypeCombo::iterator gc = genotypeCombo.begin(); gc != genotypeCombo.end(); ++gc) {
        if (!gc->genotype->homozygous) {
            lnhetscalar += multinomialCoefficientLn(gc->genotype->ploidy, gc->genotype->counts());
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
        long double diffusionPriorScalar) {

        // when we are operating on pooled samples, we will not be able to
        // ascertain the number of heterozygotes in the pool,
        // rendering P(Genotype combo | Allele frequency) meaningless
        long double priorProbabilityOfGenotypeComboG_Af = 0;
        if (!pooled) {
            priorProbabilityOfGenotypeComboG_Af = probabilityGenotypeComboGivenAlleleFrequencyln(*combo, refAllele);
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
        long double priorComboProb = priorProbabilityOfGenotypeCombo + combo->prob;

        return GenotypeComboResult(combo,
                    priorComboProb,
                    combo->prob,
                    priorProbabilityOfGenotypeCombo,
                    priorProbabilityOfGenotypeComboG_Af,
                    priorProbabilityOfGenotypeComboAf);

}

void
genotypeCombinationsPriorProbability(
        vector<GenotypeComboResult>& genotypeComboProbs,
        vector<GenotypeCombo>& bandedCombos,
        Allele& refAllele,
        long double theta,
        bool pooled,
        long double diffusionPriorScalar) {

    for (vector<GenotypeCombo>::iterator c = bandedCombos.begin(); c != bandedCombos.end(); ++c) {

        GenotypeCombo* combo = &*c;

        genotypeComboProbs.push_back(genotypeCombinationPriorProbability(combo, refAllele, theta, pooled, diffusionPriorScalar));

    }
}
