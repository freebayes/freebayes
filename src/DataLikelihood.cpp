#include "DataLikelihood.h"
#include "multichoose.h"
#include "multipermute.h"


long double
probObservedAllelesGivenGenotype(
        Sample& sample,
        Genotype& genotype,
        long double dependenceFactor,
        bool useMapQ) {

    int observationCount = sample.observationCount();
    vector<long double> alleleProbs = genotype.alleleProbabilities();
    vector<int> observationCounts = sample.alleleObservationCounts(genotype);
    int countOut = 0;
    long double prodQout = 0;  // the probability that the reads not in the genotype are all wrong

    for (Sample::iterator s = sample.begin(); s != sample.end(); ++s) {
        const string& base = s->first;
        if (!genotype.containsAllele(base)) {
            vector<Allele*>& alleles = s->second;
            if (useMapQ) {
                for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                    prodQout += (*a)->lnquality + (*a)->lnmapQuality;
                }
            } else {
                for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                    prodQout += (*a)->lnquality;
                }
            }
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

vector<pair<Genotype*, long double> >
probObservedAllelesGivenGenotypes(
        Sample& sample,
        vector<Genotype>& genotypes,
        long double dependenceFactor,
        bool useMapQ) {
    vector<pair<Genotype*, long double> > results;
    for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        results.push_back(make_pair(&*g, probObservedAllelesGivenGenotype(sample, *g, dependenceFactor, useMapQ)));
    }
    return results;
}
