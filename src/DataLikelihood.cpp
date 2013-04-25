#include "DataLikelihood.h"
#include "multichoose.h"
#include "multipermute.h"


long double
probObservedAllelesGivenGenotype(
        Sample& sample,
        Genotype& genotype,
        double dependenceFactor,
        bool useMapQ,
        Bias& observationBias,
        bool standardGLs,
        vector<Allele>& genotypeAlleles,
        double probContamination
    ) {

    int observationCount = sample.observationCount();
    vector<long double> alleleProbs = genotype.alleleProbabilities(observationBias);
    vector<int> observationCounts = genotype.alleleObservationCounts(sample);
    int countOut = 0;
    long double prodQout = 0;  // the probability that the reads not in the genotype are all wrong
    long double probObsGivenGt = 0;

    if (standardGLs) {
        for (Sample::iterator s = sample.begin(); s != sample.end(); ++s) {
            const string& base = s->first;
            if (!genotype.containsAllele(base)) {
                vector<Allele*>& alleles = s->second;
                if (useMapQ) {
                    for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                        // take the lesser of mapping quality and base quality (in log space)
                        prodQout += max((*a)->lnquality, (*a)->lnmapQuality);
                    }
                } else {
                    for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                        prodQout += (*a)->lnquality;
                    }
                }
                countOut += alleles.size();
            }
        }
    } else {
        for (Sample::iterator s = sample.begin(); s != sample.end(); ++s) { // for each read
            vector<Allele*>& alleles = s->second;
            for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                Allele& obs = **a;
                long double probi = 0;
                for (vector<Allele>::iterator b = genotypeAlleles.begin(); b != genotypeAlleles.end(); ++b) {
                    const string& base = b->currentBase;
                    double asampl = max(probContamination, min(genotype.alleleSamplingProb(base), 1-probContamination));
                    if (asampl == 0.5) {   // diploid-specific
                        if (b->isReference()) {
                            asampl += probContamination;
                        } else {
                            asampl -= probContamination;
                        }
                    }
                    if (obs.currentBase == base) {
                        probi += asampl * (1 - exp(obs.lnquality));
                    } else {
                        probi += asampl * exp(obs.lnquality);
                    }
                }
                probObsGivenGt += log(min(probi, (long double) 1.0)); // bound at 1
            }
        }
    }

    // read dependence factor, asymptotically downgrade quality values of
    // successive reads to dependenceFactor * quality
    if (standardGLs) {
        if (countOut > 1) {
            prodQout *= (1 + (countOut - 1) * dependenceFactor) / countOut;
        }

        if (sum(observationCounts) == 0) {
            return prodQout;
        } else {
            return prodQout + multinomialSamplingProbLn(alleleProbs, observationCounts);
        }
    } else {
        return probObsGivenGt;
    }

}


vector<pair<Genotype*, long double> >
probObservedAllelesGivenGenotypes(
        Sample& sample,
        vector<Genotype*>& genotypes,
        double dependenceFactor,
        bool useMapQ,
        Bias& observationBias,
        bool standardGLs,
        vector<Allele>& genotypeAlleles,
        double probContamination
    ) {
    vector<pair<Genotype*, long double> > results;
    for (vector<Genotype*>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
        results.push_back(
	    make_pair(*g,
                  probObservedAllelesGivenGenotype(
                      sample,
                      **g,
                      dependenceFactor,
                      useMapQ,
                      observationBias,
                      standardGLs,
                      genotypeAlleles,
                      probContamination)));
    }
    return results;
}
