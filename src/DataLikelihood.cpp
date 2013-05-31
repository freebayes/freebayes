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
        Contamination& contaminations,
        map<string, double>& freqs
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
            //if (!genotype.containsAllele(s->first)) {
            //    countOut += s->second.size();
            //}
            vector<Allele*>& alleles = s->second;
            for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
                Allele& obs = **a;
                long double probi = 0;
                ContaminationEstimate& contamination = contaminations.of(obs.readGroupID);
                long double qual = (1 - exp(obs.lnquality)) * (1 - exp(obs.lnmapQuality));
                for (vector<Allele>::iterator b = genotypeAlleles.begin(); b != genotypeAlleles.end(); ++b) {
                    Allele& allele = *b;
                    const string& base = b->currentBase;
                    long double q;
                    if (obs.currentBase == base) {
                        q = qual;
                    } else {
                        q = 1 - qual;
                    }
                    double asampl = genotype.alleleSamplingProb(base);
                    double freq = freqs[base];
                    if (asampl == 0.5) {
                        if (allele.isReference()) {
                            asampl = (contamination.probRefGivenHet + contamination.probRefGivenHomAlt * (1+freq));
                        } else {
                            asampl = 1 - (contamination.probRefGivenHet + contamination.probRefGivenHomAlt * (1+freq));
                        }
                    } else if (asampl == 0) {
                        asampl = contamination.probRefGivenHomAlt * (1+freq); // scale by frequency of (this) possibly contaminating allele
                    } else if (asampl == 1) {
                        asampl = 1 - (contamination.probRefGivenHomAlt * (1+freq)); // scale by frequency of possibly contaminating alleles
                    }

                    probi += asampl * q;

                }

                // bound at 1
                long double lnprobi = log(min(probi, (long double) 1.0));
                //if (countOut > 1) {
                //    lnprobi *= ((1 + (countOut - 1) * dependenceFactor) / countOut);
                //}
                probObsGivenGt += lnprobi;
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
        return isinf(probObsGivenGt) ? 0 : probObsGivenGt;// + multinomialCoefficientLn(observationCount, genotype.counts());
        //return probObsGivenGt;
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
        Contamination& contaminations,
        map<string, double>& freqs
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
                      contaminations,
                      freqs)));
    }
    return results;
}
