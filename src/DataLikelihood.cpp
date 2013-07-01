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
            vector<Allele*>& partials = sample.partials[s->first];  // mm.. cleanup
            /*
            if (!partials.empty()) {
                cerr << "partials " << partials << endl;
            }
            */
            bool onPartials = false;
            vector<Allele*>::iterator a = alleles.begin();
            for ( ; (partials.empty() && a != alleles.end()) || a != partials.end(); ++a) {
                if (a == alleles.end()) {
                    if (!partials.empty()) {
                        a = partials.begin();
                        onPartials = true;
                    } else {
                        break;
                    }
                }
                Allele& obs = **a;
                long double probi = 0;
                ContaminationEstimate& contamination = contaminations.of(obs.readGroupID);
                double scale = 1;
                long double qual = (1 - exp(obs.lnquality)) * (1 - exp(obs.lnmapQuality));
                if (onPartials) {
                    map<Allele*, set<Allele*> >::iterator r = sample.reversePartials.find(*a);
                    if (r != sample.reversePartials.end()) {
                        assert(!sample.reversePartials[*a].empty());
                        scale = (double)1/(double)sample.reversePartials[*a].size();
                    }
                }

                // TODO add partial obs, now that we have them recorded
                // how does this work?
                // each partial obs is recorded as supporting, but with observation probability scaled by the number of possible haplotypes it supports

                for (vector<Allele>::iterator b = genotypeAlleles.begin(); b != genotypeAlleles.end(); ++b) {
                    Allele& allele = *b;
                    const string& base = b->currentBase;

                    long double q;
                    if (obs.currentBase == base
                        || (onPartials && sample.observationSupports(*a, &*b))) {
                        q = qual;
                    } else {
                        q = 1 - qual;
                    }

                    if (onPartials) {
                        q *= scale; // distribute partial support evenly across supported haplotypes
                    }

                    double asampl = genotype.alleleSamplingProb(base);
                    double freq = freqs[base];

                    if (asampl == 0.5) {
                        if (allele.isReference()) {
                            asampl = (contamination.probRefGivenHet + contamination.probRefGivenHomAlt * (0.5+freq));
                        } else {
                            asampl = 1 - (contamination.probRefGivenHet + contamination.probRefGivenHomAlt * (0.5+(1-freq)));
                        }
                    } else if (asampl == 0) {
                        // scale by frequency of (this) possibly contaminating allele
                        asampl = contamination.probRefGivenHomAlt * (0.5+freq);
                    } else if (asampl == 1) {
                        // scale by frequency of (other) possibly contaminating alleles
                        asampl = 1 - (contamination.probRefGivenHomAlt * (0.5+(1-freq)));
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
