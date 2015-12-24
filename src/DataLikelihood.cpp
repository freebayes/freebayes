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

    //cerr << "P(" << genotype << " given" << endl <<  sample;

    int observationCount = sample.observationCount();
    vector<long double> alleleProbs = genotype.alleleProbabilities(observationBias);
    vector<int> observationCounts = genotype.alleleObservationCounts(sample);
    int countOut = 0;
    double countIn = 0;
    long double prodQout = 0;  // the probability that the reads not in the genotype are all wrong
    long double prodSample = 0;
    
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
        vector<Allele*> emptyA;
        vector<Allele*> emptyB;
        for (set<string>::iterator c = sample.supportedAlleles.begin();
             c != sample.supportedAlleles.end(); ++c) {

            vector<Allele*>* alleles = &emptyA;
            Sample::iterator si = sample.find(*c);
            if (si != sample.end()) alleles = &si->second;

            vector<Allele*>* partials = &emptyB;
            map<string, vector<Allele*> >::iterator pi = sample.partialSupport.find(*c);
            if (pi != sample.partialSupport.end()) partials = &pi->second;

            bool onPartials = false;
            vector<Allele*>::iterator a = alleles->begin();
            bool hasPartials = !partials->empty();
            for ( ; (!hasPartials && a != alleles->end()) || a != partials->end(); ++a) {
                if (a == alleles->end()) {
                    if (hasPartials) {
                        a = partials->begin();
                        onPartials = true;
                    } else {
                        break;
                    }
                }
                Allele& obs = **a;
                //cerr << "observation: " << obs << endl;
                long double probi = 0;
                ContaminationEstimate& contamination = contaminations.of(obs.readGroupID);
                double scale = 1;
                // note that this will underflow if we have mapping quality = 0
                // we guard against this externally, by ignoring such alignments (quality has to be > MQL0)
                long double qual = (1.0 - exp(obs.lnquality)) * (1.0 - exp(obs.lnmapQuality));
                
                if (onPartials) {
                    map<Allele*, set<Allele*> >::iterator r = sample.reversePartials.find(*a);
                    if (r != sample.reversePartials.end()) {
                        if (sample.reversePartials[*a].empty()) {
                            cerr << "partial " << *a << " has empty reverse" << endl;
                            exit(1);
                        }
                        //cerr << "partial " << *a << " supports potentially " << sample.reversePartials[*a].size() << " alleles : " << endl;
                        //for (set<Allele*>::iterator m = sample.reversePartials[*a].begin();
                        //m != sample.reversePartials[*a].end(); ++m) cerr << **m << " ";
                        //cerr << endl;
                        scale = (double)1/(double)sample.reversePartials[*a].size();
                        qual *= scale;
                    }
                }

                // TODO add partial obs, now that we have them recorded
                // how does this work?
                // each partial obs is recorded as supporting, but with observation probability scaled by the number of possible haplotypes it supports
                bool isInGenotype = false;
                long double asampl = genotype.alleleSamplingProb(obs);

                // for each of the unique genotype alleles
                for (vector<Allele>::iterator b = genotypeAlleles.begin(); b != genotypeAlleles.end(); ++b) {
                    Allele& allele = *b;
                    const string& base = allele.currentBase;
                    if (genotype.containsAllele(base)
                        && (obs.currentBase == base
                            || (onPartials && sample.observationSupports(*a, &*b)))) {
                        isInGenotype = true;
                        // use the matched allele to estimate the asampl
                        asampl = max(asampl, (long double)genotype.alleleSamplingProb(allele));
                    }
                }

                if (asampl == 0) {
                    // scale by frequency of (this) possibly contaminating allele
                    asampl = contamination.probRefGivenHomAlt;
                } else if (asampl == 1) {
                    // scale by frequency of (other) possibly contaminating alleles
                    asampl = 1 - contamination.probRefGivenHomAlt;
                } else { //if (genotype.ploidy == 2) {
                    // to deal with polyploids
                    // note that this reduces to 1 for diploid heterozygotes
                    // this term captures reference bias
                    if (obs.isReference()) {
                        asampl *= (contamination.probRefGivenHet / 0.5);
                    } else {
                        asampl *= ((1 - contamination.probRefGivenHet) / 0.5);
                    }
                }

                // distribute observation support across haplotypes
                if (!isInGenotype) {
                    prodQout += log(1-qual);
                    countOut += scale;
                } else {
                    prodSample += log(asampl*scale);
                }
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
            //cerr << "P(obs|" << genotype << ") = " << prodQout + multinomialSamplingProbLn(alleleProbs, observationCounts) << endl << endl << string(80, '@') << endl << endl;
            return prodQout + multinomialSamplingProbLn(alleleProbs, observationCounts);
            //return prodQout + samplingProbLn(alleleProbs, observationCounts);
        }
    } else {
        if (countOut > 1) {
            prodQout *= (1 + (countOut - 1) * dependenceFactor) / countOut;
        }
        long double probObsGivenGt = prodQout + prodSample;
        return isinf(probObsGivenGt) ? 0 : probObsGivenGt;
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
        Genotype& genotype = **g;
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

void
calculateSampleDataLikelihoods(
    Samples& samples,
    Results& results,
    AlleleParser* parser,
    map<int, vector<Genotype> >& genotypesByPloidy,
    Parameters& parameters,
    bool usingNull,
    Bias& observationBias,
    vector<Allele>& genotypeAlleles,
    Contamination& contaminationEstimates,
    map<string, double>& estimatedAlleleFrequencies,
    map<string, vector<vector<SampleDataLikelihood> > >& sampleDataLikelihoodsByPopulation,
    map<string, vector<vector<SampleDataLikelihood> > >& variantSampleDataLikelihoodsByPopulation,
    map<string, vector<vector<SampleDataLikelihood> > >& invariantSampleDataLikelihoodsByPopulation) {

    for (vector<string>::iterator n = parser->sampleList.begin(); n != parser->sampleList.end(); ++n) {
        //string sampleName = s->first;
        string& sampleName = *n;
        //DEBUG2("sample: " << sampleName);
        //Sample& sample = s->second;
        if (samples.find(sampleName) == samples.end()
            && !(parser->hasInputVariantAllelesAtCurrentPosition()
                 || parameters.reportMonomorphic)) {
            continue;
        }
        Sample& sample = samples[sampleName];
        vector<Genotype>& genotypes = genotypesByPloidy[parser->currentSamplePloidy(sampleName)];
        vector<Genotype*> genotypesWithObs;
        for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
            if (parameters.excludePartiallyObservedGenotypes) {
                if (g->sampleHasSupportingObservationsForAllAlleles(sample)) {
                    genotypesWithObs.push_back(&*g);
                }
            } else if (parameters.excludeUnobservedGenotypes && usingNull) {
                if (g->sampleHasSupportingObservations(sample)) {
                    //cerr << sampleName << " has suppporting obs for " << *g << endl;
                    genotypesWithObs.push_back(&*g);
                } else if (g->hasNullAllele() && g->homozygous) {
                    // this genotype will never be added if we are running in observed-only mode, but
                    // we still need it for consistency
                    genotypesWithObs.push_back(&*g);
                }
            } else {
                genotypesWithObs.push_back(&*g);
            }
        }

        // skip this sample if we have no observations supporting any of the genotypes we are going to evaluate
        if (genotypesWithObs.empty()) {
            continue;
        }

        vector<pair<Genotype*, long double> > probs
            = probObservedAllelesGivenGenotypes(sample, genotypesWithObs,
                                                parameters.RDF, parameters.useMappingQuality,
                                                observationBias, parameters.standardGLs,
                                                genotypeAlleles,
                                                contaminationEstimates,
                                                estimatedAlleleFrequencies);
            
#ifdef VERBOSE_DEBUG
        if (parameters.debug2) {
            for (vector<pair<Genotype*, long double> >::iterator p = probs.begin(); p != probs.end(); ++p) {
                cerr << parser->currentSequenceName << "," << (long unsigned int) parser->currentPosition + 1 << ","
                     << sampleName << ",likelihood," << *(p->first) << "," << p->second << endl;
            }
        }
#endif

        Result& sampleData = results[sampleName];
        sampleData.name = sampleName;
        sampleData.observations = &sample;
        for (vector<pair<Genotype*, long double> >::iterator p = probs.begin(); p != probs.end(); ++p) {
            sampleData.push_back(SampleDataLikelihood(sampleName, &sample, p->first, p->second, 0));
        }

        sortSampleDataLikelihoods(sampleData);

        string& population = parser->samplePopulation[sampleName];
        vector<vector<SampleDataLikelihood> >& sampleDataLikelihoods = sampleDataLikelihoodsByPopulation[population];
        vector<vector<SampleDataLikelihood> >& variantSampleDataLikelihoods = variantSampleDataLikelihoodsByPopulation[population];
        vector<vector<SampleDataLikelihood> >& invariantSampleDataLikelihoods = invariantSampleDataLikelihoodsByPopulation[population];

        if (parameters.genotypeVariantThreshold != 0) {
            if (sampleData.size() > 1
                && abs(sampleData.at(1).prob - sampleData.front().prob)
                < parameters.genotypeVariantThreshold) {
                variantSampleDataLikelihoods.push_back(sampleData);
            } else {
                invariantSampleDataLikelihoods.push_back(sampleData);
            }
        } else {
            variantSampleDataLikelihoods.push_back(sampleData);
        }
        sampleDataLikelihoods.push_back(sampleData);

    }
}
