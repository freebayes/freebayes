// 
// freebayes
//
// A bayesian genetic variant detector.
// 

// standard includes
//#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <float.h>

// private libraries
#include "api/BamReader.h"
#include "Fasta.h"
#include "TryCatch.h"
#include "Parameters.h"
#include "Allele.h"
#include "Sample.h"
#include "AlleleParser.h"
#include "Utility.h"

#include "multichoose.h"
#include "multipermute.h"

#include "Genotype.h"
#include "DataLikelihood.h"
#include "Marginals.h"
#include "ResultData.h"


// local helper debugging macros to improve code readability
#define DEBUG(msg) \
    if (parameters.debug) { cerr << msg << endl; }

// lower-priority messages
#ifdef VERBOSE_DEBUG
#define DEBUG2(msg) \
    if (parameters.debug2) { cerr << msg << endl; }
#else
#define DEBUG2(msg)
#endif

// must-see error messages
#define ERROR(msg) \
    cerr << msg << endl;


using namespace std; 


// freebayes main
int main (int argc, char *argv[]) {

    AlleleParser* parser = new AlleleParser(argc, argv);
    Parameters& parameters = parser->parameters;
    list<Allele*> alleles;

    Samples samples;

    ostream& out = *(parser->output);

    // this can be uncommented to force operation on a specific set of genotypes
    vector<Allele> allGenotypeAlleles;
    allGenotypeAlleles.push_back(genotypeAllele(ALLELE_GENOTYPE, "A", 1));
    allGenotypeAlleles.push_back(genotypeAllele(ALLELE_GENOTYPE, "T", 1));
    allGenotypeAlleles.push_back(genotypeAllele(ALLELE_GENOTYPE, "G", 1));
    allGenotypeAlleles.push_back(genotypeAllele(ALLELE_GENOTYPE, "C", 1));

    int allowedAlleleTypes = ALLELE_REFERENCE;
    if (parameters.allowSNPs) {
        allowedAlleleTypes |= ALLELE_SNP;
    }
    if (parameters.allowIndels) {
        allowedAlleleTypes |= ALLELE_INSERTION;
        allowedAlleleTypes |= ALLELE_DELETION;
    }
    if (parameters.allowMNPs) {
        allowedAlleleTypes |= ALLELE_MNP;
    }
    if (parameters.allowComplex) {
        allowedAlleleTypes |= ALLELE_COMPLEX;
    }

    // output VCF header
    if (parameters.output == "vcf") {
        out << parser->variantCallFile.header;
    }

    Allele nullAllele = genotypeAllele(ALLELE_NULL, "N", 1, "1N");

    unsigned long total_sites = 0;
    unsigned long processed_sites = 0;

    while (parser->getNextAlleles(samples, allowedAlleleTypes)) {

        ++total_sites;

        DEBUG2("at start of main loop");

        // don't process non-ATGC's in the reference
        string cb = parser->currentReferenceBaseString();
        if (cb != "A" && cb != "T" && cb != "C" && cb != "G") {
            DEBUG2("current reference base is N");
            continue;
        }

        if (parameters.trace) {
            for (Samples::iterator s = samples.begin(); s != samples.end(); ++s) {
                const string& name = s->first;
                for (Sample::iterator g = s->second.begin(); g != s->second.end(); ++g) {
                    vector<Allele*>& group = g->second;
                    for (vector<Allele*>::iterator a = group.begin(); a != group.end(); ++a) {
                        Allele& allele = **a;
                        parser->traceFile << parser->currentSequenceName << "," << (long unsigned int) parser->currentPosition + 1  
                            << ",allele," << name << "," << allele.readID << "," << allele.base() << ","
                            << allele.currentQuality() << "," << allele.mapQuality << endl;
                    }
                }
            }
            DEBUG2("after trace generation");
        }

        if (!parser->inTarget()) {
            DEBUG("position: " << parser->currentSequenceName << ":" << (long unsigned int) parser->currentPosition + 1
                    << " is not inside any targets, skipping");
            continue;
        }

        int coverage = countAlleles(samples);

        DEBUG("position: " << parser->currentSequenceName << ":" << (long unsigned int) parser->currentPosition + 1 << " coverage: " << coverage);

        if (!parser->hasInputVariantAllelesAtCurrentPosition()) {
            // skips 0-coverage regions
            if (coverage == 0) {
                DEBUG("no alleles left at this site after filtering");
                continue;
            } else if (coverage < parameters.minCoverage) {
                DEBUG("post-filtering coverage of " << coverage << " is less than --min-coverage of " << parameters.minCoverage);
                continue;
            }

            DEBUG2("coverage " << parser->currentSequenceName << ":" << parser->currentPosition << " == " << coverage);

            // establish a set of possible alternate alleles to evaluate at this location
            // only evaluate alleles with at least one supporting read with mapping
            // quality (MQL1) and base quality (BQL1)

            if (!sufficientAlternateObservations(samples, parameters.minAltCount, parameters.minAltFraction)) {
                DEBUG("insufficient alternate observations");
                continue;
            }
        }

        // to ensure proper ordering of output stream
        vector<string> sampleListPlusRef;

        for (vector<string>::iterator s = parser->sampleList.begin(); s != parser->sampleList.end(); ++s) {
            sampleListPlusRef.push_back(*s);
        }
        if (parameters.useRefAllele)
            sampleListPlusRef.push_back(parser->currentSequenceName);

        // establish genotype alleles using input filters
        map<string, vector<Allele*> > alleleGroups;
        groupAlleles(samples, alleleGroups);
        DEBUG2("grouped alleles by equivalence");

        vector<Allele> genotypeAlleles = parser->genotypeAlleles(alleleGroups, samples, parameters.onlyUseInputAlleles);

        // always include the reference allele as a possible genotype, even when we don't include it by default
        if (!parameters.useRefAllele) {
            vector<Allele> refAlleleVector;
            refAlleleVector.push_back(genotypeAllele(ALLELE_REFERENCE, string(1, parser->currentReferenceBase), 1, "1M"));
            genotypeAlleles = alleleUnion(genotypeAlleles, refAlleleVector);
        }

        // build haplotype alleles matching the current longest allele (often will do nothing)
        // this will adjust genotypeAlleles if changes are made
        DEBUG("building haplotype alleles, currently there are " << genotypeAlleles.size() << " genotype alleles");
        DEBUG(genotypeAlleles);
        parser->buildHaplotypeAlleles(genotypeAlleles, samples, alleleGroups, allowedAlleleTypes);
        DEBUG("built haplotype alleles, now there are " << genotypeAlleles.size() << " genotype alleles");
        DEBUG(genotypeAlleles);

        if (genotypeAlleles.size() <= 1) { // if we have only one viable allele, we don't have evidence for variation at this site
            DEBUG2("no alternate genotype alleles passed filters at " << parser->currentSequenceName << ":" << parser->currentPosition);
            continue;
        }
        DEBUG("genotype alleles: " << genotypeAlleles);

        // add the null genotype
        bool usingNull = false;
        if (parameters.excludeUnobservedGenotypes && genotypeAlleles.size() > 2) {
            genotypeAlleles.push_back(nullAllele);
            usingNull = true;
        }

        ++processed_sites;

        // generate possible genotypes

        // for each possible ploidy in the dataset, generate all possible genotypes
        vector<int> ploidies = parser->currentPloidies(samples);
        map<int, vector<Genotype> > genotypesByPloidy = getGenotypesByPloidy(ploidies, genotypeAlleles);


        DEBUG2("generated all possible genotypes:");
        if (parameters.debug2) {
            for (map<int, vector<Genotype> >::iterator s = genotypesByPloidy.begin(); s != genotypesByPloidy.end(); ++s) {
                vector<Genotype>& genotypes = s->second;
                for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
                    DEBUG2(*g);
                }
            }
        }

        Results results;
        vector<vector<SampleDataLikelihood> > sampleDataLikelihoods;
        vector<vector<SampleDataLikelihood> > variantSampleDataLikelihoods;
        vector<vector<SampleDataLikelihood> > invariantSampleDataLikelihoods;

        // TODO... segregate samples from the VCF from samples from BAMs

        DEBUG2("calculating data likelihoods");
        // calculate data likelihoods
        for (Samples::iterator s = samples.begin(); s != samples.end(); ++s) {

            string sampleName = s->first;
            //DEBUG2("sample: " << sampleName);
            Sample& sample = s->second;
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

            vector<pair<Genotype*, long double> > probs = probObservedAllelesGivenGenotypes(sample, genotypesWithObs, parameters.RDF, parameters.useMappingQuality);

            if (parameters.trace) {
                for (vector<pair<Genotype*, long double> >::iterator p = probs.begin(); p != probs.end(); ++p) {
                    parser->traceFile << parser->currentSequenceName << "," << (long unsigned int) parser->currentPosition + 1 << ","
                        << sampleName << ",likelihood," << *(p->first) << "," << p->second << endl;
                }
            }

            Result& sampleData = results[sampleName];
            sampleData.name = sampleName;
            sampleData.observations = &sample;
            for (vector<pair<Genotype*, long double> >::iterator p = probs.begin(); p != probs.end(); ++p) {
                sampleData.push_back(SampleDataLikelihood(sampleName, &sample, p->first, p->second, 0));
            }
            sortSampleDataLikelihoods(sampleData);

            //cout << exp(thisSampleDataLikelihoods.front().prob) << " - " << exp(thisSampleDataLikelihoods.at(1).prob) << endl;
            if (parameters.genotypeVariantThreshold != 0) {
                if (sampleData.size() > 1
                        && float2phred(1 - (exp(sampleData.front().prob) - exp(sampleData.at(1).prob)))
                            < parameters.genotypeVariantThreshold) {
                    //cout << "varying sample " << name << endl;
                    variantSampleDataLikelihoods.push_back(sampleData);
                } else {
                    invariantSampleDataLikelihoods.push_back(sampleData);
                }
            } else {
                variantSampleDataLikelihoods.push_back(sampleData);
            }
            sampleDataLikelihoods.push_back(sampleData);

        }
        
        DEBUG2("finished calculating data likelihoods");

        DEBUG2("obtaining genotype likelihoods input from VCF");

        // TODO, question, do these need to change wrt. NULL genotypes?
        parser->addCurrentGenotypeLikelihoods(genotypesByPloidy, sampleDataLikelihoods);
        // add these sample data likelihoods to 'invariant' likelihoods
        parser->addCurrentGenotypeLikelihoods(genotypesByPloidy, invariantSampleDataLikelihoods);

        DEBUG2("finished parsing VCF input genotype likelihoods");


        // this section is a hack to make output of trace identical to BamBayes trace
        // and also outputs the list of samples
        vector<bool> samplesWithData;
        if (parameters.trace) {
            parser->traceFile << parser->currentSequenceName << "," << (long unsigned int) parser->currentPosition + 1 << ",samples,";
            for (vector<string>::iterator s = sampleListPlusRef.begin(); s != sampleListPlusRef.end(); ++s) {
                if (parameters.trace) parser->traceFile << *s << ":";
                Results::iterator r = results.find(*s);
                if (r != results.end()) {
                    samplesWithData.push_back(true);
                } else {
                    samplesWithData.push_back(false);
                }
            }
            parser->traceFile << endl;
        }


        // calculate genotype combo likelihoods, integral over nearby genotypes
        // calculate marginals
        // and determine best genotype combination

        DEBUG2("generating banded genotype combinations from " << sampleDataLikelihoods.size() << " sample genotypes");
        list<GenotypeCombo> genotypeCombos;

        if (parameters.siteSelectionMaxIterations > 0 && !variantSampleDataLikelihoods.empty()) {
            GenotypeCombo nullCombo;
            bool addHomozygousCombos = true;
            //cerr << sampleDataLikelihoods.size() << endl;
            //cerr << variantSampleDataLikelihoods.size() << endl;
            //cerr << invariantSampleDataLikelihoods.size() << endl;
            convergentGenotypeComboSearch(
                    genotypeCombos,
                    nullCombo,  // passing an empty combo triggers use of the data likelihood max combo
                    sampleDataLikelihoods,
                    variantSampleDataLikelihoods,
                    invariantSampleDataLikelihoods,
                    samples,
                    genotypeAlleles,
                    parameters.WB,
                    parameters.TB,
                    parameters.genotypeComboStepMax,
                    parameters.TH,
                    parameters.pooled,
                    parameters.ewensPriors,
                    parameters.permute,
                    parameters.hwePriors,
                    parameters.obsBinomialPriors,
                    parameters.alleleBalancePriors,
                    parameters.diffusionPriorScalar,
                    parameters.siteSelectionMaxIterations,
                    addHomozygousCombos);
        } else {
            // handles the case where all samples have highly differentiated GLs
            if (variantSampleDataLikelihoods.empty()) {
                variantSampleDataLikelihoods = invariantSampleDataLikelihoods;
                invariantSampleDataLikelihoods.clear();
            }
            DEBUG2("generating banded genotype combinations");
            GenotypeCombo nullCombo;
            bandedGenotypeCombinationsIncludingAllHomozygousCombos(
                    genotypeCombos,
                    nullCombo, // use data likelihood max
                    sampleDataLikelihoods,
                    variantSampleDataLikelihoods,
                    invariantSampleDataLikelihoods,
                    samples,
                    genotypeAlleles,
                    parameters.WB,
                    parameters.TB,
                    parameters.genotypeComboStepMax,
                    parameters.TH,
                    parameters.pooled,
                    parameters.ewensPriors,
                    parameters.permute,
                    parameters.hwePriors,
                    parameters.obsBinomialPriors,
                    parameters.alleleBalancePriors,
                    parameters.diffusionPriorScalar);
        }

        // sort by the normalized datalikelihood + prior
        DEBUG2("sorting genotype combination likelihoods");
        GenotypeComboResultSorter gcrSorter;
        genotypeCombos.sort(gcrSorter);
        genotypeCombos.unique();

        // get posterior normalizer
        vector<long double> comboProbs;
        for (list<GenotypeCombo>::iterator gc = genotypeCombos.begin(); gc != genotypeCombos.end(); ++gc) {
            comboProbs.push_back(gc->posteriorProb);
        }
        long double posteriorNormalizer = logsumexp_probs(comboProbs);

        DEBUG2("got posterior normalizer");
        if (parameters.trace) {
            parser->traceFile << parser->currentSequenceName << "," 
                << (long unsigned int) parser->currentPosition + 1 << ",posterior_normalizer," << posteriorNormalizer << endl;
        }

        // we provide p(var|data), or the probability that the location has
        // variation between individuals relative to the probability that it
        // has no variation
        //
        // in other words:
        // p(var|d) = 1 - p(AA|d) - p(TT|d) - P(GG|d) - P(CC|d)
        //
        // the approach is go through all the homozygous combos
        // and then subtract this from 1... resolving p(var|d)

        long double pVar = 1.0;
        long double pHom = 0.0;

        long double bestComboOddsRatio = 0;

        bool hasHetCombo = false;
        bool bestOverallComboIsHet = false;
        GenotypeCombo bestCombo; // = NULL;

        string referenceBase = parser->currentReferenceBaseString();

        // calculates pvar and gets the best het combo
        for (list<GenotypeCombo>::iterator gc = genotypeCombos.begin(); gc != genotypeCombos.end(); ++gc) {
            if (gc->isHomozygous()
                    && (parameters.useRefAllele
                        || !parameters.useRefAllele && gc->alleles().front() == referenceBase)) {
                pVar -= safe_exp(gc->posteriorProb - posteriorNormalizer);
                pHom += safe_exp(gc->posteriorProb - posteriorNormalizer);
            } else if (!hasHetCombo) { // get the first het combo
                bestCombo = *gc;
                hasHetCombo = true;
                if (gc == genotypeCombos.begin()) {
                    bestOverallComboIsHet = true;
                }
            }
        }

        // if for some reason there are no het combos, use the first combo
        if (!hasHetCombo) {
            bestCombo = genotypeCombos.front();
        }
        DEBUG2("best combo: " << bestCombo);

        // odds ratio between the first and second-best combinations
        if (genotypeCombos.size() > 1) {
            bestComboOddsRatio = genotypeCombos.front().posteriorProb - (++genotypeCombos.begin())->posteriorProb;
        }

        DEBUG2("calculated pVar");

        if (parameters.trace) {
            for (list<GenotypeCombo>::iterator gc = genotypeCombos.begin(); gc != genotypeCombos.end(); ++gc) {
                vector<Genotype*> comboGenotypes;
                for (GenotypeCombo::iterator g = gc->begin(); g != gc->end(); ++g)
                    comboGenotypes.push_back((*g)->genotype);
                long double posteriorProb = gc->posteriorProb;
                long double dataLikelihoodln = gc->probObsGivenGenotypes;
                long double priorln = gc->posteriorProb;
                long double priorlnG_Af = gc->priorProbG_Af;
                long double priorlnAf = gc->priorProbAf;
                long double priorlnBin = gc->priorProbObservations;

                parser->traceFile << parser->currentSequenceName << "," << (long unsigned int) parser->currentPosition + 1 << ",genotypecombo,";

                int j = 0;
                GenotypeCombo::iterator i = gc->begin();
                for (vector<bool>::iterator d = samplesWithData.begin(); d != samplesWithData.end(); ++d) {
                    if (*d) {
                        parser->traceFile << IUPAC(*(*i)->genotype);
                        ++i;
                    } else {
                        parser->traceFile << "?";
                    }
                }
                // TODO cleanup this and above
                parser->traceFile 
                    << "," << dataLikelihoodln
                    << "," << priorln
                    << "," << priorlnG_Af
                    << "," << priorlnAf
                    << "," << priorlnBin
                    << "," << posteriorProb
                    << "," << safe_exp(posteriorProb - posteriorNormalizer)
                    << endl;
            }
        }

        DEBUG2("got bestAlleleSamplingProb");
        DEBUG("pVar = " << pVar << " " << parameters.PVL
              << " pHom = " << pHom
              << " 1 - pHom = " << 1 - pHom);

        DEBUG2("1 - " << pHom << " >= " << parameters.PVL);
        // the second clause guards against float underflow causing us not to output a position
        // practically, parameters.PVL == 0 means "report all genotypes which pass our input filters"
        if ((1 - pHom) >= parameters.PVL || parameters.PVL == 0) {
            DEBUG2("passed PVL threshold");

            GenotypeCombo bestGenotypeComboByMarginals;

            if (parameters.calculateMarginals) {

                DEBUG("calculating marginal likelihoods");

                // resample the posterior, this time without bounds on the
                // samples we vary, ensuring that we can generate marginals for
                // all sample/genotype combinations

                //SampleDataLikelihoods marginalLikelihoods = sampleDataLikelihoods;  // heavyweight copy...
                GenotypeCombo nullCombo;

                SampleDataLikelihoods nullSampleDataLikelihoods;

                genotypeCombos.clear();

                // cap the number of iterations at 2 x the number of alternate alleles
                // max it at parameters.genotypingMaxIterations (200 by default) iterations, min at 10
                int itermax = min(max(10, 2 * bestCombo.hetCount()), parameters.genotypingMaxIterations);
                //cerr << "itermax: " << itermax << endl;

                // search much longer for convergence
                convergentGenotypeComboSearch(
                        genotypeCombos,
                        nullCombo,
                        sampleDataLikelihoods, // vary everything
                        sampleDataLikelihoods,
                        nullSampleDataLikelihoods,
                        samples,
                        genotypeAlleles,
                        // XXX HACK
                        0, 0, // passing 0 for bandwidth and banddepth means "exhaustive local search"
                        // this produces properly normalized GQ's at polyallelic sites
                        //parameters.WB,
                        //parameters.TB,
                        parameters.genotypeComboStepMax,
                        parameters.TH,
                        parameters.pooled,
                        parameters.ewensPriors,
                        parameters.permute,
                        parameters.hwePriors,
                        parameters.obsBinomialPriors,
                        parameters.alleleBalancePriors,
                        parameters.diffusionPriorScalar,
                        itermax,
                        true); // add homozygous combos
                    // ^^ combo results are sorted by default

                // TODO factor out the following blocks as they are repeated from above

                // re-get posterior normalizer
                vector<long double> comboProbs;
                for (list<GenotypeCombo>::iterator gc = genotypeCombos.begin(); gc != genotypeCombos.end(); ++gc) {
                    comboProbs.push_back(gc->posteriorProb);
                }
                posteriorNormalizer = logsumexp_probs(comboProbs);

                // recalculate posterior normalizer
                pVar = 1.0;
                pHom = 0.0;
                hasHetCombo = false;
                // calculates pvar and gets the best het combo
                for (list<GenotypeCombo>::iterator gc = genotypeCombos.begin(); gc != genotypeCombos.end(); ++gc) {
                    if (gc->isHomozygous()
                            && (parameters.useRefAllele
                                || !parameters.useRefAllele && gc->alleles().front() == referenceBase)) {
                        pVar -= safe_exp(gc->posteriorProb - posteriorNormalizer);
                        pHom += safe_exp(gc->posteriorProb - posteriorNormalizer);
                    } else if (!hasHetCombo) { // get the first het combo
                        bestCombo = *gc;
                        hasHetCombo = true;
                        if (gc == genotypeCombos.begin()) {
                            bestOverallComboIsHet = true;
                        }
                    }
                }

                // if for some reason there are no het combos, use the first combo
                if (!hasHetCombo) {
                    bestCombo = genotypeCombos.front();
                }

                DEBUG2("best combo: " << bestCombo);

                // odds ratio between the first and second-best combinations
                if (genotypeCombos.size() > 1) {
                    bestComboOddsRatio = genotypeCombos.front().posteriorProb - (++genotypeCombos.begin())->posteriorProb;
                }

                marginalGenotypeLikelihoods(genotypeCombos, sampleDataLikelihoods);

                // generate the best marginal combo according to marginals, which we've sorted by
                /*
                bestCombo.clear();
                dataLikelihoodMaxGenotypeCombo(
                        bestCombo,
                        sampleDataLikelihoods,
                        parameters.TH,
                        parameters.pooled,
                        parameters.ewensPriors,
                        parameters.permute,
                        parameters.hwePriors,
                        parameters.obsBinomialPriors,
                        parameters.alleleBalancePriors,
                        parameters.diffusionPriorScalar);
                        */

                // store the marginal data likelihoods in the results, for easy parsing
                // like a vector -> map conversion...
                results.update(sampleDataLikelihoods);

            }

            map<string, int> repeats;
            if (parameters.showReferenceRepeats) {
                repeats = parser->repeatCounts(12);
            }

            vector<Allele> alts;
            if (parameters.onlyUseInputAlleles) {
                //alts = genotypeAlleles;
                for (vector<Allele>::iterator a = genotypeAlleles.begin(); a != genotypeAlleles.end(); ++a) {
                    if (!a->isReference()) {
                        alts.push_back(*a);
                    }
                }
            } else {
                // get the unique alternate alleles in this combo, sorted by frequency in the combo
                vector<pair<Allele, int> > alternates = alternateAlleles(bestCombo, referenceBase);
                for (vector<pair<Allele, int> >::iterator a = alternates.begin(); a != alternates.end(); ++a) {
                    Allele& alt = a->first;
                    if (!alt.isNull())
                        alts.push_back(alt);
                }
                // if there are no alternate alleles in the best combo, use the genotype alleles
                // XXX ...
                if (alts.empty()) {
                    for (vector<Allele>::iterator a = genotypeAlleles.begin(); a != genotypeAlleles.end(); ++a) {
                        if (!a->isReference()) {
                            alts.push_back(*a);
                        }
                    }
                }
            }

            vcf::Variant var(parser->variantCallFile);

            out << results.vcf(
                    var,
                    pHom,
                    bestComboOddsRatio,
                    samples,
                    referenceBase,
                    alts,
                    repeats,
                    parser->sampleList,
                    coverage,
                    bestCombo,
                    bestOverallComboIsHet,
                    alleleGroups,
                    genotypesByPloidy,
                    parser->sequencingTechnologies,
                    parser)
                << endl;

        } else if (!parameters.failedFile.empty()) {
            // get the unique alternate alleles in this combo, sorted by frequency in the combo
            long unsigned int position = parser->currentPosition;
            for (vector<Allele>::iterator ga =  genotypeAlleles.begin(); ga != genotypeAlleles.end(); ++ga) {
                if (ga->type == ALLELE_REFERENCE)
                    continue;
                parser->failedFile
                    << parser->currentSequenceName << "\t"
                    << position << "\t"
                    << position + ga->length << "\t"
                    << *ga << endl;
            }
            // BED format
        }
        DEBUG2("finished position");

    }

    DEBUG("total sites: " << total_sites << endl
         << "processed sites: " << processed_sites << endl
         << "ratio: " << (float) processed_sites / (float) total_sites);

    delete parser;

    return 0;

}
