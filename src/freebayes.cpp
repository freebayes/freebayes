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
#include "BamReader.h"
#include "Fasta.h"
#include "TryCatch.h"
#include "Parameters.h"
#include "Allele.h"
#include "AlleleParser.h"
#include "Utility.h"

#include "multichoose.h"
#include "multipermute.h"

#include "Genotype.h"
#include "GenotypePriors.h"
#include "DataLikelihood.h"
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

// Allele object recycling:
//
// We use the Allele freelist for performance reasons.  When an Allele object
// is destroyed, it is pushed onto this freelist.  When a new Allele object is
// created, new first checks if we have a free Allele object on the freelist.
// Because we are dynamically linked, we have to declare the freelist here,
// although it exists as a static member of the Allele class.
//
AlleleFreeList Allele::_freeList;


int main (int argc, char *argv[]) {

    AlleleParser* parser = new AlleleParser(argc, argv);
    Parameters& parameters = parser->parameters;
    list<Allele*> alleles;

    map<string, vector<Allele*> > sampleGroups;

    ostream& out = *(parser->output);

    // this can be uncommented to force operation on a specific set of genotypes
    vector<Allele> allGenotypeAlleles;
    allGenotypeAlleles.push_back(genotypeAllele(ALLELE_GENOTYPE, "A", 1));
    allGenotypeAlleles.push_back(genotypeAllele(ALLELE_GENOTYPE, "T", 1));
    allGenotypeAlleles.push_back(genotypeAllele(ALLELE_GENOTYPE, "G", 1));
    allGenotypeAlleles.push_back(genotypeAllele(ALLELE_GENOTYPE, "C", 1));
    //vector<Genotype> genotypes = allPossibleGenotypes(parameters.ploidy, genotypeAlleles);

    int allowedAlleleTypes = ALLELE_REFERENCE | ALLELE_SNP;
    if (parameters.allowIndels) {
        allowedAlleleTypes |= ALLELE_INSERTION;
        allowedAlleleTypes |= ALLELE_DELETION;
    }

    // output VCF header
    if (parameters.output == "vcf") {
        vcfHeader(out, parser->reference->filename, parser->sampleList);
    }

    while (parser->getNextAlleles(sampleGroups, allowedAlleleTypes)) {

        DEBUG2("at start of main loop");

        // don't process reference N's
        if (parser->currentReferenceBaseString() == "N") {
            continue;
        }

        //filterAlleles(alleles, allowedAlleles);
        //removeIndelMaskedAlleles(alleles, parser->currentPosition);
        
        DEBUG2("alleles filtered");

        if (parameters.trace) {
            for (map<string, vector<Allele*> >::iterator g = sampleGroups.begin(); g != sampleGroups.end(); ++g) {
                vector<Allele*>& group = g->second;
                for (vector<Allele*>::iterator a = group.begin(); a != group.end(); ++a) {
                    Allele& allele = **a;
                    parser->traceFile << parser->currentTarget->seq << "," << parser->currentPosition + 1 << ",allele,"
                        << g->first << "," << allele.readID << "," << allele.currentBase << "," 
                        << allele.currentQuality() << "," << allele.mapQuality << endl;
                }
            }
            DEBUG2("after trace generation");
        }

        int coverage = countAlleles(sampleGroups);

        DEBUG("position: " << parser->currentTarget->seq << ":" << parser->currentPosition << " coverage: " << coverage);

        // skips 0-coverage regions
        if (coverage == 0) {
            //cerr << "no alleles found at " << parser->currentTarget->seq << ":" << parser->currentPosition << endl;
            DEBUG("no alleles left at this site after filtering");
            continue;
        }

        DEBUG2("coverage " << parser->currentTarget->seq << ":" << parser->currentPosition << " == " << coverage);

        // establish a set of possible alternate alleles to evaluate at this location
        // only evaluate alleles with at least one supporting read with mapping
        // quality (MQL1) and base quality (BQL1)

        //map<string, vector<Allele*> > sampleGroups = groupAllelesBySample(alleles);
        DEBUG2("grouped alleles by sample");

        if (!sufficientAlternateObservations(sampleGroups, parameters.minAltCount, parameters.minAltFraction)) {
            DEBUG2("insufficient observations at " << parser->currentTarget->seq << ":" << parser->currentPosition);
            continue;
        }

        vector<vector<Allele*> > alleleGroups = groupAlleles(sampleGroups, &allelesEquivalent);
        DEBUG2("grouped alleles by equivalence");

        vector<string> sampleListPlusRef;
        //if (parameters.trace) {
            // figure out which samples have no data so we can print ?'s in the genotype combo trace
        for (vector<string>::iterator s = parser->sampleList.begin(); s != parser->sampleList.end(); ++s) {
            sampleListPlusRef.push_back(*s);
        }
        if (parameters.useRefAllele)
            sampleListPlusRef.push_back(parser->currentTarget->seq);
        //}

        vector<Allele> genotypeAlleles = parser->genotypeAlleles(alleleGroups, sampleGroups, allGenotypeAlleles);

        if (genotypeAlleles.size() <= 1) { // if we have only one viable alternate, we don't have evidence for variation at this site
            DEBUG2("no alternate genotype alleles passed filters at " << parser->currentTarget->seq << ":" << parser->currentPosition);
            continue;
        }
        DEBUG2("genotype alleles: " << genotypeAlleles);

        vector<Genotype> genotypes = allPossibleGenotypes(parameters.ploidy, genotypeAlleles);
        DEBUG2("generated all possible genotypes:");
        if (parameters.debug2) {
            for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
                DEBUG2(*g);
            }
        }

        Results results;

        DEBUG2("calculating data likelihoods");
        // calculate data likelihoods
        for (map<string, vector<Allele*> >::iterator sample = sampleGroups.begin();
                sample != sampleGroups.end(); ++sample) {

            string sampleName = sample->first;
            vector<Allele*>& observedAlleles = sample->second;

            vector<pair<Genotype*, long double> > probs;
            if (parameters.bamBayesDataLikelihoods) {
                 probs = bamBayesProbObservedAllelesGivenGenotypes(observedAlleles, genotypes, parameters.RDF);
            } else {
                 probs = probObservedAllelesGivenGenotypes(observedAlleles, genotypes);
            }

            map<Genotype*, long double> marginals;
            map<Genotype*, vector<long double> > rawMarginals;

            if (parameters.trace) {
                for (vector<pair<Genotype*, long double> >::iterator p = probs.begin(); p != probs.end(); ++p) {
                    parser->traceFile << parser->currentTarget->seq << "," << parser->currentPosition + 1 << ","
                        << sampleName << ",likelihood," << IUPAC2GenotypeStr(IUPAC(*(p->first))) << "," << p->second << endl;
                }
            }

            results.insert(make_pair(sampleName, ResultData(sampleName, probs, marginals, rawMarginals, observedAlleles)));

        }
        
        DEBUG2("finished calculating data likelihoods");

        // sort genotype data likelihoods and accumulate into 
        
        // XXX this section is a hack to make output of trace identical to BamBayes
        // and also outputs the list of samples
        vector<bool> samplesWithData;
        if (parameters.trace) parser->traceFile << parser->currentTarget->seq << "," << parser->currentPosition + 1 << ",samples,";
        for (vector<string>::iterator s = sampleListPlusRef.begin(); s != sampleListPlusRef.end(); ++s) {
            if (parameters.trace) parser->traceFile << *s << ":";
            Results::iterator r = results.find(*s);
            if (r != results.end()) {
                samplesWithData.push_back(true);
            } else {
                samplesWithData.push_back(false);
            }
        }
        if (parameters.trace) parser->traceFile << endl;

        vector<pair<string, vector<pair<Genotype*, long double> > > > sampleGenotypes;
        //for (Results::iterator s = results.begin(); s != results.end(); ++s) {
        for (vector<string>::iterator s = sampleListPlusRef.begin(); s != sampleListPlusRef.end(); ++s) {
            Results::iterator r = results.find(*s);
            if (r != results.end()) {
                r->second.sortDataLikelihoods();
                sampleGenotypes.push_back(make_pair(r->first, r->second.dataLikelihoods));
            }
        }

        DEBUG2("finished sorting data likelihoods");

        // calculate genotype combo likelihoods, integral over nearby genotypes
        // calculate marginals
        // and determine best genotype combination

        // TODO add in all-alternate and all-reference
        // ... or at least very add one which is consistent
        //vector<GenotypeCombo> bandedCombos = bandedGenotypeCombinations(sampleGenotypes, 2, 2);
        //vector<GenotypeCombo> bandedCombos = bandedGenotypeCombinationsIncludingBestHomozygousCombo(sampleGenotypes, 2, 2);
        DEBUG2("generating banded genotype combinations from " << genotypes.size() << " genotypes and " << sampleGenotypes.size() << " sample genotypes");
        vector<GenotypeCombo> bandedCombos = bandedGenotypeCombinationsIncludingAllHomozygousCombos(sampleGenotypes, genotypes, parameters.WB + 1, 2);
        DEBUG2("...done");

        vector<GenotypeComboResult> genotypeComboProbs;

        Allele refAllele = genotypeAllele(ALLELE_REFERENCE, string(1, parser->currentReferenceBase), 1);

        DEBUG2("calculating genotype combination likelihoods");
        for (vector<GenotypeCombo>::iterator combo = bandedCombos.begin(); combo != bandedCombos.end(); ++combo) {

            long double probabilityObservationsGivenGenotypes = 0;
            vector<Genotype*> genotypeCombo;

            for (GenotypeCombo::iterator i = combo->begin(); i != combo->end(); ++i) {
                genotypeCombo.push_back(i->second.first);
                probabilityObservationsGivenGenotypes += i->second.second;
            }

            map<int, int> frequencyCounts = countFrequencies(genotypeCombo);

            long double priorProbabilityOfGenotypeComboG_Af = 
                probabilityDiploidGenotypeComboGivenAlleleFrequencyln(*combo, refAllele);
            long double priorProbabilityOfGenotypeComboAf = 
                alleleFrequencyProbabilityln(countFrequencies(genotypeCombo), parameters.TH);
            long double priorProbabilityOfGenotypeCombo = 
                priorProbabilityOfGenotypeComboG_Af + priorProbabilityOfGenotypeComboAf;
            long double comboProb = priorProbabilityOfGenotypeCombo + probabilityObservationsGivenGenotypes;

            for (GenotypeCombo::iterator i = combo->begin(); i != combo->end(); ++i) {
                map<Genotype*, vector<long double> >& marginals = results[i->first].rawMarginals;
                Genotype* genotype = i->second.first;
                map<Genotype*, vector<long double> >::iterator marginal = marginals.find(genotype);
                if (marginal == marginals.end()) {
                    vector<long double> probs;
                    probs.push_back(comboProb);
                    marginals.insert(make_pair(genotype, probs));
                } else {
                    marginals[genotype].push_back(comboProb);
                }
            }

            genotypeComboProbs.push_back(GenotypeComboResult(*combo,
                        comboProb,
                        probabilityObservationsGivenGenotypes,
                        priorProbabilityOfGenotypeCombo,
                        priorProbabilityOfGenotypeComboG_Af,
                        priorProbabilityOfGenotypeComboAf));

        }

        DEBUG2("finished calculating genotype combination likelihoods");

        // get posterior normalizer
        vector<long double> comboProbs;
        //comboProbs.resize(genotypeComboProbs.size());
        for (vector<GenotypeComboResult>::iterator gc = genotypeComboProbs.begin(); gc != genotypeComboProbs.end(); ++gc)
            comboProbs.push_back(gc->comboProb);
        long double posteriorNormalizer = logsumexp_probs(comboProbs);
        DEBUG2("got posterior normalizer");
        if (parameters.trace) {
            parser->traceFile << parser->currentTarget->seq << "," 
                << parser->currentPosition + 1 << ",posterior_normalizer," << posteriorNormalizer << endl;
        }

        // sort by the normalized datalikelihood + prior
        sort(genotypeComboProbs.begin(), genotypeComboProbs.end(), genotypeComboResultSorter);

        DEBUG2("sorted genotype combination likelihoods");

        /*
        for (vector<pair<GenotypeCombo, long double> >::iterator c = genotypeComboProbs.begin(); c != genotypeComboProbs.end(); ++c) {
            cout << "prob:" << c->second << endl;
        }
        */
        //cout << "posteriorNormalizer = " << posteriorNormalizer << endl;
        /*
        cout << "comboProbs = [";
        for (vector<long double>::iterator i = comboProbs.begin(); i != comboProbs.end(); ++i)
            cout << *i << " ";
        cout << "]" << endl;
        */
        

        // TODO Gabor had a suggestion about this normalization, that we only normalize 'perturbed' cases
        // ....  but i think that this thought was not asserted as necessary
        // normalize marginals
        for (Results::iterator r = results.begin(); r != results.end(); ++r) {
            ResultData& d = r->second;
            //cout << r->first << endl;
            for (map<Genotype*, vector<long double> >::iterator m = d.rawMarginals.begin(); m != d.rawMarginals.end(); ++m) {
                /*
                cout << *m->first << " rawMarginals = [";
                for (vector<long double>::iterator i = m->second.begin(); i != m->second.end(); ++i)
                    cout << *i << ", ";
                cout << "]" << endl;
                cout << logsumexp_probs(m->second) << endl;
                */
                d.marginals[m->first] = logsumexp_probs(m->second) - posteriorNormalizer;
            }
        }
        DEBUG2("calculated marginal likelihoods");

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

        for (vector<GenotypeComboResult>::iterator gc = genotypeComboProbs.begin(); gc != genotypeComboProbs.end(); ++gc) {
            /*
            for (GenotypeCombo::iterator c = gc->first.begin(); c != gc->first.end(); ++c)
                cout << *c->second.first << " ";
            cout << endl;
            cout << (long double) gc->second << endl;
            */
            if (isHomozygousCombo(gc->combo)) {
                pVar -= safe_exp(gc->comboProb - posteriorNormalizer);
            }
        }
        DEBUG2("calculated pVar");

        // TODO as we only ever report heterozygous combos, we must get the
        // best non-homozygous case here even if the best case is fully
        // homozygous

        /*
        GenotypeCombo* besthc = NULL;
        //long double bestHetGenotypeComboProb;
        // these are sorted according to probability, so we just step through until we find the first het
        for (vector<GenotypeComboResult>::iterator gc = genotypeComboProbs.begin(); gc != genotypeComboProbs.end(); ++gc) {
            if (!isHomozygousCombo(gc->get<0>())) {
                besthc = &(gc->get<0>());
                //bestHetGenotypeComboProb = safe_exp(gc->get<1>() - posteriorNormalizer);
                break;
            }
        }
        if (besthc == NULL) {
            ERROR("could not find a heterozygous genotype combo at " << parser->currentTarget->seq << ":" << parser->currentPosition);
            continue; // weird situation which suggests an implementation bug, just raise error and continue
        }
        GenotypeCombo& bestHetGenotypeCombo = *besthc; // for establishing the best alternate
        */
        GenotypeCombo& bestGenotypeCombo = genotypeComboProbs.front().combo; //*besthc;
        long double bestGenotypeComboProb = genotypeComboProbs.front().comboProb;
        vector<Genotype*> bestComboGenotypes;
        for (GenotypeCombo::iterator g = bestGenotypeCombo.begin(); g != bestGenotypeCombo.end(); ++g)
            bestComboGenotypes.push_back(g->second.first);
        long double bestGenotypeComboAlleleSamplingProb = safe_exp(alleleFrequencyProbabilityln(countFrequencies(bestComboGenotypes), parameters.TH));

        if (parameters.trace) {
            for (vector<GenotypeComboResult>::iterator gc = genotypeComboProbs.begin(); gc != genotypeComboProbs.end(); ++gc) {
                vector<Genotype*> comboGenotypes;
                for (GenotypeCombo::iterator g = gc->combo.begin(); g != gc->combo.end(); ++g)
                    comboGenotypes.push_back(g->second.first);
                long double comboProb = gc->comboProb;
                long double dataLikelihoodln = gc->probObsGivenGenotypes;
                long double priorln = gc->priorProbGenotypeCombo;
                long double priorlnG_Af = gc->priorProbGenotypeComboG_Af;
                long double priorlnAf = gc->priorProbGenotypeComboAf;

                parser->traceFile << parser->currentTarget->seq << "," << parser->currentPosition + 1 << ",genotypecombo,";

                int j = 0;
                GenotypeCombo::iterator i = gc->combo.begin();
                for (vector<bool>::iterator d = samplesWithData.begin(); d != samplesWithData.end(); ++d) {
                    if (*d) {
                        parser->traceFile << IUPAC(*i->second.first);
                        ++i;
                    } else {
                        parser->traceFile << "?";
                    }
                }
                    //<< "," << gc->first
                parser->traceFile 
                    << "," << dataLikelihoodln
                    << "," << priorln
                    << "," << priorlnG_Af
                    << "," << priorlnAf
                    << "," << dataLikelihoodln + priorln
                    << "," << safe_exp(dataLikelihoodln + priorln - posteriorNormalizer)
                    << endl;
            }
        }

        DEBUG2("got bestAlleleSamplingProb");
        DEBUG2("pVar = " << pVar << " " << parameters.PVL);

        if (!parameters.suppressOutput) {
            //cerr << parser->currentPosition << " " << alleles.size() << " " << bestGenotypeComboProb << " " << genotypeComboProbs.front().second << " " <<  posteriorNormalizer << endl;

            if (parameters.output == "json") {
                out << "{ \"position\": " << parser->currentPosition + 1 // 1-based reporting, to match vcf
                    << ", \"sequence\": " << parser->currentTarget->seq
                    << ", \"best_genotype_combo\":" << bestGenotypeCombo
                    << ", \"best_genotype_combo_prob\":" << bestGenotypeComboProb 
                    << ", \"best_genotype_combo_ewens_sampling_probability\":" << bestGenotypeComboAlleleSamplingProb
                    << ", \"combos_tested\":" << bandedCombos.size()
                    << ", \"coverage\":" << coverage
                    << ", \"posterior_normalizer\":" << safe_exp(posteriorNormalizer)
                    << ", \"samples\":";
                json(out, results, parser);
                out << "}" << endl;
            }
            if (pVar >= parameters.PVL && !isHomozygousCombo(bestGenotypeCombo)) {
                if (parameters.output == "vcf") {
                    string referenceBase(1, parser->currentReferenceBase);
                    // get the unique alternate alleles in this combo, sorted by frequency in the combo
                    vector<pair<Allele, int> > alternates = alternateAlleles(bestGenotypeCombo, referenceBase);
                    Allele& bestAlt = alternates.front().first;
                    // TODO update the vcf output function to handle the reporting of multiple alternate alleles
                    out << vcf(pVar,
                            sampleGroups,
                            referenceBase,
                            bestAlt.currentBase,
                            parser->sampleList,
                            coverage,
                            bestGenotypeCombo,
                            results,
                            parser)
                        << endl;
                }
            }
        }
        DEBUG2("finished position");

    }

    delete parser;

    return 0;

}
