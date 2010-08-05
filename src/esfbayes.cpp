// 
// *Bayes
//
// A bayesian genetic variant parser.
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

// "boost" regular expression library
#include <boost/regex.hpp>

// "boost" string manipulation
#include <boost/algorithm/string/join.hpp>
// tuple
#include <boost/tuple/tuple.hpp>
// bind
#include <boost/bind.hpp>

// private libraries
#include "BamReader.h"
#include "Fasta.h"
#include "TryCatch.h"
#include "Parameters.h"
#include "Allele.h"
#include "AlleleParser.h"

#include "multichoose.h"
#include "multipermute.h"

#include "Genotype.h"
#include "GenotypePriors.h"
#include "DataLikelihood.h"
#include "ResultData.h"

#include <boost/foreach.hpp>

/*
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
*/



using namespace std; 

using boost::tuple;
using boost::make_tuple;

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

    // for now, only estimate probabilities for these genotypes
    vector<Allele> genotypeAlleles;
    genotypeAlleles.push_back(genotypeAllele(ALLELE_GENOTYPE, "A", 1));
    genotypeAlleles.push_back(genotypeAllele(ALLELE_GENOTYPE, "T", 1));
    genotypeAlleles.push_back(genotypeAllele(ALLELE_GENOTYPE, "G", 1));
    genotypeAlleles.push_back(genotypeAllele(ALLELE_GENOTYPE, "C", 1));
    vector<Genotype> genotypes = allPossibleGenotypes(2, genotypeAlleles); // generate all possible genotypes of ploidy 2


    while (parser->getNextAlleles(alleles)) {

        // skips 0-coverage regions
        if (alleles.size() == 0)
            continue;

        // TODO force calculation for samples not in this list
        map<string, vector<Allele*> > sampleGroups = groupAllelesBySample(alleles);

        Results results;

        // calculate data likelihoods
        for (map<string, vector< Allele* > >::iterator sample = sampleGroups.begin();
                sample != sampleGroups.end(); ++sample) {

            string sampleName = sample->first;
            vector<Allele*>& observedAlleles = sample->second;

            vector<pair<Genotype, long double> > probs = 
                probObservedAllelesGivenGenotypes(observedAlleles, genotypes);
            /*
            for (vector<pair<Genotype, long double> >::iterator p = probs.begin(); p != probs.end(); ++p) {
                cout << p->first << " " << p->second << endl;
            }
            */
            
            map<Genotype, long double> marginals;
            map<Genotype, vector<long double> > rawMarginals;

            results.insert(make_pair(sampleName, ResultData(sampleName, probs, marginals, rawMarginals, observedAlleles)));

        }

        // sort genotype data likelihoods and accumulate into 
        
        vector<pair<string, vector<pair<Genotype, long double> > > > sampleGenotypes;
        for (Results::iterator s = results.begin(); s != results.end(); ++s) {
            s->second.sortDataLikelihoods();
            /*
            for (vector<pair<Genotype, long double> >::iterator p = s->second.dataLikelihoods.begin(); p != s->second.dataLikelihoods.end(); ++p) {
                cout << p->first << " " << p->second << endl;
            }
            cout << endl;
            */

            sampleGenotypes.push_back(make_pair(s->first, s->second.dataLikelihoods));
        }


        // calculate genotype combo likelihoods, integral over nearby genotypes
        // calculate marginals
        // and determine best genotype combination

        vector<GenotypeCombo> bandedCombos = bandedGenotypeCombinations(sampleGenotypes, 2, 2);
        vector<pair<GenotypeCombo, long double> > genotypeComboProbs;

        for (vector<GenotypeCombo>::iterator combo = bandedCombos.begin(); combo != bandedCombos.end(); ++combo) {

            //cout << *combo << endl;

            long double probabilityObservationsGivenGenotypes = 0;
            vector<Genotype> genotypeCombo;

            for (GenotypeCombo::iterator i = combo->begin(); i != combo->end(); ++i) {
                genotypeCombo.push_back(i->second.first);
                probabilityObservationsGivenGenotypes += i->second.second;
            }

            long double priorProbabilityOfGenotypeCombo = alleleFrequencyProbabilityln(countFrequencies(genotypeCombo), parameters.TH);
            //cout << "priorProbabilityOfGenotypeCombo = " << priorProbabilityOfGenotypeCombo << endl;
            //cout << "probabilityObservationsGivenGenotypes = " << probabilityObservationsGivenGenotypes << endl;
            long double comboProb = priorProbabilityOfGenotypeCombo + probabilityObservationsGivenGenotypes;

            for (GenotypeCombo::iterator i = combo->begin(); i != combo->end(); ++i) {
                map<Genotype, vector<long double> >& marginals = results[i->first].rawMarginals;
                Genotype& genotype = i->second.first;
                //long double& prob = i->second.second;
                map<Genotype, vector<long double> >::iterator marginal = marginals.find(genotype);
                if (marginal == marginals.end()) {
                    vector<long double> probs;
                    probs.push_back(comboProb);
                    marginals.insert(make_pair(genotype, probs));
                } else {
                    marginals[genotype].push_back(comboProb);
                }
            }

            genotypeComboProbs.push_back(make_pair(*combo, comboProb));

        }
        // genotype_combo_probs = sorted(genotype_combo_probs, key=lambda c: c[1], reverse=True)
        sort(genotypeComboProbs.begin(), genotypeComboProbs.end(),
                boost::bind(&pair<GenotypeCombo, long double>::second, _1) 
                    > boost::bind(&pair<GenotypeCombo, long double>::second, _2));

        vector<long double> comboProbs;
        comboProbs.resize(genotypeComboProbs.size());
        transform(genotypeComboProbs.begin(), genotypeComboProbs.end(),
                comboProbs.begin(), boost::bind(&pair<GenotypeCombo, long double>::second, _1));

        long double posteriorNormalizer = logsumexp(comboProbs);
        //cout << "posteriorNormalizer = " << posteriorNormalizer << endl;
        /*cout << "comboProbs = [";
        for (vector<long double>::iterator i = comboProbs.begin(); i != comboProbs.end(); ++i)
            cout << *i << " ";
        cout << "]" << endl;
        cout << "posteriorNormalizer = " << posteriorNormalizer << endl;
        */

        // normalize marginals
        for (Results::iterator r = results.begin(); r != results.end(); ++r) {
            ResultData& d = r->second;
            for (map<Genotype, vector<long double> >::iterator m = d.rawMarginals.begin(); m != d.rawMarginals.end(); ++m) {
                /*
                cout << "rawMarginals = [";
                for (vector<long double>::iterator i = m->second.begin(); i != m->second.end(); ++i)
                    cout << *i << ", ";
                cout << "]" << endl;
                cout << logsumexp(m->second) << endl;
                */
                d.marginals[m->first] = logsumexp(m->second) - posteriorNormalizer;
            }
        }


        GenotypeCombo& bestGenotypeCombo = genotypeComboProbs.front().first;
        long double bestGenotypeComboProb = exp(genotypeComboProbs.front().second - posteriorNormalizer);

        if (!parameters.suppressOutput) {
            //cerr << parser->currentPosition << " " << alleles.size() << " " << bestGenotypeComboProb << " " << genotypeComboProbs.front().second << " " <<  posteriorNormalizer << endl;
            json(cout, results, *parser);
            cout << endl;
        }

    }

    delete parser;

    return 0;

}
