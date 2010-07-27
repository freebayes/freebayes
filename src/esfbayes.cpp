// 
// *Bayes
//
// A bayesian genetic variant caller.
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
#include "Caller.h"

#include "multichoose.h"
#include "multipermute.h"

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

typedef vector<
            tuple<
                string,  // sample ID
                vector<pair<Genotype, long double> >, // raw genotype probabilities
                vector<Allele*> // observations
                > 
            > Results;

int main (int argc, char *argv[]) {

    Caller* caller = new Caller(argc, argv);
    list<Allele*> alleles;

    // only estimate probabilities for these genotypes
    vector<Allele> genotypeAlleles;
    genotypeAlleles.push_back(genotypeAllele(ALLELE_SNP, "A", 1));
    genotypeAlleles.push_back(genotypeAllele(ALLELE_SNP, "T", 1));
    genotypeAlleles.push_back(genotypeAllele(ALLELE_SNP, "G", 1));
    genotypeAlleles.push_back(genotypeAllele(ALLELE_REFERENCE));
    genotypeAlleles.push_back(genotypeAllele(ALLELE_SNP, "C", 1));
    vector<vector<Allele> > genotypes = multichoose(2, genotypeAlleles);

    // set up index bands, used when generating best genotype combinations
    vector<int> indexBands;
    for (int i = 0; i < caller->parameters.WB; ++i) { indexBands.push_back(i); }

    while (caller->getNextAlleles(alleles)) {
        // skips 0-coverage regions
        if (alleles.size() == 0)
            continue;

        // TODO force calculation for samples not in this list
        map<string, vector<Allele*> > sampleGroups = groupAllelesBySample(alleles);

        // vebose; alias this with a typedef?
        Results results;

        // calculate data likelihoods
        for (map<string, vector< Allele* > >::iterator sampleAlleles = sampleGroups.begin();
                sampleAlleles != sampleGroups.end(); ++sampleAlleles) {

            vector<pair<Genotype, long double> > probs = 
                caller->probObservedAllelesGivenGenotypes(sampleAlleles->second, genotypes);
            
            normalizeGenotypeProbabilitiesln(probs);  // self-normalizes genotype probs
            // NB: if we were doing straight genotyping, this is where we would incorporate priors

            results.push_back(make_tuple(sampleAlleles->second.front()->sampleID, probs, sampleAlleles->second));

        }

        // now do ~pSnp and marginals estimation

        // ...
        // break this down into a banded estimation of the full probability
        // space (optionally do this in full to check that the banded
        // approximation approaches the full space)
        // 
        // store the terms within our bandwidth, as we'll use them several times
        // first pass; 
        // ... generate terms; sum of all probabilities for 

        // the term you want to sum over is:
        // product(caller->probObservedAllelesGivenGenotypes(...) for all individuals) * caller->priorProbGenotypes(...)
        // we'll do this for all genotype vectors within a bandwidth of changes from our best genotype vector
        // 
        // (1) best genotype vector:
        // a.) sort the vector<pair<Genotype, long double> >'s of our
        // 'results' data structure from highest to lowest probability
        for (Results::iterator sample = results.begin(); 
                sample != results.end(); ++sample) {
            vector<pair<Genotype, long double> >& probs = sample->get<1>();
            sort(probs.begin(), 
                    probs.end(), 
                    boost::bind(&pair<Genotype, long double>::second, _1) 
                        > boost::bind(&pair<Genotype, long double>::second, _2));
        }
        

        // b.) then generate a series of indexes for our genotype vector choice
        // -- The index series is all possible multiset permutations of all
        // multisets in N multichoose K where K is the number of individuals
        // and N is our bandwidth.

        // NB: we could generate this once, provided we always have data for
        // every individual or force data using priors otherwise; but for
        // simplicity and flexibility we do it here.
        //
        // the indexes of all 'best' genotype combinations
        // we will iterate through these to establish our posterior normalizer.
        // the indexes correspond to the Nth best genotype for each individual
        //
        vector<vector<int> > bandedIndexes = multichoose(results.size(), indexBands);

        map<string, map<Genotype, long double> > marginals;  // track the marginal contribution of individual genotypes, indexed by genotype

        // unclear what i meant to do here
        //pair<vector<Genotype>, long double> bestGenotypeCombo;  // track the best genotype combo and its probability

        long double posteriorNormalizer = 0;  // our posterior normalizer sum

        if (false) {

        for (vector<vector<int> >::iterator b = bandedIndexes.begin(); b != bandedIndexes.end(); ++b) {
            
            // multipermute...
            vector<vector<int> > permutations = multipermute(*b);

            // for each multiset permutation:
            // XXX looks like this is not how the optimization is done...
            //      in the previous version the banded approximation only mixed a single 2nd index at a time
            //      ... also, this could be recast as a generator to avoid memcpys

            cerr << "permutations size " << permutations.size() << endl;

            for (vector<vector<int> >::iterator p = permutations.begin(); p != permutations.end(); ++p) {
            //
            //     1) get a genotype combo from the permutation
            // 
                vector<int>& indexes = *p;
                map<Allele, int> alleleFrequencies;
                int totalCopiesAtLocus = 0; // total number of copies at this locus
                int totalGenotypePermutationsInCombo = 0;

                // the current combination
                vector<pair<string, Genotype > > combo;

                // for each individual
                Results::iterator sampleResult = results.begin();
                for (vector<int>::const_iterator i = indexes.begin(); 
                        i != indexes.end();
                        ++i, ++sampleResult) {
                    // use the index to get a reference the Nth best genotype
                    map<Allele, int> currentGenotypeAlleleFrequencies;
                    pair<Genotype, long double>& nthBest = sampleResult->get<1>().at(*i);
                    // XXX you should get the sample ID here for marginals...
                    string sampleID = sampleResult->get<0>();
                    Genotype& genotype = nthBest.first;
                    long double& prob = nthBest.second;
                    combo.push_back(make_pair(sampleID, genotype));  // store a record of this genotype/individual combination

                    // get the log probability of that genotype from results, multiply our normalizer by it
                    posteriorNormalizer += prob;
                    // sum genotype into allele frequency (af) for this combo
                    totalCopiesAtLocus += genotype.size();

                    for (vector<Allele>::iterator a = genotype.begin(); a != genotype.end(); ++a) {
                        Allele& allele = *a;

                        // increment population-wide allele frequencies
                        map<Allele, int>::iterator freq = alleleFrequencies.find(allele);
                        if (freq == alleleFrequencies.end()) {
                            alleleFrequencies.insert(pair<Allele,int>(allele, 1));
                        } else {
                            freq->second += 1;
                        }

                        // increment current genotype allele frequencies
                        freq = currentGenotypeAlleleFrequencies.find(allele);
                        if (freq == currentGenotypeAlleleFrequencies.end()) {
                            currentGenotypeAlleleFrequencies.insert(pair<Allele,int>(allele, 1));
                        } else {
                            freq->second += 1;
                        }
                    }

                    // add the number of genotype permutations to the total number of genotype permutations
                    int currentGenotypePermutations;
                    if (currentGenotypeAlleleFrequencies.size() == 1) { // one frequency accounts for all alleles in this genotype
                        currentGenotypePermutations = 1;
                    } else {
                        int factProd = 1;
                        for (map<Allele, int>::iterator c = currentGenotypeAlleleFrequencies.begin(); 
                                c != currentGenotypeAlleleFrequencies.end(); ++c) {
                            factProd *= factorial(c->second);
                        }
                        currentGenotypePermutations = factorial(genotype.size()) / factProd;
                    }
                    totalGenotypePermutationsInCombo += currentGenotypePermutations;

                }
                
            //     2) get priors, multiply by them...
            //
                // calculate priors in two components:
                //   sampling probability,
                
                int lnFrequencyFactProd = 0;
                // also, convert from allele frequencies to frequency counts
                map<int, int> alleleFrequencyCounts;
                for (map<Allele, int>::const_iterator c = alleleFrequencies.begin(); c != alleleFrequencies.end(); ++c) {
                    lnFrequencyFactProd += factorialln(c->second);
                    map<int, int>::iterator fc = alleleFrequencyCounts.find(c->second);
                    if (fc == alleleFrequencyCounts.end()) {
                        alleleFrequencyCounts.insert(fc, pair<int, int>(c->second, 1));
                    } else {
                        ++(fc->second);
                    }
                }
                long double lnGenotypeComboSamplingProb = 1 - factorialln(totalCopiesAtLocus) + lnFrequencyFactProd + log(totalGenotypePermutationsInCombo);
            
            //   ewens' sampling formula for this genotype combination
            
                long double theta = caller->parameters.TH;
                long double lnAlleleFrequencyProb;

                // apply Ewens' Sampling Formula

                // we do this in two components...
                // the 'theta scaling' component
                long double lnThetaScalar = log(theta);
                for (int h = 1; h < totalCopiesAtLocus; ++h) {
                    lnThetaScalar += log(theta + h);
                }

                // and the frequency scaling one
                long double lnFreqScalar = 0;
                for (map<int, int>::const_iterator fc = alleleFrequencyCounts.begin(); fc != alleleFrequencyCounts.end(); ++fc) {
                    int j = fc->first;
                    int alleleCount = fc->second;
                    lnFreqScalar += log(pow(theta, alleleCount)) - log(pow(j, alleleCount)) + log(factorial(alleleCount));
                }

                long double lnEwensSamplingProbability = log(factorial(totalCopiesAtLocus) - lnThetaScalar) + lnFreqScalar;

                long double lnComboProb = lnGenotypeComboSamplingProb + lnEwensSamplingProbability;
                
            //     3) sum marginals into marginals data structure
            //
                // now, for each individual genotype,
                // add the result to the marginal for the individual genotypes which are represented in this combination
                for (vector<pair<string, Genotype> >::const_iterator i = combo.begin(); i != combo.end(); ++i) {
                //map<string, map<Genotype, long double> > marginals;  // track the marginal contribution of individual genotypes, indexed by genotype
                    Genotype genotype = i->second;
                    map<Genotype, long double>& sampleMarginals = marginals[i->first];
                    map<Genotype, long double>::iterator g = sampleMarginals.find(genotype);
                    if (g == sampleMarginals.end()) {
                        sampleMarginals.insert(make_pair(genotype, exp(lnComboProb)));
                    } else {
                        g->second += exp(lnComboProb);
                    }
                }

            }
        
        }
        // normalize marginals
        // report best genotype combo
        // normalize

        }


        // report in json-formatted stream
        //
        if (!caller->parameters.suppressOutput) {
            cout << "{\"sequence\":\"" << caller->currentTarget->seq << "\","
                << "\"position\":" << caller->currentPosition + 1 << ","  /// XXX basing somehow is 1-off... 
                //<< "\"raDepth\":" << caller->registeredAlleles.size() << ","
                << "\"samples\":{";  // TODO ... quality (~pSnp)

            bool suppressComma = true; // output flag
            for (vector<tuple<string, vector<pair<Genotype, long double> >, vector<Allele*> > >::iterator 
                    sample = results.begin();
                    sample != results.end(); ++sample) {

                if (!suppressComma) { cout << ","; } else { suppressComma = false; }

                cout << "\"" << sample->get<0>() << "\":{"
                    << "\"coverage\":" << sample->get<2>().size() << ","
                    << "\"genotypes\":[";

                vector<pair<Genotype, long double> >& probs = sample->get<1>();
                sort(probs.begin(), probs.end(), genotypeCmp);
                for (vector<pair<Genotype, long double> >::iterator g = probs.begin(); 
                        g != probs.end(); ++g) {
                    if (g != probs.begin()) cout << ",";
                    cout << "[\"" << g->first << "\"," << float2phred(1 - exp(g->second)) << "]";
                }
                cout << "]";
                if (caller->parameters.outputAlleles)
                    cout << ",\"alleles\":" << json(sample->get<2>());
                cout << "}";

            }

            cout << "}}" << endl;
        }

    }

    delete caller;

    return 0;

}

// discrete elements of analysis
// 
// 1) fasta reference
// 2) bam file(s) over samples / samples
// 3) per-individual base calls (incl cigar)
// 4) priors
// 
// sets of data per individual
// and sets of data per position
// 
// for each position in the target regions (which is provided by a bed file)
// calculate the basecalls for each sample
// then, for samples for which we meet certain criteria (filters):
//     number of mismatches
//     data sufficiency (number of individual basecalls aka reads?)
//     (readmask ?)
// ... establish the probability of a snp for each possible genotype
// (which is ~ the data likelihood * the prior probablity of a snp for each sample)
// and report the top (configurable) number of possible genotypes
// 
// 
// 
// high level overview of progression:
// 
// for each region in region list
//     for each position in region
//         for each read overlapping position
//             for each base in read
//                 register base
//         evaluate prob(variation | all bases)
//         if prob(variation) >= reporting threshold
//             report variation
// 
// 
// registration of bases:
// 
//     skip clips (soft and hard)
//     presently, skip indels and only analyze aligned bases, but open development
//         to working with them in the future
//     preadmask: when we encounter a indel alignment, mask out bases in that read
//         because we are concerned about the semantics of processing it.
//     
//     we keep data on the bases, the basecall struct contains this information
// 
// 
// probability estimation
// 
//     p ( snp ) ~= ...
// 
//     p ( individual genotype | reads ) ~= ...
// 
//     p ( genotypes | basecalls ) ~= p( basecalls | genotype ) * prior( genotype ) / probability ( basecalls )
// 
// 
// 
// algorithmic core overview:
// 
// (1) individual data likelihoods
// 
// for each sample in the sample list
//     get basecalls corresponding to sample
//         for each genotype from the fixed genotype list
//             calculate the data likelihoods of p ( basecalls | genotype )   == "data likelihood" 
//                  this amounts to multiplying the quality scores from all the basecalls in that sample
// 
// (2) total genotype likelhoods for dominant genotype combinations
// 
// for each genotype combo in dominant genotype combo list
//     data likelhood p ( basecall combo | genotype combo )
// 
// 
// (3) calculate priors for dominant genotype combinations
// 
// for each genotype combo in dominant genotype combo list
//     calculate priors of that genotype combo  (well defined)
// 
// (4) calculate posterior probability of dominant genotype combinations
// 
// for each genotype combo in dominant genotype combo list
//     multiply results of corresponding (2) * (3)
// normalize  (could be a separate step)
// 
// (5) probability that of a variation given all basecalls
// 
// sum over probability of all dominant variants
// 
// (6) calculate individual sample genotype posterior marginals
// 
// for each sample
//     for each genotype
//         sum of p(genotype | reads) for fixed genotype <-- (4)
