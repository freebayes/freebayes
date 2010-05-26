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

// private libraries
#include "BamReader.h"
#include "Fasta.h"
#include "TryCatch.h"
#include "Parameters.h"
#include "Allele.h"
#include "Caller.h"

#include "multichoose.h"

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

    Caller* caller = new Caller(argc, argv);
    list<Allele*> alleles;

    // only estimate probabilities for these genotypes
    vector<Allele> genotypeAlleles;
    genotypeAlleles.push_back(genotypeAllele(ALLELE_REFERENCE));
    genotypeAlleles.push_back(genotypeAllele(ALLELE_SNP, "A", 1));
    genotypeAlleles.push_back(genotypeAllele(ALLELE_SNP, "T", 1));
    genotypeAlleles.push_back(genotypeAllele(ALLELE_SNP, "G", 1));
    genotypeAlleles.push_back(genotypeAllele(ALLELE_SNP, "C", 1));
    vector<vector<Allele> > genotypes = multichoose(2, genotypeAlleles);

    while (caller->getNextAlleles(alleles)) {
        // skips 0-coverage regions
        if (alleles.size() == 0)
            continue;

        // TODO force calculation for samples not in this list
        map<string, vector<Allele*> > sampleGroups = groupAllelesBySample(alleles);

        // vebose; alias this with a typedef?
        vector<
            tuple<
                string,  // sample ID
                vector<pair<Genotype, long double> >, // raw genotype probabilities
                vector<Allele*> // observations
                > 
            > results;

        // calculate data likelihoods
        for (map<string, vector< Allele* > >::iterator sampleAlleles = sampleGroups.begin();
                sampleAlleles != sampleGroups.end(); ++sampleAlleles) {

            vector<pair<Genotype, long double> > probs = 
                caller->probObservedAllelesGivenGenotypes(sampleAlleles->second, genotypes);
            
            normalizeGenotypeProbabilitiesln(probs);  // self-normalizes genotype probs
            // if we were doing straight genotyping, this is where we would incorporate priors

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


        // report in json-formatted stream
        //
        cout << "{\"sequence\":\"" << caller->currentTarget->seq << "\","
            << "\"position\":" << caller->currentPosition + 1 << ","  /// XXX basing somehow is 1-off... 
            //<< "\"raDepth\":" << caller->registeredAlleles.size() << ","
            << "\"samples\":{";  // TODO ... quality (~pSnp)

        bool suppressComma = true; // output flag
        for (vector<tuple<string, vector<pair<Genotype, long double> >, vector<Allele*> > >::iterator sample = results.begin();
                sample != results.end(); ++sample) {

            if (!suppressComma) { cout << ","; } else { suppressComma = false; }

            cout << "\"" << sample->get<0>() << "\":{"
                << "\"coverage\":" << sample->get<2>().size() << ","
                << "\"genotypes\":{";

            vector<pair<Genotype, long double> >& probs = sample->get<1>();
            for (vector<pair<Genotype, long double> >::iterator g = probs.begin(); 
                    g != probs.end(); ++g) {
                if (g != probs.begin()) cout << ",";
                cout << "\"" << g->first << "\":" << float2phred(1 - exp(g->second));
            }
            cout << "}";
            if (caller->parameters.outputAlleles)
                cout << ",\"alleles\":" << json(sample->get<2>());
            cout << "}";

        }

        cout << "}}" << endl;

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
