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
#include "DataLikelihood.h"

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



class ResultData
{
public:
    string name;
    vector<pair<Genotype, long double> > dataLikelihoods;
    map<Genotype, long double> marginals;
    vector<Allele*> observations;

    ResultData(string s,
        vector<pair<Genotype, long double> > d,
        map<Genotype, long double>  m,
        vector<Allele*> o)
            : name(s)
            , dataLikelihoods(d)
            , marginals(m)
            , observations(o)
    { }

    void sortDataLikelihoods(void) {
        sort(dataLikelihoods.begin(), 
                dataLikelihoods.end(), 
                boost::bind(&pair<Genotype, long double>::second, _1) 
                    > boost::bind(&pair<Genotype, long double>::second, _2));
    }

};

// maps sample names to results
typedef map<string, ResultData> Results;


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
            
            map<Genotype, long double> marginals;

            results.insert(make_pair(sampleName, ResultData(sampleName, probs, marginals, observedAlleles)));

        }

        // sort genotype data likelihoods
        
        for (Results::iterator s = results.begin(); s != results.end(); ++s) {
            s->second.sortDataLikelihoods();
        }



        // report in json-formatted stream
        //
        if (parameters.suppressOutput) {
            cout << "{\"sequence\":\"" << parser->currentTarget->seq << "\","
                << "\"position\":" << parser->currentPosition + 1 << ","  /// XXX basing somehow is 1-off... 
                //<< "\"raDepth\":" << parser->registeredAlleles.size() << ","
                << "\"samples\":{";  // TODO ... quality (~pSnp)

            bool suppressComma = true; // output flag
            for (Results::iterator p = results.begin(); p != results.end(); ++p) {

                ResultData& sample = p->second;

                if (!suppressComma) { cout << ","; } else { suppressComma = false; }

                cout << "\"" << sample.name << "\":{"
                    << "\"coverage\":" << sample.observations.size() << ","
                    << "\"genotypes\":[";

                for (map<Genotype, long double>::iterator g = sample.marginals.begin(); 
                        g != sample.marginals.end(); ++g) {
                    if (g != sample.marginals.begin()) cout << ",";
                    Genotype genotype = g->first;
                    cout << "[\"" << genotype << "\"," << float2phred(1 - exp(g->second)) << "]";
                }
                cout << "]";
                if (parser->parameters.outputAlleles)
                    cout << ",\"alleles\":" << json(sample.observations);
                cout << "}";

            }

            cout << "}}" << endl;
        }

    }

    delete parser;

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
