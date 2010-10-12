#ifndef __RESULT_DATA_H
#define __RESULT_DATA_H

#include <vector>
#include <ostream>
#include <time.h>
#include "Genotype.h"
#include "Allele.h"
#include "Utility.h"
#include "AlleleParser.h"
#include "Version.h"

using namespace std;

class ResultData;

// for sorting data likelihoods
class DataLikelihoodCompare {
public:
    bool operator()(const pair<Genotype*, long double>& a,
            const pair<Genotype*, long double>& b) {
        return a.second > b.second;
    }
};

// maps sample names to results
typedef map<string, ResultData> Results;

class ResultData
{
public:
    string name;
    vector<pair<Genotype*, long double> > dataLikelihoods;
    map<Genotype*, vector<long double> > rawMarginals;
    map<Genotype*, long double> marginals;
    vector<Allele*> observations;

    ResultData(string s,
        vector<pair<Genotype*, long double> > d,
        map<Genotype*, long double>  m,
        map<Genotype*, vector<long double> > rm,
        vector<Allele*> o)
            : name(s)
            , dataLikelihoods(d)
            , marginals(m)
            , rawMarginals(rm)
            , observations(o)
    { }

    ResultData(void) { }

    ResultData(const ResultData& r) {
        name = r.name;
        dataLikelihoods = r.dataLikelihoods;
        marginals = r.marginals;
        rawMarginals = r.rawMarginals;
        observations = r.observations;
    }

    void sortDataLikelihoods(void) {
        DataLikelihoodCompare datalikelihoodCompare;
        sort(dataLikelihoods.begin(), 
                dataLikelihoods.end(), 
                datalikelihoodCompare);
    }

    friend void json(ostream& out, Results& results, AlleleParser* parser);
    friend string vcf(
            long double comboProb,
            //long double alleleSamplingProb,
            map<string, vector<Allele*> >& sampleObservations,
            string referenceBase,
            string alternateBase,
            vector<string>& samples,
            int coverage,
            GenotypeCombo& genotypeCombo,
            Results& results,
            AlleleParser* parser);
    pair<Genotype*, long double> bestMarginalGenotype(void);

};

string dateStr(void);
void vcfHeader(ostream& out, string referenceFileName, vector<string>& samples);


#endif
