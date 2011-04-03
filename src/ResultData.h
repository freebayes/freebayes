#ifndef __RESULT_DATA_H
#define __RESULT_DATA_H

#include <vector>
#include <ostream>
#include <iomanip>
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
    map<Genotype*, long double> marginals;
    Sample* observations;

    ResultData(string s,
        vector<pair<Genotype*, long double> > d,
        map<Genotype*, long double>  m,
        Sample* o)
            : name(s)
            , dataLikelihoods(d)
            , marginals(m)
            , observations(o)
    { }

    ResultData(void) { }

    ResultData(const ResultData& r) {
        name = r.name;
        dataLikelihoods = r.dataLikelihoods;
        marginals = r.marginals;
        observations = r.observations;
    }

    void sortDataLikelihoods(void) {
        DataLikelihoodCompare datalikelihoodCompare;
        sort(dataLikelihoods.begin(), 
                dataLikelihoods.end(), 
                datalikelihoodCompare);
    }

    long double genotypeLikelihood(Genotype* g) {
        for (vector<pair<Genotype*, long double> >::iterator gl = dataLikelihoods.begin();
                gl != dataLikelihoods.end(); ++gl) {
            if (gl->first == g)
                return gl->second;
        }
    }

    friend void json(ostream& out, Results& results, AlleleParser* parser);
    friend string vcf(
            long double comboProb,
            //long double alleleSamplingProb,
            Samples& sample,
            string referenceBase,
            string alternateBase,
            Allele& altAllele,
            map<string, int> repeats,
            vector<string>& samples,
            int coverage,
            GenotypeCombo& genotypeCombo,
            bool bestOverallComboIsHet,
            map<string, vector<Allele*> >& alleleGroups,
            map<int, vector<Genotype> >& genotypesByPloidy, // pass by copy, will modify
            vector<string>& sequencingTechnologies,
            Results& results,
            AlleleParser* parser);
    pair<Genotype*, long double> bestMarginalGenotype(void);

};

string dateStr(void);
void vcfHeader(ostream& out, string referenceFileName, vector<string>& samples, Parameters& parameters, vector<string>& sequencingTechnologies);


#endif
