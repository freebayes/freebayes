#ifndef __RESULT_DATA_H
#define __RESULT_DATA_H

#include <vector>
#include <ostream>

//#include <boost/tuple/tuple.hpp>
#include <boost/bind.hpp>

#include <time.h>

#include "Genotype.h"
#include "Allele.h"
#include "Utility.h"
#include "AlleleParser.h"

using namespace std;

class ResultData;

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
        sort(dataLikelihoods.begin(), 
                dataLikelihoods.end(), 
                boost::bind(&pair<Genotype*, long double>::second, _1) 
                    > boost::bind(&pair<Genotype*, long double>::second, _2));
    }

    friend void json(ostream& out, Results& results, AlleleParser* parser);
    friend void vcf(ostream& out,
            long double comboProb,
            long double alleleSamplingProb,
            string alternateBase,
            vector<string>& samples,
            list<Allele*> observedAlleles,
            Results& results,
            AlleleParser* parser);
    pair<Genotype*, long double> bestMarginalGenotype(void);

};

string dateStr(void);
void vcfHeader(ostream& out, string referenceFileName, vector<string>& samples);



// TODO vcf output


#endif
