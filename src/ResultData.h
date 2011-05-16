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

class Result;

// for sorting data likelihoods
class DataLikelihoodCompare {
public:
    bool operator()(const pair<Genotype*, long double>& a,
            const pair<Genotype*, long double>& b) {
        return a.second > b.second;
    }
};

class Results;

class Result : public vector<SampleDataLikelihood> {

public:

    string name;
    Sample* observations;

    void sortDataLikelihoods(void) {
        SampleDataLikelihoodCompare datalikelihoodCompare;
        sort(begin(), end(), datalikelihoodCompare);
    }

    pair<Genotype*, long double> bestMarginalGenotype(void);

};

// maps sample names to results
class Results : public map<string, Result> {

public:
    void update(SampleDataLikelihoods& likelihoods) {
        for (SampleDataLikelihoods::iterator s = likelihoods.begin(); s != likelihoods.end(); ++s) {
            vector<SampleDataLikelihood>& sdls = *s;
            string& name = sdls.front().name;
            Result& result = (*this)[name];
            result.clear();
            for (vector<SampleDataLikelihood>::iterator s = sdls.begin(); s != sdls.end(); ++s) {
                result.push_back(*s);
            }
        }
    }

    string vcf(
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
            AlleleParser* parser);

};


string dateStr(void);
void vcfHeader(ostream& out, string referenceFileName, vector<string>& samples, Parameters& parameters, vector<string>& sequencingTechnologies);


#endif
