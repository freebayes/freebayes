#ifndef __RESULT_DATA_H
#define __RESULT_DATA_H

#include <vector>
#include <ostream>
#include <iomanip>
#include "Genotype.h"
#include "Allele.h"
#include "Utility.h"
#include "AlleleParser.h"
#include "../vcflib/Variant.h"
#include "Version.h"
#include "Result.h"

using namespace std;

// for sorting data likelihoods
class DataLikelihoodCompare {
public:
    bool operator()(const pair<Genotype*, long double>& a,
            const pair<Genotype*, long double>& b) {
        return a.second > b.second;
    }
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

    vcf::Variant& vcf(
        vcf::Variant& var, // variant to update
        BigFloat pHom,
        long double bestComboOddsRatio,
        //long double alleleSamplingProb,
        Samples& samples,
        string refbase,
        vector<Allele>& altAlleles,
        map<string, int> repeats,
	int genotypingIterations,
        vector<string>& sampleNames,
        int coverage,
        GenotypeCombo& genotypeCombo,
        bool bestOverallComboIsHet,
        map<string, vector<Allele*> >& alleleGroups,
        map<string, vector<Allele*> >& partialObservationGroups,
        map<Allele*, set<Allele*> >& partialSupport,
        map<int, vector<Genotype> >& genotypesByPloidy,
        vector<string>& sequencingTechnologies,
        AlleleParser* parser);

};


string dateStr(void);
void vcfHeader(ostream& out, string referenceFileName, vector<string>& samples, Parameters& parameters, vector<string>& sequencingTechnologies);


#endif
