#ifndef __SAMPLE_H
#define __SAMPLE_H

#include <string>
#include <vector>
#include <map>
#include <utility>
#include "Utility.h"
#include "Allele.h"
#include "Genotype.h"

using namespace std;

// sample tracking and allele sorting
class Sample : public map<string, vector<Allele*> > {

public:

    // the number of observations for this allele
    int observationCount(Allele& allele);

    // the number of observations for this base
    int observationCount(const string& base);

    // the total number of observations
    int observationCount(void);

    // puts alleles into the right bins if they have changed their base (as
    // occurs in the case of reference alleles)
    void sortReferenceAlleles(void);

    pair<pair<int, int>, pair<int, int> >
    baseCount(string refbase, string altbase);

    vector<Genotype*> genotypesWithEvidence(vector<Genotype>& genotypes);

    int baseCount(string base, AlleleStrand strand);

    string json(void);

    vector<int> alleleObservationCounts(Genotype& genotype);

};

class Samples : public map<string, Sample> { };



int countAlleles(Samples& samples);
// using this one...
void groupAlleles(Samples& samples, map<string, vector<Allele*> >& alleleGroups);

// filters... maybe move to its own file?
bool sufficientAlternateObservations(Samples& observations, int mincount, float minfraction);


#endif
