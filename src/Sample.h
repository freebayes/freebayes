#ifndef __SAMPLE_H
#define __SAMPLE_H

#include <string>
#include <vector>
#include <map>
#include <utility>
#include "Utility.h"
#include "Allele.h"

using namespace std;

// sample tracking and allele sorting
class Sample : public map<string, vector<Allele*> > {

    friend ostream& operator<<(ostream& out, Sample& sample);

public:

    // estimates of contamination
    double probRefGivenHomAlt;
    // and reference bias/contamination
    double probRefGivenHet;

    // the number of observations for this allele
    int observationCount(Allele& allele);

    // the number of observations for this base
    int observationCount(const string& base);

    // the total number of observations
    int observationCount(void);

    // sum of quality for the given allele
    int qualSum(Allele& allele);
    int qualSum(const string& base);

    // puts alleles into the right bins if they have changed their base (as
    // occurs in the case of reference alleles)
    void sortReferenceAlleles(void);

    pair<pair<int, int>, pair<int, int> >
    baseCount(string refbase, string altbase);

    int baseCount(string base, AlleleStrand strand);

    string json(void);

};

class Samples : public map<string, Sample> { };



int countAlleles(Samples& samples);
// using this one...
void groupAlleles(Samples& samples, map<string, vector<Allele*> >& alleleGroups);

// filters... maybe move to its own file?
bool sufficientAlternateObservations(Samples& observations, int mincount, float minfraction);


#endif
