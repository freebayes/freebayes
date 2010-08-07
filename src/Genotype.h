#ifndef __GENOTYPE_H
#define __GENOTYPE_H

#include <iostream>
#include <vector>
#include <utility> // pair
#include <algorithm>
#include <numeric>
#include <vector>
#include <iterator>
#include <math.h>
#include "Allele.h"

using namespace std;


// TODO
// Develop this stub to clean up leaky abstraction caused by direct use of
// a pair as a complex datatype.
// It's clearer to write things like 'genotypeelement->allele' and
// 'genotypeelement->count' than ..->first and ..->second
/*
class GenotypeElement : public pair<Allele, int> {

    friend ostream& operator<<(ostream& out, GenotypeElement& rhs);

public:
    Allele& allele(void) { return this->first; }
    int& count(void) { return this->second; }

};
*/


class Genotype : public vector<pair<Allele, int> > {

    friend ostream& operator<<(ostream& out, const pair<Allele, int>& rhs);
    friend ostream& operator<<(ostream& out, const Genotype& g);
    friend bool operator<(Genotype& a, Genotype& b);

public:
    
    int ploidy;
    vector<Allele> alleles;

    Genotype(vector<Allele>& ungroupedAlleles) {
        alleles = ungroupedAlleles;
        sort(alleles.begin(), alleles.end());
        vector<vector<Allele> > groups = groupAlleles_copy(alleles);
        for (vector<vector<Allele> >::const_iterator group = groups.begin(); group != groups.end(); ++group) {
            this->push_back(make_pair(group->front(), group->size()));
        }
        ploidy = getPloidy();
    }

    vector<Allele> uniqueAlleles(void);
    int getPloidy(void);
    bool containsAllele(Allele& allele);
    bool containsAlleleOtherThan(string& base);
    vector<Allele> alternateAlleles(string& refbase);
    // the probability of drawing each allele out of the genotype, ordered by allele
    vector<long double> alleleProbabilities(void);
    string str(void);
    string relativeGenotype(string& refbase);

};

// for sorting pairs of genotype, probs
/*bool genotypeCmp(pair<Genotype, long double> a, pair<Genotype, long double> b) {
    return a.second > b.second;
}
*/

vector<Genotype> allPossibleGenotypes(int ploidy, vector<Allele> potentialAlleles);

typedef vector<pair<string, pair<Genotype, long double> > > GenotypeCombo;

vector<GenotypeCombo>
bandedGenotypeCombinations(
        vector<pair<string, vector<pair<Genotype, long double> > > >& sampleGenotypes,
        int bandwidth, int banddepth);

ostream& operator<<(ostream& out, vector<GenotypeCombo>& combo);
ostream& operator<<(ostream& out, GenotypeCombo& g);


#endif
