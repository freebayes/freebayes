#ifndef __GENOTYPE_H
#define __GENOTYPE_H

#include <iostream>
#include <vector>
#include <utility> // pair
#include <algorithm>
#include <numeric>
#include <vector>
#include <iterator>
#include "Allele.h"

using namespace std;


class Genotype : public vector<pair<int, Allele> > {

    friend ostream& operator<<(ostream& out, pair<int, Allele>& rhs);
    friend ostream& operator<<(ostream& out, Genotype& g);

public:
    
    int ploidy;
    vector<Allele> alleles;

    Genotype(vector<Allele>& ungroupedAlleles) {
        alleles = ungroupedAlleles;
        vector<vector<Allele> > groups = groupAlleles_copy(alleles);
        for (vector<vector<Allele> >::const_iterator group = groups.begin(); group != groups.end(); ++group) {
            this->push_back(make_pair(group->size(), group->front()));
        }
        ploidy = getPloidy();
    }

    int getPloidy(void);
    bool containsAllele(Allele& allele);
    // the probability of drawing each allele out of the genotype, ordered by allele
    vector<long double> alleleProbabilities(void);

};

// for sorting pairs of genotype, probs
/*bool genotypeCmp(pair<Genotype, long double> a, pair<Genotype, long double> b) {
    return a.second > b.second;
}
*/

vector<Genotype> allPossibleGenotypes(int ploidy, vector<Allele> potentialAlleles);

vector<vector<pair<string, Genotype> > >
bandedGenotypeCombinations(vector<pair<string, vector<Genotype> > > sampleGenotypes,
        int bandwidth, int banddepth);


#endif
