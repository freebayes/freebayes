#include "Genotype.h"
#include "multichoose.h"
#include "multipermute.h"


int Genotype::getPloidy(void) {
    int result = 0;
    for (Genotype::const_iterator i = this->begin(); i != this->end(); ++i) {
        result += i->first;
    }
    return result;
}

bool Genotype::containsAllele(Allele& allele) {
    for (Genotype::iterator i = this->begin(); i != this->end(); ++i) {
        Allele& b = i->second;
        if (allele == b)
            return true;
    }
    return false;
}

// the probability of drawing each allele out of the genotype, ordered by allele
vector<long double> Genotype::alleleProbabilities(void) {
    vector<long double> probs;
    for (vector<pair<int, Allele> >::const_iterator a = this->begin(); a != this->end(); ++a) {
        probs.push_back((long double) a->first / (long double) ploidy);
    }
    return probs;
}

ostream& operator<<(ostream& out, pair<int, Allele>& rhs) {
    out << rhs.first << rhs.second.alternateSequence;
    return out;
}

ostream& operator<<(ostream& out, Genotype& g) {
    Genotype::iterator i = g.begin(); ++i;
    out << g.front();
    for (;i != g.end(); ++i) {
        out << " " << *i;
    }
    return out;
}

vector<Genotype> allPossibleGenotypes(int ploidy, vector<Allele> potentialAlleles) {
    vector<Genotype> genotypes;
    vector<vector<Allele> > alleleCombinations = multichoose(ploidy, potentialAlleles);
    for (vector<vector<Allele> >::iterator combo = alleleCombinations.begin(); combo != alleleCombinations.end(); ++combo) {
        genotypes.push_back(Genotype(*combo));
    }
    return genotypes;
}

/*
def banded_genotype_combinations(sample_genotypes, bandwidth, band_depth):
    # always provide the 'best' case
    yield [(sample, genotypes[0]) for sample, genotypes in sample_genotypes]
    for i in range(1, bandwidth):
        for j in range(1, band_depth):  # band_depth is the depth to which we explore the bandwith... TODO explain better
            indexes = j * [i] + (len(sample_genotypes) - j) * [0]
            for index_permutation in multiset.permutations(indexes):
                yield [(sample, genotypes[index]) for index, (sample, genotypes) in zip(index_permutation, sample_genotypes)] 
*/

vector<vector<pair<string, Genotype> > >
bandedGenotypeCombinations(vector<pair<string, vector<Genotype> > > sampleGenotypes,
        int bandwidth, int banddepth) {

    int nsamples = sampleGenotypes.size();
    vector<vector<pair<string, Genotype> > > combos;
    vector<pair<string, Genotype> > firstcombo;
    for (vector<pair<string, vector<Genotype> > >::const_iterator s = sampleGenotypes.begin();
            s != sampleGenotypes.end(); ++s) {
        firstcombo.push_back(make_pair(s->first, s->second.at(0)));
    }
    combos.push_back(firstcombo);
    for (int i = 0; i < bandwidth; ++i) {
        for (int j = 1; j < banddepth; ++j) {
            vector<int> indexes;
            for (int h = 0; h < j; ++h)
                indexes.push_back(i);
            for (int h = 0; h < (nsamples - j); ++h)
                indexes.push_back(0);
            vector<vector<int> > indexPermutations = multipermute(indexes);
            for (vector<vector<int> >::const_iterator p = indexPermutations.begin(); p != indexPermutations.end(); ++p) {
                vector<pair<string, Genotype> > combo;
                vector<int>::const_iterator n = p->begin();
                for (vector<pair<string, vector<Genotype> > >::const_iterator s = sampleGenotypes.begin();
                        s != sampleGenotypes.end(); ++s, ++n) {
                    combo.push_back(make_pair(s->first, s->second.at(*n)));
                }
                combos.push_back(combo);
            }
        }
    }
    return combos;
}


