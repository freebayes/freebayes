#include "Genotype.h"
#include "multichoose.h"
#include "multipermute.h"


vector<Allele> Genotype::uniqueAlleles(void) {
    vector<Allele> uniques;
    for (Genotype::iterator g = this->begin(); g != this->end(); ++g)
        uniques.push_back(g->first);
    return uniques;
}

int Genotype::getPloidy(void) {
    int result = 0;
    for (Genotype::const_iterator i = this->begin(); i != this->end(); ++i) {
        result += i->second;
    }
    return result;
}

vector<Allele> Genotype::alternateAlleles(string& base) {
    vector<Allele> alleles;
    for (Genotype::iterator i = this->begin(); i != this->end(); ++i) {
        Allele& b = i->first;
        if (base != b.base())
            alleles.push_back(b);
    }
    return alleles;
}

string Genotype::relativeGenotype(string& refbase) {
    string rg;
    for (Genotype::iterator i = this->begin(); i != this->end(); ++i) {
        Allele& b = i->first;
        if (refbase != b.base()) {
            for (int j = 0; j < i->second; ++j)
                rg += "1/";
        } else {
            for (int j = 0; j < i->second; ++j)
                rg += "0/";
        }
    }
    return rg.substr(0, rg.size() - 1); // chop trailing '/'
}

bool Genotype::containsAlleleOtherThan(string& base) {
    for (Genotype::iterator i = this->begin(); i != this->end(); ++i) {
        Allele& b = i->first;
        if (base != b.base())
            return true;
    }
    return false;
}

bool Genotype::containsAllele(Allele& allele) {
    for (Genotype::iterator i = this->begin(); i != this->end(); ++i) {
        Allele& b = i->first;
        if (allele == b)
            return true;
    }
    return false;
}

bool Genotype::homozygous(void) {
    return this->size() == 1;
}

// the probability of drawing each allele out of the genotype, ordered by allele
vector<long double> Genotype::alleleProbabilities(void) {
    vector<long double> probs;
    for (vector<pair<Allele, int> >::const_iterator a = this->begin(); a != this->end(); ++a) {
        probs.push_back((long double) a->second / (long double) ploidy);
    }
    return probs;
}

string Genotype::str(void) {
    string s;
    for (Genotype::const_iterator allele = this->begin(); allele != this->end(); ++allele) {
        for (int i = 0; i < allele->second; ++i)
            s += allele->first.alternateSequence;
    }
    return s;
}

ostream& operator<<(ostream& out, const pair<Allele, int>& rhs) {
    out << rhs.first.alternateSequence << rhs.second;
    return out;
}

ostream& operator<<(ostream& out, const Genotype& g) {
    Genotype::const_iterator i = g.begin(); ++i;
    out << g.front();
    for (;i != g.end(); ++i) {
        out << " " << *i;
    }
    return out;
}

ostream& operator<<(ostream& out, vector<GenotypeCombo>& g) {
    for (vector<GenotypeCombo>::iterator i = g.begin(); i != g.end(); ++i) {
        out << *i << endl;
    }
    return out;
}

ostream& operator<<(ostream& out, GenotypeCombo& g) {
    GenotypeCombo::iterator i = g.begin(); ++i;
    out << "{\"" << g.front().first << "\":[\"" << g.front().second.first << "\"," << exp(g.front().second.second) << "]";
    for (;i != g.end(); ++i) {
        out << ", \"" << i->first << "\":[\"" << i->second.first << "\"," << exp(i->second.second) << "]";
    }
    out << "}";
    return out;
}


bool operator<(Genotype& a, Genotype& b) {
    // genotypes of different ploidy are evaluated according to their relative ploidy
    if (a.ploidy != b.ploidy)
        return a.ploidy < b.ploidy;
    // because our constructor sorts each Genotype.alleles, we assume that we
    // have two equivalently sorted vectors to work with
    Genotype::iterator ai = a.begin();
    Genotype::iterator bi = b.begin();
    // step through each genotype, and if we find a difference between either
    // their allele or count return a<b
    for (; ai != a.end() && bi != b.end(); ++ai, ++bi) {
        if (ai->first != bi->first)
            return ai->first < bi->first;
        else if (ai->second != bi->second)
            return ai->second < bi->second;
    }
    return false; // if the two are equal, then we return false per C++ convention
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

vector<GenotypeCombo>
bandedGenotypeCombinations(
        vector<pair<string, vector<pair<Genotype, long double> > > >& sampleGenotypes,
        int bandwidth, int banddepth) {

    int nsamples = sampleGenotypes.size();
    vector<GenotypeCombo> combos;
    /*
    GenotypeCombo firstcombo;
    for (vector<pair<string, vector<pair<Genotype, long double> > > >::const_iterator s = sampleGenotypes.begin();
            s != sampleGenotypes.end(); ++s) {
        firstcombo.push_back(make_pair(s->first, s->second.at(0)));
    }
    combos.push_back(firstcombo);
    */
    for (int i = 0; i < bandwidth; ++i) {
        for (int j = 1; j < banddepth; ++j) {
            vector<int> indexes;
            for (int h = 0; h < j; ++h)
                indexes.push_back(i);
            for (int h = 0; h < (nsamples - j); ++h)
                indexes.push_back(0);
            vector<vector<int> > indexPermutations = multipermute(indexes);
            for (vector<vector<int> >::const_iterator p = indexPermutations.begin(); p != indexPermutations.end(); ++p) {
                GenotypeCombo combo;
                vector<int>::const_iterator n = p->begin();
                for (vector<pair<string, vector<pair<Genotype, long double> > > >::const_iterator s = sampleGenotypes.begin();
                        s != sampleGenotypes.end(); ++s, ++n) {
                    combo.push_back(make_pair(s->first, s->second.at(*n)));
                }
                combos.push_back(combo);
            }
        }
    }
    return combos;
}


vector<GenotypeCombo>
bandedGenotypeCombinationsIncludingBestHomozygousCombo(
        vector<pair<string, vector<pair<Genotype, long double> > > >& sampleGenotypes,
        int bandwidth, int banddepth) {

    vector<GenotypeCombo> combos = bandedGenotypeCombinations(sampleGenotypes, bandwidth, banddepth);
    // is there already a homozygous combo?
    bool hasHomozygousCombo = false;
    for (vector<GenotypeCombo>::iterator c = combos.begin(); c != combos.end(); ++c) {
        bool allhomozygous = true;
        for (GenotypeCombo::iterator gc = c->begin(); gc != c->end(); ++gc) {
            if (!gc->second.first.homozygous()) {
                allhomozygous = false;
                break;
            }
        }
        if (allhomozygous) {
            hasHomozygousCombo = true;
            break;
        }
    }
    if (!hasHomozygousCombo) {
        GenotypeCombo homozygousCombo;
        // push back the best homozygous combo
        for(vector<pair<string, vector<pair<Genotype, long double> > > >::iterator s = sampleGenotypes.begin();
                s != sampleGenotypes.end(); ++s) {
            for (vector<pair<Genotype, long double> >::iterator g = s->second.begin(); g != s->second.end(); ++g) {
                if (g->first.homozygous()) {
                    homozygousCombo.push_back(make_pair(s->first, *g));
                    break;
                }
            }
        }
        combos.push_back(homozygousCombo);
    }

    return combos;

}
