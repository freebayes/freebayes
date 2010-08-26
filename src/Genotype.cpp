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

vector<int> Genotype::alleleCountsInObservations(vector<Allele*> observations) {
    vector<int> counts;
    for (Genotype::iterator i = this->begin(); i != this->end(); ++i) {
        int count = 0;
        Allele& b = i->first;
        for (vector<Allele*>::iterator o = observations.begin(); o != observations.end(); ++o) {
            Allele& obs = **o;
            if (obs.base() == b.base())
                ++count;
        }
        counts.push_back(count);
    }
    return counts;
}

/*
void inOutObservationCounts(vector<Allele*> observations,
        map<Allele, int>& inCounts, 
        map<Allele, int>& outCounts) {
    vector<Allele> uniques = uniqueAlleles();
    for (vector<Allele>::iterator a = uniques.begin(); a != uniques.end(); ++a) {
        inCounts.insert(make_pair(*a, 0));
    }
    map<Allele, int>::iterator count;
    for (vector<Allele*>::iterator o = observations.begin(); o != observations.end(); ++o) {
        Allele& a = **o;
        count = inCounts.find(a);
        if (count != inCounts.end()) {
            count->second += 1;
        } else {
            count = outCounts.find(a);
            if (count != outCounts.end()) {
                count->second += 1;
            } else {
                outCounts.insert(make_pair(a, 1));
            }
        }
    }
}
*/

vector<Allele> Genotype::alternateAlleles(string& base) {
    vector<Allele> alleles;
    for (Genotype::iterator i = this->begin(); i != this->end(); ++i) {
        Allele& b = i->first;
        if (base != b.base())
            alleles.push_back(b);
    }
    return alleles;
}

string Genotype::relativeGenotype(string& refbase, string& altbase) {
    vector<string> rg;
    for (Genotype::iterator i = this->begin(); i != this->end(); ++i) {
        Allele& b = i->first;
        if (b.base() == altbase && refbase != b.base()) {
            for (int j = 0; j < i->second; ++j)
                rg.push_back("1/");
        } else {
            for (int j = 0; j < i->second; ++j)
                rg.push_back("0/");
        }
    }
    sort(rg.begin(), rg.end()); // enforces the same ordering for all genotypes
    string result = accumulate(rg.begin(), rg.end(), string(""));
    return result.substr(0, result.size() - 1); // chop trailing '/'
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

string IUPAC(Genotype& genotype) {
    const string g = genotype.str();
    if (g == "AA") return "A";
    if (g == "AC") return "M";
    if (g == "AG") return "R";
    if (g == "AT") return "W";
    if (g == "CA") return "M";
    if (g == "CC") return "C";
    if (g == "CG") return "S";
    if (g == "CT") return "Y";
    if (g == "GA") return "R";
    if (g == "GC") return "S";
    if (g == "GG") return "G";
    if (g == "GT") return "K";
    if (g == "TA") return "W";
    if (g == "TC") return "Y";
    if (g == "TG") return "K";
    if (g == "TT") return "T";
    return g;
}

ostream& operator<<(ostream& out, const pair<Allele, int>& rhs) {
    out << rhs.first.alternateSequence << rhs.second;
    //for (int i = 0; i < rhs.second; ++i)
    //    out << rhs.first.alternateSequence;
    return out;
}

ostream& operator<<(ostream& out, const Genotype& g) {
    Genotype::const_iterator i = g.begin(); ++i;
    out << g.front();
    for (;i != g.end(); ++i) {
        out << *i;
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
    out << "{\"" << g.front().first << "\":[\"" << *(g.front().second.first) << "\"," << exp(g.front().second.second) << "]";
    for (;i != g.end(); ++i) {
        out << ", \"" << i->first << "\":[\"" << *(i->second.first) << "\"," << exp(i->second.second) << "]";
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
        vector<pair<string, vector<pair<Genotype*, long double> > > >& sampleGenotypes,
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
                for (vector<pair<string, vector<pair<Genotype*, long double> > > >::const_iterator s = sampleGenotypes.begin();
                        s != sampleGenotypes.end(); ++s, ++n) {
                    combo.push_back(make_pair(s->first, s->second.at(*n)));
                }
                combos.push_back(combo);
            }
        }
    }
    return combos;
}

// TODO we have to evaluate all homozygous combos

vector<GenotypeCombo>
bandedGenotypeCombinationsIncludingAllHomozygousCombos(
        vector<pair<string, vector<pair<Genotype*, long double> > > >& sampleGenotypes,
        vector<Genotype>& genotypes,
        int bandwidth, int banddepth) {

    // obtain the combos

    vector<GenotypeCombo> combos = bandedGenotypeCombinations(sampleGenotypes, bandwidth, banddepth);

    // determine which homozygous combos we already have

    map<Genotype*, bool> genotypesWithHomozygousCombos;

    for (vector<GenotypeCombo>::iterator c = combos.begin(); c != combos.end(); ++c) {
        bool allSameAndHomozygous = true;
        GenotypeCombo::iterator gc = c->begin();
        Genotype* genotype;
        if (gc->second.first->homozygous()) {
            genotype = gc->second.first;
        } else {
            continue;
        }
        for (; gc != c->end(); ++gc) {
            if (! (gc->second.first == genotype) ) {
                allSameAndHomozygous = false;
                break;
            }
        }
        if (allSameAndHomozygous) {
            genotypesWithHomozygousCombos[genotype] == true;
        }
    }

    // accumulate the needed homozygous combos

    map<Genotype*, GenotypeCombo> homozygousCombos;

    for (vector<Genotype>::iterator genotype = genotypes.begin(); genotype != genotypes.end(); ++genotype) {
        if (!genotype->homozygous())
            continue;
        map<Genotype*, bool>::iterator g = genotypesWithHomozygousCombos.find(&*genotype);
        if (g == genotypesWithHomozygousCombos.end()) {
            for (vector<pair<string, vector<pair<Genotype*, long double> > > >::const_iterator s = sampleGenotypes.begin();
                    s != sampleGenotypes.end(); ++s) {
                for (vector<pair<Genotype*, long double> >::const_iterator d = s->second.begin(); d != s->second.end(); ++d) {
                    if (d->first == &*genotype) {
                        homozygousCombos[d->first].push_back(make_pair(s->first, make_pair(d->first, d->second)));
                    }
                }
            }
        }
    }

    for (map<Genotype*, GenotypeCombo>::iterator c = homozygousCombos.begin(); c != homozygousCombos.end(); ++c) {
        combos.push_back(c->second);
    }
    
    return combos;

}

vector<GenotypeCombo>
bandedGenotypeCombinationsIncludingBestHomozygousCombo(
        vector<pair<string, vector<pair<Genotype*, long double> > > >& sampleGenotypes,
        int bandwidth, int banddepth) {

    vector<GenotypeCombo> combos = bandedGenotypeCombinations(sampleGenotypes, bandwidth, banddepth);
    // is there already a homozygous combo?
    bool hasHomozygousCombo = false;
    for (vector<GenotypeCombo>::iterator c = combos.begin(); c != combos.end(); ++c) {
        bool allhomozygous = true;
        for (GenotypeCombo::iterator gc = c->begin(); gc != c->end(); ++gc) {
            if (!gc->second.first->homozygous()) {
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
        for(vector<pair<string, vector<pair<Genotype*, long double> > > >::iterator s = sampleGenotypes.begin();
                s != sampleGenotypes.end(); ++s) {
            for (vector<pair<Genotype*, long double> >::iterator g = s->second.begin(); g != s->second.end(); ++g) {
                if (g->first->homozygous()) {
                    homozygousCombo.push_back(make_pair(s->first, *g));
                    break;
                }
            }
        }
        combos.push_back(homozygousCombo);
    }

    return combos;

}

bool isHomozygousCombo(GenotypeCombo& combo) {
    GenotypeCombo::iterator g = combo.begin();
    Genotype* genotype = g->second.first;
    if (!genotype->homozygous())
        return false;
    for (; g != combo.end(); ++g) {
        if (g->second.first != genotype)
            return false;
    }
    return true;
}

pair<int, int> alternateAndReferenceCount(vector<Allele*>& observations, string& refbase, string altbase) {
    int altcount = 0;
    int refcount = 0;
    for (vector<Allele*>::iterator allele = observations.begin(); allele != observations.end(); ++allele) {
        if ((*allele)->base() == refbase)
            ++refcount;
        else if ((*allele)->base() == altbase)
            ++altcount;
    }
    return make_pair(altcount, refcount);
}


GenotypeComboMap genotypeCombo2Map(GenotypeCombo& gc) {
    GenotypeComboMap gcm;
    for (GenotypeCombo::iterator g = gc.begin(); g != gc.end(); ++g) {
        gcm[g->first] = g->second;
    }
    return gcm;
}
