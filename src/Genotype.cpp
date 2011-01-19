#include "Genotype.h"
#include "multichoose.h"
#include "multipermute.h"


vector<Allele> Genotype::uniqueAlleles(void) {
    vector<Allele> uniques;
    for (Genotype::iterator g = this->begin(); g != this->end(); ++g)
        uniques.push_back(g->allele);
    return uniques;
}

int Genotype::getPloidy(void) {
    int result = 0;
    for (Genotype::const_iterator i = this->begin(); i != this->end(); ++i) {
        result += i->count;
    }
    return result;
}

vector<int> Genotype::alleleCountsInObservations(Sample& sample) {
    vector<int> counts;
    for (Genotype::iterator i = this->begin(); i != this->end(); ++i) {
        int count = 0;
        Allele& b = i->allele;
        counts.push_back(sample.observationCount(b));
    }
    return counts;
}

vector<int> Genotype::counts(void) {
    vector<int> counts;
    for (Genotype::iterator i = this->begin(); i != this->end(); ++i) {
        counts.push_back(i->count);
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
        Allele& b = i->allele;
        if (base != b.currentBase)
            alleles.push_back(b);
    }
    return alleles;
}

int Genotype::alleleFrequency(const string& base) {
    map<string, int>::iterator ge = alleleCounts.find(base);
    if (ge == alleleCounts.end()) {
        return 0;
    } else {
        return ge->second;
    }
}

int Genotype::alleleFrequency(Allele& allele) {
    map<string, int>::iterator ge = alleleCounts.find(allele.currentBase);
    if (ge == alleleCounts.end()) {
        return 0;
    } else {
        return ge->second;
    }
}

string Genotype::relativeGenotype(string& refbase, string& altbase) {
    vector<string> rg;
    for (Genotype::iterator i = this->begin(); i != this->end(); ++i) {
        Allele& b = i->allele;
        if (b.currentBase == altbase && refbase != b.currentBase) {
            for (int j = 0; j < i->count; ++j)
                rg.push_back("1/");
        } else {
            for (int j = 0; j < i->count; ++j)
                rg.push_back("0/");
        }
    }
    sort(rg.begin(), rg.end()); // enforces the same ordering for all genotypes
    reverse(rg.begin(), rg.end()); // 1/0 ordering, or 1/1/0 etc.
    string result = accumulate(rg.begin(), rg.end(), string(""));
    return result.substr(0, result.size() - 1); // chop trailing '/'
}

bool Genotype::containsAllele(const string& base) {
    map<string, int>::iterator ge = alleleCounts.find(base);
    if (ge == alleleCounts.end()) {
        return false;
    } else {
        return true;
    }
}

bool Genotype::containsAllele(Allele& allele) {
    map<string, int>::iterator ge = alleleCounts.find(allele.currentBase);
    if (ge == alleleCounts.end()) {
        return false;
    } else {
        return true;
    }
}

bool Genotype::isHomozygous(void) {
    return this->size() == 1;
}

// the probability of drawing each allele out of the genotype, ordered by allele
vector<long double> Genotype::alleleProbabilities(void) {
    vector<long double> probs;
    for (vector<GenotypeElement>::const_iterator a = this->begin(); a != this->end(); ++a) {
        probs.push_back((long double) a->count / (long double) ploidy);
    }
    return probs;
}

string Genotype::str(void) {
    string s;
    for (Genotype::const_iterator ge = this->begin(); ge != this->end(); ++ge) {
        for (int i = 0; i < ge->count; ++i)
            s += ge->allele.alternateSequence;
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

string IUPAC2GenotypeStr(string iupac) {
    if (iupac == "A") return "AA";
    if (iupac == "M") return "AC";
    if (iupac == "R") return "AG";
    if (iupac == "W") return "AT";
    if (iupac == "C") return "CC";
    if (iupac == "S") return "CG";
    if (iupac == "Y") return "CT";
    if (iupac == "G") return "GG";
    if (iupac == "K") return "GT";
    if (iupac == "T") return "TT";
    return iupac;
}

ostream& operator<<(ostream& out, const GenotypeElement& rhs) {
    out << rhs.allele.alternateSequence << rhs.count;
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
    out << "{\"" << g.front().sampleName << "\":[\"" << *(g.front().genotype) << "\"," << exp(g.front().prob) << "]";
    for (;i != g.end(); ++i) {
        out << ", \"" << i->sampleName << "\":[\"" << *(i->genotype) << "\"," << exp(i->prob) << "]";
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
        if (ai->allele != bi->allele)
            return ai->allele < bi->allele;
        else if (ai->count != bi->count)
            return ai->count < bi->count;
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


int GenotypeCombo::numberOfAlleles(void) {
    int count = 0;
    for (GenotypeCombo::iterator g = this->begin(); g != this->end(); ++g) {
        count += g->genotype->ploidy;
    }
    return count;
}

void GenotypeCombo::initAlleleFrequencies(void) {
    alleleFrequencies = countAlleles();
}

int GenotypeCombo::alleleFrequency(Allele& allele) {
    map<string, int>::iterator f = alleleFrequencies.find(allele.currentBase);
    if (f == alleleFrequencies.end()) {
        return 0;
    } else {
        return f->second;
    }
}


void GenotypeCombo::updateAlleleFrequencies(Genotype* oldGenotype, Genotype* newGenotype) {
    // remove allele frequency information for old genotype
    for (Genotype::iterator g = oldGenotype->begin(); g != oldGenotype->end(); ++g) {
        GenotypeElement& ge = *g;
        alleleFrequencies[ge.allele.currentBase] -= ge.count;
    }
    // add allele frequency information for new genotype
    for (Genotype::iterator g = newGenotype->begin(); g != newGenotype->end(); ++g) {
        GenotypeElement& ge = *g;
        alleleFrequencies[ge.allele.currentBase] += ge.count;
    }
    // remove allele frequencies which are now 0 or below
    for (map<string, int>::iterator af = alleleFrequencies.begin(); af != alleleFrequencies.end(); ++af) {
        if (af->second <= 0) {
            alleleFrequencies.erase(af);
        }
    }
}

map<string, int> GenotypeCombo::countAlleles(void) {
    map<string, int> alleleCounts;
    for (GenotypeCombo::iterator g = this->begin(); g != this->end(); ++g) {
        for (Genotype::iterator a = g->genotype->begin(); a != g->genotype->end(); ++a) {
            map<string, int>::iterator c = alleleCounts.find(a->allele.currentBase);
            if (c != alleleCounts.end()) {
                c->second += a->count;
            } else {
                alleleCounts.insert(make_pair(a->allele.currentBase, a->count));
            }
        }
    }
    return alleleCounts;
}

map<int, int> GenotypeCombo::countFrequencies(void) {
    map<int, int> frequencyCounts;
    for (map<string, int>::iterator a = alleleFrequencies.begin(); a != alleleFrequencies.end(); ++a) {
        map<int, int>::iterator c = frequencyCounts.find(a->second);
        if (c != frequencyCounts.end()) {
            c->second += 1;
        } else {
            frequencyCounts[a->second] = 1;
        }
    }
    return frequencyCounts;
}

vector<int> GenotypeCombo::counts(void) {
    map<string, int> alleleCounts = countAlleles();
    vector<int> counts;
    for (map<string, int>::iterator a = alleleCounts.begin(); a != alleleCounts.end(); ++a) {
        counts.push_back(a->second);
    }
    return counts;
}

// returns true if the combination is 100% homozygous
bool GenotypeCombo::isHomozygous(void) {
    GenotypeCombo::iterator g = begin();
    Genotype* genotype = g->genotype;
    if (!genotype->homozygous) {
        return false;
    } else {
        Allele& allele = g->genotype->front().allele;
        for (; g != end(); ++g) {
            if (!g->genotype->homozygous || g->genotype->front().allele != allele)
                return false;
        }
        return true;
    }
}

void
bandedGenotypeCombinations(
    vector<GenotypeCombo>& combos,
    SampleGenotypesAndProbs& sampleGenotypes,
    int bandwidth, int banddepth) {

    int nsamples = sampleGenotypes.size();

    // generate the best genotype combination according to data likelihoods
    GenotypeCombo comboKing;
    for (SampleGenotypesAndProbs::const_iterator s = sampleGenotypes.begin();
            s != sampleGenotypes.end(); ++s) {
        const pair<Genotype*, long double>& p = s->second.at(0);
        comboKing.push_back(SampleGenotypeProb(s->first, p.first, p.second));
        comboKing.prob += p.second;
    }
    comboKing.initAlleleFrequencies();

    // overview:
    //
    // For each order of indexes in the bandwidth and banddepth, Obtain all
    // multiset permutations of a set of indexes.  Then use these indexes to
    // get the nth-best genotype from each individual's set of genotypes for
    // which we have data likelihoods (sampleGenotypes), and turn this set into
    // a genotype combination.  Update the combination probability inline here
    // so we don't incur O(N^2) penalty calculating the probability within our
    // genotypeCombinationPriors calculation loop, where we estimate the
    // posterior probability of each genotype combination given its data
    // likelihood and the prior probability of the distribution of alleles it
    // represents.
    // 
    // example (bandwidth = 2, banddepth = 1)
    //  indexes:      0 0 0 0 1, 0 0 0 1 1
    //
    //  permutations: 0 0 0 0 1 
    //                0 0 0 1 0
    //                0 0 1 0 0
    //                0 1 0 0 0
    //                1 0 0 0 0
    //                1 1 0 0 0 
    //                0 1 1 0 0 
    //                1 0 1 0 0 
    //                0 1 0 1 0 
    //                0 0 1 1 0 
    //                1 0 0 1 0 
    //                0 1 0 0 1 
    //                0 0 1 0 1 
    //                0 0 0 1 1 
    //                1 0 0 0 1 
    //
    // We then convert these permutation to genotype combinations by using the
    // index to pick the nth-best genotype according to sorted individual
    // genotype data likelihoods.
    //
    // In addition to this simple case, We can flexibly extend this to larger
    // search spaces by changing the depth and width of the deviations from the
    // data likelihood maximizer (aka 'king').
    //
    for (int i = 0; i <= bandwidth; ++i) {
        for (int j = 1; j <= banddepth; ++j) {
            vector<int> indexes;
            for (int h = 0; h < j; ++h)
                indexes.push_back(i);
            for (int h = 0; h < (nsamples - j); ++h)
                indexes.push_back(0);
            vector<vector<int> > indexPermutations = multipermute(indexes);
            for (vector<vector<int> >::const_iterator p = indexPermutations.begin(); p != indexPermutations.end(); ++p) {
                combos.push_back(comboKing); // copy the king, and then we'll modify it according to the indicies
                GenotypeCombo& combo = combos.back();
                bool useCombo = true;
                GenotypeCombo::iterator currentSampleGenotype = combo.begin();
                vector<int>::const_iterator n = p->begin();
                for (SampleGenotypesAndProbs::const_iterator s = sampleGenotypes.begin();
                        s != sampleGenotypes.end(); ++s, ++n, ++currentSampleGenotype) {
                    int offset = *n;
                    if (offset > 0) {
                        // ignore this combo if it's beyond the bounds of the individual's set of genotypes
                        if (offset >= s->second.size()) {
                            useCombo = false;
                            break;
                        }
                        const pair<Genotype*, long double>& p = s->second.at(offset);
                        // get the old and new genotypes, which we compare to
                        // change the cached counts and probability of the
                        // combo
                        Genotype* oldGenotype = currentSampleGenotype->genotype;
                        Genotype* newGenotype = p.first;
                        combo.updateAlleleFrequencies(oldGenotype, newGenotype);
                        // replace genotype with new genotype
                        currentSampleGenotype->genotype = newGenotype;
                        // update data likelihood sum for combo
                        long double dataLikelihood = p.second;
                        // find difference
                        long double diff = s->second.at(0).second - dataLikelihood;
                        // adjust combination total data likelihood
                        combo.prob -= diff;
                    }
                    //combo.push_back(SampleGenotypeProb(s->first, p.first, dataLikelihood));
                }
                if (!useCombo) {
                    combos.erase(combos.end() - 1);
                }
            }
        }
    }

}


void
bandedGenotypeCombinationsIncludingAllHomozygousCombos(
    vector<GenotypeCombo>& combos,
    SampleGenotypesAndProbs& sampleGenotypes,
    map<int, vector<Genotype> >& genotypesByPloidy,
    int bandwidth, int banddepth) {

    // obtain the combos

    bandedGenotypeCombinations(combos, sampleGenotypes, bandwidth, banddepth);

    // determine which homozygous combos we already have

    map<Allele, bool> allelesWithHomozygousCombos;

    for (vector<GenotypeCombo>::iterator c = combos.begin(); c != combos.end(); ++c) {
        bool allSameAndHomozygous = true;
        GenotypeCombo::iterator gc = c->begin();
        Genotype* genotype;
        if (gc->genotype->homozygous) {
            genotype = gc->genotype;
        } else {
            continue;
        }
        for (; gc != c->end(); ++gc) {
            if (! (gc->genotype == genotype) ) {
                allSameAndHomozygous = false;
                break;
            }
        }
        if (allSameAndHomozygous) {
            allelesWithHomozygousCombos[genotype->front().allele] == true;
        }
    }

    // accumulate the needed homozygous combos

    map<Allele, GenotypeCombo> homozygousCombos;

    for (map<int, vector<Genotype> >::iterator pg = genotypesByPloidy.begin(); pg != genotypesByPloidy.end(); ++pg) {
        vector<Genotype>& genotypes = pg->second;
        for (vector<Genotype>::iterator g = genotypes.begin(); g != genotypes.end(); ++g) {
            Genotype* genotype = &*g;
            if (!genotype->homozygous)
                continue;
            Allele& allele = genotype->front().allele;
            map<Allele, bool>::iterator g = allelesWithHomozygousCombos.find(allele);
            if (g == allelesWithHomozygousCombos.end()) {
                // we need to make a new combo
                // iterate through the sample genotype vector
                GenotypeCombo& combo = homozygousCombos[allele];
                for (vector<pair<string, vector<pair<Genotype*, long double> > > >::const_iterator s = sampleGenotypes.begin();
                        s != sampleGenotypes.end(); ++s) {
                    // for each sample genotype, if the genotype is the same as our currently needed genotype, push it back onto a new combo
                    for (vector<pair<Genotype*, long double> >::const_iterator d = s->second.begin(); d != s->second.end(); ++d) {
                        // this check is ploidy-independent
                        if (d->first->homozygous && d->first->front().allele == allele) {
                            combo.push_back(SampleGenotypeProb(s->first, d->first, d->second));
                        }
                    }
                }
            }
        }
    }

    // accumulate homozygous combos and set their combo data probabilities
    for (map<Allele, GenotypeCombo>::iterator c = homozygousCombos.begin(); c != homozygousCombos.end(); ++c) {
        GenotypeCombo& gc = c->second;
        gc.prob = 0;
        for (GenotypeCombo::iterator sgp = gc.begin(); sgp != gc.end(); ++sgp) {
            gc.prob += sgp->prob; // set up data likelihood for combo
        }
        gc.initAlleleFrequencies();  // cache allele frequency information
        combos.push_back(gc);
    }

}

void
bandedGenotypeCombinationsIncludingBestHomozygousCombo(
    vector<GenotypeCombo>& combos,
    SampleGenotypesAndProbs& sampleGenotypes,
    int bandwidth, int banddepth) {

    bandedGenotypeCombinations(combos, sampleGenotypes, bandwidth, banddepth);
    // is there already a homozygous combo?
    bool hasHomozygousCombo = false;
    for (vector<GenotypeCombo>::iterator c = combos.begin(); c != combos.end(); ++c) {
        bool allhomozygous = true;
        for (GenotypeCombo::iterator gc = c->begin(); gc != c->end(); ++gc) {
            if (!gc->genotype->homozygous) {
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
                if (g->first->homozygous) {
                    homozygousCombo.push_back(SampleGenotypeProb(s->first, g->first, g->second));
                    break;
                }
            }
        }
        combos.push_back(homozygousCombo);
    }

}

pair<int, int> alternateAndReferenceCount(vector<Allele*>& observations, string& refbase, string altbase) {
    int altcount = 0;
    int refcount = 0;
    for (vector<Allele*>::iterator allele = observations.begin(); allele != observations.end(); ++allele) {
        if ((*allele)->currentBase == refbase)
            ++refcount;
        else if ((*allele)->currentBase == altbase)
            ++altcount;
    }
    return make_pair(altcount, refcount);
}


void genotypeCombo2Map(GenotypeCombo& gc, GenotypeComboMap& gcm) {
    for (GenotypeCombo::iterator g = gc.begin(); g != gc.end(); ++g) {
        gcm[g->sampleName] = make_pair(g->genotype, g->prob);
    }
}


// returns a list of the alternate alleles represented by the given genotype
// combo sorted by frequency
vector<pair<Allele, int> > alternateAlleles(GenotypeCombo& combo, string referenceBase) {

    map<Allele, int> alternates;

    for (GenotypeCombo::iterator g = combo.begin(); g != combo.end(); ++g) {
        vector<Allele> alts = g->genotype->alternateAlleles(referenceBase);
        for (vector<Allele>::iterator a = alts.begin(); a != alts.end(); ++a) {
            if (alternates.find(*a) == alternates.end()) {
                alternates[*a] = 1;
            } else {
                alternates[*a] += 1;
            }
        }
    }

    vector<pair<Allele, int> > sortedAlternates;

    for (map<Allele, int>::iterator a = alternates.begin(); a != alternates.end(); ++a) {
        sortedAlternates.push_back(make_pair(a->first, a->second));
    }

    AllelePairIntCompare alleleCountCompare;
    sort(sortedAlternates.begin(), sortedAlternates.end(), alleleCountCompare);

    return sortedAlternates;

}

int Genotype::containedAlleleTypes(void) {
    int t = 0;
    for (Genotype::iterator g = begin(); g != end(); ++g) {
        t |= g->allele.type;
    }
    return t;
}
