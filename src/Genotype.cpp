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

vector<int> Genotype::counts(void) {
    vector<int> counts;
    for (Genotype::iterator i = this->begin(); i != this->end(); ++i) {
        counts.push_back(i->count);
    }
    return counts;
}

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
        } else if (b.currentBase != altbase && refbase != b.currentBase) {
            for (int j = 0; j < i->count; ++j)
                rg.push_back("./");
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

string IUPAC2GenotypeStr(string iupac, int ploidy) {
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
    for (int i = 0; i < rhs.count; ++i)
        out << rhs.allele.alternateSequence;
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
    out << "{\"" << g.front()->name << "\":[\"" << *(g.front()->genotype) << "\"," << exp(g.front()->prob) << "]";
    for (;i != g.end(); ++i) {
        out << ", \"" << (*i)->name << "\":[\"" << *((*i)->genotype) << "\"," << exp((*i)->prob) << "]";
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
    for (map<string, int>::iterator f = alleleFrequencies.begin(); f != alleleFrequencies.end(); ++f) {
        count += f->second;
    }
    return count;
}

// initializes cached counts associated with each GenotypeCombo
void GenotypeCombo::init(bool useObsExpectations) {
    for (GenotypeCombo::iterator s = begin(); s != end(); ++s) {
        const SampleDataLikelihood& sdl = **s;
        const Sample& sample = *sdl.sample;
        for (Genotype::iterator a = sdl.genotype->begin(); a != sdl.genotype->end(); ++a) {
            const string& alleleBase = a->allele.currentBase;

            // allele frequencies in selected genotypes in combo
            map<string, int>::iterator c = alleleFrequencies.find(alleleBase);
            if (c != alleleFrequencies.end()) {
                c->second += a->count;
            } else {
                alleleFrequencies.insert(make_pair(alleleBase, a->count));
            }

            if (useObsExpectations) {
                // observational frequencies for binomial priors
                Sample::const_iterator as = sample.find(alleleBase);
                if (as != sample.end()) {
                    vector<Allele*> alleles = as->second;
                    for (vector<Allele*>::iterator o = alleles.begin(); o != alleles.end(); ++o) {
                        const Allele& allele = **o;
                        if (allele.basesLeft() >= allele.basesRight()) {
                            alleleReadPlacementCounts[alleleBase].first += 1;
                        } else {
                            alleleReadPlacementCounts[alleleBase].second += 1;
                        }
                        if (allele.strand == STRAND_FORWARD) {
                            alleleStrandCounts[alleleBase].first += 1;
                        } else {
                            alleleStrandCounts[alleleBase].second += 1;
                        }
                    }
                }
            }
        }
    }
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

void GenotypeCombo::updateCachedCounts(
        Sample* sample,
        Genotype* oldGenotype,
        Genotype* newGenotype,
        bool useObsExpectations) {

    // TODO can we improve efficiency by only adjusting for bases which are actually changed

    // remove allele frequency information for old genotype
    for (Genotype::iterator g = oldGenotype->begin(); g != oldGenotype->end(); ++g) {
        GenotypeElement& ge = *g;
        const string& base = ge.allele.currentBase;
        alleleFrequencies[base] -= ge.count;
        if (useObsExpectations) {
            Sample::iterator s = sample->find(base);
            if (s != sample->end()) {
                const vector<Allele*>& alleles = s->second;
                int forward_strand = 0;
                int reverse_strand = 0;
                int placed_left = 0;
                int placed_right = 0;
                for (vector<Allele*>::const_iterator a = alleles.begin(); a != alleles.end(); ++a) {
                    const Allele& allele = **a;
                    if (allele.strand == STRAND_FORWARD) {
                        ++forward_strand;
                    } else {
                        ++reverse_strand;
                    }
                    if (allele.basesLeft() >= allele.basesRight()) {
                        ++placed_left;
                    } else {
                        ++placed_right;
                    }
                }
                alleleStrandCounts[base].first -= forward_strand;
                alleleStrandCounts[base].second -= forward_strand;
                alleleReadPlacementCounts[base].first -= placed_left;
                alleleReadPlacementCounts[base].second -= placed_right;
            }
        }
    }

    // add allele frequency information for new genotype
    for (Genotype::iterator g = newGenotype->begin(); g != newGenotype->end(); ++g) {
        GenotypeElement& ge = *g;
        const string& base = ge.allele.currentBase;
        alleleFrequencies[base] += ge.count;
        if (useObsExpectations) {
            Sample::iterator s = sample->find(base);
            if (s != sample->end()) {
                const vector<Allele*>& alleles = s->second;
                int forward_strand = 0;
                int reverse_strand = 0;
                int placed_left = 0;
                int placed_right = 0;
                for (vector<Allele*>::const_iterator a = alleles.begin(); a != alleles.end(); ++a) {
                    const Allele& allele = **a;
                    if (allele.strand == STRAND_FORWARD) {
                        ++forward_strand;
                    } else {
                        ++reverse_strand;
                    }
                    if (allele.basesLeft() >= allele.basesRight()) {
                        ++placed_left;
                    } else {
                        ++placed_right;
                    }
                }
                alleleStrandCounts[base].first += forward_strand;
                alleleStrandCounts[base].second += forward_strand;
                alleleReadPlacementCounts[base].first += placed_left;
                alleleReadPlacementCounts[base].second += placed_right;
            }
        }
    }

    // remove allele frequencies which are now 0 or below
    for (map<string, int>::iterator af = alleleFrequencies.begin();
            af != alleleFrequencies.end(); ++af) {
        if (af->second <= 0) {
            alleleFrequencies.erase(af);
        }
    }

    if (useObsExpectations) {

        for (map<string, pair<int, int> >::iterator as = alleleStrandCounts.begin();
                as != alleleStrandCounts.end(); ++as) {
            if (as->second.first <= 0) {
                as->second.first = 0;
            }
            if (as->second.second <= 0) {
                as->second.second = 0;
            }
            if (as->second.first == 0 && as->second.second == 0) {
                alleleStrandCounts.erase(as);
            }
        }

        for (map<string, pair<int, int> >::iterator ap = alleleReadPlacementCounts.begin();
                ap != alleleReadPlacementCounts.end(); ++ap) {
            if (ap->second.first <= 0) {
                ap->second.first = 0;
            }
            if (ap->second.second <= 0) {
                ap->second.second = 0;
            }
            if (ap->second.first == 0 && ap->second.second == 0) {
                alleleReadPlacementCounts.erase(ap);
            }
        }

    }

}

map<string, int> GenotypeCombo::countAlleles(void) {
    map<string, int> alleleCounts;
    for (GenotypeCombo::iterator g = this->begin(); g != this->end(); ++g) {
        SampleDataLikelihood& sdl = **g;
        for (Genotype::iterator a = sdl.genotype->begin(); a != sdl.genotype->end(); ++a) {
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
    //map<string, int> alleleCounts = countAlleles();
    vector<int> counts;
    for (map<string, int>::iterator a = alleleFrequencies.begin(); a != alleleFrequencies.end(); ++a) {
        counts.push_back(a->second);
    }
    return counts;
}

// how many copies of the locus are in the whole genotype combination?
int GenotypeCombo::ploidy(void) {
    int copies;
    for (map<string, int>::iterator a = alleleFrequencies.begin(); a != alleleFrequencies.end(); ++a) {
        copies += a->second;
    }
    return copies;
}

vector<long double> GenotypeCombo::alleleProbs(void) {
    vector<long double> probs;
    long double copies = ploidy();
    for (map<string,int>::iterator a = alleleFrequencies.begin(); a != alleleFrequencies.end(); ++a) {
        probs.push_back(a->second / copies);
    }
    return probs;
}

vector<string> GenotypeCombo::alleles(void) {
    vector<string> bases;
    for (map<string, int>::iterator a = alleleFrequencies.begin(); a != alleleFrequencies.end(); ++a) {
        bases.push_back(a->first);
    }
    return bases;
}

// returns true if the combination is 100% homozygous
bool GenotypeCombo::isHomozygous(void) {
    GenotypeCombo::iterator g = begin();
    Genotype* genotype = (*g)->genotype;
    if (!genotype->homozygous) {
        return false;
    } else {
        Allele& allele = (*g)->genotype->front().allele;
        for (; g != end(); ++g) {
            if (!(*g)->genotype->homozygous || (*g)->genotype->front().allele != allele)
                return false;
        }
        return true;
    }
}

void
bandedGenotypeCombinations(
    vector<GenotypeCombo>& combos,
    SampleDataLikelihoods& sampleDataLikelihoods,
    Samples& samples,
    bool useObsExpectations,
    int bandwidth, int banddepth,
    float logStepMax) {

    int nsamples = sampleDataLikelihoods.size();

    // generate the best genotype combination according to data likelihoods
    GenotypeCombo comboKing;
    for (SampleDataLikelihoods::iterator s = sampleDataLikelihoods.begin();
            s != sampleDataLikelihoods.end(); ++s) {
        SampleDataLikelihood* sdl = &s->front();
        comboKing.push_back(sdl);
        comboKing.prob += sdl->prob;
    }
    comboKing.init(useObsExpectations);

    // overview:
    //
    // For each order of indexes in the bandwidth and banddepth, Obtain all
    // multiset permutations of a set of indexes.  Then use these indexes to
    // get the nth-best genotype from each individual's set of genotypes for
    // which we have data likelihoods (sampleDataLikelihoods), and turn this set into
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
            bool reuseLastCombo = false;
            for (vector<vector<int> >::const_iterator p = indexPermutations.begin(); p != indexPermutations.end(); ++p) {
                if (reuseLastCombo) {  // reuse the last combo if we've skipped it, saving a few % runtime copying combos
                    reuseLastCombo = false;
                } else {
                    combos.push_back(comboKing); // copy the king, and then we'll modify it according to the indicies
                }
                GenotypeCombo& combo = combos.back();
                GenotypeCombo::iterator currentSampleGenotypeItr = combo.begin();
                vector<int>::const_iterator n = p->begin();
                for (SampleDataLikelihoods::iterator s = sampleDataLikelihoods.begin();
                        s != sampleDataLikelihoods.end(); ++s, ++n, ++currentSampleGenotypeItr) {
                    SampleDataLikelihood& oldsdl = **currentSampleGenotypeItr;
                    SampleDataLikelihood*& oldsdl_ptr = *currentSampleGenotypeItr;
                    vector<SampleDataLikelihood>& sdls = *s;
                    int offset = *n;
                    if (offset > 0) {
                        // ignore this combo if it's beyond the bounds of the individual's set of genotypes
                        if (offset >= s->size()) {
                            reuseLastCombo = true;
                            break;
                        }
                        // ignore this combo if the swapped genotype has a data likelihood more than logStepMax
                        // from the best data likelihood.  logStepMax == -1 indicates no filtering
                        if (offset > 0 && logStepMax >= 0 && sdls.front().prob - sdls.at(offset).prob > logStepMax) {
                            reuseLastCombo = true;
                            break;
                        }
                        SampleDataLikelihood* newsdl = &sdls.at(offset);
                        // get the old and new genotypes, which we compare to
                        // change the cached counts and probability of the
                        // combo
                        combo.updateCachedCounts(oldsdl.sample, oldsdl.genotype, newsdl->genotype, useObsExpectations);
                        // replace genotype with new genotype
                        oldsdl_ptr = newsdl;
                        // find data likelihood difference from ComboKing
                        long double diff = oldsdl.prob - newsdl->prob;
                        // adjust combination total data likelihood
                        combo.prob -= diff;
                    }
                }
            }
            if (reuseLastCombo) {
                combos.erase(combos.end() - 1);
            }
        }
    }
}


void
bandedGenotypeCombinationsIncludingAllHomozygousCombos(
    vector<GenotypeCombo>& combos,
    SampleDataLikelihoods& sampleDataLikelihoods,
    Samples& samples,
    bool useObsExpectations,
    map<int, vector<Genotype> >& genotypesByPloidy,
    vector<Allele>& genotypeAlleles,
    int bandwidth, int banddepth,
    float logStepMax) {

    // obtain the combos

    bandedGenotypeCombinations(combos, sampleDataLikelihoods, samples, useObsExpectations, bandwidth, banddepth, logStepMax);

    // determine which homozygous combos we already have

    map<Allele, bool> allelesWithHomozygousCombos;

    for (vector<GenotypeCombo>::iterator c = combos.begin(); c != combos.end(); ++c) {
        bool allSameAndHomozygous = true;
        GenotypeCombo::iterator gc = c->begin();
        Genotype* genotype;
        if ((*gc)->genotype->homozygous) {
            genotype = (*gc)->genotype;
        } else {
            continue;
        }
        for (; gc != c->end(); ++gc) {
            if (! ((*gc)->genotype == genotype) ) {
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

    for (vector<Allele>::iterator a = genotypeAlleles.begin(); a != genotypeAlleles.end(); ++a) {
        Allele& allele = *a;
        map<Allele, bool>::iterator g = allelesWithHomozygousCombos.find(allele);
        if (g == allelesWithHomozygousCombos.end()) {
            // we need to make a new combo
            // iterate through the sample genotype vector
            GenotypeCombo& combo = homozygousCombos[allele];
            for (SampleDataLikelihoods::iterator s = sampleDataLikelihoods.begin();
                    s != sampleDataLikelihoods.end(); ++s) {
                // for each sample genotype, if the genotype is the same as our currently needed genotype, push it back onto a new combo
                for (vector<SampleDataLikelihood>::iterator d = s->begin(); d != s->end(); ++d) {
                    SampleDataLikelihood& sdl = *d;
                    // this check is ploidy-independent
                    if (sdl.genotype->homozygous && sdl.genotype->front().allele == allele) {
                        combo.push_back(&sdl);
                        break;
                    }
                }
            }
        }
    }

    // accumulate homozygous combos and set their combo data probabilities
    for (map<Allele, GenotypeCombo>::iterator c = homozygousCombos.begin(); c != homozygousCombos.end(); ++c) {
        GenotypeCombo& gc = c->second;
        gc.prob = 0;
        for (GenotypeCombo::iterator sdl = gc.begin(); sdl != gc.end(); ++sdl) {
            gc.prob += (*sdl)->prob; // set up data likelihood for combo
        }
        gc.init(useObsExpectations);  // cache allele frequency information
        combos.push_back(gc);
    }

}

void
bandedGenotypeCombinationsIncludingBestHomozygousCombo(
    vector<GenotypeCombo>& combos,
    SampleDataLikelihoods& sampleDataLikelihoods,
    Samples& samples,
    bool useObsExpectations,
    int bandwidth, int banddepth, float logStepMax) {

    bandedGenotypeCombinations(combos, sampleDataLikelihoods, samples, useObsExpectations, bandwidth, banddepth, logStepMax);
    // is there already a homozygous combo?
    bool hasHomozygousCombo = false;
    for (vector<GenotypeCombo>::iterator c = combos.begin(); c != combos.end(); ++c) {
        bool allhomozygous = true;
        for (GenotypeCombo::iterator gc = c->begin(); gc != c->end(); ++gc) {
            if (!(*gc)->genotype->homozygous) {
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
        for(SampleDataLikelihoods::iterator s = sampleDataLikelihoods.begin();
                s != sampleDataLikelihoods.end(); ++s) {
            for (vector<SampleDataLikelihood>::iterator g = s->begin(); g != s->end(); ++g) {
                if (g->genotype->homozygous) {
                    homozygousCombo.push_back(&*g);
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
        gcm[(*g)->name] = *g;;
    }
}


// returns a list of the alternate alleles represented by the given genotype
// combo sorted by frequency
vector<pair<Allele, int> > alternateAlleles(GenotypeCombo& combo, string referenceBase) {

    map<Allele, int> alternates;

    for (GenotypeCombo::iterator g = combo.begin(); g != combo.end(); ++g) {
        vector<Allele> alts = (*g)->genotype->alternateAlleles(referenceBase);
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


vector<int> Genotype::alleleObservationCounts(Sample& sample) {
    vector<int> counts;
    for (Genotype::iterator i = begin(); i != end(); ++i) {
        int count = 0;
        Allele& b = i->allele;
        counts.push_back(sample.observationCount(b));
    }
    return counts;
}
