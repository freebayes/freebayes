#include "Allele.h"
#include "multichoose.h"
#include "TryCatch.h"


int Allele::referenceOffset(void) const {
    /*cout << readID << " offset checked " << referencePosition - position << " against position " << position 
        << " allele length " << length << " str length " << referenceSequence.size() << " qstr size " << qualityString.size() << endl;
        */
    return *currentReferencePosition - position;
}

// called prior to using the allele in analysis
void Allele::update(void) {
    currentBase = base();
    quality = currentQuality();
    lnquality = phred2ln(quality);
    if (type == ALLELE_REFERENCE) {
        basesLeft = bpLeft + referenceOffset();
        basesRight = bpRight - referenceOffset();
    }
}

// quality of subsequence of allele
const int Allele::subquality(int startpos, int len) const {
    int start = startpos - position;
    int sum = 0;
    for (int i = start; i < len; ++i) {
        sum += baseQualities.at(i);
    }
    return sum;
}

// quality of subsequence of allele
const long double Allele::lnsubquality(int startpos, int len) const {
    return phred2ln(subquality(startpos, len));
}

const int Allele::subquality(const Allele &a) const {
    int sum = 0;
    int rp = a.position - position;
    int l = a.length;
    int L = l;
    int spanstart = 0;
    int spanend = 1;
    //int l = min((int) a.length, (int) baseQualities.size() - start);
    if (a.type == ALLELE_INSERTION) {
        L = l + 2;
        if (L > baseQualities.size()) {
            L = baseQualities.size();
            spanstart = 0;
        } else {
            // set lower bound to 0
            if (rp < (L / 2)) {
                spanstart = 0;
            } else {
                spanstart = rp - (L / 2);
            }
            // set upper bound to the string length
            if (spanstart + L > baseQualities.size()) {
                spanstart = baseQualities.size() - L;
            }
        }
        //string qualstr = baseQualities.substr(spanstart, L);
        spanend = spanstart + L;
    } else if (a.type == ALLELE_DELETION) {
        L = l + 2;
        if (L > baseQualities.size()) {
            L = baseQualities.size();
            spanstart = 0;
        } else {
            // set lower bound to 0
            if (rp < 1) {
                spanstart = 0;
            } else {
                spanstart = rp - 1;
            }
            // set upper bound to the string length
            if (spanstart + L > baseQualities.size()) {
                spanstart = baseQualities.size() - L;
            }
        }
        spanend = spanstart + L;
    } else if (a.type == ALLELE_MNP) {
        L = l;
        if (L > baseQualities.size()) {
            L = baseQualities.size();
            spanstart = 0;
        } else {
            if (rp < 1) {
                spanstart = 0;
            } else {
                spanstart = rp;
            }
            // impossible
            if (spanstart + L > baseQualities.size()) {
                spanstart = baseQualities.size() - L;
            }
        }
        spanend = spanstart + L;
    }
    // TODO ALLELE_COMPLEX
    for (int i = spanstart; i < spanend; ++i) {
        sum += baseQualities.at(i);
    }
    return sum * (l / L);
}

const long double Allele::lnsubquality(const Allele& a) const {
    return phred2ln(subquality(a));
}

void updateAllelesCachedData(vector<Allele*>& alleles) {
    for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        (*a)->update();
    }
}

/*
const int Allele::basesLeft(void) const {
    if (type == ALLELE_REFERENCE) {
        return bpLeft + referenceOffset();
    } else {
        return bpLeft;
    }
}

const int Allele::basesRight(void) const {
    if (type == ALLELE_REFERENCE) {
        return bpRight - referenceOffset();
    } else {
        return bpRight;
    }
}
*/

// quality at a given reference position
const short Allele::currentQuality(void) const {
    switch (this->type) {
        case ALLELE_REFERENCE:
            TRY {
            return baseQualities.at(referenceOffset());
            } CATCH;
            break;
        case ALLELE_INSERTION:
        case ALLELE_DELETION:
        case ALLELE_SNP:
        case ALLELE_MNP:
        case ALLELE_COMPLEX:
            return quality;
            break;
    }
}

const long double Allele::lncurrentQuality(void) const {
    return phred2ln(currentQuality());
}

string Allele::typeStr(void) {

    string t;

    switch (this->type) {
        case ALLELE_GENOTYPE:
            t = "genotype";
            break;
        case ALLELE_REFERENCE:
            t = "reference";
            break;
        case ALLELE_MNP:
            t = "mnp";
            break;
        case ALLELE_SNP:
            t = "snp";
            break;
        case ALLELE_INSERTION:
            t = "insertion";
            break;
        case ALLELE_DELETION:
            t = "deletion";
            break;
        case ALLELE_COMPLEX:
            t = "complex";
            break;
        default:
            t = "unknown";
            break;
    }

    return t;

}

bool Allele::isReference(void) {
    return type == ALLELE_REFERENCE;
}

bool Allele::isSNP(void) {
    return type == ALLELE_SNP;
}

bool Allele::isInsertion(void) {
    return type == ALLELE_INSERTION;
}

bool Allele::isDeletion(void) {
    return type == ALLELE_DELETION;
}

bool Allele::isMNP(void) {
    return type == ALLELE_MNP;
}

bool Allele::isComplex(void) {
    return type == ALLELE_COMPLEX;
}

const string Allele::base(void) const { // the base of this allele

    switch (this->type) {
        case ALLELE_REFERENCE:
            if (genotypeAllele)
                return alternateSequence;
            else
                return currentReferenceBase;
            break;
        case ALLELE_GENOTYPE:
            return alternateSequence;
            break;
        case ALLELE_SNP:
            return alternateSequence;
            break;
        case ALLELE_MNP:
            return alternateSequence;
            break;
        case ALLELE_INSERTION:
            return "I" + alternateSequence;
            break;
        case ALLELE_DELETION:
            return "D" + convert(length);
            break;
        case ALLELE_COMPLEX:
            return "C" + alternateSequence;
            break;
        default:
            break;
    }

}

const bool Allele::masked(void) const {

    // guard against uninitialized indelMask
    if (indelMask.size() == 0)
        return false;

    if (genotypeAllele)
        return false;

    switch (this->type) {
        case ALLELE_GENOTYPE:
            return false;
            break;
        case ALLELE_REFERENCE:
            return indelMask.at(referenceOffset());
            break;
        case ALLELE_SNP:
            return indelMask.at(0);
            break;
        case ALLELE_INSERTION: // XXX presently these are masked by default...
        case ALLELE_DELETION:
        case ALLELE_MNP:
        case ALLELE_COMPLEX:
            return true;
            break;
        default:
            break;
    }

}

string stringForAllele(Allele &allele) {

    stringstream out;
    if (!allele.genotypeAllele) {
        out.precision(1);
        out 
            << allele.sampleID << "\t"
            << allele.readID << "\t"
            << allele.typeStr() << "\t" 
            << scientific << fixed << allele.position << "\t"
            << allele.length << "\t"
            << (allele.strand == STRAND_FORWARD ? "+" : "-") << "\t"
            << allele.referenceSequence << "\t"
            << allele.alternateSequence << "\t"
            << allele.quality << "\t" << endl;
    } else {
        out << allele.typeStr() << "\t"
            << allele.length << "\t"
            << allele.alternateSequence
            << endl;
    }

    return out.str();
}

string stringForAlleles(vector<Allele> &alleles) {
    stringstream out;
    for (vector<Allele>::iterator allele = alleles.begin(); allele != alleles.end(); allele++) {
        out << stringForAllele(*allele) << endl;
    }
    return out.str();
}

string json(vector<Allele*> &alleles) {
    stringstream out;
    vector<Allele*>::iterator a = alleles.begin();
    out << "[" << (*a)->json(); ++a;
    for (; a != alleles.end(); ++a)
        out << "," << (*a)->json();
    out << "]";
    return out.str();
}

string json(Allele*& allele) { return allele->json(); }

string Allele::json(void) {
    stringstream out;
    if (!genotypeAllele) {
        out << "{\"id\":\"" << readID << "\""
            << ",\"type\":\"" << typeStr() << "\""
            << ",\"length\":" << ((type == ALLELE_REFERENCE) ? 1 : length)
            << ",\"position\":" << position 
            << ",\"strand\":\"" << (strand == STRAND_FORWARD ? "+" : "-") << "\"";
        if (type == ALLELE_REFERENCE ) {
            out << ",\"base\":\"" << alternateSequence.at(referenceOffset()) << "\""
                //<< ",\"reference\":\"" << allele.referenceSequence.at(referenceOffset) << "\""
                << ",\"quality\":" << currentQuality();
        } else {
            out << ",\"base\":\"" << alternateSequence << "\""
                //<< ",\"reference\":\"" << allele.referenceSequence << "\""
                << ",\"quality\":" << quality;
        }
        out << "}";

    } else {
        out << "{\"type\":\"" << typeStr() << "\"";
        switch (type) {
            case ALLELE_REFERENCE:
                out << "}";
                break;
            default:
                out << "\",\"length\":" << length 
                    << ",\"alt\":\"" << alternateSequence << "\"}";
                break;
        }
    }
    return out.str();
}

ostream &operator<<(ostream &out, vector<Allele*> &alleles) {
    vector<Allele*>::iterator a = alleles.begin();
    out << **a++;
    while (a != alleles.end())
        out << "|" << **a++;
    return out;
}

ostream &operator<<(ostream &out, vector<Allele> &alleles) {
    vector<Allele>::iterator a = alleles.begin();
    out << *a++;
    while (a != alleles.end())
        out << "|" << *a++;
    return out;
}

ostream &operator<<(ostream &out, list<Allele*> &alleles) {
    list<Allele*>::iterator a = alleles.begin();
    out << **a++;
    while (a != alleles.end())
        out << "|" << **a++;
    return out;
}

ostream &operator<<(ostream &out, Allele* &allele) {
    out << *allele;
    return out;
}

ostream &operator<<(ostream &out, Allele &allele) {

    if (!allele.genotypeAllele) {
        // << &allele << ":" 
        out.precision(1);
        out << allele.sampleID
            << ":" << allele.readID 
            << ":" << allele.typeStr() 
            << ":" << allele.length 
            << ":" << scientific << fixed << allele.position 
            << ":" << (allele.strand == STRAND_FORWARD ? "+" : "-")
            << ":" << allele.alternateSequence
            //<< ":" << allele.referenceSequence
            << ":" << allele.quality;
    } else {
        out << allele.typeStr() 
            << ":" << allele.length 
            << ":" << (string) allele.alternateSequence;
    }
    return out;
}

bool operator<(const Allele &a, const Allele &b) {
    //cerr << "allele<" << endl;
    return a.currentBase < b.currentBase;
}

// TODO fixme??
// alleles are equal if they represent the same reference-relative variation or
// sequence, which we encode as a string and compare here
bool operator==(const Allele &a, const Allele &b) {
    //cerr << "allele==" << endl;
    return a.currentBase == b.currentBase;
}

bool operator!=(const Allele& a, const Allele& b) {
    return ! (a == b);
}

bool Allele::equivalent(Allele &b) {

    if (type != b.type) {
        return false;
    } else {
        switch (type) {
            case ALLELE_REFERENCE:
                // reference alleles are, by definition, always equivalent
                return true;
                break;
            case ALLELE_SNP:
            case ALLELE_MNP:
                if (alternateSequence == b.alternateSequence)
                    return true;
                break;
            case ALLELE_DELETION:
                if (length == b.length)
                    return true;
                break;
            case ALLELE_INSERTION:
                if (length == b.length 
                    && alternateSequence == b.alternateSequence)
                    return true;
                break;
            case ALLELE_COMPLEX:
                if (length == b.length
                    && alternateSequence == b.alternateSequence)
                    return true;
                break;
            default:
                break;
        }
    }

    return false;
}

bool areHomozygous(vector<Allele*>& alleles) {
    Allele* prev = alleles.front();
    for (vector<Allele*>::iterator allele = alleles.begin() + 1; allele != alleles.end(); ++allele) {
        if (**allele != *prev) {
            return false;
        }
    }
    return true;
}

// counts alleles which satisfy operator==
map<Allele, int> countAlleles(vector<Allele*>& alleles) {
    map<Allele, int> counts;
    for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        Allele& allele = **a;
        map<Allele, int>::iterator f = counts.find(allele);
        if (f == counts.end()) {
            counts[allele] = 1;
        } else {
            counts[allele] += 1;
        }
    }
    return counts;
}

// counts alleles which satisfy operator==
map<string, int> countAllelesString(vector<Allele*>& alleles) {
    map<string, int> counts;
    for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        Allele& thisAllele = **a;
        const string& allele = thisAllele.currentBase;
        map<string, int>::iterator f = counts.find(allele);
        if (f == counts.end()) {
            counts[allele] = 1;
        } else {
            counts[allele] += 1;
        }
    }
    return counts;
}

map<string, int> countAllelesString(vector<Allele>& alleles) {
    map<string, int> counts;
    for (vector<Allele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        Allele& thisAllele = *a;
        const string& allele = thisAllele.currentBase;
        map<string, int>::iterator f = counts.find(allele);
        if (f == counts.end()) {
            counts[allele] = 1;
        } else {
            counts[allele] += 1;
        }
    }
    return counts;
}


map<Allele, int> countAlleles(vector<Allele>& alleles) {
    map<Allele, int> counts;
    for (vector<Allele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        Allele& allele = *a;
        map<Allele, int>::iterator f = counts.find(allele);
        if (f == counts.end()) {
            counts[allele] = 1;
        } else {
            counts[allele] += 1;
        }
    }
    return counts;
}

map<Allele, int> countAlleles(list<Allele*>& alleles) {
    map<Allele, int> counts;
    for (list<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        Allele& allele = **a;
        map<Allele, int>::iterator f = counts.find(allele);
        if (f == counts.end()) {
            counts[allele] = 1;
        } else {
            counts[allele] += 1;
        }
    }
    return counts;
}


map<string, vector<Allele*> > groupAllelesBySample(list<Allele*>& alleles) {
    map<string, vector<Allele*> > groups;
    for (list<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        Allele*& allele = *a;
        groups[allele->sampleID].push_back(allele);
    }
    return groups;
}

vector<Allele> uniqueAlleles(list<Allele*>& alleles) {
    vector<Allele> uniques;
    map<Allele, int> counts = countAlleles(alleles);
    for (map<Allele, int>::iterator c = counts.begin(); c != counts.end(); ++c) {
        uniques.push_back(c->first);
    }
    return uniques;
}

void groupAllelesBySample(list<Allele*>& alleles, map<string, vector<Allele*> >& groups) {
    for (list<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        Allele*& allele = *a;
        groups[allele->sampleID].push_back(allele);
    }
}

void groupAlleles(map<string, vector<Allele*> >& sampleGroups, map<string, vector<Allele*> >& alleleGroups) {
    for (map<string, vector<Allele*> >::iterator sample = sampleGroups.begin(); sample != sampleGroups.end(); ++sample) {
        for (vector<Allele*>::iterator allele = sample->second.begin(); allele != sample->second.end(); ++allele) {
            alleleGroups[(*allele)->currentBase].push_back(*allele);
        }
    }
}

vector<vector<Allele*> > groupAlleles(map<string, vector<Allele*> > &sampleGroups, bool (*fncompare)(Allele &a, Allele &b)) {
    vector<vector<Allele*> > groups;
    for (map<string, vector<Allele*> >::iterator sg = sampleGroups.begin(); sg != sampleGroups.end(); ++sg) {
        vector<Allele*>& alleles = sg->second;
        for (vector<Allele*>::iterator oa = alleles.begin(); oa != alleles.end(); ++oa) {
            bool unique = true;
            for (vector<vector<Allele*> >::iterator ag = groups.begin(); ag != groups.end(); ++ag) {
                if ((*fncompare)(**oa, *ag->front())) {
                    ag->push_back(*oa);
                    unique = false;
                    break;
                }
            }
            if (unique) {
                vector<Allele*> trueAlleleGroup;
                trueAlleleGroup.push_back(*oa);
                groups.push_back(trueAlleleGroup);
            }
        }
    }
    return groups;
}

vector<vector<Allele*> >  groupAlleles(list<Allele*> &alleles, bool (*fncompare)(Allele* &a, Allele* &b)) {
    vector<vector<Allele*> > groups;
    for (list<Allele*>::iterator oa = alleles.begin(); oa != alleles.end(); ++oa) {
        bool unique = true;
        for (vector<vector<Allele*> >::iterator ag = groups.begin(); ag != groups.end(); ++ag) {
            if ((*fncompare)(*oa, ag->front())) {
                ag->push_back(*oa);
                unique = false;
                break;
            }
        }
        if (unique) {
            vector<Allele*> trueAlleleGroup;
            trueAlleleGroup.push_back(*oa);
            groups.push_back(trueAlleleGroup);
        }
    }
    return groups;
}

vector<vector<Allele*> >  groupAlleles(list<Allele> &alleles, bool (*fncompare)(Allele &a, Allele &b)) {
    vector<vector<Allele*> > groups;
    for (list<Allele>::iterator oa = alleles.begin(); oa != alleles.end(); ++oa) {
        bool unique = true;
        for (vector<vector<Allele*> >::iterator ag = groups.begin(); ag != groups.end(); ++ag) {
            if ((*fncompare)(*oa, *ag->front())) {
                ag->push_back(&*oa);
                unique = false;
                break;
            }
        }
        if (unique) {
            vector<Allele*> trueAlleleGroup;
            trueAlleleGroup.push_back(&*oa);
            groups.push_back(trueAlleleGroup);
        }
    }
    return groups;
}

vector<vector<Allele*> >  groupAlleles(vector<Allele*> &alleles, bool (*fncompare)(Allele &a, Allele &b)) {
    vector<vector<Allele*> > groups;
    for (vector<Allele*>::iterator oa = alleles.begin(); oa != alleles.end(); ++oa) {
        bool unique = true;
        for (vector<vector<Allele*> >::iterator ag = groups.begin(); ag != groups.end(); ++ag) {
            if ((*fncompare)(**oa, *ag->front())) {
                ag->push_back(*oa);
                unique = false;
                break;
            }
        }
        if (unique) {
            vector<Allele*> trueAlleleGroup;
            trueAlleleGroup.push_back(*oa);
            groups.push_back(trueAlleleGroup);
        }
    }
    return groups;
}

vector<vector<Allele*> >  groupAlleles(vector<Allele> &alleles, bool (*fncompare)(Allele &a, Allele &b)) {
    vector<vector<Allele*> > groups;
    for (vector<Allele>::iterator oa = alleles.begin(); oa != alleles.end(); ++oa) {
        bool unique = true;
        for (vector<vector<Allele*> >::iterator ag = groups.begin(); ag != groups.end(); ++ag) {
            if ((*fncompare)(*oa, *ag->front())) {
                ag->push_back(&*oa);
                unique = false;
                break;
            }
        }
        if (unique) {
            vector<Allele*> trueAlleleGroup;
            trueAlleleGroup.push_back(&*oa);
            groups.push_back(trueAlleleGroup);
        }
    }
    return groups;
}

vector<vector<Allele> >  groupAlleles_copy(list<Allele> &alleles, bool (*fncompare)(Allele &a, Allele &b)) {
    vector<vector<Allele> > groups;
    for (list<Allele>::iterator oa = alleles.begin(); oa != alleles.end(); ++oa) {
        bool unique = true;
        for (vector<vector<Allele> >::iterator ag = groups.begin(); ag != groups.end(); ++ag) {
            if ((*fncompare)(*oa, ag->front())) {
                ag->push_back(*oa);
                unique = false;
                break;
            }
        }
        if (unique) {
            vector<Allele> trueAlleleGroup;
            trueAlleleGroup.push_back(*oa);
            groups.push_back(trueAlleleGroup);
        }
    }
    return groups;
}

vector<vector<Allele> >  groupAlleles_copy(vector<Allele> &alleles, bool (*fncompare)(Allele &a, Allele &b)) {
    vector<vector<Allele> > groups;
    for (vector<Allele>::iterator oa = alleles.begin(); oa != alleles.end(); ++oa) {
        bool unique = true;
        for (vector<vector<Allele> >::iterator ag = groups.begin(); ag != groups.end(); ++ag) {
            if ((*fncompare)(*oa, ag->front())) {
                ag->push_back(*oa);
                unique = false;
                break;
            }
        }
        if (unique) {
            vector<Allele> trueAlleleGroup;
            trueAlleleGroup.push_back(*oa);
            groups.push_back(trueAlleleGroup);
        }
    }
    return groups;
}

vector<vector<Allele> >  groupAlleles_copy(vector<Allele> &alleles) {
    vector<vector<Allele> > groups;
    for (vector<Allele>::iterator oa = alleles.begin(); oa != alleles.end(); ++oa) {
        bool unique = true;
        for (vector<vector<Allele> >::iterator ag = groups.begin(); ag != groups.end(); ++ag) {
            if (*oa == ag->front()) {
                ag->push_back(*oa);
                unique = false;
                break;
            }
        }
        if (unique) {
            vector<Allele> trueAlleleGroup;
            trueAlleleGroup.push_back(*oa);
            groups.push_back(trueAlleleGroup);
        }
    }
    return groups;
}

bool Allele::sameSample(Allele &other) { return this->sampleID == other.sampleID; }

bool allelesSameType(Allele* &a, Allele* &b) { return a->type == b->type; }

bool allelesEquivalent(Allele* &a, Allele* &b) { return a->equivalent(*b); }

bool allelesSameSample(Allele* &a, Allele* &b) { return a->sampleID == b->sampleID; }

bool allelesSameType(Allele &a, Allele &b) { return a.type == b.type; }

bool allelesEquivalent(Allele &a, Allele &b) { return a.equivalent(b); }

bool allelesSameSample(Allele &a, Allele &b) { return a.sampleID == b.sampleID; }

bool allelesEqual(Allele &a, Allele &b) { return a == b; }

vector<Allele> genotypeAllelesFromAlleleGroups(vector<vector<Allele*> > &groups) {

    vector<Allele> results;
    for (vector<vector<Allele*> >::iterator g = groups.begin(); g != groups.end(); ++g)
        results.push_back(genotypeAllele(*g->front()));
    return results;

}

vector<Allele> genotypeAllelesFromAlleles(vector<Allele*> &alleles) {

    vector<Allele> results;
    for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a)
        results.push_back(genotypeAllele(**a));
    return results;

}

vector<Allele> genotypeAllelesFromAlleleGroups(vector<vector<Allele> > &groups) {

    vector<Allele> results;
    for (vector<vector<Allele> >::iterator g = groups.begin(); g != groups.end(); ++g)
        results.push_back(genotypeAllele(g->front()));
    return results;

}

vector<Allele> genotypeAllelesFromAlleles(vector<Allele> &alleles) {

    vector<Allele> results;
    for (vector<Allele>::iterator a = alleles.begin(); a != alleles.end(); ++a)
        results.push_back(genotypeAllele(*a));
    return results;

}

Allele genotypeAllele(Allele &a) {
    return Allele(a.type, a.alternateSequence, a.length, a.referenceLength, a.position);
}

Allele genotypeAllele(AlleleType type, string alt, unsigned int len, unsigned int reflen, long double pos) {
    return Allele(type, alt, len, reflen, pos);
}

/*
vector<Allele> uniqueAlleleObservations(vector<vector<Allele> > &alleleObservations) {
}
*/


/*
void* AlleleFreeList::NewAllele() {
    if (_allocs > _tick_allocs) {
        _allocs = 0; // reset allocation counter
        // cerr << "_size = " << _size
        //  << " _min_size = " << _min_size
        //  << " _max_size = " << _max_size << endl;
        Resize(_size - _min_size);  // resize to max flux in freelist size
        _min_size = _size; // reset size counters
        _max_size = _size;
    } else {
        ++_allocs;
    }
    if (_p != NULL) {
        void* mem = _p;
        _p = _p->pNext();
        --_size;
        if (_size < _min_size) {
            _min_size = _size;
        }
        return mem;
    } else {
        return ::new char [sizeof (Allele)];
    }
}

void AlleleFreeList::Recycle(void* mem) {
    Allele* allele = static_cast<Allele*> (mem);
    allele->_pNext = _p;
    _p = allele;
    ++_size;
    ++_allocs;
    if (_size > _max_size) {
        _max_size = _size;
    }
}

void AlleleFreeList::Resize(int new_size) {
    //cerr << "resizing free list from " << _size << " to " << new_size << endl;
    while (_size > new_size && _p != NULL) {
        char * mem = reinterpret_cast<char *> (_p);
        _p = _p->pNext();
        ::delete [] mem;
        --_size;
    }
    //cerr << "new size " << _size << endl;
}

AlleleFreeList::~AlleleFreeList() {
    Purge();
}

void AlleleFreeList::Purge() {
    while (_p != NULL) {
        char * mem = reinterpret_cast<char *> (_p);
        _p = _p->pNext();
        ::delete [] mem;
        --_size;
    }
}
*/

int allowedAlleleTypes(vector<AlleleType>& allowedEnumeratedTypes) {
    int allowedTypes = 0;// (numberOfPossibleAlleleTypes, false);
    for (vector<AlleleType>::iterator t = allowedEnumeratedTypes.begin(); t != allowedEnumeratedTypes.end(); ++t) {
        allowedTypes |= *t;
    }
    return allowedTypes;
}

void filterAlleles(list<Allele*>& alleles, int allowedTypes) {

    for (list<Allele*>::iterator allele = alleles.begin(); allele != alleles.end(); ++allele) {
        bool allowed = false;
        if (!(allowedTypes & (*allele)->type))
            *allele = NULL;
    }
    alleles.erase(remove(alleles.begin(), alleles.end(), (Allele*)NULL), alleles.end());

}

// removes alleles which are indelmasked at position
void removeIndelMaskedAlleles(list<Allele*>& alleles, long double position) {

    for (list<Allele*>::iterator allele = alleles.begin(); allele != alleles.end(); ++allele) {
        //cerr << *allele << " " << (*allele)->indelMask.size() << " " << (*allele)->referenceOffset() << endl;
        if ((*allele)->masked()) {
            *allele = NULL;
        }
    }
    alleles.erase(remove(alleles.begin(), alleles.end(), (Allele*)NULL), alleles.end());

}

int countAlleles(map<string, vector<Allele*> >& sampleGroups) {

    int count = 0;
    for (map<string, vector<Allele*> >::iterator sg = sampleGroups.begin(); sg != sampleGroups.end(); ++sg) {
        count += sg->second.size();
    }
    return count;

}

int countAllelesWithBase(vector<Allele*>& alleles, string base) {

    int count = 0;
    for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        if ((*a)->currentBase == base)
            ++count;
    }
    return count;

}

int baseCount(vector<Allele*>& alleles, string base, AlleleStrand strand) {

    int count = 0;
    for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        if ((*a)->currentBase == base && (*a)->strand == strand)
            ++count;
    }
    return count;

}

pair<pair<int, int>, pair<int, int> >
baseCount(vector<Allele*>& alleles, string refbase, string altbase) {
    
    int forwardRef = 0;
    int reverseRef = 0;
    int forwardAlt = 0;
    int reverseAlt = 0;

    for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        string base = (*a)->currentBase;
        AlleleStrand strand = (*a)->strand;
        if (base == refbase) {
            if (strand == STRAND_FORWARD)
                ++forwardRef;
            else if (strand == STRAND_REVERSE)
                ++reverseRef;
        } else if (base == altbase) {
            if (strand == STRAND_FORWARD)
                ++forwardAlt;
            else if (strand == STRAND_REVERSE)
                ++reverseAlt;
        }
    }

    return make_pair(make_pair(forwardRef, forwardAlt), make_pair(reverseRef, reverseAlt));

}

// combines the two alleles into a complex variant, updates important data
void Allele::mergeAllele(const Allele& newAllele) {
    type = ALLELE_COMPLEX;
    alternateSequence += newAllele.alternateSequence;
    length += newAllele.length; // hmmm
    referenceLength += newAllele.referenceLength;
    basesRight = newAllele.basesRight;
    bpRight = newAllele.bpRight;
    currentBase = base();
    quality += newAllele.quality;
    lnquality += newAllele.lnquality;
}

unsigned int Allele::getLengthOnReference(void) {

    switch (type) {
        case ALLELE_SNP:
        case ALLELE_MNP:
        case ALLELE_REFERENCE:
            return length;
            break;
        case ALLELE_INSERTION:
            return 0;
            break;
        case ALLELE_DELETION:
            return length;
            break;
        default:
            return length;
            break;
    }

}

vector<Allele> alleleUnion(vector<Allele>& a1, vector<Allele>& a2) {
    map<string, Allele> alleleSet;
    vector<Allele> results;
    for (vector<Allele>::iterator a = a1.begin(); a != a1.end(); ++a) {
        alleleSet.insert(make_pair(a->base(), *a));
    }
    for (vector<Allele>::iterator a = a2.begin(); a != a2.end(); ++a) {
        alleleSet.insert(make_pair(a->base(), *a));
    }
    for (map<string, Allele>::iterator a = alleleSet.begin(); a != alleleSet.end(); ++a) {
        results.push_back(a->second);
    }
    return results;
}
