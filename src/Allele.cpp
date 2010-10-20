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
}

void updateAllelesCachedData(vector<Allele*>& alleles) {
    for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        (*a)->update();
    }
}

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
        case ALLELE_MISMATCH:
            t = "mismatch";
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
        case ALLELE_CNV:
            t = "cnv";
            break;
        default:
            t = "unknown";
            break;
    }

    return t;

}

const string Allele::base(void) const { // the base of this allele

    if (genotypeAllele)
        return alternateSequence;

    switch (this->type) {
        case ALLELE_REFERENCE:
            return string(1, *currentReferenceBase);
            break;
        case ALLELE_GENOTYPE:
            return alternateSequence;
            break;
        //case ALLELE_MISMATCH:
            //break;
        case ALLELE_SNP:
            return alternateSequence;
            break;
        case ALLELE_INSERTION:
            return "I" + length; // unclear what to do here
            break;
        case ALLELE_DELETION:
            return "D" + length;
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
            return true;
            break;
        case ALLELE_DELETION:
            return true;
            break;
        default:
            break;
    }

}

string stringForAllele(Allele &allele) {

    stringstream out;
    if (!allele.genotypeAllele) {
        out 
            << allele.sampleID << "\t"
            << allele.readID << "\t"
            << allele.typeStr() << "\t" 
            << allele.position << "\t"
            << allele.length << "\t"
            << (allele.strand == STRAND_FORWARD ? "+" : "-") << "\t"
            << allele.referenceSequence << "\t"
            << allele.alternateSequence << "\t"
            << allele.qualityString << "\t"
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

string json(map<string, vector<Allele*> > &alleles) {
    stringstream out;
    out << "[";
    bool first = true;
    for (map<string, vector<Allele*> >::iterator g = alleles.begin(); g != alleles.end(); ++g) {
        vector<Allele*>& alleles = g->second;
        for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
            if (!first) { out << ","; } else { first = false; }
            out << json(**a);
        }
    }
    out << "]";
    return out.str();
}

string json(vector<Allele*> &alleles) {
    stringstream out;
    vector<Allele*>::iterator a = alleles.begin();
    out << "[" << json(**a); ++a;
    for (; a != alleles.end(); ++a)
        out << "," << json(**a);
    out << "]";
    return out.str();
}

string json(Allele*& allele) { return json(*allele); }

string json(Allele& allele) {
    int referenceOffset = allele.referenceOffset();
    stringstream out;
    if (!allele.genotypeAllele) {
        out << "{\"id\":\"" << allele.readID << "\""
            << ",\"type\":\"" << allele.typeStr() << "\""
            << ",\"length\":" << ((allele.type == ALLELE_REFERENCE) ? 1 : allele.length)
            << ",\"position\":" << allele.position 
            << ",\"strand\":\"" << (allele.strand == STRAND_FORWARD ? "+" : "-") << "\"";
        if (allele.type == ALLELE_REFERENCE ) {
            out << ",\"base\":\"" << allele.alternateSequence.at(referenceOffset) << "\""
                //<< ",\"reference\":\"" << allele.referenceSequence.at(referenceOffset) << "\""
                << ",\"quality\":" << allele.currentQuality();
        } else {
            out << ",\"base\":\"" << allele.alternateSequence << "\""
                //<< ",\"reference\":\"" << allele.referenceSequence << "\""
                << ",\"quality\":" << allele.quality;
        }
        out << "}";

    } else {
        out << "{\"type\":\"" << allele.typeStr() << "\"";
        switch (allele.type) {
            case ALLELE_REFERENCE:
                out << "}";
                break;
            default:
                out << "\",\"length\":" << allele.length 
                    << ",\"alt\":\"" << allele.alternateSequence << "\"}";
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
        out << allele.readID 
            << ":" << allele.typeStr() 
            << ":" << allele.length 
            << ":" << allele.position 
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
                if (currentBase == b.currentBase)
                    return true;
                break;
            case ALLELE_SNP:
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
            alleleGroups[(*allele)->base()].push_back(*allele);
        }
    }
}

void groupAlleles(Samples& samples, map<string, vector<Allele*> >& alleleGroups) {
    for (Samples::iterator s = samples.begin(); s != samples.end(); ++s) {
        Sample& sample = s->second;
        for (Sample::iterator g = sample.begin(); g != sample.end(); ++g) {
            const string& base = g->first;
            const vector<Allele*>& alleles = g->second;
            vector<Allele*>& group = alleleGroups[base];
            group.reserve(group.size() + distance(alleles.begin(), alleles.end()));
            group.insert(group.end(), alleles.begin(), alleles.end());
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
    return Allele(a.type, a.alternateSequence, a.length, a.position);
}

Allele genotypeAllele(AlleleType type, string alt, unsigned int len) {
    return Allele(type, alt, len);
}

/*
vector<Allele> uniqueAlleleObservations(vector<vector<Allele> > &alleleObservations) {
}
*/


void* AlleleFreeList::NewAllele() {
    if (_p != NULL) {
        void* mem = _p;
        _p = _p->pNext();
        return mem;
    } else {
        return ::new char [sizeof (Allele)];
    }
}

void AlleleFreeList::Recycle(void* mem) {
    Allele* allele = static_cast<Allele*> (mem);
    allele->_pNext = _p;
    _p = allele;
    //++_allocs;
}

AlleleFreeList::~AlleleFreeList() {
    Purge();
}

void AlleleFreeList::Purge() {
    while (_p != NULL) {
        char * mem = reinterpret_cast<char *> (_p);
        _p = _p->pNext();
        ::delete [] mem;
    }
}

bool sufficientAlternateObservations(Samples& samples, int mincount, float minfraction) {

    int totalAlternateCount = 0;
    int totalReferenceCount = 0;

    for (Samples::iterator s = samples.begin(); s != samples.end(); ++s) {

        //cerr << s->first << endl;
        Sample& sample = s->second;
        int alternateCount = 0;
        int observationCount = 0;

        for (Sample::iterator group = sample.begin(); group != sample.end(); ++group) {
            const string& base = group->first;
            //cerr << base << endl;
            vector<Allele*>& alleles = group->second;
            //cerr << alleles.size() << endl;
            if (alleles.size() == 0)
                continue;
            if (alleles.front()->type != ALLELE_REFERENCE) {
                alternateCount += alleles.size();
            } else {
                totalReferenceCount += alleles.size();
            }
            observationCount += alleles.size();
        }

        if (alternateCount >= mincount && ((float) alternateCount / (float) observationCount) >= minfraction)
            return true;
        totalAlternateCount += alternateCount;
    
    }

    // always analyze if we have more alternate observations than reference observations
    // this is meant to catch the case in which the reference is the rare allele
    // it will probably also catch cases in which we have very low coverage
    if (totalReferenceCount < totalAlternateCount) {
        return true;
    }

    return false;

}

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
void removeIndelMaskedAlleles(list<Allele*>& alleles, long unsigned int position) {

    for (list<Allele*>::iterator allele = alleles.begin(); allele != alleles.end(); ++allele) {
        //cerr << *allele << " " << (*allele)->indelMask.size() << " " << (*allele)->referenceOffset() << endl;
        if ((*allele)->masked()) {
            *allele = NULL;
        }
    }
    alleles.erase(remove(alleles.begin(), alleles.end(), (Allele*)NULL), alleles.end());

}

int countAlleles(Samples& samples) {

    int count = 0;
    for (Samples::iterator s = samples.begin(); s != samples.end(); ++s) {
        Sample& sample = s->second;
        for (Sample::iterator sg = sample.begin(); sg != sample.end(); ++sg) {
            count += sg->second.size();
        }
    }
    return count;

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
