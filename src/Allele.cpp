#include "Allele.h"
#include "multichoose.h"
#include "TryCatch.h"


int Allele::referenceOffset(Position referencePosition) {
    /*cout << readID << " offset checked " << referencePosition - position << " against position " << position 
        << " allele length " << length << " str length " << referenceSequence.size() << " qstr size " << qualityString.size() << endl;
        */
    return referencePosition - position;
}

// quality at a given reference position
short Allele::Quality(Position referencePosition) {
    switch (this->type) {
        case ALLELE_REFERENCE:
            TRY {
            return qualityChar2ShortInt(qualityString.at(referenceOffset(referencePosition)));
            } CATCH;
            break;
        case ALLELE_INSERTION:
        case ALLELE_DELETION:
        case ALLELE_SNP:
            return quality;
            break;
    }
}

string Allele::Type(void) {

    string t;

    switch (this->type) {
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

string stringForAllele(Allele &allele) {

    stringstream out;
    if (!allele.genotypeAllele) {
        out 
            << allele.sampleID << "\t"
            << allele.readID << "\t"
            << allele.Type() << "\t" 
            << allele.position << "\t"
            << allele.length << "\t"
            << (allele.strand == STRAND_FORWARD ? "+" : "-") << "\t"
            << allele.referenceSequence << "\t"
            << allele.alternateSequence << "\t"
            << allele.qualityString << "\t"
            << allele.quality << "\t" << endl;
    } else {
        out << allele.Type() << "\t";
        switch (allele.type) {
            case ALLELE_REFERENCE:
                break;
            default:
                out << allele.length << "\t"
                    << allele.alternateSequence;
                break;
        }
        out << endl;
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
    out << "[" << json(**a++);
    while (a != alleles.end())
        out << "," << json(**a++);
    out << "]";
    return out.str();
}

string json(Allele*& allele) { return json(*allele); }

string json(Allele& allele) {
    stringstream out;
    if (!allele.genotypeAllele) {
            // << &allele << ":" 
            out << "{\"id\":\"" << allele.readID 
                << "\",\"type\":\"" << allele.Type() 
                << "\",\"length\":" << allele.length 
                << ",\"position\":" << allele.position 
                << ",\"strand\":\"" << (allele.strand == STRAND_FORWARD ? "+" : "-")
                << "\",\"alt\":\"" << allele.alternateSequence
                << "\",\"reference\":\"" << allele.referenceSequence
                << "\",\"quality\":" << allele.quality << "}";
    } else {
        out << "{\"type\":\"" << allele.Type() << "\"";
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
                << ":" << allele.Type() 
                << ":" << allele.length 
                << ":" << allele.position 
                << ":" << (allele.strand == STRAND_FORWARD ? "+" : "-")
                << ":" << allele.alternateSequence
                << ":" << allele.referenceSequence
                << ":" << allele.quality;
    } else {
        out << allele.Type();
        switch (allele.type) {
            case ALLELE_REFERENCE:
                break;
            default:
                out << ":" << allele.length << ":" << (string) allele.alternateSequence;
                break;
        }
    }
    return out;
}

// for sorting alleles by type
bool operator<(Allele &a, Allele &b) {
    return a.type < b.type;
}

// for sorting pairs of genotype, probs
bool genotypeCmp(pair<Genotype, long double> a, pair<Genotype, long double> b) {
    return a.second > b.second;
}

// impossible:
/*
bool operator<(Allele* &a, Allele* &b) {
    return a->type < b->type;
}

// possible...
bool allelePtrCmp(Allele* &a, Allele* &b) {
    return a->type < b->type;
}
*/

bool operator==(Allele &a, Allele &b) {

    return a.type == b.type 
        && a.referenceName == b.referenceName 
        && a.position == b.position 
        && a.length == b.length
        && a.alternateSequence == b.alternateSequence;

}

bool Allele::equivalent(Allele &b) {

    if (type != b.type) {
        return false;
    } else {
        switch (type) {
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
            case ALLELE_REFERENCE:
                return true;
                break;
            default:
                break;
        }
    }

    return false;
}

map<string, vector<Allele*> > groupAllelesBySample(list<Allele*>& alleles) {
    map<string, vector<Allele*> > groups;
    groupAllelesBySample(alleles, groups);
    return groups;
}

void groupAllelesBySample(list<Allele*>& alleles, map<string, vector<Allele*> >& groups) {
    for (list<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
        Allele*& allele = *a;
        groups[allele->sampleID].push_back(allele);
    }
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

bool Allele::sameSample(Allele &other) { return this->sampleID == other.sampleID; }

bool allelesSameType(Allele* &a, Allele* &b) { return a->type == b->type; }

bool allelesEquivalent(Allele* &a, Allele* &b) { return a->equivalent(*b); }

bool allelesSameSample(Allele* &a, Allele* &b) { return a->sampleID == b->sampleID; }

bool allelesSameType(Allele &a, Allele &b) { return a.type == b.type; }

bool allelesEquivalent(Allele &a, Allele &b) { return a.equivalent(b); }

bool allelesSameSample(Allele &a, Allele &b) { return a.sampleID == b.sampleID; }

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
    ++_allocs;
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
