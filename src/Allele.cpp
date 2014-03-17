#include "Allele.h"
#include "multichoose.h"
#include "TryCatch.h"


int Allele::referenceOffset(void) const {
    /*cout << readID << " offset checked " << referencePosition - position << " against position " << position 
        << " allele length " << length << " str length " << referenceSequence.size() << " qstr size " << qualityString.size() << endl;
        */
    return *currentReferencePosition - position;
}

void Allele::setQuality(void) {
    quality = currentQuality();
    lnquality = phred2ln(quality);
}

int Allele::bpLeft(void) {
    return position - alignmentStart;
}

int Allele::bpRight(void) {
    return alignmentEnd - (position + referenceLength);
}

// called prior to using the allele in analysis
// called again when haplotype alleles are built, in which case the "currentBase" is set to the alternate sequence of the allele
void Allele::update(int haplotypeLength) {
    if (haplotypeLength == 1) {
        if (type == ALLELE_REFERENCE) {
            currentBase = string(1, *currentReferenceBase);
        } else {
            currentBase = base();
        }
    } else {
        currentBase = base();
    }
    // should be done after setting currentBase to haplotypeLength
    if (isReference()) setQuality();
    basesLeft = bpLeft();
    basesRight = bpRight();
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
    //cerr << readID << " " << position << "-" << position + length << " " << alternateSequence.size() << " vs " << baseQualities.size() << endl;
    switch (this->type) {
        case ALLELE_REFERENCE:
            // should check a different way... this is wrong
            // it will catch it all the time,
            if (currentBase.size() > 1) {
                return averageQuality(baseQualities);
            } else {
                int off = referenceOffset();
                if (off < 0 || off > baseQualities.size()) {
                    return 0;
                } else {
                    return baseQualities.at(off);
                }
            }
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

string Allele::typeStr(void) const {

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
        case ALLELE_NULL:
            t = "null";
            break;
        default:
            t = "unknown";
            break;
    }

    return t;

}

bool Allele::isReference(void) const {
    return type == ALLELE_REFERENCE;
}

bool Allele::isSNP(void) const {
    return type == ALLELE_SNP;
}

bool Allele::isInsertion(void) const {
    return type == ALLELE_INSERTION;
}

bool Allele::isDeletion(void) const {
    return type == ALLELE_DELETION;
}

bool Allele::isMNP(void) const {
    return type == ALLELE_MNP;
}

bool Allele::isComplex(void) const {
    return type == ALLELE_COMPLEX;
}

bool Allele::isNull(void) const {
    return type == ALLELE_NULL;
}

const string Allele::base(void) const { // the base of this allele

    switch (this->type) {
    case ALLELE_REFERENCE:
        if (genotypeAllele)
            return alternateSequence;
        else
            return currentBase;
        break;
    case ALLELE_GENOTYPE:
        return alternateSequence;
        break;
        /*
    case ALLELE_GENOTYPE:
        return alternateSequence; // todo fix
        break;
    case ALLELE_REFERENCE:
        return "R:" + convert(position) + ":" + cigar + ":" + alternateSequence;
        break;
        */
    case ALLELE_SNP:
        return "S:" + convert(position) + ":" + cigar + ":" + alternateSequence;
        break;
    case ALLELE_MNP:
        return "M:" + convert(position) + ":" + cigar + ":" + alternateSequence;
        break;
    case ALLELE_INSERTION:
        return "I:" + convert(position) + ":" + cigar + ":" + alternateSequence;
        break;
    case ALLELE_DELETION:
        return "D:" + convert(position) + ":" + cigar;
        break;
    case ALLELE_COMPLEX:
        return "C:" + convert(position) + ":" + cigar + ":" + alternateSequence;
        break;
    case ALLELE_NULL:
        return "N:" + convert(position) + ":" + alternateSequence;
    default:
        break;
    }

}

string stringForAllele(const Allele &allele) {

    stringstream out;
    if (!allele.genotypeAllele) {
        out.precision(1);
        out 
            << allele.sampleID << ":"
            << allele.readID << ":"
            << allele.typeStr() << ":"
            << allele.cigar << ":"
            << scientific << fixed << allele.position << ":"
            << allele.length << ":"
            << (allele.strand == STRAND_FORWARD ? "+" : "-") << ":"
            << allele.referenceSequence << ":"
            << allele.alternateSequence << ":"
            << allele.quality << ":"
            << allele.basesLeft << ":"
            << allele.basesRight;
    } else {
        out << allele.typeStr() << ":"
            << allele.cigar << ":"
            << scientific << fixed << allele.position << ":"
            << allele.length << ":"
            << allele.alternateSequence;
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
        int prec = out.precision();
        // << &allele << ":" 
        out.precision(1);
        out << allele.sampleID
            << ":" << allele.readID 
            << ":" << allele.typeStr() 
            << ":" << allele.length 
            << ":" << allele.referenceLength
            << ":" << scientific << fixed << allele.position 
            << ":" << (allele.strand == STRAND_FORWARD ? "+" : "-")
            << ":" << allele.alternateSequence
            //<< ":" << allele.referenceSequence
            << ":" << allele.repeatRightBoundary
            << ":" << allele.cigar
            << ":" << allele.lnmapQuality
            << ":" << allele.lnquality;
        out.precision(prec);
    } else {
        out << allele.typeStr() 
            << ":" << allele.cigar
            << ":" << scientific << fixed << allele.position
            << ":" << allele.length 
            << ":" << (string) allele.alternateSequence;
    }
    out.precision(5);
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
                    && alternateSequence == b.alternateSequence
                    && cigar == b.cigar)
                    return true;
                break;
            case ALLELE_NULL:
                return alternateSequence == b.alternateSequence;
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

// in haplotype calling, if alleles have the same alternate sequence, the
// should have the same cigar and position.  this function picks the most
// common allele observation per alternate sequence and homogenizes the rest to
// the same if they are not reference alleles
void homogenizeAlleles(map<string, vector<Allele*> >& alleleGroups, string& refseq, Allele& refallele) {
    map<string, map<string, int> > equivs;
    map<string, Allele*> homogenizeTo;
    // find equivalencies between alleles
    // base equivalency is self
    for (map<string, vector<Allele*> >::iterator g = alleleGroups.begin(); g != alleleGroups.end(); ++g) {
        Allele& allele = *g->second.front();
        if (allele.isReference()) {
            continue;
        }
        equivs[allele.alternateSequence][g->first]++;
    }
    // 
    for (map<string, map<string, int> >::iterator e = equivs.begin(); e != equivs.end(); ++e) {
        string altseq = e->first;
        map<string, int>& group = e->second;
        map<int, string> ordered;
        for (map<string, int>::iterator f = group.begin(); f != group.end(); ++f) {
            // pick the best by count
            ordered[f->second] = f->first;
        }

        // choose the most common group
        string& altbase = ordered.rbegin()->second;
        if (altseq == refseq) {
            homogenizeTo[altseq] = &refallele;
        } else {
            homogenizeTo[altseq] = alleleGroups[altbase].front();
        }
    }
    for (map<string, vector<Allele*> >::iterator g = alleleGroups.begin(); g != alleleGroups.end(); ++g) {
        vector<Allele*>& alleles = g->second;
        if (alleles.front()->isReference()) {
            continue;
        }
        string& altseq = alleles.front()->alternateSequence;
        Allele* toallele = homogenizeTo[altseq];
        string& cigar = toallele->cigar;
        AlleleType type = toallele->type;
        long int position = toallele->position;
        for (vector<Allele*>::iterator a = alleles.begin(); a != alleles.end(); ++a) {
            (*a)->cigar = cigar;
            (*a)->type = type;
            (*a)->position = position;
            (*a)->update();
        }
    }
}

void resetProcessedFlag(map<string, vector<Allele*> >& alleleGroups) {
    for (map<string, vector<Allele*> >::iterator g = alleleGroups.begin(); g != alleleGroups.end(); ++g) {
        for (vector<Allele*>::iterator a = g->second.begin(); a != g->second.end(); ++a) {
            (*a)->processed = false;
        }
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
    return Allele(a.type, a.alternateSequence, a.length, a.referenceLength, a.cigar, a.position, a.repeatRightBoundary);
}

Allele genotypeAllele(AlleleType type, string alt, unsigned int len, string cigar, unsigned int reflen, long int pos, long int rrbound) {
    return Allele(type, alt, len, reflen, cigar, pos, rrbound);
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

string Allele::readSeq(void) {
    string r;
    for (vector<Allele>::iterator a = alignmentAlleles->begin(); a != alignmentAlleles->end(); ++a) {
        r.append(a->alternateSequence);
    }
    return r;
}

string Allele::read5p(void) {
    string r;
    vector<Allele>::const_reverse_iterator a = alignmentAlleles->rbegin();
    while (&*a != this) {
        ++a;
    }
    if ((a+1) != alignmentAlleles->rend()) ++a;
    while (a != alignmentAlleles->rend()) {
        r = a->alternateSequence + r;
        ++a;
    }
    r.append(alternateSequence);
    return r;
}

string Allele::read3p(void) {
    string r = alternateSequence;
    vector<Allele>::const_iterator a = alignmentAlleles->begin();
    while (&*a != this) {
        ++a;
    }
    if ((a+1) != alignmentAlleles->end()) ++a;
    while (a != alignmentAlleles->end()) {
        r.append(a->alternateSequence);
        ++a;
    }
    return r;
}

string Allele::read5pNonNull(void) {
    string r = alternateSequence;
    vector<Allele>::const_reverse_iterator a = alignmentAlleles->rbegin();
    while (&*a != this) {
        ++a;
    }
    while (a != alignmentAlleles->rend() && !a->isNull()) {
        if (&*a != this) {
            r = a->alternateSequence + r;
        }
        ++a;
    }
    return r;
}

string Allele::read3pNonNull(void) {
    string r = alternateSequence;
    vector<Allele>::const_iterator a = alignmentAlleles->begin();
    while (&*a != this) {
        ++a;
    }
    while (a != alignmentAlleles->end() && !a->isNull()) {
        if (&*a != this) {
            r.append(a->alternateSequence);
        }
        ++a;
    }
    return r;
}

int Allele::read5pNonNullBases(void) {
    int bp = 0;
    vector<Allele>::const_reverse_iterator a = alignmentAlleles->rbegin();
    while (&*a != this) {
        ++a;
    }
    while (a != alignmentAlleles->rend() && !a->isNull()) {
        if (&*a != this) {
            //cerr << "5p bp = " << bp << " adding " << stringForAllele(*a) << " to " << stringForAllele(*this) << endl;
            bp += a->alternateSequence.size();
        }
        ++a;
    }
    return bp;
}

int Allele::read3pNonNullBases(void) {
    int bp = 0;
    vector<Allele>::const_iterator a = alignmentAlleles->begin();
    while (&*a != this) {
        ++a;
    }
    while (a != alignmentAlleles->end() && !a->isNull()) {
        if (&*a != this) {
            //cerr << "3p bp = " << bp << " adding " << stringForAllele(*a) << " to " << stringForAllele(*this) << endl;
            bp += a->alternateSequence.size();
        }
        ++a;
    }
    return bp;
}

// adjusts the allele to have a new start
// returns the ref/alt sequence obtained by subtracting length from the left end of the allele
void Allele::subtract(
        int subtractFromRefStart,
        int subtractFromRefEnd,
        string& substart,
        string& subend,
        vector<pair<int, string> >& cigarStart,
        vector<pair<int, string> >& cigarEnd,
        vector<short>& qsubstart,
        vector<short>& qsubend
    ) {

    substart.clear();
    subend.clear();
    cigarStart.clear();
    cigarEnd.clear();
    qsubstart.clear();
    qsubend.clear();

    // prepare to adjust cigar
    list<pair<int, string> > cigarL = splitCigarList(cigar);

    // walk the cigar string to determine where to make the left cut in the alternate sequence
    int subtractFromAltStart = 0;
    if (subtractFromRefStart) {
        int refbpstart = subtractFromRefStart;
        pair<int, string> c;
        while (!cigarL.empty()) {
            c = cigarL.front();
            cigarL.pop_front();
            char op = c.second[0];
            switch (op) {
                case 'M':
                case 'X':
                case 'N':
                    refbpstart -= c.first;
                    subtractFromAltStart += c.first;
                    break;
                case 'I':
                    subtractFromAltStart += c.first;
                    break;
                case 'D':
                    refbpstart -= c.first;
                    break;
                default:
                    break;
            }

            cigarStart.push_back(c);

            if (refbpstart < 0) {
                // split/adjust the last cigar element
                cigarL.push_front(c);
                cigarL.front().first = -refbpstart;
                cigarStart.back().first += refbpstart;
                switch (op) {
                    case 'M':
                    case 'X':
                    case 'N':
                    case 'I':
                        subtractFromAltStart += refbpstart;
                        break;
                    case 'D':
                    default:
                        break;
                }
                break; // we're done
            }
        }
    }


    int subtractFromAltEnd = 0;
    // walk the cigar string to determine where to make the right cut in the alternate sequence
    if (subtractFromRefEnd) {
        int refbpend = subtractFromRefEnd;
        pair<int, string> c;
        while (!cigarL.empty() && refbpend > 0) {
            c = cigarL.back();
            cigarL.pop_back();
            char op = c.second[0];
            switch (op) {
                case 'M':
                case 'X':
                case 'N':
                    subtractFromAltEnd += c.first;
                    refbpend -= c.first;
                    break;
                case 'I':
                    subtractFromAltEnd += c.first;
                    break;
                case 'D':
                    refbpend -= c.first;
                    break;
                default:
                    break;
            }

            cigarEnd.insert(cigarEnd.begin(), c);

            if (refbpend < 0) {
                // split/adjust the last cigar element
                cigarL.push_back(c);
                cigarL.back().first = -refbpend;
                cigarEnd.front().first += refbpend;
                switch (op) {
                    case 'M':
                    case 'X':
                    case 'I':
                    case 'N':
                        subtractFromAltEnd += refbpend;
                        break;
                    case 'D':
                    default:
                        break;
                }
                break; // drop out of loop, we're done
            }
        }
    }

    // adjust the alternateSequence
    substart = alternateSequence.substr(0, subtractFromAltStart);
    subend = alternateSequence.substr(alternateSequence.size() - subtractFromAltEnd, subtractFromAltEnd);
    alternateSequence.erase(0, subtractFromAltStart);
    alternateSequence.erase(alternateSequence.size() - subtractFromAltEnd, subtractFromAltEnd);

    // adjust the quality string
    qsubstart.insert(qsubstart.begin(), baseQualities.begin(), baseQualities.begin() + subtractFromAltStart);
    qsubend.insert(qsubend.begin(), baseQualities.begin() + baseQualities.size() - subtractFromAltEnd, baseQualities.end());
    baseQualities.erase(baseQualities.begin(), baseQualities.begin() + subtractFromAltStart);
    baseQualities.erase(baseQualities.begin() + baseQualities.size() - subtractFromAltEnd, baseQualities.end());

    // reset the cigar
    cigarL.erase(remove_if(cigarL.begin(), cigarL.end(), isEmptyCigarElement), cigarL.end());
    cigar = joinCigarList(cigarL);

    // reset the length
    length = alternateSequence.size();

    // update the type specification
    updateTypeAndLengthFromCigar();

    // adjust the position
    position += subtractFromRefStart; // assumes the first-base of the alleles is reference==, not ins

    //referenceLength -= subtractFromRefStart;
    //referenceLength -= subtractFromRefEnd;

    referenceLength = referenceLengthFromCigar();

}

void Allele::subtractFromStart(int bp, string& seq, vector<pair<int, string> >& cig, vector<short>& quals) {
    string emptystr;
    vector<pair<int, string> > emptycigar;
    vector<short> emptyquals;
    subtract(bp, 0, seq, emptystr, cig, emptycigar, quals, emptyquals);
}

void Allele::subtractFromEnd(int bp, string& seq, vector<pair<int, string> >& cig, vector<short>& quals) {
    string emptystr;
    vector<pair<int, string> > emptycigar;
    vector<short> emptyquals;
    subtract(0, bp, emptystr, seq, emptycigar, cig, emptyquals, quals);
}

void Allele::addToStart(string& seq, vector<pair<int, string> >& cig, vector<short>& quals) {
    string emptystr;
    vector<pair<int, string> > emptycigar;
    vector<short> emptyquals;
    add(seq, emptystr, cig, emptycigar, quals, emptyquals);
}

void Allele::addToEnd(string& seq, vector<pair<int, string> >& cig, vector<short>& quals) {
    string emptystr;
    vector<pair<int, string> > emptycigar;
    vector<short> emptyquals;
    add(emptystr, seq, emptycigar, cig, emptyquals, quals);
}

void Allele::add(
        string& addToStart,
        string& addToEnd,
        vector<pair<int, string> >& cigarStart,
        vector<pair<int, string> >& cigarEnd,
        vector<short>& qaddToStart,
        vector<short>& qaddToEnd
    ) {

    // adjust the position
    for (vector<pair<int, string> >::iterator c = cigarStart.begin(); c != cigarStart.end(); ++c) {
        switch (c->second[0]) {
            case 'M':
            case 'X':
            case 'D':
            case 'N':
                position -= c->first;
                break;
            case 'I':
            default:
                break;
        }
    }

    // prepare to adjust cigar
    vector<pair<int, string> > cigarV = splitCigar(cigar);

    // adjust the cigar
    if (!cigarStart.empty()) {
        if (cigarStart.back().second == cigarV.front().second) {
            // merge
            cigarV.front().first += cigarStart.back().first;
            cigarStart.pop_back();
        }
    }
    cigarV.insert(cigarV.begin(), cigarStart.begin(), cigarStart.end());

    if (!cigarEnd.empty()) {
        if (cigarEnd.front().second == cigarV.back().second) {
            // merge
            cigarV.back().first += cigarEnd.front().first;
            cigarEnd.pop_back();
        } else {
            cigarV.insert(cigarV.end(), cigarEnd.begin(), cigarEnd.end());
        }
    }

    // adjust the alternateSequence
    alternateSequence.insert(0, addToStart);
    alternateSequence.append(addToEnd);

    // adjust the quality string
    baseQualities.insert(baseQualities.begin(), qaddToStart.begin(), qaddToStart.end());
    baseQualities.insert(baseQualities.end(), qaddToEnd.begin(), qaddToEnd.end());

    // reset the cigar
    cigarV.erase(remove_if(cigarV.begin(), cigarV.end(), isEmptyCigarElement), cigarV.end());
    cigar = joinCigar(cigarV);

    updateTypeAndLengthFromCigar();

    // reset referenceLength
    referenceLength = referenceLengthFromCigar();

}

void Allele::updateTypeAndLengthFromCigar(void) {

    vector<pair<int, string> > cigarV = splitCigar(cigar);

    map<char, int> cigarTypes;
    map<char, int> cigarLengths;
    for (vector<pair<int, string> >::iterator c = cigarV.begin(); c != cigarV.end(); ++c) {
        ++cigarTypes[c->second[0]];
        cigarLengths[c->second[0]] += c->first;
    }
    if (cigarTypes.size() == 1) {
        switch (cigarTypes.begin()->first) {
            case 'M':
                type = ALLELE_REFERENCE;
                break;
            case 'I':
                type = ALLELE_INSERTION;
                break;
            case 'D':
                type = ALLELE_DELETION;
                break;
            case 'X':
                if (cigarLengths['X'] > 1) {
                    type = ALLELE_MNP;
                } else {
                    type = ALLELE_SNP;
                }
                break;
            case 'N':
                type = ALLELE_NULL;
                break;
            default:
                break;
        }
    } else if (cigarTypes.size() == 2) {
        if (cigarTypes['M'] > 0) {
            if (cigarTypes['I'] == 1) {
                type = ALLELE_INSERTION;
            } else if (cigarTypes['D'] == 1) {
                type = ALLELE_DELETION;
            } else if (cigarTypes['X'] == 1) {
                if (cigarLengths['X'] > 1) {
                    type = ALLELE_MNP;
                } else {
                    type = ALLELE_SNP;
                }
            } else {
                type = ALLELE_COMPLEX;
            }
        } else {
            type = ALLELE_COMPLEX;
        }
    } else {
        type = ALLELE_COMPLEX;
    }

    // recalculate allele length and quality, based on type
    switch (type) {
        case ALLELE_REFERENCE:
            length = alternateSequence.size();
            break;
        case ALLELE_SNP:
        case ALLELE_MNP:
            length = cigarLengths['X'];
            break;
        case ALLELE_INSERTION:
            length = cigarLengths['I'];
            break;
        case ALLELE_DELETION:
            length = cigarLengths['D'];
            break;
        case ALLELE_COMPLEX:
            length = alternateSequence.size();
            break;
        case ALLELE_NULL:
            length = alternateSequence.size();
            break;
        default:
            break;
    }

}

int referenceLengthFromCigar(string& cigar) {
    int r = 0;
    vector<pair<int, string> > cigarV = splitCigar(cigar);
    for (vector<pair<int, string> >::iterator c = cigarV.begin(); c != cigarV.end(); ++c) {
        switch (c->second[0]) {
        case 'M':
        case 'X':
        case 'D':
        case 'N':
            r += c->first;
            break;
        case 'I':
        default:
            break;
        }
    }
    return r;
}

int Allele::referenceLengthFromCigar(void) {
    int r = 0;
    vector<pair<int, string> > cigarV = splitCigar(cigar);
    for (vector<pair<int, string> >::iterator c = cigarV.begin(); c != cigarV.end(); ++c) {
        switch (c->second[0]) {
        case 'M':
        case 'X':
        case 'D':
        case 'N':
            r += c->first;
            break;
        case 'I':
        default:
            break;
        }
    }
    return r;
}


// combines the two alleles into a complex variant, updates important data
void Allele::mergeAllele(const Allele& newAllele, AlleleType newType) {
    //cout << stringForAllele(*this) << endl << stringForAllele(newAllele) << endl;
    type = newType;
    alternateSequence += newAllele.alternateSequence;
    length += newAllele.length; // hmmm
    basesRight = newAllele.basesRight;
    baseQualities.insert(baseQualities.end(), newAllele.baseQualities.begin(), newAllele.baseQualities.end());
    currentBase = base();
    // XXX note that we don't add Q values for intermingled gaps in combined alleles
    if (newAllele.type != ALLELE_REFERENCE) {
        quality = min(newAllele.quality, quality);
        lnquality = max(newAllele.lnquality, lnquality);
        //quality = minQuality(baseQualities);
        //lnquality = log(quality);
    } else {
        quality = averageQuality(baseQualities);
        lnquality = log(quality);
        basesRight += newAllele.referenceLength;
    }
    if (newAllele.type != ALLELE_REFERENCE) {
        repeatRightBoundary = newAllele.repeatRightBoundary;
    }
    cigar = mergeCigar(cigar, newAllele.cigar);
    referenceLength = referenceLengthFromCigar();
    //cout << stringForAllele(*this) << endl << endl;
}

void Allele::squash(void) {
    // will trigger destruction of this allele in the AlleleParser
    length = 0;
    position = 0;
}

unsigned int Allele::getLengthOnReference(void) {
    return referenceLengthFromCigar();
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

bool isEmptyAllele(const Allele& allele) {
    return allele.length == 0;
}

bool isDividedIndel(const Allele& allele) {
    vector<pair<int, string> > cigarV = splitCigar(allele.cigar);
    if (cigarV.front().second == "D" || cigarV.front().second == "I") {
        return true;
    } else {
        return false;
    }
}

// returns true if this indel is not properly flanked by reference-matching sequence
bool isUnflankedIndel(const Allele& allele) {
    if (allele.isReference() || allele.isSNP() || allele.isMNP()) {
        return false;
    } else {
        vector<pair<int, string> > cigarV = splitCigar(allele.cigar);
        if (cigarV.back().second == "D"
            || cigarV.back().second == "I"
            || cigarV.front().second == "D"
            || cigarV.front().second == "I") {
            return true;
        } else {
            return false;
        }
    }
}

bool isEmptyAlleleOrIsDividedIndel(const Allele& allele) {
    return isEmptyAllele(allele) || isDividedIndel(allele);
}
