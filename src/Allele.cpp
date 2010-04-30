#include "Allele.h"
#include "multichoose.h"

Allele::Allele(AlleleType t, 
            string refname, 
            Position pos, 
            unsigned int len, 
            string refallele, 
            string alt, 
            string sampleid,
            string readid,
            bool strnd, 
            short qual,
            string qstr, 
            short mapqual)
    : type(t)
    , referenceName(refname)
    , position(pos)
    , length(len)
    , referenceSequence(refallele)
    , alternateSequence(alt)
    , sampleID(sampleid)
    , readID(readid)
    , strand(strnd ? STRAND_FORWARD : STRAND_REVERSE)
    , quality((qual == -1) ? averageQuality(qstr) : qual) // passing -1 as quality triggers this calculation
    , qualityString(qstr)
    , mapQuality(mapqual) 
    , genotypeAllele(false)
{ }

// for constructing genotype alleles
Allele::Allele(AlleleType t,
        string alt,
        unsigned int len,
        Position pos,
        bool gallele) 
    : type(t)
    , alternateSequence(alt)
    , length(len)
    , quality(0)
    , position(pos)
    , genotypeAllele(true)
{ }


// quality at a given reference position
short Allele::Quality(Position referencePosition) {
    switch (this->type) {
        case ALLELE_REFERENCE:
            // XXX TODO!!!!
            //return qualityChar2ShortInt(*qualityString.substr(referencePosition - position, 1).c_str());
            //break;
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

ostream &operator<<(ostream &out, Allele &allele) {

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

    return out;
}

ostream &operator<<(ostream &out, vector<Allele> &alleles) {
    for (vector<Allele>::iterator allele = alleles.begin(); allele != alleles.end(); allele++) {
        out << *allele;
    }
    return out;
}

string stringForAlleles(vector<Allele> &av) {
    stringstream o;
    for (vector<Allele>::iterator a = av.begin(); a != av.end(); a++) {
        o << stringForAllele(*a) << "|";
    }
    string s = o.str();
    return s.substr(0,s.size()-1);
}

string stringForAllele(Allele &allele) {

    stringstream o;
    if (!allele.genotypeAllele) {
            o << allele.sampleID << ":" << allele.Type() << ":" 
                << allele.length << (allele.strand == STRAND_FORWARD ? "+" : "-")
                << ":" << allele.alternateSequence;
    } else {
        o << allele.Type();
        switch (allele.type) {
            case ALLELE_REFERENCE:
                break;
            default:
                o << ":" << allele.length << ":" << (string) allele.alternateSequence;
                break;
        }
    }
    return o.str();
}

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
                if (alternateSequence == b.alternateSequence && position == b.position)
                    return true;
                break;
            case ALLELE_DELETION:
                if (length == b.length && position == b.position)
                    return true;
                break;
            case ALLELE_INSERTION:
                if (length == b.length 
                    && position == b.position 
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

vector<vector<Allele> >  groupAlleles(vector<Allele> &alleles, bool (*fncompare)(Allele &a, Allele &b)) {
    vector<vector<Allele> > groups;
    for (vector<Allele>::iterator oa = alleles.begin(); oa != alleles.end(); ++oa) {
        bool unique = true;
        for (vector<vector<Allele> >::iterator ag = groups.begin(); ag != groups.end(); ++ag) {
            // is this just comparing allele pointers, or am i dereferencing properly?
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

bool Allele::sameSample(Allele &other) {
    return this->sampleID == other.sampleID;
}

bool allelesSameType(Allele &a, Allele &b) {
    return a.type == b.type;
}

bool allelesEquivalent(Allele &a, Allele &b) {
    return a.equivalent(b);
}

bool allelesSameSample(Allele &a, Allele &b) {
    return a.sampleID == b.sampleID;
}

vector<Allele> genotypeAllelesFromAlleleGroups(vector<vector<Allele> > &groups) {

    vector<Allele> results;
    for (vector<vector<Allele> >::iterator g = groups.begin(); g != groups.end(); g++)
        results.push_back(genotypeAllele(g->front()));
    return results;

}

Allele genotypeAllele(Allele &a) {
    return Allele(a.type, a.alternateSequence, a.length, a.position);
}

/*
vector<Allele> uniqueAlleleObservations(vector<vector<Allele> > &alleleObservations) {
}
*/
