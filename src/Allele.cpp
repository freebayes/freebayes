#include "Allele.h"

Allele::Allele(AlleleType t, 
        string refname, 
        Position pos, 
        unsigned int len, 
        string refallele,
        string alt, 
        string sample, 
        bool strnd, 
        short qual,
        string qstr,
        short mapqual) 
    : type(t)
    , referenceName(refname)
    , position(pos)
    , length(len)
    , referenceAllele(refallele)
    , alternateAllele(alt)
    , sampleID(sample)
    , strand(strnd ? STRAND_FORWARD : STRAND_REVERSE)
    , quality((qual == -1) ? averageQuality(qstr) : qual) // passing -1 as quality triggers this calculation
    , qualityString(qstr)
    , mapQuality(mapqual) 
{ }

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

    out << "type: " << allele.Type() << endl
        << "sequence name: " << allele.referenceName << endl
        << "position: " << allele.position << endl
        << "strand: " << (allele.strand == STRAND_FORWARD ? "+" : "-") << endl
        << "length: " << allele.length << endl
        << "reference allele: " << allele.referenceAllele << endl
        << "alternate allele: " << allele.alternateAllele << endl
        << "sample id: " << allele.sampleID << endl
        << "quality: " << allele.quality << endl
        << "quality string: " << allele.qualityString << endl
        << "read mapping quality: " << allele.mapQuality;

    return out;
}
