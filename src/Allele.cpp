#include "Allele.h"

Allele::Allele(AlleleType t, 
        ReferenceID refID, 
        string refname, 
        Position pos, 
        int len, 
        string refallele,
        string alt, 
        SampleID sample, 
        bool strnd, 
        short qual,
        short mapqual) 
    : type(t)
    , referenceID(refID)
    , referenceName(refname)
    , position(pos)
    , length(len)
    , referenceAllele(refallele)
    , alternateAllele(alt)
    , sampleID(sample)
    , strand(strnd ? STRAND_FORWARD : STRAND_REVERSE)
    , quality(qual) 
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
        << "sequence id: " << allele.referenceID << endl
        << "sequence name: " << allele.referenceName << endl
        << "position: " << allele.position << endl
        << "strand: " << (allele.strand == STRAND_FORWARD ? "+" : "-") << endl
        << "length: " << allele.length << endl
        << "reference allele: " << allele.referenceAllele << endl
        << "alternate allele: " << allele.alternateAllele << endl
        << "sample id: " << allele.sampleID << endl
        << "quality: " << allele.quality << endl
        << "read mapping quality: " << allele.mapQuality;

    return out;
}
