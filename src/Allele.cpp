#include "Allele.h"

Allele::Allele(AlleleType t, 
        ReferenceID refID, 
        Position pos, 
        int len, 
        string alt, 
        SampleID sample, 
        Strand strnd, 
        short qual) 
    : type(t)
    , referenceID(refID)
    , position(pos)
    , length(len)
    , alternate(alt)
    , sampleID(sample)
    , strand(strnd)
    , quality(qual) 
{ }

string Allele::Type(void) {

    string t;

    switch (this->type) {
        case ALLELE_REFERENCE:
            t = "reference";
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
        << "reference name: " << allele.referenceName << endl
        << "position: " << allele.position << endl
        << "sequence id: " << allele.referenceID << endl
        << "alternate sequence: " << allele.alternate << endl
        << "strand: " << (allele.strand == STRAND_FORWARD ? "+" : "-") << endl
        << "sample id: " << allele.sampleID << endl;

    return out;
}
