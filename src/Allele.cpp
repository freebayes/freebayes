#include "Allele.h"


Allele::Allele(AlleleType type, 
           ReferenceID referenceID,
           Position position, 
           int length,
           string alternate, 
           SampleID sampleid,
           Strand strand,
           short quality);
    : type(type)
    , referenceID(referenceID)
    , position(position)
    , length(length)
    , alternate(alternate)
    , sampleID(sampleID)
    , strand(strand)
    , quality(quality)
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
        << "sequence id: " << allele.sequenceID << endl
        << "alternate sequence: " << allele.alternate << endl
        << "strand: " << (allele.strand == STRAND_PLUS ? "plus" : "minus") << endl
        << "sample id: " << allele.sampleID << endl;

    return out;
}
