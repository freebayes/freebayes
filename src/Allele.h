#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#ifndef _ALLELE_H
#define _ALLELE_H

using namespace std;

// a structure describing an allele
//

enum AlleleType {
    ALLELE_REFERENCE,
    ALLELE_SNP, 
    ALLELE_INSERTION, 
    ALLELE_DELETION,
    ALLELE_CNV
};

enum Strand {
    STRAND_FORWARD,
    STRAND_REVERSE
};

typedef string SampleID;
typedef unsigned int ReferenceID;
typedef long unsigned int Position;

class Allele {

    friend ostream &operator<<(ostream &out, Allele &a);

public:

    // type
    AlleleType type;
    // reference
    string referenceName;
    // sequence / chromosome
    ReferenceID referenceID;
    // alternate sequence
    string alternate;
    // position 1-based against reference
    Position position;
    // and event length (deletion implies 0, snp implies 1, insertion >1)
    int length;
    // strand
    Strand strand;
    // representative sample ID
    SampleID sampleID;
    // supporting reads
    //vector<BamAlignment*> supportingAlignments;
    short quality; // quality score associated with this allele

    Allele(AlleleType type, 
           ReferenceID referenceID,
           Position position, 
           int length,
           string alternate, 
           SampleID sampleid,
           Strand strand,
           short quality);

    // return a string representation of the type, for output
    string Type(void);

};

#endif
