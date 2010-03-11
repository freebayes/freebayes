#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

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
    STRAND_PLUS,
    STRAND_MINUS
};

typedef unsigned int SampleID;
typedef unsigned int SequenceID;
typedef long unsigned int Position;

class Allele {

    friend ostream &operator<<(ostream &out, Allele &a);

public:

    // type
    AlleleType type;
    // reference
    string referenceName;
    // sequence / chromosome
    SequenceID sequenceID;
    // alternate sequence
    string alternate;
    // position
    Position position;
    // strand
    Strand strand;
    // representative sample ID
    SampleID sampleID;

    Allele(AlleleType type, 
           string referenceName, 
           SequenceID sequenceid,
           Position position, 
           string alternate, 
           Strand strand, 
           SampleID sampleid);

    // return a string representation of the type, for output
    string Type(void);

};
