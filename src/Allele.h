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
    ALLELE_MISMATCH, 
    ALLELE_SNP, 
    ALLELE_INSERTION, 
    ALLELE_DELETION,
    ALLELE_CNV
};

enum Strand {
    STRAND_FORWARD,
    STRAND_REVERSE
};

typedef long unsigned int Position;

class Allele {

    friend ostream &operator<<(ostream &out, Allele &a);

public:

    // type
    AlleleType type;
    // reference
    string referenceName; // TODO, set me
    // reference sequence at this allele
    string referenceAllele;
    // alternate sequence
    string alternateAllele;
    // position 1-based against reference
    Position position;
    // and event length (deletion implies 0, snp implies 1, insertion >1)
    int length;
    // strand
    Strand strand;
    // representative sample ID
    string sampleID;
    // supporting reads
    //vector<BamAlignment*> supportingAlignments;
    short quality; // base quality score associated with this allele
    short mapQuality; // map quality for the originating read

    Allele(AlleleType type, 
           string refname,
           Position position, 
           int length,
           string refseq,
           string alt, 
           string sampleid,
           bool strand,
           short quality,
           short mapQual);

    // return a string representation of the type, for output
    string Type(void);

};

#endif
