#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <limits>
#include <sstream>
#include "Utility.h"

#ifndef _ALLELE_H
#define _ALLELE_H

using namespace std;

// a structure describing an allele

enum AlleleType {
    ALLELE_NONE,
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

    friend string stringForAllele(Allele &a);
    friend string stringForAlleles(vector<Allele> &av);

    friend ostream &operator<<(ostream &out, vector<Allele> &a);

    friend ostream &operator<<(ostream &out, Allele &a);
    friend bool operator==(Allele &a, Allele &b);

public:

    AlleleType type;        // type of the allele, enumerated above
    string referenceName;   // reference name, for sanity checking
    string referenceSequence; // reference sequence or "" (in case of insertions)
    string alternateSequence; // alternate sequence or "" (in case of deletions and reference alleles)
    Position position;      // position 1-based against reference
    unsigned int length;    // and event length (deletion implies 0, snp implies 1, insertion >1)
    Strand strand;          // strand, true = +, false = -
    string sampleID;        // representative sample ID
    string readID;          // id of the read which the allele is drawn from
    string qualityString;   // quality string drawn from sequencer
    short quality;          // base quality score associated with this allele
    short mapQuality;       // map quality for the originating read
    bool genotypeAllele;    // if this is an abstract 'genotype' allele

    // default constructor, for converting alignments into allele observations
    Allele(AlleleType type, 
           string refname,
           Position position, 
           unsigned int length,
           string refseq,
           string alt, 
           string sampleid,
           string readid,
           bool strand,
           short qual,
           string qualityStr,
           short mapQual);

    // for constructing abstract 'genotype' alleles
    Allele(AlleleType t, string alt, unsigned int len, Position pos, bool gallele = true);

    bool equivalent(Allele &a);  // heuristic 'equivalency' between two alleles, which depends on their type
    string Type(void); // return a string representation of the allele type, for output
    short Quality(Position referencePosition);  // for getting the quality of a given position in multi-bp alleles
    bool sameSample(Allele &other);  // if the other allele has the same sample as this one

};

//typedef vector<Allele> Genotype; // unused

vector<vector<Allele> > groupAlleles(vector<Allele> &alleles, bool (*fncompare)(Allele &a, Allele &b));
bool allelesSameType(Allele &a, Allele &b);
bool allelesEquivalent(Allele &a, Allele &b);
bool allelesSameSample(Allele &a, Allele &b);
vector<Allele> genotypeAllelesFromAlleleGroups(vector<vector<Allele> > &groups);
Allele genotypeAllele(Allele& a);

#endif
