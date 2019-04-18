#ifndef __INDEL_ALLELE_H
#define __INDEL_ALLELE_H

#include <string>
#include <iostream>
#include <sstream>

using namespace std;

class FBIndelAllele {
    friend ostream& operator<<(ostream&, const FBIndelAllele&);
    friend bool operator==(const FBIndelAllele&, const FBIndelAllele&);
    friend bool operator!=(const FBIndelAllele&, const FBIndelAllele&);
    friend bool operator<(const FBIndelAllele&, const FBIndelAllele&);
public:
    bool insertion;
    int length;
    int position;
    int readPosition;
    string sequence;
    bool splice;

    bool homopolymer(void);

    FBIndelAllele(bool i, int l, int p, int rp, string s, bool n)
    : insertion(i), length(l), position(p), readPosition(rp), sequence(s), splice(n)
    { }
};

bool FBhomopolymer(string sequence);
ostream& operator<<(ostream& out, const FBIndelAllele& indel);
bool operator==(const FBIndelAllele& a, const FBIndelAllele& b);
bool operator!=(const FBIndelAllele& a, const FBIndelAllele& b);
bool operator<(const FBIndelAllele& a, const FBIndelAllele& b);

#endif
