#ifndef __INDEL_ALLELE_H
#define __INDEL_ALLELE_H

#include <string>
#include <iostream>
#include <sstream>

using namespace std;

class IndelAllele {
    friend ostream& operator<<(ostream&, const IndelAllele&);
    friend bool operator==(const IndelAllele&, const IndelAllele&);
    friend bool operator!=(const IndelAllele&, const IndelAllele&);
    friend bool operator<(const IndelAllele&, const IndelAllele&);
public:
    bool insertion;
    int length;
    int position;
    string sequence;

    bool homopolymer(void);

    IndelAllele(bool i, int l, int p, string s)
        : insertion(i), length(l), position(p), sequence(s)
    { }
};

bool homopolymer(string sequence);
ostream& operator<<(ostream& out, const IndelAllele& indel);
bool operator==(const IndelAllele& a, const IndelAllele& b);
bool operator!=(const IndelAllele& a, const IndelAllele& b);
bool operator<(const IndelAllele& a, const IndelAllele& b);

#endif
