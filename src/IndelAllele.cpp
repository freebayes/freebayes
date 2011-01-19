#include "IndelAllele.h"

using namespace std;


bool IndelAllele::homopolymer(void) {
    string::iterator s = sequence.begin();
    char c = *s++;
    while (s != sequence.end()) {
        if (c != *s++) return false;
    }
    return true;
}

bool homopolymer(string sequence) {
    string::iterator s = sequence.begin();
    char c = *s++;
    while (s != sequence.end()) {
        if (c != *s++) return false;
    }
    return true;
}

ostream& operator<<(ostream& out, const IndelAllele& indel) {
    string t = indel.insertion ? "i" : "d";
    out << t <<  ":" << indel.position << ":" << indel.readPosition << ":" << indel.sequence;
    return out;
}

bool operator==(const IndelAllele& a, const IndelAllele& b) {
    return (a.insertion == b.insertion
            && a.length == b.length
            && a.position == b.position
            && a.sequence == b.sequence);
}

bool operator!=(const IndelAllele& a, const IndelAllele& b) {
    return !(a==b);
}

bool operator<(const IndelAllele& a, const IndelAllele& b) {
    ostringstream as, bs;
    as << a;
    bs << b;
    return as.str() < bs.str();
}
