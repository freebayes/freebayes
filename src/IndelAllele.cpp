#include "IndelAllele.h"

using namespace std;


bool FBIndelAllele::homopolymer(void) {
    string::iterator s = sequence.begin();
    char c = *s++;
    while (s != sequence.end()) {
        if (c != *s++) return false;
    }
    return true;
}

bool FBhomopolymer(string sequence) {
    string::iterator s = sequence.begin();
    char c = *s++;
    while (s != sequence.end()) {
        if (c != *s++) return false;
    }
    return true;
}

ostream& operator<<(ostream& out, const FBIndelAllele& indel) {
    string t = indel.insertion ? "i" : "d";
    out << t <<  ":" << indel.position << ":" << indel.readPosition
        << ":" << indel.sequence << ":" << (indel.splice?"splice":"");
    return out;
}

bool operator==(const FBIndelAllele& a, const FBIndelAllele& b) {
    return (a.insertion == b.insertion
            && a.length == b.length
            && a.position == b.position
            && a.sequence == b.sequence
            && a.splice == b.splice);
}

bool operator!=(const FBIndelAllele& a, const FBIndelAllele& b) {
    return !(a==b);
}

bool operator<(const FBIndelAllele& a, const FBIndelAllele& b) {
    ostringstream as, bs;
    as << a;
    bs << b;
    return as.str() < bs.str();
}
