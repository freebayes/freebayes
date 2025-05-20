#include "vcflib/allele.hpp"

#include <iostream>
#include <tuple> // for std::tie

namespace vcflib {

using namespace std;

ostream& operator<<(ostream& out, const VariantAllele& var) {
    out << var.position << " " << var.ref << " -> " << var.alt;
    return out;
}

VariantAllele operator+(const VariantAllele& a, const VariantAllele& b) {
    return VariantAllele(a.ref + b.ref, a.alt + b.alt, a.position);
}

bool operator<(const VariantAllele& a, const VariantAllele& b) {
    return std::tie(a.position, a.ref, a.alt) < std::tie(b.position, b.ref, b.alt);

}

bool operator==(const VariantAllele& a, const VariantAllele& b) {
    return a.ref == b.ref && a.alt == b.alt && a.position == b.position;
}

bool VariantAllele::is_pure_indel(void) {
    return !ref.empty() && alt.empty() || !alt.empty() && ref.empty();
}

// shift 1bp in between the two variants to be in the "left" (first) allele
void shift_mid_left(VariantAllele& a, VariantAllele& b) {
    if (!b.is_pure_indel()) {
        a.alt.append(b.alt.substr(0,1));
        a.ref.append(b.ref.substr(0,1));
        b.alt = b.alt.substr(1);
        b.ref = b.ref.substr(1);
        ++b.position;
    } else {
        a.alt.append(b.alt);
        a.ref.append(b.ref);
        b.alt.clear();
        b.ref.clear();
        b.position = 0;
    }
}

// shift 1bp in between the two variants to be in the "right" (second) allele
void shift_mid_right(VariantAllele& a, VariantAllele& b) {
    if (!a.is_pure_indel()) {
        b.alt = a.alt.substr(a.alt.size()-1,1) + b.alt;
        b.ref = a.ref.substr(a.ref.size()-1,1) + b.ref;
        a.alt = a.alt.substr(0,a.alt.size()-1);
        a.ref = a.ref.substr(0,a.alt.size()-1);
        --b.position;
    } else {
        // a is pure indel
        // if del, the position will shift when merging
        if (!a.ref.empty() && a.alt.empty()) {
            b.position = a.position;
        }
        // else if pure ins, no change in position
        // but in any case we will combine
        b.alt = a.alt + b.alt;
        b.ref = a.ref + b.ref;
        a.alt.clear();
        a.ref.clear();
        a.position = 0;
    }
}

}
