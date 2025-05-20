#pragma once

#include "allele.hpp"

#include <vector>
#include <list>

namespace vcflib {

using namespace std;

typedef vector<pair<int, char> > Cigar;
string varCigar(const vector<VariantAllele>& vav, bool xForMismatch = false);
string mergeCigar(const string& c1, const string& c2);
vector<pair<int, char> > splitUnpackedCigar(const string& cigarStr);
vector<pair<int, string> > old_splitCigar(const string& cigarStr);
vector<pair<int, char> > splitCigar(const string& cigarStr);
list<pair<int, char> > splitCigarList(const string& cigarStr);
int cigarRefLen(const vector<pair<int, char> >& cigar);
vector<pair<int, char> > cleanCigar(const vector<pair<int, char> >& cigar);
string old_joinCigar(const vector<pair<int, string> >& cigar);
string joinCigar(const vector<pair<int, char> >& cigar);
string joinCigarList(const list<pair<int, char> >& cigar);
bool isEmptyCigarElement(const pair<int, char>& elem);

}
