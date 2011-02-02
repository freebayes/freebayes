#ifndef BEDREADER_H
#define BEDREADER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <map>
#include <iterator>
#include <algorithm>

using namespace std;

// stores the posiitional information of a bed target entry
class BedTarget {

public:

    string seq;  // sequence name
    int left;    // left position
    int right;   // right position, adjusted to 0-base
    string desc; // descriptive information, target name typically

    BedTarget(string s, int l, int r, string d = "")
        : seq(s)
        , left(l)
        , right(r)
        , desc(d)
    { }

};


class BedReader : public ifstream {

public:
    BedReader(string& fname) {
        open(fname.c_str());
    }
    vector<BedTarget> entries(void);

};

#endif

