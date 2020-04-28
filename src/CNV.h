#ifndef __CNV_H
#define __CNV_H

#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <stdlib.h>
#include <algorithm>
#include <tuple>
#include "split.h"

using namespace std;

typedef map<string, map<string, vector<tuple<long int, long int, int> > > > SampleSeqCNVMap;

class CNVMap {

public:
    CNVMap(void) : defaultPloidy(2) { }
    void setDefaultPloidy(int defploidy);
    void setSamplePloidy(const string& sample, int ploidy);
    bool load(string const& filename);
    int ploidy(string const& sample, string const& seq, long int position);
    void setPloidy(string const& sample, string const& seq, long int start, long int end, int ploidy);

private:
    // note: this map is stored as 0-based, end position exclusive
    SampleSeqCNVMap sampleSeqCNV;
    int defaultPloidy;
    map<string, int> samplePloidy;

};

#endif
