#ifndef BEDREADER_CPP
#define BEDREADER_CPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <map>
#include <iterator>
#include <algorithm>
#include <cmath>
#include "split.h"
#include "BedReader.h"

using namespace std;

vector<BedTarget> BedReader::entries(void) {

    vector<BedTarget> entries;

    if (!is_open()) {
        cerr << "bed targets file is not open" << endl;
        exit(1);
    }

    string line;
    while (std::getline(*this, line)) {
        vector<string> fields = split(line, '\t');
        BedTarget entry(fields[0], atoi(fields[1].c_str()), atoi(fields[2].c_str()), fields[3]);
        entries.push_back(entry);
    }

    return entries;

}

#endif
