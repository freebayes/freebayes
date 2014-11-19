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
#include "Utility.h"
#include "BedReader.h"
#include "../intervaltree/IntervalTree.h"

using namespace std;

vector<BedTarget> BedReader::entries(void) {

    vector<BedTarget> entries;

    if (!is_open()) {
        cerr << "bed targets file is not open" << endl;
        exit(1);
    }

    string line;
    while (std::getline(*this, line)) {
        // BED is base-numbered, 0-origin, half-open.  This parse turns that
        // into base-numbered, 0-origin, fully-closed for internal use.  All
        // coordinates used internally should be in the latter, and coordinates
        // from the user in the former should be converted immediately to the
        // internal format.
        vector<string> fields = split(line, " \t");
        BedTarget entry(strip(fields[0]),
                        atoi(strip(fields[1]).c_str()),
                        atoi(strip(fields[2]).c_str()) - 1, // use inclusive format internally
                        (fields.size() >= 4) ? strip(fields[3]) : "");
        entries.push_back(entry);
    }

    return entries;

}

bool BedReader::targetsContained(string& seq, long left, long right) {
    vector<Interval<BedTarget*> > results;
    intervals[seq].findContained(left, right, results);
    return !results.empty();
}

bool BedReader::targetsOverlap(string& seq, long left, long right) {
    vector<Interval<BedTarget*> > results;
    intervals[seq].findOverlapping(left, right, results);
    return !results.empty();
}

vector<BedTarget*> BedReader::targetsContaining(BedTarget& target) {
    vector<Interval<BedTarget*> > results;
    intervals[target.seq].findContained(target.left, target.right, results);
    vector<BedTarget*> contained;
    for (vector<Interval<BedTarget*> >::iterator r = results.begin(); r != results.end(); ++r) {
        contained.push_back(r->value);
    }
    return contained;
}

vector<BedTarget*> BedReader::targetsOverlapping(BedTarget& target) {
    vector<Interval<BedTarget*> > results;
    intervals[target.seq].findOverlapping(target.left, target.right, results);
    vector<BedTarget*> overlapping;
    for (vector<Interval<BedTarget*> >::iterator r = results.begin(); r != results.end(); ++r) {
        overlapping.push_back(r->value);
    }
    return overlapping;
}

#endif
