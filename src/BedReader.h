#ifndef FREEBAYES_BEDREADER_H
#define FREEBAYES_BEDREADER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <map>
#include <iterator>
#include <algorithm>
#include <intervaltree/IntervalTree.h>
#include "split.h"

using namespace std;

// stores the posiitional information of a bed target entry
class BedTarget {

public:

    string seq;  // sequence name
    int left;    // left position
    int right;   // right position, adjusted to 0-base inclusive
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
    vector<BedTarget> targets;
    map<string, IntervalTree<int, BedTarget*> > intervals; // intervals by reference sequence

    vector<BedTarget> entries(void);

    bool targetsContained(string& seq, long left, long right);
    bool targetsOverlap(string& seq, long left, long right);
    vector<BedTarget*> targetsContaining(BedTarget& target);
    vector<BedTarget*> targetsOverlapping(BedTarget& target);

    BedReader(void) { }

    BedReader(string& fname) {
        openFile(fname);
        buildIntervals();
    }

    void openFile(string& fname) {
        open(fname.c_str());
        targets = entries();
    }

    void buildIntervals(void) {
        map<string, IntervalTree<int, BedTarget*>::interval_vector> intervalsBySeq;
        for (vector<BedTarget>::iterator t = targets.begin(); t != targets.end(); ++t) {
            intervalsBySeq[t->seq].push_back(Interval<int, BedTarget*>(t->left, t->right, &*t));
        }
        for (map<string, IntervalTree<int, BedTarget*>::interval_vector>::const_iterator s = intervalsBySeq.begin(); s != intervalsBySeq.end(); ++s) {
            intervals[s->first] = IntervalTree<int, BedTarget*>((IntervalTree<int, BedTarget*>::interval_vector&&)s->second);
        }
    }

};

#endif
