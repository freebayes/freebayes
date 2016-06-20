#ifndef __NONCALL_H
#define __NONCALL_H

#include <string>
#include <vector>
#include <set>
#include <map>
#include <utility>
#include "Utility.h"
#include "Allele.h"
#include "Sample.h"

using namespace std;

class NonCall {
public:
    NonCall(void)
        : refCount(0)
        , reflnQ(0)
        , altCount(0)
        , altlnQ(0)
        , minDepth(0)
        , nCount(0)

    { }
    NonCall(int rc, long double rq, int ac, long double aq, int mdp)
        : refCount(rc)
        , reflnQ(rq)
        , altCount(ac)
        , altlnQ(aq)
        , minDepth(mdp)
    { }
    int refCount;
    int altCount;
    int minDepth;
    int nCount  ; 
    long double reflnQ;
    long double altlnQ;
};

class NonCalls : public map<string, map<long, map<string, NonCall> > > {
public:
    void record(const string& seqName, long pos, const Samples& samples);
    NonCall aggregateAll(void);
    void aggregatePerSample(map<string, NonCall>& perSite);
    pair<string, long> firstPos(void);
    pair<string, long> lastPos(void);
};

#endif
