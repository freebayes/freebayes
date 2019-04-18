#include "NonCall.h"

NonCall NonCalls::aggregateAll(void) {
    NonCall aggregate;
    bool first = true;
    for (NonCalls::const_iterator nc = this->begin(); nc != this->end(); ++nc) {
        for (map<long, map<string, NonCall> >::const_iterator p
            = nc->second.begin(); p != nc->second.end(); ++p) {
            for (map<string, NonCall>::const_iterator s = p->second.begin();
                 s != p->second.end(); ++s) {
                const NonCall& nonCall = s->second;
                aggregate.refCount += nonCall.refCount;
                aggregate.altCount += nonCall.altCount;
                aggregate.reflnQ += nonCall.reflnQ;
                aggregate.altlnQ += nonCall.altlnQ;
                aggregate.nCount += 1;
                if (first) {
                      aggregate.minDepth = nonCall.refCount + nonCall.altCount;
                      first = false;
                } else {
                    aggregate.minDepth = min(aggregate.minDepth,
                      nonCall.refCount + nonCall.altCount);
                  }
            }
        }
    }
    return aggregate;
}

void NonCalls::aggregatePerSample(map<string, NonCall>& perSample) {
    set<string> seen;
    for (NonCalls::const_iterator nc = this->begin(); nc != this->end(); ++nc) {
        for (map<long, map<string, NonCall> >::const_iterator p = nc->second.begin();
             p != nc->second.end(); ++p) {
            for (map<string, NonCall>::const_iterator s = p->second.begin();
                 s != p->second.end(); ++s) {
                const string& name = s->first;
                const NonCall& nonCall = s->second;
                NonCall& aggregate = perSample[name];
                aggregate.refCount += nonCall.refCount;
                aggregate.altCount += nonCall.altCount;
                aggregate.reflnQ += nonCall.reflnQ;
                aggregate.altlnQ += nonCall.altlnQ;
                aggregate.nCount += 1;
                if (!seen.count(name)) {
                    aggregate.minDepth = nonCall.refCount + nonCall.altCount;
                    seen.insert(name);
                } else {
                    aggregate.minDepth = min(aggregate.minDepth, nonCall.refCount + nonCall.altCount);
                }
            }
        }
    }
}

void NonCalls::record(const string& seqName, long pos, const Samples& samples) {
    map<string, NonCall>& site = (*this)[seqName][pos];
    for (Samples::const_iterator s = samples.begin(); s != samples.end(); ++s) {
        // tally ref and non-ref alleles
        const string& name = s->first;
        const Sample& sample = s->second;
        NonCall& noncall = site[name];
        for (Sample::const_iterator a = sample.begin(); a != sample.end(); ++a) {
            const vector<Allele*>& alleles = a->second;
            for (vector<Allele*>::const_iterator o = alleles.begin(); o != alleles.end(); ++o) {
                Allele& allele = **o;
                if (allele.isReference()) {
                    ++noncall.refCount;
                    noncall.reflnQ += allele.lnquality;
                } else {
                    ++noncall.altCount;
                    noncall.altlnQ += allele.lnquality;
                }
            }
        }
    }
}

pair<string, long> NonCalls::firstPos(void) {
    const string& startChrom = begin()->first;
    long startPos = begin()->second.begin()->first;
    return make_pair(startChrom, startPos);
}

pair<string, long> NonCalls::lastPos(void) {
    const string& endChrom = rbegin()->first;
    long endPos = rbegin()->second.rbegin()->first;
    return make_pair(endChrom, endPos);
}
