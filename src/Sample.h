#ifndef __SAMPLE_H
#define __SAMPLE_H

#include <string>
#include <vector>
#include <map>
#include <utility>
#include "Utility.h"
#include "Allele.h"

using namespace std;

class StrandBaseCounts {

public:
    int forwardRef;
    int forwardAlt;
    int reverseRef;
    int reverseAlt;

StrandBaseCounts(void)
    : forwardRef(0), forwardAlt(0), reverseRef(0), reverseAlt(0)
    { }

StrandBaseCounts(int fr,
                 int fa,
                 int rr,
                 int ra)
    : forwardRef(fr)
        , forwardAlt(fa)
        , reverseRef(rr)
        , reverseAlt(ra) { }

};

// sample tracking and allele sorting
class Sample : public map<string, vector<Allele*> > {

    friend ostream& operator<<(ostream& out, Sample& sample);

public:

    // includes both fully and partially-supported observations after adding partial obs
    set<string> supportedAlleles;
    void setSupportedAlleles(void);

    // partial support for alleles, such as for observations that partially overlap the calling window
    map<string, vector<Allele*> > partialSupport;

    // for fast scaling of qualities for partial supports
    map<Allele*, set<Allele*> > reversePartials;

    // clear the above
    void clearPartialObservations(void);

    // set of partial observations (keys of the above map) cached for faster GL calculation
    //vector<Allele*> partialObservations;

    // if the observation (partial or otherwise) supports the allele
    bool observationSupports(Allele* obs, Allele* allele);

    // the number of observations for this allele
    int observationCount(Allele& allele);
    double observationCountInclPartials(Allele& allele);
    double partialObservationCount(Allele& allele);

    // the number of observations for this base
    int observationCount(const string& base);
    double observationCountInclPartials(const string& base);
    double partialObservationCount(const string& base);
    int partialObservationCount(void);

    // the total number of observations
    int observationCount(void);
    int observationCountInclPartials(void);

    // sum of quality for the given allele
    // (includes partial support)
    int qualSum(Allele& allele);
    int qualSum(const string& base);
    double partialQualSum(Allele& allele);
    double partialQualSum(const string& base);

    // puts alleles into the right bins if they have changed their base (as
    // occurs in the case of reference alleles)
    void sortReferenceAlleles(void);

    StrandBaseCounts strandBaseCount(string refbase, string altbase);

    int baseCount(string base, AlleleStrand strand);

    string json(void);

};

class Samples : public map<string, Sample> {
public:
    map<string, double> estimatedAlleleFrequencies(void);
    void assignPartialSupport(vector<Allele>& alleles,
                              vector<Allele*>& partialObservations,
                              map<string, vector<Allele*> >& partialObservationGroups,
                              map<Allele*, set<Allele*> >& partialObservationSupport,
                              unsigned long haplotypeStart,
                              int haplotypeLength);

    int observationCount(void);
    double observationCountInclPartials(void);

    int observationCount(Allele& allele);
    double observationCountInclPartials(Allele& allele);
    double partialObservationCount(Allele& allele);

    int observationCount(const string& base);
    double observationCountInclPartials(const string& base);
    double partialObservationCount(const string& base);

    int qualSum(Allele& allele);
    int qualSum(const string& base);
    double partialQualSum(Allele& allele);
    double partialQualSum(const string& base);

    void clearFullObservations(void);
    void clearPartialObservations(void);
    void setSupportedAlleles(void);
};

int countAlleles(Samples& samples);
// using this one...
void groupAlleles(Samples& samples, map<string, vector<Allele*> >& alleleGroups);

// filters... maybe move to its own file?
bool sufficientAlternateObservations(Samples& observations, int mincount, float minfraction);


#endif
