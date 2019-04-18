#ifndef CONTAMINATION_H
#define CONTAMINATION_H

#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "split.h"

using namespace std;

class ContaminationEstimate {
public:
    double probRefGivenHet;
    double probRefGivenHomAlt;
    double refBias;
ContaminationEstimate(void) : probRefGivenHet(0.5), probRefGivenHomAlt(0), refBias(0) { }
ContaminationEstimate(double ra, double aa) : probRefGivenHet(ra), probRefGivenHomAlt(aa)
    {
        refBias = probRefGivenHet * 2 - 1;
    }
};

class Contamination : public map<string, ContaminationEstimate> {
public:
    ContaminationEstimate defaultEstimate;
    void open(string& file);
    double probRefGivenHet(string& sample);
    double probRefGivenHomAlt(string& sample);
    double refBias(string& sample);
    ContaminationEstimate& of(string& sample);
Contamination(void) : defaultEstimate(ContaminationEstimate(0.5, 0)) { }
Contamination(double ra, double aa) : defaultEstimate(ContaminationEstimate(ra, aa)) { }
};

#endif
