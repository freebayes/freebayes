#ifndef CONTAMINATION_CPP
#define CONTAMINATION_CPP

#include "Contamination.h"
#include "convert.h"

void Contamination::open(string& file) {

    ifstream input;
    input.open(file.c_str());
    if (!input.is_open()) {
        cerr << "contamination estimates file " << file << " is not open" << endl;
        exit(1);
    }

    string line;
    int last;
    while (std::getline(input, line)) {
        vector<string> fields = split(line, " \t");
        if (fields.size() != 3) {
            cerr << "could not parse contamination estimate:" << endl
                 << line << endl
                 << "should be of the form:" << endl
                 << "sample p(read=R|genotype=AR) p(read=A|genotype=AA)" << endl;
            exit(1);
        }
        string sample = fields[0];
        ContaminationEstimate c;
        convert(fields[1], c.probRefGivenHet);
        convert(fields[2], c.probRefGivenHomAlt);
        if (sample == "*") { // default
            defaultEstimate = c;
        } else {
            insert(make_pair(sample, c));
        }
    }
    input.close();
}

double Contamination::probRefGivenHet(string& sample) {
    Contamination::iterator s = find(sample);
    if (s != end()) {
        return s->second.probRefGivenHet;
    } else {
        return defaultEstimate.probRefGivenHet;
    }
}

double Contamination::probRefGivenHomAlt(string& sample) {
    Contamination::iterator s = find(sample);
    if (s != end()) {
        return s->second.probRefGivenHomAlt;
    } else {
        return defaultEstimate.probRefGivenHomAlt;
    }
}

double Contamination::refBias(string& sample) {
    Contamination::iterator s = find(sample);
    if (s != end()) {
        return s->second.refBias;
    } else {
        return defaultEstimate.refBias;
    }
}

ContaminationEstimate& Contamination::of(string& sample) {
    Contamination::iterator s = find(sample);
    if (s != end()) {
        return s->second;
    } else {
        return defaultEstimate;
    }
}

#endif
