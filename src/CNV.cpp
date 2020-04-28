#include "CNV.h"

bool CNVMap::load(string const& filename) {
    string line;
    ifstream cnvFile(filename.c_str(), ios::in);
    if (cnvFile.is_open()) {
        while (getline (cnvFile, line)) {
            vector<string> fields = split(line, " \t");
            if (fields.size() == 2) {
                const string& sample = fields.at(0);
                int ploidy = atoi(fields.at(1).c_str());
                setSamplePloidy(sample, ploidy);
            } else if (fields.size() == 5) {
                // note conversion between 1 and 0 based
                const string& sample = fields.at(3);
                const string& seq = fields.at(0);
                long int start = atol(fields.at(1).c_str());
                long int end = atol(fields.at(2).c_str());
                int ploidy = atoi(fields.at(4).c_str());
                setPloidy(sample, seq, start, end, ploidy);
            } else {
                cerr << "error [freebayes::CNVMap] could not parse CNVmap line " << line << endl;
                exit(1);
            }
        }
    } else {
        return false;
    }
    return true;
}

void CNVMap::setDefaultPloidy(int defploidy) {
    defaultPloidy = defploidy;
}

void CNVMap::setSamplePloidy(const string& sample, int ploidy) {
    samplePloidy[sample] = ploidy;
}

void CNVMap::setPloidy(string const& sample, string const& seq, long int start, long int end, int ploidy) {
    sampleSeqCNV[sample][seq].push_back(make_tuple(start, end, ploidy));
}

int CNVMap::ploidy(string const& sample, string const& seq, long int position) {

    int basePloidy = (samplePloidy.find(sample) != samplePloidy.end()) ? samplePloidy[sample] : defaultPloidy;

    if (sampleSeqCNV.empty()) {
        return basePloidy;
    }

    SampleSeqCNVMap::iterator scnv = sampleSeqCNV.find(sample);

    if (scnv == sampleSeqCNV.end()) {
        return basePloidy;
    } else {
        map<string, vector<tuple<long int, long int, int> > >::iterator c = scnv->second.find(seq);
        if (c == scnv->second.end()) {
            return basePloidy;
        } else {
            vector<tuple<long int, long int, int> >& cnvs = c->second;

            vector<tuple<long int, long int, int> >::iterator i = upper_bound(cnvs.begin(), cnvs.end(),
                position, [](long int position, tuple<long int, long int, int> const& element) {
                    return position < get<1>(element);
                });
            if (i == cnvs.end()) {
                return basePloidy;
            }

            long int start = get<0>(*i);
            int copyNumber = get<2>(*i);
            if (start <= position) {
                return copyNumber;
            } else {
                return basePloidy;
            }
        }
    }

}
