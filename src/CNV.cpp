#include "CNV.h"

bool CNVMap::load(string const& filename) {
    string line;
    ifstream cnvFile(filename.c_str(), ios::in);
    if (cnvFile.is_open()) {
        while (getline (cnvFile, line)) {
            vector<string> fields = split(line, " \t");
            // note conversion between 1 and 0 based
            setPloidy(fields.at(3), fields.at(0), atol(fields.at(1).c_str()), atol(fields.at(2).c_str()), atoi(fields.at(4).c_str()));
        }
    } else {
        return false;
    }
    return true;
}

void CNVMap::setDefaultPloidy(int defploidy) {
    defaultPloidy = defploidy;
}

void CNVMap::setPloidy(string const& sample, string const& seq, long int start, long int end, int ploidy) {
    sampleSeqCNV[sample][seq][make_pair(start, end)] = ploidy;
}

int CNVMap::ploidy(string const& sample, string const& seq, long int position) {

    if (sampleSeqCNV.empty()) {
        return defaultPloidy;
    }

    SampleSeqCNVMap::iterator scnv = sampleSeqCNV.find(sample);

    if (scnv == sampleSeqCNV.end()) {
        return defaultPloidy;
    } else {
        map<string, map<pair<long int, long int>, int> >::iterator c = scnv->second.find(seq);
        if (c == scnv->second.end()) {
            return defaultPloidy;
        } else {
            map<pair<long int, long int>, int>& cnvs = c->second;
            for (map<pair<long int, long int>, int>::iterator i = cnvs.begin(); i != cnvs.end(); ++i) {
                pair<long int, long int> range = i->first;
                int copyNumber = i->second;
                if (range.first <= position && range.second > position) {
                    return copyNumber;
                } else if (position > range.first && position > range.second) {
                    // we've passed any potential matches in this sequence, and the map
                    // is sorted by pair, so we don't have any matching ranges
                    break;
                }
            }
            return defaultPloidy;
        }
    }

}
