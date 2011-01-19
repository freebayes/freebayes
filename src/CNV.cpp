#include "CNV.h"

bool CNVMap::load(int defploidy, string const& filename) {
    defaultPloidy = defploidy;
    string line;
    ifstream cnvFile(filename.c_str(), ios::in);
    if (cnvFile.is_open()) {
        while (getline (cnvFile, line)) {
            vector<string> fields = split(line, " \t");
            setPloidy(fields.at(3), fields.at(0), atol(fields.at(1).c_str()), atol(fields.at(2).c_str()), atoi(fields.at(4).c_str()));
        }
    } else {
        return false;
    }
    return true;
}

void CNVMap::setPloidy(string const& sample, string const& seq, long int start, long int end, int ploidy) {
    sampleSeqCNV[sample][seq][make_pair(start, end)] = ploidy;
}

int CNVMap::ploidy(string const& sample, string const& seq, long int position) {
    map<pair<long int, long int>, int>& cnvs = sampleSeqCNV[sample][seq];
    for (map<pair<long int, long int>, int>::iterator i = cnvs.begin(); i != cnvs.end(); ++i) {
        pair<long int, long int> range = i->first;
        int copyNumber = i->second;
        if (range.first <= position && range.second > position) {
            return copyNumber;
        } else if (position < range.first) {
            // we've passed any potential matches in this sequence, and the map
            // is sorted by pair, so we don't have any matching ranges
            break;
        }
    }
    return defaultPloidy;
}
