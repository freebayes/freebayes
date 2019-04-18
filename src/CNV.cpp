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
    sampleSeqCNV[sample][seq][make_pair(start, end)] = ploidy;
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
        map<string, map<pair<long int, long int>, int> >::iterator c = scnv->second.find(seq);
        if (c == scnv->second.end()) {
            return basePloidy;
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
            return basePloidy;
        }
    }

}
