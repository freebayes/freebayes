#ifndef BIAS_CPP
#define BIAS_CPP

#include "Bias.h"
#include "convert.h"

void Bias::open(string& file) {

    ifstream input;
    input.open(file.c_str());
    if (!input.is_open()) {
        cerr << "allele reference bias description " << file << " is not open" << endl;
        exit(1);
    }

    string line;
    bool firstrecord = true;
    int last;
    while (std::getline(input, line)) {
        vector<string> fields = split(line, " \t");
        if (firstrecord) { 
            convert(fields[0], minLength);
            last = minLength - 1;
        }
        convert(fields[0], maxLength);
        if (maxLength != last + 1) {
            cerr << "gap or out-of-order bias list in " << file << endl;
            cerr << line << endl;
            exit(1);
        } else {
            last = maxLength;
        }
        long double dbias;
        convert(fields[1], dbias);
        biases.push_back(dbias);
    }
    input.close();
}

long double Bias::bias(int length) {
    if (biases.empty()) return 1; // no bias
    if (length < minLength) {
        return biases.front();
    } else if (length > maxLength) {
        return biases.back();
    } else {
        return biases.at(length - minLength);
    }
}

bool Bias::empty(void) {
    return biases.empty();
}


#endif
