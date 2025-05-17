#ifndef FREEBAYES_BIAS_H
#define FREEBAYES_BIAS_H

#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <fastahack/split.h>

using namespace std;

class Bias {
    
    int minLength;
    int maxLength;
    vector<long double> biases;

public:

    Bias(void) : minLength(0), maxLength(0) { }
    void open(string& file);
    long double bias(int length);
    bool empty(void);

};


#endif
