// utility functions
//
#ifndef _UTILITY_H
#define _UTILITY_H

#include <math.h>
#include <vector>
#include <string>

short qualityChar2ShortInt(char c);
//long double phred2float(int qual);
long double phred2float(int qual);
short float2phred(long double prob);
// here 'joint' means 'probability that we have a vector entirely composed of true bases'
short jointQuality(const std::vector<short>& quals);
short jointQuality(const std::string& qualstr);
// 
short averageQuality(const std::string& qualstr);
//unsigned int factorial(int n);

#endif
