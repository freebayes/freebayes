// utility functions
//
#ifndef _UTILITY_H
#define _UTILITY_H

#include <math.h>
#include <vector>

short qualityChar2ShortInt(char c);
double phred2float(int qual);
short float2phred(double prob);
short jointQuality(const std::vector<short>& quals);

#endif
