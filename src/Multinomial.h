#ifndef __MULTINOMIAL_H
#define __MULTINOMIAL_H

#include "Utility.h"
#include <vector>

long double multinomial(vector<long double> probs, vector<int> obs);
long double multinomialln(vector<long double> probs, vector<int> obs);

#endif
