// utility functions
//
#ifndef _UTILITY_H
#define _UTILITY_H

#include <math.h>
#include <vector>
#include <string>
#include <algorithm>
#include <string>
#include <float.h>
#include <iostream>

using namespace std;

long double factorial(int);
short qualityChar2ShortInt(char c);
char qualityInt2Char(short i);
//long double phred2float(int qual);
long double phred2ln(int qual);
int ln2phred(long double prob);
long double phred2float(int qual);
int float2phred(long double prob);
long double powln(long double m, int n);
// here 'joint' means 'probability that we have a vector entirely composed of true bases'
short jointQuality(const std::vector<short>& quals);
short jointQuality(const std::string& qualstr);
std::vector<short> qualities(const std::string& qualstr);
// 
short averageQuality(const std::string& qualstr);
//unsigned int factorial(int n);
bool stringInVector(string item, vector<string> items);

long double gammaln( long double x);
long double factorial( int n);
long double factorialln( int n);
long double cofactor( int n, int i);
long double cofactorln( int n, int i);

long double safe_exp(long double ln);
long double logsumexp_probs(const vector<long double>& lnv);
long double logsumexp(const vector<long double>& lnv);

long double betaln(const vector<long double>& alphas);
long double beta(const vector<long double>& alphas);



#endif
