// utility functions
//
#ifndef _UTILITY_H
#define _UTILITY_H

#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <string>
#include <float.h>
#include <iostream>
#include <map>
#include <time.h>
#include "convert.h"

using namespace std;

long double factorial(int);
short qualityChar2ShortInt(char c);
long double qualityChar2LongDouble(char c);
long double lnqualityChar2ShortInt(char c);
char qualityInt2Char(short i);
//long double phred2float(int qual);
long double phred2ln(int qual);
long double ln2phred(long double prob);
long double ln2log10(long double prob);
long double log102ln(long double prob);
long double phred2float(int qual);
long double float2phred(long double prob);
long double powln(long double m, int n);
// here 'joint' means 'probability that we have a vector entirely composed of true bases'
long double jointQuality(const std::vector<short>& quals);
long double jointQuality(const std::string& qualstr);
std::vector<short> qualities(const std::string& qualstr);
// 
long double sumQuality(const std::string& qualstr);
long double averageQuality(const std::string& qualstr);
//unsigned int factorial(int n);
bool stringInVector(string item, vector<string> items);
int upper(int c); // helper to below, wraps toupper
string uppercase(string s);
bool allATGC(string& s);
string strip(string const& str, char const* separators = " \t");

int binomialCoefficient(int n, int k);
long double binomialProb(int k, int n, long double p);
long double binomialProbln(int k, int n, long double p);

long double poissonpln(int observed, int expected);
long double poissonp(int observed, int expected);
long double poissonPvalLn(int a, int b);

long double gammaln( long double x);
long double factorial( int n);
long double factorialln( int n);
long double __factorialln( int n);

class FactorialCache : public map<int, long double> {
public:
    long double factorialln(int n) {
        map<int, long double>::iterator f = find(n);
        if (f == end()) {
            long double fln = __factorialln(n);
            insert(make_pair(n, fln));
            return fln;
        } else {
            return f->second;
        }
    }
};

long double cofactor( int n, int i);
long double cofactorln( int n, int i);

long double harmonicSum(int n);

long double safedivide(long double a, long double b);

long double safe_exp(long double ln);
long double logsumexp_probs(const vector<long double>& lnv);
long double logsumexp(const vector<long double>& lnv);

long double betaln(const vector<long double>& alphas);
long double beta(const vector<long double>& alphas);

long double hoeffding(double successes, double trials, double prob);
long double hoeffdingln(double successes, double trials, double prob);

int levenshteinDistance(const std::string source, const std::string target);
bool isTransition(string& ref, string& alt);

string dateStr(void);

long double string2float(const string& s);
long double log10string2ln(const string& s);

#endif
