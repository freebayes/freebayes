// utility functions
//
#include "Utility.h"
#include "Sum.h"
#include "Product.h"

#define PHRED_MAX 50000.0 // max Phred seems to be about 43015 (?), could be an underflow bug...

using namespace std;

short qualityChar2ShortInt(char c) {
    return static_cast<short>(c) - 33;
}

long double qualityChar2LongDouble(char c) {
    return static_cast<long double>(c) - 33;
}

long double lnqualityChar2ShortInt(char c) {
    return log(static_cast<short>(c) - 33);
}

char qualityInt2Char(short i) {
    return static_cast<char>(i + 33);
}

long double phred2ln(int qual) {
    return M_LN10 * qual * -.1;
}

long double ln2phred(long double prob) {
    return -10 * M_LOG10E * prob;
}

long double phred2float(int qual) {
    return pow(10, qual * -.1);
}

long double float2phred(long double prob) {
    if (prob == 1)
        return PHRED_MAX;
    long double p = -10 * (long double) log10(prob);
    if (p < 0 || p > PHRED_MAX) // int overflow guard
        return PHRED_MAX;
    else
        return p;
}

long double powln(long double m, int n) {
    return m * n;
}

// the probability that we have a completely true vector of qualities
long double jointQuality(const std::vector<short>& quals) {
    std::vector<long double> probs;
    for (int i = 0; i<quals.size(); ++i) {
        probs.push_back(phred2float(quals[i]));
    }
    // product of probability we don't have a true event for each element
    long double prod = 1 - probs.front();
    for (int i = 1; i<probs.size(); ++i) {
        prod *= 1 - probs.at(i);
    }
    // and then invert it again to get probability of an event
    return 1 - prod;
}

long double jointQuality(const std::string& qualstr) {

    long double jq = 1;
    // product of probability we don't have a true event for each element
    for (string::const_iterator q = qualstr.begin(); q != qualstr.end(); ++q) {
        jq *= 1 - phred2float(qualityChar2ShortInt(*q));
    }

    // and then invert it again to get probability of an event
    return 1 - jq;

}

std::vector<short> qualities(const std::string& qualstr) {

    std::vector<short> quals;
    for (int i=0; i<qualstr.size(); i++)
        quals.push_back(qualityChar2ShortInt(qualstr.at(i)));

    return quals;

}

long double sumQuality(const std::string& qualstr) {
    long double qual = 0;
    for (string::const_iterator q = qualstr.begin(); q != qualstr.end(); ++q)
            qual += qualityChar2LongDouble(*q);
    return qual;
}

// crudely averages quality scores in phred space
long double averageQuality(const std::string& qualstr) {

    long double qual = 0; //(long double) *max_element(quals.begin(), quals.end());
    for (string::const_iterator q = qualstr.begin(); q != qualstr.end(); ++q)
            qual += qualityChar2LongDouble(*q);
    return qual /= qualstr.size();

}

bool stringInVector(string item, vector<string> items) {
    for (vector<string>::iterator i = items.begin(); i != items.end(); ++i) {
        if (item == *i) {
            return true;
        }
    }
    return false;
}

bool allATGC(string& s) {
    for (string::iterator c = s.begin(); c != s.end(); ++c) {
        char b = *c;
        if (b != 'A' && b != 'T' && b != 'G' && b != 'C') {
            return false;
        }
    }
    return true;
}

int upper(int c) {
    return toupper((unsigned char) c);
}

string uppercase(string s) {
    transform(s.begin(), s.end(), s.begin(), upper);
    return s;
}

string strip(string const& str, char const* separators) {
    string::size_type const first = str.find_first_not_of(separators);
    return (first == string::npos) ? string()
        : str.substr(first, str.find_last_not_of(separators) - first + 1);
}


int binomialCoefficient(int n, int k) {
    int i = 1;
    int result = n - k + i++;
    while (i <= k) {
        result *= (n - k + i) / i;
        ++i;
    }
    return result;
}

// k successes in n trials with prob of success p
long double binomialProb(int k, int n, long double p) {
    return factorial(n) / (factorial(k) * factorial(n - k)) * pow(p, k) * pow(1 - p, n - k);
}

long double binomialProbln(int k, int n, long double p) {
    return factorialln(n) - (factorialln(k) + factorialln(n - k)) + powln(log(p), k) + powln(log(1 - p), n - k);
}

long double poissonpln(int observed, int expected) {
    return ((log(expected) * observed) - expected) - factorialln(observed);
}

long double poissonp(int observed, int expected) {
    return (double) pow((double) expected, (double) observed) * (double) pow(M_E, (double) -expected) / factorial(observed);
}


// given the expected number of events is the max of a and b
// what is the probability that we might observe less than the observed?
long double poissonPvalLn(int a, int b) {

    int expected, observed;
    if (a > b) {
        expected = a; observed = b;
    } else {
        expected = b; observed = a;
    }

    vector<long double> probs;
    for (int i = 0; i < observed; ++i) {
        probs.push_back(poissonpln(i, expected));
    }

    return logsumexp_probs(probs);

}


long double gammaln(
		     long double x
		     ) {

  long double cofactors[] = { 76.18009173, 
                              -86.50532033,
                              24.01409822,
                              -1.231739516,
                              0.120858003E-2,
                              -0.536382E-5 };    

  long double x1 = x - 1.0;
  long double tmp = x1 + 5.5;
  tmp -= (x1 + 0.5) * log(tmp);
  long double ser = 1.0;
  for (int j=0; j<=5; j++) {
    x1 += 1.0;
    ser += cofactors[j]/x1;
  }
  long double y =  (-1.0 * tmp + log(2.50662827465 * ser));

  return y;
}

long double factorial(
		      int n
		      ) {
  if (n < 0) {
    return (long double)0.0;
  }
  else if (n == 0) {
    return (long double)1.0;
  }
  else {
    return exp(gammaln(n + 1.0));
  }
}

long double factorialln(
			int n
			) {
  if (n < 0) {
    return (long double)-1.0;
  }
  else if (n == 0) {
    return (long double)0.0;
  }
  else {
    return gammaln(n + 1.0);
  }
}

long double cofactor(
		     int n, 
		     int i
		     ) {
  if ((n < 0) || (i < 0) || (n < i)) {
    return (long double)0.0;
  }
  else if (n == i) {
    return (long double)1.0;
  }
  else {
    return exp(gammaln(n + 1.0) - gammaln(i + 1.0) - gammaln(n-i + 1.0));
  }
}

long double cofactorln(
		       int n, 
		       int i
		       ) {
  if ((n < 0) || (i < 0) || (n < i)) {
    return (long double)-1.0;
  }
  else if (n == i) {
    return (long double)0.0;
  }
  else {
    return gammaln(n + 1.0) - gammaln(i + 1.0) - gammaln(n-i + 1.0);
  }
}

// prevent underflows by returning exp(LDBL_MIN_EXP) if exponentiation will produce an underflow
long double safe_exp(long double ln) {
    if (ln < LDBL_MIN_EXP) {  // -16381
        return LDBL_MIN;      // 3.3621e-4932
    } else {
        return exp(ln);
    }
}

// 'safe' log summation for probabilities
long double logsumexp_probs(const vector<long double>& lnv) {
    vector<long double>::const_iterator i = lnv.begin();
    long double maxN = *i;
    ++i;
    for (; i != lnv.end(); ++i) {
        if (*i > maxN)
            maxN = *i;
    }
    long double sum = 0;
    for (vector<long double>::const_iterator i = lnv.begin(); i != lnv.end(); ++i) {
        sum += safe_exp(*i - maxN);
    }
    return maxN + log(sum);
}

// unsafe, kept for potential future use
long double logsumexp(const vector<long double>& lnv) {
    long double maxAbs, minN, maxN, c;
    vector<long double>::const_iterator i = lnv.begin();
    long double n = *i;
    maxAbs = n; maxN = n; minN = n;
    ++i;
    for (; i != lnv.end(); ++i) {
        n = *i;
        if (n > maxN)
            maxN = n;
        if (fabs(n) > maxAbs)
            maxAbs = fabs(n);
        if (n < minN)
            minN = n;
    }
    if (maxAbs > maxN) {
        c = minN;
    } else {
        c = maxN;
    }
    long double sum = 0;
    for (vector<long double>::const_iterator i = lnv.begin(); i != lnv.end(); ++i) {
        sum += exp(*i - c);
    }
    return c + log(sum);
}

long double betaln(const vector<long double>& alphas) {
    vector<long double> gammalnAlphas;
    gammalnAlphas.resize(alphas.size());
    transform(alphas.begin(), alphas.end(), gammalnAlphas.begin(), gammaln);
    return sum(gammalnAlphas) - gammaln(sum(alphas));
}

long double beta(const vector<long double>& alphas) {
    return exp(betaln(alphas));
}

long double hoeffding(double successes, double trials, double prob) {
    return 0.5 * exp(-2 * pow(trials * prob - successes, 2) / trials);
}

// the sum of the harmonic series 1, n
long double harmonicSum(int n) {
    long double r = 0;
    long double i = 1;
    while (i <= n) {
        r += 1 / i;
        ++i;
    }
    return r;
}

bool isTransition(string& ref, string& alt) {
    if (((ref == "A" && alt == "G") || (ref == "G" && alt == "A")) ||
        ((ref == "C" && alt == "T") || (ref == "T" && alt == "C"))) {
        return true;
    } else {
        return false;
    }
}

// Levenshtein Distance Algorithm: C++ Implementation
// by Anders Sewerin Johansen
// http://www.merriampark.com/ldcpp.htm

int levenshteinDistance(const std::string source, const std::string target) {

  // Step 1

  const int n = source.length();
  const int m = target.length();
  if (n == 0) {
    return m;
  }
  if (m == 0) {
    return n;
  }

  // Good form to declare a TYPEDEF

  typedef std::vector< std::vector<int> > Tmatrix; 

  Tmatrix matrix(n+1);

  // Size the vectors in the 2.nd dimension. Unfortunately C++ doesn't
  // allow for allocation on declaration of 2.nd dimension of vec of vec

  for (int i = 0; i <= n; i++) {
    matrix[i].resize(m+1);
  }

  // Step 2

  for (int i = 0; i <= n; i++) {
    matrix[i][0]=i;
  }

  for (int j = 0; j <= m; j++) {
    matrix[0][j]=j;
  }

  // Step 3

  for (int i = 1; i <= n; i++) {

    const char s_i = source[i-1];

    // Step 4

    for (int j = 1; j <= m; j++) {

      const char t_j = target[j-1];

      // Step 5

      int cost;
      if (s_i == t_j) {
        cost = 0;
      }
      else {
        cost = 1;
      }

      // Step 6

      const int above = matrix[i-1][j];
      const int left = matrix[i][j-1];
      const int diag = matrix[i-1][j-1];
      int cell = min( above + 1, min(left + 1, diag + cost));

      // Step 6A: Cover transposition, in addition to deletion,
      // insertion and substitution. This step is taken from:
      // Berghel, Hal ; Roach, David : "An Extension of Ukkonen's 
      // Enhanced Dynamic Programming ASM Algorithm"
      // (http://www.acm.org/~hlb/publications/asm/asm.html)

      if (i>2 && j>2) {
        int trans=matrix[i-2][j-2]+1;
        if (source[i-2]!=t_j) trans++;
        if (s_i!=target[j-2]) trans++;
        if (cell>trans) cell=trans;
      }

      matrix[i][j]=cell;
    }
  }

  // Step 7

  return matrix[n][m];
}
