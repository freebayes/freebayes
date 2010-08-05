// utility functions
//
#include "Utility.h"
#include "Sum.h"
#include "Product.h"

using namespace std;

short qualityChar2ShortInt(char c) {
    return static_cast<short>(c) - 33;
}

long double phred2ln(int qual) {
    return M_LN10 * qual * -.1;
}

int ln2phred(long double prob) {
    return -10 * M_LOG10E * prob;
}

long double phred2float(int qual) {
    return pow(10,qual * -.1);
}

int float2phred(long double prob) {
    return std::min(-10 * (long double) log10(prob), (long double) 99);
}

long double powln(long double m, int n) {
    long double r = 0;
    for (int i = 0; i < n; ++i) {
        r += m;
    }
    return r;
}

// TODO do the following in phred (log) space for perf boost

// the probability that we have a completely true vector of qualities
short jointQuality(const std::vector<short>& quals) {
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
    return float2phred(1 - prod);
}

short jointQuality(const std::string& qualstr) {

    std::vector<short> quals;
    for (int i=0; i<qualstr.size(); i++)
        quals.push_back(qualityChar2ShortInt(qualstr.at(i)));

    return jointQuality(quals);

}

short averageQuality(const std::string& qualstr) {

    long double q = 0;
    for (int i=0; i<qualstr.size(); i++)
        q += phred2float(qualityChar2ShortInt(qualstr.at(i)));
    return float2phred(q / qualstr.size());

}

bool stringInVector(string item, vector<string> items) {
    for (vector<string>::iterator i = items.begin(); i != items.end(); ++i) {
        if (item == *i) {
            return true;
        }
    }
    return false;
}

// from Function-Math.cpp

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
    transform(alphas.begin(), alphas.end(), gammalnAlphas.begin(), gammaln);
    return sum(gammalnAlphas) - gammaln(sum(alphas));
}

long double beta(const vector<long double>& alphas) {
    return exp(betaln(alphas));
}
