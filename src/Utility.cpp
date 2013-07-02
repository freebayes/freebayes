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

 long double ln2log10(long double prob) {
     return M_LOG10E * prob;
 }

 long double log102ln(long double prob) {
     return M_LN10 * prob;
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
         return PHRED_MAX;  // guards against "-0"
     long double p = -10 * (long double) log10(prob);
     if (p < 0 || p > PHRED_MAX) // int overflow guard
         return PHRED_MAX;
     else
         return p;
 }

 long double big2phred(const BigFloat& prob) {
     return -10 * (long double) (ttmath::Log(prob, (BigFloat)10)).ToDouble();
 }

 long double nan2zero(long double x) {
     if (x != x) {
     return 0;
     } else {
     return x;
     }
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

 long double minQuality(const std::string& qualstr) {
     long double qual = 0;
     for (string::const_iterator q = qualstr.begin(); q != qualstr.end(); ++q) {
         long double nq = qualityChar2LongDouble(*q);
         if (qual == 0) {
             qual = nq;
         } else if (nq < qual) {
             qual = nq;
         }
     }
     return qual;
 }

 // crudely averages quality scores in phred space
 long double averageQuality(const std::string& qualstr) {

     long double qual = 0; //(long double) *max_element(quals.begin(), quals.end());
     for (string::const_iterator q = qualstr.begin(); q != qualstr.end(); ++q)
             qual += qualityChar2LongDouble(*q);
     return qual / qualstr.size();

 }

 long double averageQuality(const vector<short>& qualities) {

     long double qual = 0;
     for (vector<short>::const_iterator q = qualities.begin(); q != qualities.end(); ++q) {
         qual += *q;
     }
     return qual / qualities.size();
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

long double __binomialProbln(int k, int n, long double p) {
    return factorialln(n) - (factorialln(k) + factorialln(n - k)) + powln(log(p), k) + powln(log(1 - p), n - k);
}

long double binomialCoefficientLn(int k, int n) {
    return factorialln(n) - (factorialln(k) + factorialln(n - k));
}

BinomialCache binomialCache;

long double binomialProbln(int k, int n, long double p) {
    return binomialCache.binomialProbln(k, n, p);
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

FactorialCache factorialCache;

long double factorialln(int n) {
    return factorialCache.factorialln(n);
}

long double __factorialln(
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

BigFloat big_exp(long double ln) {
    BigFloat x, result;
    x.FromDouble(ln);
    result = ttmath::Exp(x);
    return result;
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
    BigFloat sum = 0;
    for (vector<long double>::const_iterator i = lnv.begin(); i != lnv.end(); ++i) {
        sum += big_exp(*i - maxN);
    }
    BigFloat maxNb; maxNb.FromDouble(maxN);
    BigFloat bigResult = maxNb + ttmath::Ln(sum);
    long double result;
    return bigResult.ToDouble();
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

long double hoeffdingln(double successes, double trials, double prob) {
    return log(0.5) + (-2 * pow(trials * prob - successes, 2) / trials);
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

// current date string in YYYYMMDD format
string dateStr(void) {

    time_t rawtime;
    struct tm* timeinfo;
    char buffer[80];

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer, 80, "%Y%m%d", timeinfo);

    return string(buffer);

}

long double string2float(const string& s) {
    long double r;
    convert(s, r);
    return r;
}

long double log10string2ln(const string& s) {
    long double r;
    convert(s, r);
    return log102ln(r);
}

long double safedivide(long double a, long double b) {
    if (b == 0) {
        if (a == 0) {
            return 1;
        } else {
            return 0;
        }
    } else {
        return a / b;
    }
}

string mergeCigar(const string& c1, const string& c2) {
    vector<pair<int, string> > cigar1 = splitCigar(c1);
    vector<pair<int, string> > cigar2 = splitCigar(c2);
    // check if the middle elements are the same
    if (cigar1.back().second == cigar2.front().second) {
        cigar1.back().first += cigar2.front().first;
        cigar2.erase(cigar2.begin());
    }
    for (vector<pair<int, string> >::iterator c = cigar2.begin(); c != cigar2.end(); ++c) {
        cigar1.push_back(*c);
    }
    return joinCigar(cigar1);
}

vector<pair<int, string> > splitCigar(const string& cigarStr) {
    vector<pair<int, string> > cigar;
    string number;
    string type;
    // strings go [Number][Type] ...
    for (string::const_iterator s = cigarStr.begin(); s != cigarStr.end(); ++s) {
        char c = *s;
        if (isdigit(c)) {
            if (type.empty()) {
                number += c;
            } else {
                // signal for next token, push back the last pair, clean up
                cigar.push_back(make_pair(atoi(number.c_str()), type));
                number.clear();
                type.clear();
                number += c;
            }
        } else {
            type += c;
        }
    }
    if (!number.empty() && !type.empty()) {
        cigar.push_back(make_pair(atoi(number.c_str()), type));
    }
    return cigar;
}

list<pair<int, string> > splitCigarList(const string& cigarStr) {
    list<pair<int, string> > cigar;
    string number;
    string type;
    // strings go [Number][Type] ...
    for (string::const_iterator s = cigarStr.begin(); s != cigarStr.end(); ++s) {
        char c = *s;
        if (isdigit(c)) {
            if (type.empty()) {
                number += c;
            } else {
                // signal for next token, push back the last pair, clean up
                cigar.push_back(make_pair(atoi(number.c_str()), type));
                number.clear();
                type.clear();
                number += c;
            }
        } else {
            type += c;
        }
    }
    if (!number.empty() && !type.empty()) {
        cigar.push_back(make_pair(atoi(number.c_str()), type));
    }
    return cigar;
}

string joinCigar(const vector<pair<int, string> >& cigar) {
    string cigarStr;
    for (vector<pair<int, string> >::const_iterator c = cigar.begin(); c != cigar.end(); ++c) {
        if (c->first) {
            cigarStr += convert(c->first) + c->second;
        }
    }
    return cigarStr;
}

string joinCigarList(const list<pair<int, string> >& cigar) {
    string cigarStr;
    for (list<pair<int, string> >::const_iterator c = cigar.begin(); c != cigar.end(); ++c) {
        cigarStr += convert(c->first) + c->second;
    }
    return cigarStr;
}

bool isEmptyCigarElement(const pair<int, string>& elem) {
    return elem.first == 0;
}


// string * overload
// from http://stackoverflow.com/a/5145880
std::string operator*(std::string const &s, size_t n)
{
    std::string r;  // empty string
    r.reserve(n * s.size());
    for (size_t i=0; i<n; i++)
        r += s;
    return r;
}

// normalize vector sum to 1
void normalizeSumToOne(vector<long double>& v) {
    long double sum = 0;
    for (vector<long double>::iterator i = v.begin(); i != v.end(); ++i) {
	sum += *i;
    }
    for (vector<long double>::iterator i = v.begin(); i != v.end(); ++i) {
	*i /= sum;
    }
}

// splits the file on '\n', adds the resulting values to v
void addLinesFromFile(vector<string>& v, const string& f) {
    ifstream ifs;
    ifs.open(f.c_str(), ifstream::in);
    if (!ifs.is_open()) {
        cerr << "could not open " << f << endl;
        exit(1);
    }
    string line;
    while (std::getline(ifs, line)) {
        v.push_back(line);
    }
}
