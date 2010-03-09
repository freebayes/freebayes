//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Function-Math
// Math methods
// Copyright 2006 Gabor T. Marth, Boston College
// All rights reserved
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef FUNCTION_MATH_CPP
#define FUNCTION_MATH_CPP

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <cmath>

using namespace std;

#include "Type-Hash.h"
#include "Function-Math.h"
#include "TryCatch.h"

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// static
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

map<string, string > InitIUPAC(void) {

  map<string, string > IUPAC;

  IUPAC["A"] = "A"; 
  IUPAC["C"] = "C"; 
  IUPAC["G"] = "G"; 
  IUPAC["T"] = "T"; 
  IUPAC["AA"] = "A"; 
  IUPAC["AC"] = "M"; 
  IUPAC["AG"] = "R"; 
  IUPAC["AT"] = "W"; 
  IUPAC["CA"] = "M"; 
  IUPAC["CC"] = "C"; 
  IUPAC["CG"] = "S"; 
  IUPAC["CT"] = "Y"; 
  IUPAC["GA"] = "R"; 
  IUPAC["GC"] = "S"; 
  IUPAC["GG"] = "G"; 
  IUPAC["GT"] = "K"; 
  IUPAC["TA"] = "W"; 
  IUPAC["TC"] = "Y"; 
  IUPAC["TG"] = "K"; 
  IUPAC["TT"] = "T"; 
  
  return IUPAC;
}

// IUPAC ambiguity coding
static map<string, string > IUPAC = InitIUPAC();

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// templates
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// sorts keys of a hash in order of associated value
//------------------------------------------------------------------------------
template< typename keyType, typename valueType >
vector<keyType> sortKeysByValue(map<keyType, valueType, std::less<keyType> > hash, bool descend) {

  // instantiate inverse hash as a multimap
  multimap<valueType, keyType, std::less<valueType> > inverseHash;

  // load elements of hash into inverseHash
  for (typename map<keyType, valueType, std::less<keyType> >::const_iterator iter = hash.begin(); iter != hash.end(); iter++) {
    keyType key = iter->first;
    valueType value = iter->second;
    inverseHash.insert(typename multimap<valueType, keyType, std::less<valueType> >::value_type(value, key));
  }

  // compose vector of original keys sorted by original values
  vector<keyType> sortedKeys;
  for(typename multimap<valueType, keyType, std::less<valueType> >::const_iterator iter = inverseHash.begin(); 
      iter != inverseHash.end(); iter++) {
    keyType key = iter->second;
    sortedKeys.push_back(key);
  }

  // reverse if descending order was required
  if (descend) {
    reverse(sortedKeys.begin(), sortedKeys.end());
  }

  // return
  return sortedKeys;
}

//------------------------------------------------------------------------------
// prints a vector
//------------------------------------------------------------------------------
template < class T >
void printVector(const vector< T > &vectorRef, const bool separator, const bool bracket, ostream &output2) {
  if ( vectorRef.empty() ) {
    if (bracket) {
      output2 << "()";
    }
  }
  else {
    if (bracket) {
      output2 << "(";
    }
    if (separator) {
      std::ostream_iterator< T > output( output2, " " );
      std::copy( vectorRef.begin(), vectorRef.end(), output );
    }
    else {
      std::ostream_iterator< T > output( output2 );
      std::copy( vectorRef.begin(), vectorRef.end(), output );
    }
    if (bracket) {
      output2 << ")";
    }
  }
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// functions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// gammaln
//------------------------------------------------------------------------------
long double gammaln(
		     long double x
		     ) {

  //----------------------------------------------------------------------------
  // constants
  //----------------------------------------------------------------------------
  // XXX TODO these cofactors are drawn from population genetics, but their
  // genesis is not noted anywhere in this source.
  vector<long double> cofactors;
  cofactors.push_back(76.18009173);
  cofactors.push_back(-86.50532033);
  cofactors.push_back(24.01409822);
  cofactors.push_back(-1.231739516);
  cofactors.push_back(0.120858003E-2);
  cofactors.push_back(-0.536382E-5);    

  //----------------------------------------------------------------------------
  // compute
  //----------------------------------------------------------------------------
  long double x1 = x - 1.0;
  long double tmp = x1 + 5.5;
  tmp -= (x1 + 0.5) * log(tmp);
  long double ser = 1.0;
  for (int j=0; j<=5; j++) {
    x1 += 1.0;
    ser += cofactors[j]/x1;
  }
  long double y =  (-1.0 * tmp + log(2.50662827465 * ser));

  // return
  return y;
}

//------------------------------------------------------------------------------
// factorial
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
// factorialln
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
// cofactor
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
// cofactorln
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
// logGenotypeLikelihoods (IUPAC genotype code output)
//------------------------------------------------------------------------------
map<string, long double > logGenotypeLikelihoods(
							       vector<Basecall> basecalls,
							       bool diploid,
							       long double dependenceFactor,
							       bool debug
							       ) {

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Notations:
  //     Gi: genotype of a given individual (e.g. AA, AG, GG)
  //     Bi: set of base calls for the individual (e.g. A,A,A,G,A,A,G)
  //     Ai: set of true bases estimated by the base calls (e.g. A,A,G,A,G,G)
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // variables
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // output
  //----------------------------------------------------------------------------
  map<string, long double > gl;

  //----------------------------------------------------------------------------
  // constants
  //----------------------------------------------------------------------------
  long double LOGFACTOR = log((long double)10.0) / ((long double)-10.0); 
  long double ln3 = log((long double)3.0);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // calculate genotype likelihoods
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // Diploid individual:
  //----------------------------------------------------------------------------
  if (diploid) {

    // variants
    vector<string> genotypes;
    string AC = "AC"; genotypes.push_back(AC);
    string AG = "AG"; genotypes.push_back(AG);
    string AT = "AT"; genotypes.push_back(AT);
    string CG = "CG"; genotypes.push_back(CG);
    string CT = "CT"; genotypes.push_back(CT);
    string GT = "GT"; genotypes.push_back(GT);

    //--------------------------------------------------------------------------
    // cycle through each variant
    //--------------------------------------------------------------------------
    for (vector<string>::const_iterator genIter = genotypes.begin();
         genIter != genotypes.end(); genIter++) {
      string genotype = *genIter;

      //------------------------------------------------------------------------
      // get variant alleles
      //------------------------------------------------------------------------
      string a1, a2;
      TRY { a1 = genotype.substr(0, 1); } CATCH;
      TRY { a2 = genotype.substr(1, 1); } CATCH;
    
      //------------------------------------------------------------------------
      // define and initialize individual data likelihoods
      //------------------------------------------------------------------------
      long double lnProbBiGivenGi11 = 0;
      long double lnProbBiGivenGi22 = 0;
      long double lnProbBiGivenGi12 = 0;

      //------------------------------------------------------------------------
      // allele counts
      //------------------------------------------------------------------------
      unsigned int c1 = 0; // count of a1 alleles
      unsigned int c2 = 0; // count of a2 alleles
      unsigned int c3 = 0; // count of other (not a1 or a2) alleles

     // aggregate allele specific quality values
      long double sumQa1 = 0;
      long double sumQa2 = 0;
      long double sumQa3 = 0;

      //------------------------------------------------------------------------
      // cycle through each base call, and update sum allele quality values
      //   and allele counts
      //------------------------------------------------------------------------
      for(vector<Basecall>::const_iterator iter = basecalls.begin();
          iter != basecalls.end(); iter++) {

        //----------------------------------------------------------------------
        // retrieve base call and call data
        //----------------------------------------------------------------------
        Basecall bc = *iter;

        string b = bc.base;
        short q = bc.qual;
        string s = bc.strand;

        //----------------------------------------------------------------------
        // calculate update terms for lnProbBiGivenGi terms
        // tally number of a1, a2, and other base calls
        //----------------------------------------------------------------------
        if (b == a1) {

          // update terms
          sumQa1 += (long double)q;

         // increment a1 allele count
          c1++;
        } 
        else if (b == a2) {

          // update terms
          sumQa2 += (long double)q;

          // increment a2 allele count
          c2++;
        }
        else {

          // update terms
          sumQa3 += (long double)q;

          // increment other allele count
          c3++;
        }
      }

      //------------------------------------------------------------------------
      // downgrade aggregate quality values for a1 and a2 using dependence factor
      //------------------------------------------------------------------------
      if (c1 > 1) {
        sumQa1 *= ((long double)1.0 + ((long double)c1 - (long double)1.0) * dependenceFactor) / (long double)c1;
      }
      if (c2 > 1) {
        sumQa2 *= ((long double)1.0 + ((long double)c2 - (long double)1.0) * dependenceFactor) / (long double)c2;
      }

      //------------------------------------------------------------------------
      // assign ln data probabilities
      //------------------------------------------------------------------------
      lnProbBiGivenGi11 = LOGFACTOR * (sumQa2 + sumQa3) - (c2+c3) * ln3;
      lnProbBiGivenGi22 = LOGFACTOR * (sumQa1 + sumQa3) - (c1+c3) * ln3;
      lnProbBiGivenGi12 = LOGFACTOR * sumQa3 - c3 * ln3;

      //------------------------------------------------------------------------
      // calculate (c1+c2+c3 over c1)
      //------------------------------------------------------------------------
      //      long double lnBinomialfactor = cofactorln(c1+c2+c3, c1) + (c1+c2+c3) * log(0.5);
      long double lnBinomialfactor = (c1+c2+c3) * log(0.5);
      if (debug) {
        cout << "lnBinomialfactor=" << lnBinomialfactor << " binomialFactor=" << exp(lnBinomialfactor) << endl;
      }

      //------------------------------------------------------------------------
      // update lnProbBiGivenGi12 with binomial factor
      //------------------------------------------------------------------------
      lnProbBiGivenGi12 += lnBinomialfactor;

      //------------------------------------------------------------------------
      // Register genotype log likelihoods.
      //------------------------------------------------------------------------
      gl[IUPAC[a1+a1]] = lnProbBiGivenGi11;
      gl[IUPAC[a2+a2]] = lnProbBiGivenGi22;
      gl[IUPAC[a1+a2]] = lnProbBiGivenGi12;
      gl[IUPAC[a2+a1]] = lnProbBiGivenGi12;

      if (debug) {
        cout << " a1=" << a1 << " a2=" << a2;
        cout << " lnP(Bi|Gi=11)=" << lnProbBiGivenGi11;
        cout << " lnP(Bi|Gi=22)=" << lnProbBiGivenGi22;
        cout << " lnP(Bi|Gi=12)=" << lnProbBiGivenGi12;
        cout << endl;
      }
    }
  }

  //----------------------------------------------------------------------------
  // Haploid individual
  //----------------------------------------------------------------------------
  else {

    // genotypes
    vector<string> genotypes;
    string A = "A"; genotypes.push_back(A);
    string C = "C"; genotypes.push_back(C);
    string G = "G"; genotypes.push_back(G);
    string T = "T"; genotypes.push_back(T);

    //--------------------------------------------------------------------------
    // cycle through each variant
    //--------------------------------------------------------------------------
    for (vector<string>::const_iterator genIter = genotypes.begin();
         genIter != genotypes.end(); genIter++) {

      //------------------------------------------------------------------------
      // get variant allele
      //------------------------------------------------------------------------
      string a = *genIter;

      // allele counts
      unsigned int c = 0;
      unsigned int c3 = 0;

      // aggregate allele specific quality value
      long double sumQa = 0;
      long double sumQa3 = 0;

      //------------------------------------------------------------------------
      // define and initialize individual data likelihoods
      //------------------------------------------------------------------------
      long double lnProbBiGivenGi = 0;

      //------------------------------------------------------------------------
      // cycle through each base call, and update data likelihoods
      //------------------------------------------------------------------------
      for(vector<Basecall>::const_iterator iter = basecalls.begin();
          iter != basecalls.end(); iter++) {

        // retrieve base call
        Basecall bc = *iter;

        // retrieve base call data
        string b = bc.base;
        short q = bc.qual;
        string s = bc.strand;

        if (b == a) {

          // update terms
          sumQa += (long double)q;

          // increment a allele count
          c++;
        } 
        else {

          // update terms
          sumQa3 += (long double)q;

          // increment other allele count
          c3++;
        }
      }
     
      //------------------------------------------------------------------------
      // assign ln data probabilities
      //------------------------------------------------------------------------
      lnProbBiGivenGi = LOGFACTOR * sumQa3 - c3 * ln3;
 
      //------------------------------------------------------------------------
      // Register genotype log likelihoods.
      //------------------------------------------------------------------------
      gl[IUPAC[a]] = lnProbBiGivenGi;
    }
  }

  //----------------------------------------------------------------------------
  // return
  //----------------------------------------------------------------------------
  return gl;
}

//------------------------------------------------------------------------------
// varBayesProb
//------------------------------------------------------------------------------
Variation varBayesProb(
		       map<string, vector<Basecall> > indBasecalls,
		       map<string, SampleInfo > indSampleInfo,
		       bool diploid,
		       long double THETA,
		       long double MU,
		       bool banded,
		       int WB,
		       int TB,
		       int TR,
		       string allele1,
		       string allele2,
		       bool cgp,
		       bool debug2
		       ) {

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Notations:
  //     Gi: genotype of a given individual (e.g. AA, AG, GG)
  //     Bi: set of base calls for the individual (e.g. A,A,A,G,A,A,G)
  //     Ai: set of true bases estimated by the base calls (e.g. A,A,G,A,G,G)
  //
  //     Go: genotype vector of all individuals
  //     Bo: base calls for all individuals
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // constants
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  long double LOGFACTOR = log(10) / (-10.0); 
  string Gi1; Gi1 += "1";
  string Gi2; Gi2 += "2";
  string Gi11; Gi11 += "11";
  string Gi22; Gi22 += "22";
  string Gi12; Gi12 += "12";

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // 1. calculate ln P(Bi|Gi) for each individual where Bi is the set of 
  //      base calls, and Gi is the genotype of that individual
  // 2. make special genotype combinations (i.e. two homomorphic combinations
  //      and the nominally "best" genotype combo (without considering priors)
  // 3. determine allele frequencies for nominally "best" genotype combo
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  
  //----------------------------------------------------------------------------
  // Initializations.
  //----------------------------------------------------------------------------

  // Define individual list.
  vector<string> individuals;

  // Define ln Prob(R|G) for each individual
  //
  // XXX turn this into an unordered map from sample name to a struct holding Gi -> prob mapping
  map<string, map<string, long double > > LnProbBiGivenGi;

  //----------------------------------------------------------------------------
  // Process through each individual with base calls.
  //----------------------------------------------------------------------------
  if (debug2) {
    clog << "Individual genotype likelihoods P(Bi|Ai), P(Ai|Gi), and P(Bi|Gi):"<< endl;    
  }

  for (map<string, vector<Basecall> >::const_iterator indIter = indBasecalls.begin();
       indIter != indBasecalls.end(); indIter++) {

    // Retreive individual.
    string ind = indIter->first;

    // Store individual in individual list.
    individuals.push_back(ind);

    // Filter irrelevant (not allele1 or allele 2) and low quality base calls, and
    //   count remaining alleles.
    vector<Basecall> basecallsScreened;
    int count1 = 0; int count2 = 0;
    for(vector<Basecall>::const_iterator iter = indBasecalls[ind].begin();
	iter != indBasecalls[ind].end(); iter++) {
      Basecall bc = *iter;
      
      // Filter irrelevant base calls.
      if ((bc.base != allele1) && (bc.base != allele2)) {
	continue;
      }

      // Register base call.
      basecallsScreened.push_back(bc);
      
      // Update allele counts.
      if (bc.base == allele1) {count1++;}
      else if (bc.base == allele2) {count2++;}
    }
    
    // total number of sequences
    int count = count1 + count2;
    
    //--------------------------------------------------------------------------
    // Calculate P(Bi|Gi) for this individual.
    //--------------------------------------------------------------------------
    if (debug2) {
      clog << "  ind=" << ind;    
    }

    // Diploid individual:
    if (diploid) {

      //------------------------------------------------------------------------
      // Calculate ln P(Bi|Ai) probabilities. 
      //------------------------------------------------------------------------
      long double lnProbBiGivenAi11 = 0;
      long double lnProbBiGivenAi22 = 0;
      long double lnProbBiGivenAi12 = 0;
      for(vector<Basecall>::const_iterator iter = basecallsScreened.begin();
	  iter != basecallsScreened.end(); iter++) {
	Basecall bc = *iter;
	short q = bc.qual;

	// Update those probabilities where this base represents an error.
	long double lnProbCorrect = log(1 - pow(10, q/(-10.0)));
	if (bc.base == allele1) {
	  lnProbBiGivenAi22 += LOGFACTOR * q;
	  lnProbBiGivenAi11 += lnProbCorrect;
	  lnProbBiGivenAi12 += lnProbCorrect;
	}
	else if (bc.base == allele2) {
	  lnProbBiGivenAi11 += LOGFACTOR * q;
	  lnProbBiGivenAi22 += lnProbCorrect;
	  lnProbBiGivenAi12 += lnProbCorrect;
	}
      }
      
      //------------------------------------------------------------------------
      // Calculate ln P(Ai|Gi) probabilities.
      //------------------------------------------------------------------------
      long double lnProbAi11GivenGi11 = 0;
      long double lnProbAi22GivenGi22 = 0;
      long double lnProbAi11GivenGi12 = count * log(0.5);
      long double lnProbAi22GivenGi12 = count * log(0.5);

      // maybe this should also be just count * log(0.5); 
      // this is because we are only considering a single base combination, not
      // all base combinations with (count, count1) allele frequency
      long double lnProbAi12GivenGi12 = cofactorln(count, count1) + (count * log(0.5));
      //      long double lnProbAi12GivenGi12 = count * log(0.5);

      //------------------------------------------------------------------------
      // Calculate ln P(Bi|Gi) = SUM(all Ai) ln P(Bi|Ai) + ln P(Ai|Gi).
      // Shortcut: we are only considering one variant allele combination A i.e.
      //   the one suggested by the base calls B... if the individual read 
      //   coverage is not too high (<20), and the base quality values are high
      //   (>=30) than the error of the approx should be < 1-2%.
      // We could possibly do a "banded" improvement by considering allele combos
      //   one edit distance away... this would probably raise the probability.
      //   of the het genotype: P(Bi|Gi12) somewhat.
      // We could also implement the "recursive" algorithm for this case, keeping
      //   only the dominant terms in each round.
      //------------------------------------------------------------------------
      long double lnProbBiGivenGi11 = lnProbBiGivenAi11 + lnProbAi11GivenGi11;
      long double lnProbBiGivenGi22 = lnProbBiGivenAi22 + lnProbAi22GivenGi22;
      long double lnProbBiGivenGi12;

      // If the most likely base combination is Gi11 or Gi22 then ProbS12GivenG12 is
      //   already accounted for.
      if ((count1 == 0) || (count2 == 0)) {
	lnProbBiGivenGi12 = log(
			      exp(lnProbBiGivenAi11 + lnProbAi11GivenGi12)
			      + exp(lnProbBiGivenAi22 + lnProbAi22GivenGi12)
			      );
      }
      else {
	lnProbBiGivenGi12 = log(
			      exp(lnProbBiGivenAi11 + lnProbAi11GivenGi12)
			      + exp(lnProbBiGivenAi22 + lnProbAi22GivenGi12)
			      + exp(lnProbBiGivenAi12 + lnProbAi12GivenGi12)
			      );
      }
      if (debug2) {
	clog << " lnP(Bi|Ai=11)=" << lnProbBiGivenAi11;
	clog << " lnP(Bi|Ai=22)=" << lnProbBiGivenAi22;
	clog << " lnP(Bi|Ai12)=" << lnProbBiGivenAi12;
	clog << " lnP(Ai=11|Gi=11)=" << lnProbAi11GivenGi11;
	clog << " lnP(Ai=22|Gi=22)=" << lnProbAi22GivenGi22;
	clog << " lnP(Ai=11|Gi=12)=" << lnProbAi11GivenGi12;
	clog << " lnP(Ai=22|Gi=12)=" << lnProbAi22GivenGi12;
	clog << " lnP(Ai=12|Gi=12)=" << lnProbAi12GivenGi12;
	clog << " lnP(Bi|Gi=11)=" << lnProbBiGivenGi11;
	clog << " lnP(Bi|Gi=22)=" << lnProbBiGivenGi22;
	clog << " lnP(Bi|Gi=12)=" << lnProbBiGivenGi12;


	clog << " lnProbBiGivenAi11 + lnProbAi11GivenGi12=" << lnProbBiGivenAi11 + lnProbAi11GivenGi12;
	clog << " lnProbBiGivenAi22 + lnProbAi22GivenGi12=" << lnProbBiGivenAi22 + lnProbAi22GivenGi12;
	clog << " lnProbBiGivenAi12 + lnProbAi12GivenGi12=" << lnProbBiGivenAi12 + lnProbAi12GivenGi12;


	clog << " P(Bi|Ai=11)=" << exp(lnProbBiGivenAi11);
	clog << " P(Bi|Ai=22)=" << exp(lnProbBiGivenAi22);
	clog << " P(Bi|Ai=12)=" << exp(lnProbBiGivenAi12);
	clog << " P(Ai=11|Gi=11)=" << exp(lnProbAi11GivenGi11);
	clog << " P(Ai=22|Gi=22)=" << exp(lnProbAi22GivenGi22);
	clog << " P(Ai=11|Gi=12)=" << exp(lnProbAi11GivenGi12);
	clog << " P(Ai=22|Gi=12)=" << exp(lnProbAi22GivenGi12);
	clog << " P(Ai=12|Gi=12)=" << exp(lnProbAi12GivenGi12);
	clog << " P(Bi|Gi=11)=" << exp(lnProbBiGivenGi11);
	clog << " P(Bi|Gi=22)=" << exp(lnProbBiGivenGi22);
	clog << " P(Bi|Gi=12)=" << exp(lnProbBiGivenGi12)<< endl;
      }

      // Register genotype log likelihoods.
      LnProbBiGivenGi[ind][Gi11] = lnProbBiGivenGi11;
      LnProbBiGivenGi[ind][Gi22] = lnProbBiGivenGi22;
      LnProbBiGivenGi[ind][Gi12] = lnProbBiGivenGi12;
    }

    // Haploid individual:
    else {

      //------------------------------------------------------------------------
      // calculate ln P(Bi|Ai) probabilities
      //------------------------------------------------------------------------
      long double lnProbBiGivenAi11 = 0;
      long double lnProbBiGivenAi22 = 0;
      for(vector<Basecall>::const_iterator iter = basecallsScreened.begin();
	  iter != basecallsScreened.end(); iter++) {
	Basecall bc = *iter;
      
	// update those probabilities where this base represents and error
	if (bc.base == allele1) {
	  lnProbBiGivenAi22 += LOGFACTOR * bc.qual;
	}
	else if (bc.base == allele2) {
	  lnProbBiGivenAi11 += LOGFACTOR * bc.qual;
	}
      }
      
      //------------------------------------------------------------------------
      // calculate ln P(Ai|Gi) probabilities
      //------------------------------------------------------------------------
      long double lnProbAi11GivenGi1 = 0;
      long double lnProbAi22GivenGi2 = 0;

      //------------------------------------------------------------------------
      // calculate ln P(Bi|Gi) = ln P(Bi|Ai) + ln P(Ai|Gi) probabilities
      //------------------------------------------------------------------------
      long double lnProbBiGivenGi1 = lnProbBiGivenAi11 + lnProbAi11GivenGi1;
      long double lnProbBiGivenGi2 = lnProbBiGivenAi22 + lnProbAi22GivenGi2;

      // register
      LnProbBiGivenGi[ind][Gi1] = lnProbBiGivenGi1;
      LnProbBiGivenGi[ind][Gi2] = lnProbBiGivenGi2;
    }
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // do necessary management of unrelateds and children to deal with trio data
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // assign pedigree status: child, mother, father, unrelated
  //----------------------------------------------------------------------------
  map<string, string > pedigreeStatus, fatherId, motherId;

  // assign initial pedigree status based on SampleInfo
  for (map<string, SampleInfo >::const_iterator indIter = indSampleInfo.begin();
       indIter != indSampleInfo.end(); indIter++) {

    // retreive individual and samle info
    string ind = indIter->first;
    SampleInfo si = indIter->second;

    // initalize pedigree status as U (unrelated)
    pedigreeStatus[ind] = "U";
    fatherId[ind] = "";
    motherId[ind] = "";
    if (si.fatherId != "" && si.motherId != "") {
      pedigreeStatus[ind] = "C";
      fatherId[ind] = si.fatherId;
      motherId[ind] = si.motherId;
    }
  }

  // check if children have at least one parent also present
  // if not, reassign child as an unrelated
  // if there is only one parent present, make list of dummy parents to add to 
  //   individual list
  map<string, bool > isDummyParent;
  for (map<string, string >::const_iterator indIter = pedigreeStatus.begin();
       indIter != pedigreeStatus.end(); indIter++) {
    string ind = indIter->first;
    string ps = indIter->second;

    if (ps == "C") {
      // reads from neither mother nor father --> reclassify child as unrelated
      if (pedigreeStatus.count(fatherId[ind]) == 0 && pedigreeStatus.count(motherId[ind]) == 0) {

	// reclassify as unrelated
	pedigreeStatus[ind] = "U";
      }
      // no reads from father --> make dummy father
      else if (pedigreeStatus.count(fatherId[ind]) == 0) {
	
	// create entries for dummy parent and add to list
	string dummyParent = fatherId[ind];
	pedigreeStatus[dummyParent] = "U";
	fatherId[dummyParent] = "";
	motherId[dummyParent] = "";
	isDummyParent[dummyParent] = true;
      }
      // no reads from mother --> make dummy mother
      else if (pedigreeStatus.count(motherId[ind]) == 0) {
	
	// create entries for dummy parent and add to list
	string dummyParent = motherId[ind];
	pedigreeStatus[dummyParent] = "U";
	fatherId[dummyParent] = "";
	motherId[dummyParent] = "";
	isDummyParent[dummyParent] = true;
      }
    }
  }

  // add dummy parents to individual list, and add P(Bi|Gi) values
  for (map<string, bool >::const_iterator indIter = isDummyParent.begin();
       indIter != isDummyParent.end(); indIter++) {
    string ind = indIter->first;

    // add to list of individuals
    individuals.push_back(ind);

    // assign probabilities
    if (diploid) {
      LnProbBiGivenGi[ind][Gi11] = 0;
      LnProbBiGivenGi[ind][Gi22] = 0;
      LnProbBiGivenGi[ind][Gi12] = 0;
    }
    else {
      LnProbBiGivenGi[ind][Gi1] = 0;
      LnProbBiGivenGi[ind][Gi2] = 0;
    }
  }

  //----------------------------------------------------------------------------
  // sort unrelateds (including parents) and children into separate lists
  //----------------------------------------------------------------------------
  vector<string> unrelateds, children;
  for (vector<string>::const_iterator indIter = individuals.begin();
       indIter != individuals.end(); indIter++) {

    // get individual
    string ind = *indIter;
    
    if (pedigreeStatus[ind] == "U") {
      unrelateds.push_back(ind);
    }
    else {
      children.push_back(ind);
    }
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Build special genotype combinations: mono allele 1, mono allele 2.
  // Build highest data likelihood genotype combination if banded algorithm 
  // Register corresponding ln P(Bo|Go) data likelihoods
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  
  // genotype combination vector with no elements
  GenotypeCombo GoNull;

  // special genotype combinations
  GenotypeCombo GoAllGi1, GoAllGi2, GoAllGi11, GoAllGi22, GoMixGi;

  // P(Bo|Go) for special genotype combinations
  long double lnProbBoGivenGoAllGi1 = 0;
  long double lnProbBoGivenGoAllGi2 = 0;
  long double lnProbBoGivenGoAllGi11 = 0;
  long double lnProbBoGivenGoAllGi22 = 0;
  long double lnProbBoGivenGoMixGi = 0;
  
  // build special genotype combinations
  for (vector<string>::const_iterator indIter = individuals.begin();
       indIter != individuals.end(); indIter++) {
    
    // get individual
    string ind = *indIter;
    
    // diploid DNA
    if (diploid) {
      
      //------------------------------------------------------------------------
      // add Gi11 to all-Gi11 genotype combo and update probability
      //------------------------------------------------------------------------
      GoAllGi11[ind] = Gi11;
      lnProbBoGivenGoAllGi11 += LnProbBiGivenGi[ind][Gi11];
      
      //------------------------------------------------------------------------
      // add Gi22 to all-Gi22 genotype combo and update probability
      //------------------------------------------------------------------------
      GoAllGi22[ind] = Gi22;
      lnProbBoGivenGoAllGi22 += LnProbBiGivenGi[ind][Gi22];
      
      //------------------------------------------------------------------------
      // determine highest likelihood individual genotype (Gi) and add
      // to highest data likelihood genotype combo if banded alg.
      //------------------------------------------------------------------------
      if (banded) {
	// initialize as Gi11
	string GiBest = Gi11; long double lnProbBiGivenGiBest = LnProbBiGivenGi[ind][Gi11];
	
	// test Gi22
	if (LnProbBiGivenGi[ind][Gi22] > lnProbBiGivenGiBest) {
	  GiBest = Gi22; lnProbBiGivenGiBest = LnProbBiGivenGi[ind][Gi22];
	}
	
	// test Gi12
	if (LnProbBiGivenGi[ind][Gi12] > lnProbBiGivenGiBest) {
	  GiBest = Gi12; lnProbBiGivenGiBest = LnProbBiGivenGi[ind][Gi12];
	}
      
	// add best Gi to highest data likelihood genotype combo
	GoMixGi[ind] = GiBest;
	lnProbBoGivenGoMixGi += lnProbBiGivenGiBest;
      }
    }
    
    // haploid DNA
    else {
      //------------------------------------------------------------------------
      // add Gi1 to all-Gi1 Go combo and update probability
      //------------------------------------------------------------------------
      GoAllGi1[ind] = Gi1;
      lnProbBoGivenGoAllGi1 += LnProbBiGivenGi[ind][Gi1];
      
      //------------------------------------------------------------------------
      // add Gi2 to all-Gi2 Go combo and update probability
      //------------------------------------------------------------------------
      GoAllGi2[ind] = Gi2;
      lnProbBoGivenGoAllGi2 += LnProbBiGivenGi[ind][Gi2];
      
      //------------------------------------------------------------------------
      // Determine highest likelihood individual genotype (Gi) and add
      //   to highest data likelihood genotype combo if banded alg.
      //------------------------------------------------------------------------
      if (banded) {
	// initialize as Gi11
	string GiBest = Gi1; long double lnProbBiGivenGiBest = LnProbBiGivenGi[ind][Gi1];
      
	// test Gi2
	if (LnProbBiGivenGi[ind][Gi2] > lnProbBiGivenGiBest) {
	  GiBest = Gi2; lnProbBiGivenGiBest = LnProbBiGivenGi[ind][Gi2];
	}
      
	// add best Gi to highest data likelihood genotype combo
	GoMixGi[ind] = GiBest;
	lnProbBoGivenGoMixGi += lnProbBiGivenGiBest;
      }
    }
  }
  
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Create the highest-contributing subset of Go genotype combinations with
  //   either recursive method or banded approximation, & calculate ln P(Bo|Go).
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  
  //----------------------------------------------------------------------------
  // Define ln P(Bo|Go) for all genotype combinations.
  // 
  //----------------------------------------------------------------------------
  GenotypeComboProbability LnProbBoGivenGo;

  //----------------------------------------------------------------------------
  // banded approximation
  //----------------------------------------------------------------------------
  if (banded) {

    //--------------------------------------------------------------------------
    // Ceate a list of all genotype combinations <= WB edit distance away
    //   from highest data likelihood genotype combo (GoMixGi).
    //--------------------------------------------------------------------------
      
    // Register nominal "best" genotype combination
    LnProbBoGivenGo[GoMixGi] = lnProbBoGivenGoMixGi;

    // Iterate through edit distances (limit edit distances to 
    //   the number of individuals.
    for (int w=1; w<=min(WB, (int)individuals.size()) - 1; w++) {
      
      // make new LnProbBoGivenGo map
      GenotypeComboProbability LnProbBoGivenGoNew;

      // Iterate through every element of current list
      for (GenotypeComboProbability::const_iterator GoIter = LnProbBoGivenGo.begin();
	   GoIter != LnProbBoGivenGo.end(); GoIter++) {
	
	// Retreive genotype.and probability
	GenotypeCombo Go = GoIter->first;
	long double probGo = GoIter->second;

	//----------------------------------------------------------------------
	// Make every new 1-edit distance genotype combo from current list.
	//----------------------------------------------------------------------

	// Iterate through every index in Go genotype combination.
	for (GenotypeCombo::const_iterator indIter = Go.begin();
	     indIter != Go.end(); indIter++) {
	  string ind = indIter->first;
	  string Gi = indIter->second;

	  long double probGi = LnProbBiGivenGi[ind][Gi];

	  if (diploid) {
	    // Make new combo with Gi=Gi11 for individual=ind and add to list if new.
	    GenotypeCombo GoReplaceWithGi11 = Go;
	    GoReplaceWithGi11[ind] = Gi11;
	    LnProbBoGivenGoNew[GoReplaceWithGi11] = probGo - probGi + LnProbBiGivenGi[ind][Gi11];
	    
	    // Make new combo with Gi22 at GiIndex and add to list if new.
	    GenotypeCombo GoReplaceWithGi22 = Go;
	    GoReplaceWithGi22[ind] = Gi22;
	    LnProbBoGivenGoNew[GoReplaceWithGi22] = probGo - probGi + LnProbBiGivenGi[ind][Gi22];
	    
	    // Make new combo with Gi12 at GiIndex and add to list if new.
	    map<string, string > GoReplaceWithGi12 = Go;
	    GoReplaceWithGi12[ind] = Gi12;
	    LnProbBoGivenGoNew[GoReplaceWithGi12] = probGo - probGi + LnProbBiGivenGi[ind][Gi12];
	  }
	  else {
	    // Make new combo with Gi1 at GiIndex and add to list if new.
	    GenotypeCombo GoReplaceWithGi1 = Go;
	    GoReplaceWithGi1[ind] = Gi1;
	    LnProbBoGivenGoNew[GoReplaceWithGi1] = probGo - probGi + LnProbBiGivenGi[ind][Gi1];
	    
	    // Make new combo with Gi2 at GiIndex and add to list if new.
	    GenotypeCombo GoReplaceWithGi2 = Go;
	    GoReplaceWithGi2[ind] = Gi2;
	    LnProbBoGivenGoNew[GoReplaceWithGi2] = probGo - probGi + LnProbBiGivenGi[ind][Gi2];
	  }
	}
      }

      //------------------------------------------------------------------------
      // update LnProbBoGivenGo with LnProbBoGivenGoNew (or with dominant terms
      //   thereof)
      //------------------------------------------------------------------------
      
      // If TB == 0 (means do not cap number of terms) or number of terms is
      //   lower than allowed, simply replace hashes
      if ((TB == 0) || (LnProbBoGivenGoNew.size() <= TR)) {
	LnProbBoGivenGo = LnProbBoGivenGoNew;
      }
      
      // otherwise only take the best TB entries (after including priors) 
      else {
	
	// rank genotype combinations Go according to decreasing P(Bo|Go)
	vector<GenotypeCombo> GoListSorted = sortKeysByValue(LnProbBoGivenGoNew, true);
	
	// initialize term counter
	int t = 0;
	for(vector<GenotypeCombo >::const_iterator GoIter = GoListSorted.begin();
	    GoIter != GoListSorted.end(); GoIter++) {
	  
	  // increment term counter
	  t++;
	  if (t > TB) {
	    break;
	  }
	  
	  // retreive combo
	  GenotypeCombo Go = *GoIter;
	  
	  // copy element
	  LnProbBoGivenGo[Go] = LnProbBoGivenGoNew[Go];
	}
      }
    }
    
    //--------------------------------------------------------------------------
    // always add back two monomorphic genotype combinations to list
    //--------------------------------------------------------------------------    

    // diploid DNA
    if (diploid) {
      LnProbBoGivenGo[GoAllGi11] = lnProbBoGivenGoAllGi11;
      LnProbBoGivenGo[GoAllGi22] = lnProbBoGivenGoAllGi22;
    }

    // haploid DNA
    else {
      LnProbBoGivenGo[GoAllGi1] = lnProbBoGivenGoAllGi1;
      LnProbBoGivenGo[GoAllGi2] = lnProbBoGivenGoAllGi2;
    }
  }

  //----------------------------------------------------------------------------
  // recursive method
  //----------------------------------------------------------------------------
  else {
    
    // initializa probability for null genotype combination
    LnProbBoGivenGo[GoNull] = 0;
  
    // cycle through each individual
    for (vector<string>::const_iterator indIter = individuals.begin();
	 indIter != individuals.end(); indIter++) {

      // get individual
      string ind = *indIter;

      // new hash containing values updated for this individual
      GenotypeComboProbability LnProbBoGivenGoNew;

      // iterate through all possible genotypes for this invididual
      for (map<string, long double >::const_iterator GiIter = LnProbBiGivenGi[ind].begin();
	   GiIter !=  LnProbBiGivenGi[ind].end(); GiIter++) {
	
	// retreive genotype and probability
	string Gi = GiIter->first;
	long double lnProb = GiIter->second;
	
	// iterate through all current genotype combinations
	for (GenotypeComboProbability::const_iterator GoIter = LnProbBoGivenGo.begin();
	     GoIter != LnProbBoGivenGo.end(); GoIter++) {

	  // retreive genotype combination and corresponding prob ln value
	  GenotypeCombo Go = GoIter->first;
	  long double LnProb = GoIter->second;
	  
	  // make new genotype combo
	  GenotypeCombo GoNew = Go;
	  GoNew[ind] = Gi;
	  
	  // update LnProbBoGivenGo with this individual
	  LnProbBoGivenGoNew[GoNew] = LnProb + lnProb;
	}
      }
      
      //------------------------------------------------------------------------
      // update LnProbBoGivenGo with LnProbBoGivenGoNew (or with dominant terms
      //   thereof)
      //------------------------------------------------------------------------
      
      // If TR == 0 (means do not cap number of terms) or number of terms is
      // lower than allowed, simply replace hashes
      if ((TR == 0) || (LnProbBoGivenGoNew.size() <= TR)) {
	LnProbBoGivenGo = LnProbBoGivenGoNew;
      }
      
      // otherwise only take the best TR entries (after including priors) 
      else {
	
	// rank genotype combinations Go according to decreasing P(Bo|Go)
	vector<GenotypeCombo> GoListSorted = sortKeysByValue(LnProbBoGivenGoNew, true);
	
	// initialize term counter
	int t = 0;
	for(vector<GenotypeCombo >::const_iterator GoIter = GoListSorted.begin();
	    GoIter != GoListSorted.end(); GoIter++) {
	  
	  // increment term counter
	  t++;
	  if (t > TR) {
	    break;
	  }
	  
	  // retreive combo
	  GenotypeCombo Go = *GoIter;
	  
	  // copy element
	  LnProbBoGivenGo[Go] = LnProbBoGivenGoNew[Go];
	}
      }
    }

    //--------------------------------------------------------------------------
    // add probabilities of two monomorphic genotype combinations
    //--------------------------------------------------------------------------

    // diploid DNA
    if (diploid) {
      LnProbBoGivenGo[GoAllGi11] = lnProbBoGivenGoAllGi11;
      LnProbBoGivenGo[GoAllGi22] = lnProbBoGivenGoAllGi22;
    }

    // haploid DNA
    else {
      LnProbBoGivenGo[GoAllGi1] = lnProbBoGivenGoAllGi1;
      LnProbBoGivenGo[GoAllGi2] = lnProbBoGivenGoAllGi2;
    }
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Calculate ln P(Go|Bo) for set of genotype combinations produced
  //   by either the banded or the recursive method
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // define ln P(G|B) and ln Prior(G) for all genotype combinations
  //----------------------------------------------------------------------------
  GenotypeComboProbability LnProbGoGivenBo, LnProbGo;

  
  // iterate through every genotype combination considered
  for (GenotypeComboProbability::const_iterator GoIter = LnProbBoGivenGo.begin();
       GoIter != LnProbBoGivenGo.end(); GoIter++) {
    
    // retreive genotype
    GenotypeCombo Go = GoIter->first;
    
    //--------------------------------------------------------------------------
    // calculate prior
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // unrelateds: calculate the part of the prior that comes from the unrelateds
    //--------------------------------------------------------------------------
    int n=0, f1=0, h=0;
    for (vector<string>::const_iterator indIter = unrelateds.begin();
	 indIter != unrelateds.end(); indIter++) {
      
      // retreive individual and genotype
      string ind = *indIter;
      string Gi = Go[ind];
      
      // update sample size, allele frequency, and number of hets
      if (Gi == Gi11) {
	n += 2;
	f1 += 2;
      }
      else if (Gi == Gi22) {
	n += 2;
      }
      else if (Gi == Gi12) {
	n += 2;
	f1 += 1;
	h += 1;
      }
      else if (Gi == Gi1) {
	n += 1;
	f1 += 1;
      }
      else if (Gi == Gi2) {
	n += 1;
      }
    }
    
    // calculate that part of the prior that comes from the unrelateds
    long double lnPriorUnrelated = 0;
    if (f1 == 0 || f1 == n) {
      // monomorphic case
      lnPriorUnrelated += log (1 - THETA * sumOneOverI(n-1)) - log(2.0);
    }
    else {
      // polymorphic case
      lnPriorUnrelated += log(THETA) + h * log(2.0) + log(n) - cofactorln(n, f1) - log(f1) - log(n-f1);
    }
    
    //--------------------------------------------------------------------------
    // children: calculate that part of the prior that comes from the children
    //--------------------------------------------------------------------------
    long double lnPriorChildren = 0;
    for (vector<string>::const_iterator indIter = children.begin();
	 indIter != children.end(); indIter++) {
      
      // retreive individual and genotype
      string ind = *indIter;
      string Gi = Go[ind];
      
      // retreive parent genotypes
      string GiMother = Go[motherId[ind]];
      string GiFather = Go[fatherId[ind]];
      
      // calculate Pr(GiChild | GiMother, GiFather) from transmission 
      // probilities, parental genotypes and somatic mutation likelihood)
      lnPriorChildren += log(childGivenParentsGenotypeProbability(Gi, GiMother, GiFather, MU));
    }
    
    //------------------------------------------------------------------------
    // calculate and register total prior for this genotype combo
    //------------------------------------------------------------------------
    long double lnPrior = lnPriorUnrelated + lnPriorChildren;
    LnProbGo[Go] = lnPrior;
    
    //------------------------------------------------------------------------
    // calculate and register P(Go|Bo) probability
    //------------------------------------------------------------------------
    LnProbGoGivenBo[Go] = LnProbBoGivenGo[Go] + lnPrior;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Calculate P(Go|Bo) values from raw ln P(Go|Bo) values
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // normalize ln P(Go|Bo) by the largest value to avoid 0 probabilities after
  // exponentiation
  //----------------------------------------------------------------------------

  // sort genotype combinations according to descending ln P(G|B)
  vector<GenotypeCombo> GoListSorted = sortKeysByValue(LnProbGoGivenBo, true);
  GenotypeCombo GoKing = GoListSorted[0];
  long double lnProbGoGivenBoKing = LnProbGoGivenBo[GoKing];

  // 1. bring up log quantities to the largest probLn value (to avoid all 0 probs)
  // 2. calculate probability normalization factor
  long double sumProb = 0;
  for(GenotypeComboProbability::const_iterator GoIter = LnProbGoGivenBo.begin();
      GoIter != LnProbGoGivenBo.end(); GoIter++) {
    GenotypeCombo Go = GoIter->first;
    long double lnProb =  LnProbGoGivenBo[Go] - lnProbGoGivenBoKing;

    sumProb += exp(lnProb);
  }

  //----------------------------------------------------------------------------
  // 1. normalize probabilities by prob normalization factor
  // 2. calculate final P(Go|Bo) for each genotype combination
  // 3. calculate individual genotype probabilities
  //----------------------------------------------------------------------------

  // define final P(Go|Bo) hash
  GenotypeComboProbability ProbGoGivenBo;

  // define individual genotype probability distribution hash
  map<string, map<string, long double > > ProbGiGivenBo;

  // 1. initialize individual genotype probability hash
  // 2. generate monomorphic genotype combinations
  GenotypeCombo GoMonoGi1, GoMonoGi2;
  for (vector<string>::const_iterator indIter = individuals.begin();
       indIter != individuals.end(); indIter++) {
    string ind = *indIter;
    if (diploid) {

      // initialize genotype probability
      ProbGiGivenBo[ind][Gi11] = 0;
      ProbGiGivenBo[ind][Gi22] = 0;
      ProbGiGivenBo[ind][Gi12] = 0;

      // update monomorphic genotype combinations
      GoMonoGi1[ind] = Gi11;
      GoMonoGi2[ind] = Gi22;
    }
    else {

      // initialize genotype probability
      ProbGiGivenBo[ind][Gi1] = 0;
      ProbGiGivenBo[ind][Gi2] = 0;

      // update monomorphic genotype combinations
      GoMonoGi1[ind] = Gi1;
      GoMonoGi2[ind] = Gi2;
    }
  }

  // process every genotype combination for which probability was calculated
  for(GenotypeComboProbability::const_iterator GoIter = LnProbGoGivenBo.begin();
      GoIter != LnProbGoGivenBo.end(); GoIter++) {
    GenotypeCombo Go = GoIter->first;
    long double lnProbGoGivenBoWork = GoIter->second;
    
    // calculate final P(G|B) for this genotype combination
    ProbGoGivenBo[Go] = exp(lnProbGoGivenBoWork - lnProbGoGivenBoKing) / sumProb;
    
    // update individual genotype probabilities if regured
    if (cgp) {
      for (int i=0; i<individuals.size(); i++) {
	
	// retreive individual and genotype
	string ind = individuals[i];
	string Gi = Go[ind];
	
	// update genotype probability
	//      ProbGiGivenBo[ind][Gi] += lnProbGoGivenBoWork;
	ProbGiGivenBo[ind][Gi] += ProbGoGivenBo[Go];
      }
    }
  }

  // calculate P(SNP)
  long double pSnp = 1.0;
  if (ProbGoGivenBo.count(GoMonoGi1) > 0) {
    pSnp -= ProbGoGivenBo[GoMonoGi1];
  }
  if (ProbGoGivenBo.count(GoMonoGi2) > 0) {
    pSnp -= ProbGoGivenBo[GoMonoGi2];
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // debugging printouts if necessary
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  if (debug2) {

    //--------------------------------------------------------------------------
    // print probabilities for genotype combinations
    //--------------------------------------------------------------------------
    clog << "Data likelihood, prior, and posterior P(Bo|Go), P(Go), P(Go|Bo):" << endl;
    for(GenotypeComboProbability::const_iterator GoIter = ProbGoGivenBo.begin();
      GoIter != ProbGoGivenBo.end(); GoIter++) {
      GenotypeCombo Go = GoIter->first;
      long double probGoGivenBo = GoIter->second;

      long double lnProbBoGivenGo = LnProbBoGivenGo[Go];
      long double probBoGivenGo = exp(lnProbBoGivenGo);

      long double lnProbGo = LnProbGo[Go];
      long double probGo = exp(lnProbGo);

      long double lnProbGoGivenBo = LnProbGoGivenBo[Go];

      clog << "  Go=(";
      for (int i=0; i<individuals.size(); i++) {
	string ind = individuals[i];
	string Gi = Go[ind];
	//	clog << "ind=" << ind << " Gi=" << Gi << ":";
	clog << ind << ":" << Gi << " ";
      }
      clog << ")";
      clog << " lnP(Bo|Go)=" << lnProbBoGivenGo;
      clog << " P(Bo|Go)=" << probBoGivenGo;
      clog << " lnP(Go)=" << lnProbGo;
      clog << " P(Go)=" << probGo;
      clog << " lnP(Go|Bo)=" << lnProbGoGivenBo;
      clog << " P(Go|Bo)=" << probGoGivenBo << endl;
    }

    //--------------------------------------------------------------------------
    // print individual genotype probabilities
    //--------------------------------------------------------------------------
    clog << "Individual genotype probabilities P(Gi|Bo):" << endl;

    // cycle through each individual
    for (vector<string>::const_iterator indIter = individuals.begin();
	 indIter != individuals.end(); indIter++) {
      string ind = *indIter;
      clog << "  ind=" << ind;
      for (map<string, long double >::const_iterator GiIter = ProbGiGivenBo[ind].begin();
	   GiIter !=  ProbGiGivenBo[ind].end(); GiIter++) {
	string Gi = GiIter->first;
	long double probGiGivenBo = GiIter->second;
	clog << " P(Gi=" << Gi << "|Bo)=" << probGiGivenBo;
      }
      clog << endl;
    }
    
    //--------------------------------------------------------------------------
    // print P(SNP)
    //--------------------------------------------------------------------------
    clog << "SNP probability:" << endl;
    clog << "  P(SNP)=" << pSnp << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // assign output Variation data structure
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // define
  //----------------------------------------------------------------------------
  Variation var;

  //----------------------------------------------------------------------------
  // assign simple fields
  //----------------------------------------------------------------------------

  // assign ploidy
  var.diploid = diploid;

  // assign alleles
  var.alleles.push_back(allele1);
  var.alleles.push_back(allele2);

  // assign P(SNP)
  var.pSnp = pSnp;

  //----------------------------------------------------------------------------
  // assign invidual genotype log data likelihoods
  //----------------------------------------------------------------------------
  var.individualGenotypeLogDataLikelihood = LnProbBiGivenGi;

  //----------------------------------------------------------------------------
  // assign invidual genotype probabilities
  //----------------------------------------------------------------------------
  map<string, map<string, long double > > individualGenotypeProbability;
  for (map<string, map<string, long double > >::const_iterator 
	 indIter = ProbGiGivenBo.begin(); indIter != ProbGiGivenBo.end(); indIter++) {
    string ind = indIter->first;

    // skip if this individual is a dummy parent
    if (isDummyParent[ind]) {
      continue;
    }

    map<string, long double > genotypeProb = indIter->second;
    for (map<string, long double >::const_iterator 
	   GiIter = genotypeProb.begin(); GiIter != genotypeProb.end(); GiIter++) {
      string Gi = GiIter->first;
      long double p = GiIter->second;

      // recode genotype according to expected output format
      string t;
      if (Gi == Gi1) {t = allele1;}
      else if (Gi == Gi2) {t = allele2;}
      else if (Gi == Gi11) {t = allele1 + allele1;}
      else if (Gi == Gi22) {t = allele2 + allele2;}
      else if (Gi == Gi12) {t = allele1 + allele2;}

      // register
      individualGenotypeProbability[ind][t] = p;
    }
  }

  // assign individualGenotypeProbability
  var.individualGenotypeProbability = individualGenotypeProbability;

  //----------------------------------------------------------------------------
  // return Variation data structure
  //----------------------------------------------------------------------------
  return var;
}

//------------------------------------------------------------------------------
// posteriorProb2 -- genotype probs based always on at least 3 values
//------------------------------------------------------------------------------
Variation posteriorProb2(
        ///                   vvv this could become a class, as there are only so many genotypes (3 different 0/0, 0/1, 1/1)
        //       and the individuals could be a vector
			map<string, map<string, long double > > LnProbBiGivenGi,
			map<string, bool > DataSufficientGi,
			vector<string> individuals,
			bool diploid,
			long double THETA,
			bool banded,
			int WB,
			int TB,
			bool includeMonoB,
			int TR,
			string allele1,
			string allele2,
			bool debug2
			) {

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Notations:
  //     Gi: genotype of a given individual (e.g. AA, AG, GG)
  //     Bi: set of base calls for the individual (e.g. A,A,A,G,A,A,G)
  //     Ai: set of true bases estimated by the base calls (e.g. A,A,G,A,G,G)
  //
  //     Go: genotype vector of all individuals
  //     Bo: base calls for all individuals
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // constants
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // constants
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  long double LOGFACTOR = log(10) / (-10.0); 

  // haploid genotypes
  string Gi1 = IUPAC[allele1];
  string Gi2 = IUPAC[allele2];

  // diploid genotypes
  string Gi11 = IUPAC[allele1+allele1];
  string Gi22 = IUPAC[allele2+allele2];
  string Gi12 = IUPAC[allele1+allele2];

  // undetermined genotype
  string GiN = "?";

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Build special genotype combinations: mono allele 1, mono allele 2.
  // Build highest data likelihood genotype combination if banded algorithm 
  // Register corresponding ln P(Bo|Go) data likelihoods
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  
  // genotype combination vector with no elements
  string GoNull;

  // special genotype combinations
  string GoAllGi1, GoAllGi2, GoAllGi11, GoAllGi22, GoMixGi;

  // P(Bo|Go) for special genotype combinations
  long double lnProbBoGivenGoAllGi1 = 0;
  long double lnProbBoGivenGoAllGi2 = 0;
  long double lnProbBoGivenGoAllGi11 = 0;
  long double lnProbBoGivenGoAllGi22 = 0;
  long double lnProbBoGivenGoMixGi = 0;
  
  // build special genotype combinations
  for (vector<string>::const_iterator indIter = individuals.begin();
       indIter != individuals.end(); indIter++) {
    
    // get individual
    string ind = *indIter;
    
    // diploid DNA
    if (diploid) {
      if (DataSufficientGi[ind]) {
	//----------------------------------------------------------------------
	// add Gi11 to all-Gi11 genotype combo and update probability
	//----------------------------------------------------------------------
	GoAllGi11 += Gi11;
	lnProbBoGivenGoAllGi11 += LnProbBiGivenGi[ind][Gi11];
      
	//----------------------------------------------------------------------
	// add Gi22 to all-Gi22 genotype combo and update probability
	//----------------------------------------------------------------------
	GoAllGi22 += Gi22;
	lnProbBoGivenGoAllGi22 += LnProbBiGivenGi[ind][Gi22];
      
	//----------------------------------------------------------------------
	// determine highest likelihood individual genotype (Gi) and add
	// to highest data likelihood genotype combo
	//----------------------------------------------------------------------

	// initialize as Gi12
	string GiBest = Gi12; long double lnProbBiGivenGiBest = LnProbBiGivenGi[ind][Gi12];
      
	// test Gi11
	if (LnProbBiGivenGi[ind][Gi11] > lnProbBiGivenGiBest) {
	  GiBest = Gi11; lnProbBiGivenGiBest = LnProbBiGivenGi[ind][Gi11];
	}
	
	// test Gi22
	if (LnProbBiGivenGi[ind][Gi22] > lnProbBiGivenGiBest) {
	  GiBest = Gi22; lnProbBiGivenGiBest = LnProbBiGivenGi[ind][Gi22];
	}
	
	// add best Gi to highest data likelihood genotype combo
	GoMixGi += GiBest;
	lnProbBoGivenGoMixGi += lnProbBiGivenGiBest;
      }
      else {
	GoAllGi11 += GiN;
	GoAllGi22 += GiN;
	GoMixGi += GiN;
      }
    }
    
    // haploid DNA
    else {
      if (DataSufficientGi[ind]) {
	//----------------------------------------------------------------------
	// add Gi1 to all-Gi1 Go combo and update probability
	//----------------------------------------------------------------------
	GoAllGi1 += Gi1;
	lnProbBoGivenGoAllGi1 += LnProbBiGivenGi[ind][Gi1];
      
	//----------------------------------------------------------------------
	// add Gi2 to all-Gi2 Go combo and update probability
	//----------------------------------------------------------------------
	GoAllGi2 += Gi2;
	lnProbBoGivenGoAllGi2 += LnProbBiGivenGi[ind][Gi2];
      
	//----------------------------------------------------------------------
	// Determine highest likelihood individual genotype (Gi) and add
	//   to highest data likelihood genotype combo
	//----------------------------------------------------------------------

	// initialize as Gi11
	string GiBest = Gi1; long double lnProbBiGivenGiBest = LnProbBiGivenGi[ind][Gi1];
	
	// test Gi2
	if (LnProbBiGivenGi[ind][Gi2] > lnProbBiGivenGiBest) {
	  GiBest = Gi2; lnProbBiGivenGiBest = LnProbBiGivenGi[ind][Gi2];
	}
	
	// add best Gi to highest data likelihood genotype combo
	GoMixGi += GiBest;
	lnProbBoGivenGoMixGi += lnProbBiGivenGiBest;
      }
      else {
	GoAllGi11 += GiN;
	GoAllGi22 += GiN;
	GoMixGi += GiN;
      }
    }
  }
  
  if (debug2) {
    cerr << "BEST Go=" << GoMixGi << " lnProbBoGivenGoMixGi=" << lnProbBoGivenGoMixGi << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Create the highest-contributing subset of Go genotype combinations with
  //   either recursive method or banded approximation, & calculate ln P(Bo|Go).
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  
  //----------------------------------------------------------------------------
  // Define ln P(Bo|Go) for all genotype combinations.
  //----------------------------------------------------------------------------
  map<string, long double > LnProbBoGivenGo;

  //----------------------------------------------------------------------------
  // banded approximation
  //----------------------------------------------------------------------------
  if (banded) {

    //--------------------------------------------------------------------------
    // Create a list of all genotype combinations <= WB edit distance away
    //   from highest data likelihood genotype combo (GoMixGi) and from
    //   monomorphics
    //--------------------------------------------------------------------------
      
    // Register nominal "best" genotype combination in initial list
    LnProbBoGivenGo[GoMixGi] = lnProbBoGivenGoMixGi;

    // if monomorphic genotype combinations are to be included in the initial list
    if (includeMonoB) {

      // Register monomorphics in initial list
      if (diploid) {
	LnProbBoGivenGo[GoAllGi11] = lnProbBoGivenGoAllGi11;
	LnProbBoGivenGo[GoAllGi22] = lnProbBoGivenGoAllGi22;
      }
      else {
	LnProbBoGivenGo[GoAllGi1] = lnProbBoGivenGoAllGi1;
	LnProbBoGivenGo[GoAllGi2] = lnProbBoGivenGoAllGi2;
      }
    }

    // Iterate through edit distances. Limit edit distances to 
    //   the number of individuals, 
    for (unsigned int w=1; w<=min(WB, (int)individuals.size()); w++) {
      
      // make new LnProbBoGivenGo
      map<string, long double > LnProbBoGivenGoNew;

      // Iterate through every element of current list
      for (map<string, long double >::const_iterator GoIter = LnProbBoGivenGo.begin();
	   GoIter != LnProbBoGivenGo.end(); GoIter++) {
	
	// Retreive genotype.and probability
	string Go = GoIter->first;
	long double lnProbGo = GoIter->second;

	//----------------------------------------------------------------------
	// Make every new 1-edit distance genotype combo from current list if
	//   there is data sufficience. Otherwise skip
	//----------------------------------------------------------------------

	// Iterate through every individual in Go genotype combination.
	for (unsigned int indIndex = 0; indIndex < individuals.size(); ++indIndex) {

	  // get individual
	  string ind = individuals[indIndex];
	  
	  // only continue if there is sufficient data for this sample
	  if (DataSufficientGi[ind]) {
	    
	    // get individual genotype

        string Gi;
	    TRY { Gi = Go.substr(indIndex, 1); } CATCH;
	    
	    // get individual log genotype probability
	    long double lnProbGi = LnProbBiGivenGi[ind][Gi];
	    
	    if (diploid) {
	      // Make new combo with Gi=Gi11 for individual=ind and add to list if new.
	      string GoReplaceWithGi11 = Go;
	      GoReplaceWithGi11[indIndex] = Gi11[0]; // HACK
	      LnProbBoGivenGoNew[GoReplaceWithGi11] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi11];
	      
	      // Make new combo with Gi22 at GiIndex and add to list if new.
	      string GoReplaceWithGi22 = Go;
	      GoReplaceWithGi22[indIndex] = Gi22[0]; // HACK
	      LnProbBoGivenGoNew[GoReplaceWithGi22] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi22];
	      
	      // Make new combo with Gi12 at GiIndex and add to list if new.
	      string GoReplaceWithGi12 = Go;
	      GoReplaceWithGi12[indIndex] = Gi12[0];
	      LnProbBoGivenGoNew[GoReplaceWithGi12] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi12];
	    }
	    else {
	      // Make new combo with Gi1 at GiIndex and add to list if new.
	      string GoReplaceWithGi1 = Go;
	      GoReplaceWithGi1[indIndex] = Gi1[0]; // HACK
	      LnProbBoGivenGoNew[GoReplaceWithGi1] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi1];
	      
	      // Make new combo with Gi2 at GiIndex and add to list if new.
	      string GoReplaceWithGi2 = Go;
	      GoReplaceWithGi2[indIndex] = Gi2[0]; // HACK
	      LnProbBoGivenGoNew[GoReplaceWithGi2] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi2];
	    }
	  }
	}
      }

      //------------------------------------------------------------------------
      // update LnProbBoGivenGo with LnProbBoGivenGoNew (or with dominant terms
      //   thereof)
      //------------------------------------------------------------------------
      
      // If TB == 0 (means do not cap number of terms) or number of terms is
      //   lower than allowed, simply replace hashes
      if ((TB == 0) || (LnProbBoGivenGoNew.size() <= TR)) {
	LnProbBoGivenGo = LnProbBoGivenGoNew;
      }
      
      // otherwise only take the best TB entries (after including priors) 
      else {
	
	// rank genotype combinations Go according to decreasing P(Bo|Go)
	vector<string> GoListSorted = sortKeysByValue(LnProbBoGivenGoNew, true);
	
	// initialize term counter
	int t = 0;
	for(vector<string>::const_iterator GoIter = GoListSorted.begin();
	    GoIter != GoListSorted.end(); GoIter++) {
	  
	  // increment term counter
	  t++;
	  if (t > TB) {
	    break;
	  }
	  
	  // retreive combo
	  string Go = *GoIter;
	  
	  // copy element
	  LnProbBoGivenGo[Go] = LnProbBoGivenGoNew[Go];
	}
      }
    }
    
    //--------------------------------------------------------------------------
    // always add back two monomorphic genotype combinations to list
    //--------------------------------------------------------------------------    

    // diploid DNA
    if (diploid) {
      LnProbBoGivenGo[GoAllGi11] = lnProbBoGivenGoAllGi11;
      LnProbBoGivenGo[GoAllGi22] = lnProbBoGivenGoAllGi22;
    }

    // haploid DNA
    else {
      LnProbBoGivenGo[GoAllGi1] = lnProbBoGivenGoAllGi1;
      LnProbBoGivenGo[GoAllGi2] = lnProbBoGivenGoAllGi2;
    }
  }

  //----------------------------------------------------------------------------
  // recursive method
  // !!!FIX FOR undetermined genotypes!!!
  //----------------------------------------------------------------------------
  else {
    
    // initializ probability for null genotype combination
    LnProbBoGivenGo[GoNull] = 0;
  
    // cycle through each individual
    for (vector<string>::const_iterator indIter = individuals.begin();
	 indIter != individuals.end(); indIter++) {

      // get individual
      string ind = *indIter;

      // new hash containing values updated for this individual
      map<string, long double > LnProbBoGivenGoNew;

      // iterate through all possible genotypes for this invididual
      for (map<string, long double >::const_iterator GiIter = LnProbBiGivenGi[ind].begin();
	   GiIter !=  LnProbBiGivenGi[ind].end(); GiIter++) {
	
	// retreive genotype and probability
	string Gi = GiIter->first;
	long double lnProb = GiIter->second;
	
	// iterate through all current genotype combinations
	for (map<string, long double >::const_iterator GoIter = LnProbBoGivenGo.begin();
	     GoIter != LnProbBoGivenGo.end(); GoIter++) {

	  // retreive genotype combination and corresponding prob ln value
	  string Go = GoIter->first;
	  long double LnProb = GoIter->second;
	  
	  // make new genotype combo
	  string GoNew = Go + Gi;
	  
	  // update LnProbBoGivenGo with this individual
	  LnProbBoGivenGoNew[GoNew] = LnProb + lnProb;
	}
      }
      
      //------------------------------------------------------------------------
      // update LnProbBoGivenGo with LnProbBoGivenGoNew (or with dominant terms
      //   thereof)
      //------------------------------------------------------------------------
      
      // If TR == 0 (means do not cap number of terms) or number of terms is
      // lower than allowed, simply replace hashes
      if ((TR == 0) || (LnProbBoGivenGoNew.size() <= TR)) {
	LnProbBoGivenGo = LnProbBoGivenGoNew;
      }
      
      // otherwise only take the best TR entries (after including priors) 
      else {
	
	// rank genotype combinations Go according to decreasing P(Bo|Go)
	vector<string> GoListSorted = sortKeysByValue(LnProbBoGivenGoNew, true);
	
	// initialize term counter
	int t = 0;
	for(vector<string>::const_iterator GoIter = GoListSorted.begin();
	    GoIter != GoListSorted.end(); GoIter++) {
	  
	  // increment term counter
	  t++;
	  if (t > TR) {
	    break;
	  }
	  
	  // retreive combo
	  string Go = *GoIter;
	  
	  // copy element
	  LnProbBoGivenGo[Go] = LnProbBoGivenGoNew[Go];
	}
      }
    }

    //--------------------------------------------------------------------------
    // add probabilities of two monomorphic genotype combinations
    //--------------------------------------------------------------------------

    // diploid DNA
    if (diploid) {
      LnProbBoGivenGo[GoAllGi11] = lnProbBoGivenGoAllGi11;
      LnProbBoGivenGo[GoAllGi22] = lnProbBoGivenGoAllGi22;
    }

    // haploid DNA
    else {
      LnProbBoGivenGo[GoAllGi1] = lnProbBoGivenGoAllGi1;
      LnProbBoGivenGo[GoAllGi2] = lnProbBoGivenGoAllGi2;
    }
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Calculate ln P(Go|Bo) for set of genotype combinations produced
  //   by either the banded or the recursive method
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // define ln P(G|B) and ln Prior(G) for all genotype combinations
  //----------------------------------------------------------------------------
  map<string, long double > LnProbGoGivenBo, LnProbGo;

  // iterate through every genotype combination considered
  for (map<string, long double >::const_iterator GoIter = LnProbBoGivenGo.begin();
       GoIter != LnProbBoGivenGo.end(); GoIter++) {
    
    // retreive genotype
    string Go = GoIter->first;
    
    // XXX factor out for legibility?
    //--------------------------------------------------------------------------
    // calculate prior
    //--------------------------------------------------------------------------
    int n=0, f1=0, h=0;

    // Iterate through every individual in Go genotype combination.
    for (unsigned int indIndex = 0; indIndex < individuals.size(); ++indIndex) {
      
      // get individual
      string ind = individuals[indIndex];
      
      // only count if there is sufficient data for this sample
      if (DataSufficientGi[ind]) {
	    
	// get individual genotype
    string Gi;
	TRY { Gi = Go.substr(indIndex, 1); } CATCH;
	
	// update sample size, allele frequency, and number of hets
	if (diploid) {
	  if (Gi == Gi11) {
	    n += 2;
	    f1 += 2;
	  }
	  else if (Gi == Gi22) {
	    n += 2;
	  }
	  else if (Gi == Gi12) {
	    n += 2;
	    f1 += 1;
	    h += 1;
	  }
	}
	else {
	  if (Gi == Gi1) {
	    n += 1;
	    f1 += 1;
	  }
	  else if (Gi == Gi2) {
	    n += 1;
	  }
	}
      }
    }
    
    //------------------------------------------------------------------------
    // calculate and register total prior for this genotype combo
    //------------------------------------------------------------------------
    long double lnPrior = 0;
    if (n == 0) {
      // no reads
      lnPrior = log (1 - THETA) - log(2.0);
    }
    else if (f1 == 0 || f1 == n) {
      // monomorphic case, reads > 0
      lnPrior = log (1 - THETA * sumOneOverI(n-1)) - log(2.0);
    }
    else {
      // polymorphic case, reads > 0
      lnPrior = log(THETA) + h * log(2.0) + log(n) - cofactorln(n, f1) - log(f1) - log(n-f1);
    }
    LnProbGo[Go] = lnPrior;
    
    //------------------------------------------------------------------------
    // calculate and register raw P(Go|Bo) probability
    //------------------------------------------------------------------------
    LnProbGoGivenBo[Go] = LnProbBoGivenGo[Go] + lnPrior;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Calculate P(Go|Bo) values from raw ln P(Go|Bo) values
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // normalize ln P(Go|Bo) by the largest value to avoid 0 probabilities after
  // exponentiation
  //----------------------------------------------------------------------------

  // sort genotype combinations according to descending ln P(G|B)
  vector<string> GoListSorted = sortKeysByValue(LnProbGoGivenBo, true);

  // get genotype combination with the highest probability 
  string GoKing = GoListSorted[0];

  // get the highest probability corresponding to that genotype combination
  long double lnProbBoGivenGoKing = LnProbBoGivenGo[GoKing];
  long double lnProbGoGivenBoKing = LnProbGoGivenBo[GoKing];

  // 1. bring up log quantities to the largest probLn value (to avoid all 0 probs)
  // 2. calculate probability normalization factor
  long double sumProb = 0;
  for(map<string, long double >::const_iterator GoIter = LnProbGoGivenBo.begin();
      GoIter != LnProbGoGivenBo.end(); GoIter++) {

    // retreive genotype combination and corresponding log probability
    string Go = GoIter->first;
    long double lnProbGoGivenBo = GoIter->second;

    // normalize by greatest genotype probability in set
    long double lnProb =  lnProbGoGivenBo - lnProbGoGivenBoKing;

    // update total probability
    sumProb += exp(lnProb);
  }

  //----------------------------------------------------------------------------
  // 1. normalize probabilities by prob normalization factor
  // 2. calculate final P(Go|Bo) for each genotype combination
  //----------------------------------------------------------------------------

  // define final P(Go|Bo) hash
  map<string, long double > ProbGoGivenBo;

  // process every genotype combination for which probability was calculated
  for(map<string, long double >::const_iterator GoIter = LnProbGoGivenBo.begin();
      GoIter != LnProbGoGivenBo.end(); GoIter++) {

    // retrieve genotype combination and its log probability
    string Go = GoIter->first;
    long double lnProbGoGivenBo = GoIter->second;

    // calculate final P(G|B) for this genotype combination
    //   (since sumProb was calculated using values normalized by lnProbGoGivenBoKing,
    //    we need to do the same here).
    ProbGoGivenBo[Go] = exp(lnProbGoGivenBo - lnProbGoGivenBoKing) / sumProb;

    // print debugging info if required
    if (debug2) {
      // print Go, lnProbGoGivenBo and ProbGoGivenBo
      cerr << "Go=" << Go << " LnProbBoGivenGo=" << LnProbBoGivenGo[Go] << " LnProbGo[Go]=" << LnProbGo[Go] << " lnProbGoGivenBo=" << lnProbGoGivenBo << " lnProbGoGivenBo-lnProbGoGivenBoKing=" << lnProbGoGivenBo - lnProbGoGivenBoKing << " ProbGoGivenBo=" << ProbGoGivenBo[Go] << endl;
    }    
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // calculate P(SNP)
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // SNP probability
  long double pSnp = 1.0;

  // subtract probabilities of monomorphic genotype combinations
  if (diploid) {
    if (ProbGoGivenBo.count(GoAllGi11) > 0) {
      pSnp -= ProbGoGivenBo[GoAllGi11];
    }
    if (ProbGoGivenBo.count(GoAllGi22) > 0) {
      pSnp -= ProbGoGivenBo[GoAllGi22];
    }
  }
  else {
    if (ProbGoGivenBo.count(GoAllGi1) > 0) {
      pSnp -= ProbGoGivenBo[GoAllGi1];
    }
    if (ProbGoGivenBo.count(GoAllGi2) > 0) {
      pSnp -= ProbGoGivenBo[GoAllGi2];
    }
  }


  if (debug2) {
    cerr << "pSnp=" << pSnp << " ProbGoGivenBo[GoAllGi11]=" << ProbGoGivenBo[GoAllGi11] << " ProbGoGivenBo[GoAllGi22]=" << ProbGoGivenBo[GoAllGi22] << endl;
  }

  // fix SNP probability if needed
  if (pSnp < 0) {pSnp = 0;}

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Calculate individual posterior genotype probabilities
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // Define individual genotype probability distribution hash
  //----------------------------------------------------------------------------
  map<string, map<string, long double > > ProbGiGivenBo;

  //----------------------------------------------------------------------------
  // Estimate genotype probabilities for each individual
  //----------------------------------------------------------------------------
  for (unsigned int indIndex=0; indIndex<individuals.size(); indIndex++) {
	
    // Retrieve individual
    string ind = individuals[indIndex];

    /*
    cout << "indIndex=" << indIndex << " ind=" << ind << endl;
    */

    // Make new map of log genotype likelihoods for genotype probabilty estimation
    map<string, long double > LnProbBoGivenGoGenotype;

    // Register best posterior probability genotype combination in initial list
    LnProbBoGivenGoGenotype[GoKing] = lnProbBoGivenGoKing;

    /*
    // Register monomorphics in initial list
    if (diploid) {
      LnProbBoGivenGoGenotype[GoAllGi11] = lnProbBoGivenGoAllGi11;
      LnProbBoGivenGoGenotype[GoAllGi22] = lnProbBoGivenGoAllGi22;
    }
    else {
      LnProbBoGivenGoGenotype[GoAllGi1] = lnProbBoGivenGoAllGi1;
      LnProbBoGivenGoGenotype[GoAllGi2] = lnProbBoGivenGoAllGi2;
    }
    */

    //--------------------------------------------------------------------------
    // Make 1-edit distance genotype combinations and calculate data likelihoods
    //--------------------------------------------------------------------------

    // Temp log genotype likelihood map
    map<string, long double > LnProbBoGivenGoGenotypeNew;

    // Iterate through every element of current list
    for (map<string, long double >::const_iterator GoIter = LnProbBoGivenGoGenotype.begin();
	 GoIter != LnProbBoGivenGoGenotype.end(); GoIter++) {
      
      // Retreive genotype.and probability
      string Go = GoIter->first;
      long double lnProbGo = GoIter->second;

      //------------------------------------------------------------------------
      // Make every new 1-edit distance genotype combo from current list for
      //   this individual
      //------------------------------------------------------------------------

      // get individual genotype
      string Gi;
      TRY { Gi = Go.substr(indIndex, 1); } CATCH;

      // get individual log genotype probability
      long double lnProbGi = 0;

      // if data suffient, there was an actual genotype with a corresponding log probability
      if (DataSufficientGi[ind]) {
	lnProbGi = LnProbBiGivenGi[ind][Gi];
      }

      if (diploid) {
	// Make new combo with Gi=Gi11 for individual=ind and add to list if new.
	string GoReplaceWithGi11 = Go;
	GoReplaceWithGi11[indIndex] = Gi11[0]; // HACK
	LnProbBoGivenGoGenotypeNew[GoReplaceWithGi11] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi11];

	// Make new combo with Gi22 at GiIndex and add to list if new.
	string GoReplaceWithGi22 = Go;
	GoReplaceWithGi22[indIndex] = Gi22[0]; // HACK
	LnProbBoGivenGoGenotypeNew[GoReplaceWithGi22] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi22];
	    
	// Make new combo with Gi12 at GiIndex and add to list if new.
	string GoReplaceWithGi12 = Go;
	GoReplaceWithGi12[indIndex] = Gi12[0];
	LnProbBoGivenGoGenotypeNew[GoReplaceWithGi12] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi12];
      }
      else {
	// Make new combo with Gi1 at GiIndex and add to list if new.
	string GoReplaceWithGi1 = Go;
	GoReplaceWithGi1[indIndex] = Gi1[0]; // HACK
	LnProbBoGivenGoGenotypeNew[GoReplaceWithGi1] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi1];
	
	// Make new combo with Gi2 at GiIndex and add to list if new.
	string GoReplaceWithGi2 = Go;
	GoReplaceWithGi2[indIndex] = Gi2[0]; // HACK
	LnProbBoGivenGoGenotypeNew[GoReplaceWithGi2] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi2];
      }
    }

    // Update genotype likeihood map with all the new values
    LnProbBoGivenGoGenotype = LnProbBoGivenGoGenotypeNew;

    //--------------------------------------------------------------------------
    // Make new map of log genotype likelihoods for genotype probabilty estimation
    //--------------------------------------------------------------------------
    map<string, long double > LnProbGoGivenBoGenotype;

    //--------------------------------------------------------------------------
    // Estimate genotype likelihoods
    //--------------------------------------------------------------------------

    // initialize genotype probability
    if (diploid) {
      ProbGiGivenBo[ind][Gi11] = 0;
      ProbGiGivenBo[ind][Gi22] = 0;
      ProbGiGivenBo[ind][Gi12] = 0;
    }
    else {
      ProbGiGivenBo[ind][Gi1] = 0;
      ProbGiGivenBo[ind][Gi2] = 0;
    }

    // initialize Gi-specific posterior sum terms
    long double sumProbGi1GivenBo = 0;
    long double sumProbGi2GivenBo = 0;
    long double sumProbGi11GivenBo = 0;
    long double sumProbGi22GivenBo = 0;
    long double sumProbGi12GivenBo = 0;

    // iterate through every genotype combination considered
    for (map<string, long double >::const_iterator GoIter = LnProbBoGivenGoGenotype.begin();
	 GoIter != LnProbBoGivenGoGenotype.end(); GoIter++) {
    
      // retreive genotype
      string Go = GoIter->first;
    
      /*
      cout << "  Go=" << Go << " " << Go.substr(0, indIndex) << "{"<< Go[indIndex] << "}" << Go.substr(indIndex+1, Go.size() - indIndex - 1);
      */

      // initialize probability value
      long double lnProbGoGivenBo = 0;

      //------------------------------------------------------------------------
      // retrieve or calculate anew posterior term for this genotype combination
      //------------------------------------------------------------------------

      // if P(Go|Bo) already calculated, use it
      if (LnProbGoGivenBo.count(Go) > 0) {
	lnProbGoGivenBo = LnProbGoGivenBo[Go];
      }

      // otherwise calculate anew
      else {

	//----------------------------------------------------------------------
	// calculate prior
	//----------------------------------------------------------------------
	int n=0, f1=0, h=0;

	// Iterate through every individual in Go genotype combination.
	for (unsigned int indIndex = 0; indIndex < individuals.size(); ++indIndex) {
      
	  // get individual
	  string ind = individuals[indIndex];
      
	  // if data suffient, count it
	  if (DataSufficientGi[ind]) {
	    // get individual genotype
        string Gi;
	    TRY { Gi = Go.substr(indIndex, 1); } CATCH;
	    
	    // get individual log genotype probability
	    long double lnProbGi = LnProbBiGivenGi[ind][Gi];
	    
	    // update sample size, allele frequency, and number of hets
	    if (diploid) {
	      if (Gi == Gi11) {
		n += 2;
		f1 += 2;
	      }
	      else if (Gi == Gi22) {
		n += 2;
	      }
	      else if (Gi == Gi12) {
		n += 2;
		f1 += 1;
		h += 1;
	      }
	    }
	    else {
	      if (Gi == Gi1) {
		n += 1;
		f1 += 1;
	      }
	      else if (Gi == Gi2) {
		n += 1;
	      }
	    }
	  }
	}
	 
	// calculate and register total prior for this genotype combo
	long double lnPrior = 0;
	if (n == 0) {
	  // no reads
	  lnPrior = log (1 - THETA) - log(2.0);
	}
	else if (f1 == 0 || f1 == n) {
	  // monomorphic case, reads > 0
	  lnPrior = log (1 - THETA * sumOneOverI(n-1)) - log(2.0);
	}
	else {
	  // polymorphic case, reads > 0
	  lnPrior = log(THETA) + h * log(2.0) + log(n) - cofactorln(n, f1) - log(f1) - log(n-f1);
	}
	
	//--------------------------------------------------------------------
	// calculate and register raw P(Go|Bo) probability
	//--------------------------------------------------------------------
	lnProbGoGivenBo = LnProbBoGivenGoGenotype[Go] + lnPrior;
      }

      // Register genotype combo probability
      LnProbGoGivenBoGenotype[Go] = lnProbGoGivenBo;

      /*
      cout << " lnProbGoGivenBo=" << lnProbGoGivenBo << endl;
      */
    }

    //------------------------------------------------------------------------
    // determine highest genotype combo probability
    //------------------------------------------------------------------------

    // sort genotype combinations according to descending ln P(G|B)
    vector<string> GoListSorted = sortKeysByValue(LnProbGoGivenBoGenotype, true);
    
    // get genotype combination with the highest probability 
    string GoQueen = GoListSorted[0];
    
    // get the highest probability corresponding to that genotype combination
    long double lnProbBoGivenGoQueen = LnProbBoGivenGoGenotype[GoQueen];
    long double lnProbGoGivenBoQueen = LnProbGoGivenBoGenotype[GoQueen];

    /*
    cout << "  GoQueen=" << GoQueen << " lnProbBoGivenGoQueen=" << lnProbBoGivenGoQueen << " lnProbGoGivenBoQueen=" << lnProbGoGivenBoQueen << endl;
    */
    //------------------------------------------------------------------------
    // calculate final posterior genotype probabilities
    //------------------------------------------------------------------------
    for (map<string, long double >::const_iterator GoIter = LnProbGoGivenBoGenotype.begin();
	 GoIter != LnProbGoGivenBoGenotype.end(); GoIter++) {
    
      // retreive genotype and probability
      string Go = GoIter->first;
      long double lnProbGoGivenBo = GoIter->second;

      //----------------------------------------------------------------------
      // aggregate this probability with corresponding Gi for this individual
      //----------------------------------------------------------------------

      // get individual genotype
      string Gi;
      TRY { Gi = Go.substr(indIndex, 1); } CATCH;
      
      if (diploid) {
	if (Gi == Gi11) {
	  sumProbGi11GivenBo += exp(lnProbGoGivenBo - lnProbGoGivenBoQueen);
	}
	else if (Gi == Gi22) {
	  sumProbGi22GivenBo += exp(lnProbGoGivenBo - lnProbGoGivenBoQueen);
	}
	else if (Gi == Gi12) {
	  sumProbGi12GivenBo += exp(lnProbGoGivenBo - lnProbGoGivenBoQueen);
	}
      }
      else {
	if (Gi == Gi1) {
	  sumProbGi1GivenBo += exp(lnProbGoGivenBo - lnProbGoGivenBoQueen);
	}
	else if (Gi == Gi2) {
	  sumProbGi2GivenBo += exp(lnProbGoGivenBo - lnProbGoGivenBoQueen);
	}
      }

      /*
      cout << "  Go=" << Go << " " << Go.substr(0, indIndex) << "{"<< Go[indIndex] << "}" << Go.substr(indIndex+1, Go.size() - indIndex - 1);
      cout << "    lnProbGoGivenBo - lnProbGoGivenBoQueen=" << lnProbGoGivenBo - lnProbGoGivenBoQueen << " exp()=" << exp(lnProbGoGivenBo - lnProbGoGivenBoQueen) << endl;
      cout << "    sumProbGi11GivenBo=" << sumProbGi11GivenBo << endl;
      cout << "    sumProbGi22GivenBo=" << sumProbGi22GivenBo << endl;
      cout << "    sumProbGi12GivenBo=" << sumProbGi12GivenBo << endl;
      */
    }

    // calculate final genotype likelihoods
    if (diploid) {
      long double sumProb = sumProbGi11GivenBo + sumProbGi22GivenBo + sumProbGi12GivenBo;
      ProbGiGivenBo[ind][Gi11] = sumProbGi11GivenBo / sumProb;
      ProbGiGivenBo[ind][Gi22] = sumProbGi22GivenBo / sumProb;
      ProbGiGivenBo[ind][Gi12] = sumProbGi12GivenBo / sumProb;

      /*
      cout << "  subProb=" << sumProb << " ProbGiGivenBo[ind][Gi11]=" << ProbGiGivenBo[ind][Gi11] << " ProbGiGivenBo[ind][Gi22]=" << ProbGiGivenBo[ind][Gi22] << " ProbGiGivenBo[ind][Gi12]=" << ProbGiGivenBo[ind][Gi12] << endl;
      */
    }
    else {
      long double sumProb = sumProbGi1GivenBo + sumProbGi2GivenBo;
      ProbGiGivenBo[ind][Gi1] = sumProbGi1GivenBo / sumProb;
      ProbGiGivenBo[ind][Gi2] = sumProbGi2GivenBo / sumProb;
    }
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // assign output Variation data structure
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // define
  //----------------------------------------------------------------------------
  Variation var;

  //----------------------------------------------------------------------------
  // assign simple fields
  //----------------------------------------------------------------------------

  // assign ploidy
  var.diploid = diploid;

  // assign alleles
  var.alleles.push_back(allele1);
  var.alleles.push_back(allele2);

  // assign P(SNP)
  var.pSnp = pSnp;

  //----------------------------------------------------------------------------
  // assign invidual genotype log data likelihoods
  //----------------------------------------------------------------------------
  var.individualGenotypeLogDataLikelihood = LnProbBiGivenGi;

  //----------------------------------------------------------------------------
  // assign invidual genotype probabilities
  //----------------------------------------------------------------------------
  map<string, map<string, long double > > individualGenotypeProbability;
  for (map<string, map<string, long double > >::const_iterator 
	 indIter = ProbGiGivenBo.begin(); indIter != ProbGiGivenBo.end(); indIter++) {
    string ind = indIter->first;
    map<string, long double > genotypeProb = indIter->second;
    for (map<string, long double >::const_iterator 
	   GiIter = genotypeProb.begin(); GiIter != genotypeProb.end(); GiIter++) {
      string Gi = GiIter->first;
      long double p = GiIter->second;

      // recode genotype according to expected output format
      string t;
      if (diploid) {
	if (Gi == Gi11) {t = allele1 + allele1;}
	else if (Gi == Gi22) {t = allele2 + allele2;}
	else if (Gi == Gi12) {t = allele1 + allele2;}
      }
      else {
	if (Gi == Gi1) {t = allele1;}
	else if (Gi == Gi2) {t = allele2;}
      }

      // register
      individualGenotypeProbability[ind][t] = p;
    }
  }

  // assign individualGenotypeProbability
  var.individualGenotypeProbability = individualGenotypeProbability;

  //----------------------------------------------------------------------------
  // return Variation data structure
  //----------------------------------------------------------------------------
  return var;
}

//------------------------------------------------------------------------------
// posteriorProb -- genotype probs based always on at least 3 values
//------------------------------------------------------------------------------
Variation posteriorProb(
			map<string, map<string, long double > > LnProbBiGivenGi,
			vector<string> individuals,
			bool diploid,
			long double THETA,
			bool banded,
			int WB,
			int TB,
			bool includeMonoB,
			int TR,
			string allele1,
			string allele2,
			bool debug2
			) {

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Notations:
  //     Gi: genotype of a given individual (e.g. AA, AG, GG)
  //     Bi: set of base calls for the individual (e.g. A,A,A,G,A,A,G)
  //     Ai: set of true bases estimated by the base calls (e.g. A,A,G,A,G,G)
  //
  //     Go: genotype vector of all individuals
  //     Bo: base calls for all individuals
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // constants
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // constants
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  long double LOGFACTOR = log(10) / (-10.0); 

  // haploid genotypes
  string Gi1 = IUPAC[allele1];
  string Gi2 = IUPAC[allele2];

  // diploid genotypes
  string Gi11 = IUPAC[allele1+allele1];
  string Gi22 = IUPAC[allele2+allele2];
  string Gi12 = IUPAC[allele1+allele2];

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Build special genotype combinations: mono allele 1, mono allele 2.
  // Build highest data likelihood genotype combination if banded algorithm 
  // Register corresponding ln P(Bo|Go) data likelihoods
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  
  // genotype combination vector with no elements
  string GoNull;

  // special genotype combinations
  string GoAllGi1, GoAllGi2, GoAllGi11, GoAllGi22, GoMixGi;

  // P(Bo|Go) for special genotype combinations
  long double lnProbBoGivenGoAllGi1 = 0;
  long double lnProbBoGivenGoAllGi2 = 0;
  long double lnProbBoGivenGoAllGi11 = 0;
  long double lnProbBoGivenGoAllGi22 = 0;
  long double lnProbBoGivenGoMixGi = 0;
  
  // build special genotype combinations
  for (vector<string>::const_iterator indIter = individuals.begin();
       indIter != individuals.end(); indIter++) {
    
    // get individual
    string ind = *indIter;
    
    // diploid DNA
    if (diploid) {
      
      //------------------------------------------------------------------------
      // add Gi11 to all-Gi11 genotype combo and update probability
      //------------------------------------------------------------------------
      GoAllGi11 += Gi11;
      lnProbBoGivenGoAllGi11 += LnProbBiGivenGi[ind][Gi11];
      
      //------------------------------------------------------------------------
      // add Gi22 to all-Gi22 genotype combo and update probability
      //------------------------------------------------------------------------
      GoAllGi22 += Gi22;
      lnProbBoGivenGoAllGi22 += LnProbBiGivenGi[ind][Gi22];
      
      //------------------------------------------------------------------------
      // determine highest likelihood individual genotype (Gi) and add
      // to highest data likelihood genotype combo
      //------------------------------------------------------------------------

      // initialize as Gi12
      string GiBest = Gi12; long double lnProbBiGivenGiBest = LnProbBiGivenGi[ind][Gi12];
      
      // test Gi11
      if (LnProbBiGivenGi[ind][Gi11] > lnProbBiGivenGiBest) {
	GiBest = Gi11; lnProbBiGivenGiBest = LnProbBiGivenGi[ind][Gi11];
      }
      
      // test Gi22
      if (LnProbBiGivenGi[ind][Gi22] > lnProbBiGivenGiBest) {
	GiBest = Gi22; lnProbBiGivenGiBest = LnProbBiGivenGi[ind][Gi22];
      }
      
      // add best Gi to highest data likelihood genotype combo
      GoMixGi += GiBest;
      lnProbBoGivenGoMixGi += lnProbBiGivenGiBest;
    }
    
    // haploid DNA
    else {
      //------------------------------------------------------------------------
      // add Gi1 to all-Gi1 Go combo and update probability
      //------------------------------------------------------------------------
      GoAllGi1 += Gi1;
      lnProbBoGivenGoAllGi1 += LnProbBiGivenGi[ind][Gi1];
      
      //------------------------------------------------------------------------
      // add Gi2 to all-Gi2 Go combo and update probability
      //------------------------------------------------------------------------
      GoAllGi2 += Gi2;
      lnProbBoGivenGoAllGi2 += LnProbBiGivenGi[ind][Gi2];
      
      //------------------------------------------------------------------------
      // Determine highest likelihood individual genotype (Gi) and add
      //   to highest data likelihood genotype combo
      //------------------------------------------------------------------------

      // initialize as Gi11
      string GiBest = Gi1; long double lnProbBiGivenGiBest = LnProbBiGivenGi[ind][Gi1];
      
      // test Gi2
      if (LnProbBiGivenGi[ind][Gi2] > lnProbBiGivenGiBest) {
	GiBest = Gi2; lnProbBiGivenGiBest = LnProbBiGivenGi[ind][Gi2];
      }
      
      // add best Gi to highest data likelihood genotype combo
      GoMixGi += GiBest;
      lnProbBoGivenGoMixGi += lnProbBiGivenGiBest;
    }
  }
  
  if (debug2) {
    cerr << "BEST Go=" << GoMixGi << " lnProbBoGivenGoMixGi=" << lnProbBoGivenGoMixGi << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Create the highest-contributing subset of Go genotype combinations with
  //   either recursive method or banded approximation, & calculate ln P(Bo|Go).
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  
  //----------------------------------------------------------------------------
  // Define ln P(Bo|Go) for all genotype combinations.
  //----------------------------------------------------------------------------
  map<string, long double > LnProbBoGivenGo;

  //----------------------------------------------------------------------------
  // banded approximation
  //----------------------------------------------------------------------------
  if (banded) {

    //--------------------------------------------------------------------------
    // Create a list of all genotype combinations <= WB edit distance away
    //   from highest data likelihood genotype combo (GoMixGi) and from
    //   monomorphics
    //--------------------------------------------------------------------------
      
    // Register nominal "best" genotype combination in initial list
    LnProbBoGivenGo[GoMixGi] = lnProbBoGivenGoMixGi;

    // if monomorphic genotype combinations are to be included in the initial list
    if (includeMonoB) {

      // Register monomorphics in initial list
      if (diploid) {
	LnProbBoGivenGo[GoAllGi11] = lnProbBoGivenGoAllGi11;
	LnProbBoGivenGo[GoAllGi22] = lnProbBoGivenGoAllGi22;
      }
      else {
	LnProbBoGivenGo[GoAllGi1] = lnProbBoGivenGoAllGi1;
	LnProbBoGivenGo[GoAllGi2] = lnProbBoGivenGoAllGi2;
      }
    }

    // Iterate through edit distances. Limit edit distances to 
    //   the number of individuals, 
    for (unsigned int w=1; w<=min(WB, (int)individuals.size()); w++) {
      
      // make new LnProbBoGivenGo map
      map<string, long double > LnProbBoGivenGoNew;

      // Iterate through every element of current list
      for (map<string, long double >::const_iterator GoIter = LnProbBoGivenGo.begin();
	   GoIter != LnProbBoGivenGo.end(); GoIter++) {
	
	// Retreive genotype.and probability
	string Go = GoIter->first;
	long double lnProbGo = GoIter->second;

	//----------------------------------------------------------------------
	// Make every new 1-edit distance genotype combo from current list.
	//----------------------------------------------------------------------

	// Iterate through every individual in Go genotype combination.
	for (unsigned int indIndex = 0; indIndex < individuals.size(); ++indIndex) {

	  // get individual
	  string ind = individuals[indIndex];
	  
	  // get individual genotype
      string Gi;
	  TRY { Gi = Go.substr(indIndex, 1); } CATCH;

	  // get individual log genotype probability
	  long double lnProbGi = LnProbBiGivenGi[ind][Gi];

	  if (diploid) {
	    // Make new combo with Gi=Gi11 for individual=ind and add to list if new.
	    string GoReplaceWithGi11 = Go;
	    GoReplaceWithGi11[indIndex] = Gi11[0]; // HACK
	    LnProbBoGivenGoNew[GoReplaceWithGi11] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi11];

	    // Make new combo with Gi22 at GiIndex and add to list if new.
	    string GoReplaceWithGi22 = Go;
	    GoReplaceWithGi22[indIndex] = Gi22[0]; // HACK
	    LnProbBoGivenGoNew[GoReplaceWithGi22] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi22];
	    
	    // Make new combo with Gi12 at GiIndex and add to list if new.
	    string GoReplaceWithGi12 = Go;
	    GoReplaceWithGi12[indIndex] = Gi12[0];
	    LnProbBoGivenGoNew[GoReplaceWithGi12] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi12];
	  }
	  else {
	    // Make new combo with Gi1 at GiIndex and add to list if new.
	    string GoReplaceWithGi1 = Go;
	    GoReplaceWithGi1[indIndex] = Gi1[0]; // HACK
	    LnProbBoGivenGoNew[GoReplaceWithGi1] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi1];
	    
	    // Make new combo with Gi2 at GiIndex and add to list if new.
	    string GoReplaceWithGi2 = Go;
	    GoReplaceWithGi2[indIndex] = Gi2[0]; // HACK
	    LnProbBoGivenGoNew[GoReplaceWithGi2] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi2];
	  }
	}
      }

      //------------------------------------------------------------------------
      // update LnProbBoGivenGo with LnProbBoGivenGoNew (or with dominant terms
      //   thereof)
      //------------------------------------------------------------------------
      
      // If TB == 0 (means do not cap number of terms) or number of terms is
      //   lower than allowed, simply replace hashes
      if ((TB == 0) || (LnProbBoGivenGoNew.size() <= TR)) {
	LnProbBoGivenGo = LnProbBoGivenGoNew;
      }
      
      // otherwise only take the best TB entries (after including priors) 
      else {
	
	// rank genotype combinations Go according to decreasing P(Bo|Go)
	vector<string> GoListSorted = sortKeysByValue(LnProbBoGivenGoNew, true);
	
	// initialize term counter
	int t = 0;
	for(vector<string>::const_iterator GoIter = GoListSorted.begin();
	    GoIter != GoListSorted.end(); GoIter++) {
	  
	  // increment term counter
	  t++;
	  if (t > TB) {
	    break;
	  }
	  
	  // retreive combo
	  string Go = *GoIter;
	  
	  // copy element
	  LnProbBoGivenGo[Go] = LnProbBoGivenGoNew[Go];
	}
      }
    }
    
    //--------------------------------------------------------------------------
    // always add back two monomorphic genotype combinations to list
    //--------------------------------------------------------------------------    

    // diploid DNA
    if (diploid) {
      LnProbBoGivenGo[GoAllGi11] = lnProbBoGivenGoAllGi11;
      LnProbBoGivenGo[GoAllGi22] = lnProbBoGivenGoAllGi22;
    }

    // haploid DNA
    else {
      LnProbBoGivenGo[GoAllGi1] = lnProbBoGivenGoAllGi1;
      LnProbBoGivenGo[GoAllGi2] = lnProbBoGivenGoAllGi2;
    }
  }

  //----------------------------------------------------------------------------
  // recursive method
  //----------------------------------------------------------------------------
  else {
    
    // initializa probability for null genotype combination
    LnProbBoGivenGo[GoNull] = 0;
  
    // cycle through each individual
    for (vector<string>::const_iterator indIter = individuals.begin();
	 indIter != individuals.end(); indIter++) {

      // get individual
      string ind = *indIter;

      // new hash containing values updated for this individual
      map<string, long double > LnProbBoGivenGoNew;

      // iterate through all possible genotypes for this invididual
      for (map<string, long double >::const_iterator GiIter = LnProbBiGivenGi[ind].begin();
	   GiIter !=  LnProbBiGivenGi[ind].end(); GiIter++) {
	
	// retreive genotype and probability
	string Gi = GiIter->first;
	long double lnProb = GiIter->second;
	
	// iterate through all current genotype combinations
	for (map<string, long double >::const_iterator GoIter = LnProbBoGivenGo.begin();
	     GoIter != LnProbBoGivenGo.end(); GoIter++) {

	  // retreive genotype combination and corresponding prob ln value
	  string Go = GoIter->first;
	  long double LnProb = GoIter->second;
	  
	  // make new genotype combo
	  string GoNew = Go + Gi;
	  
	  // update LnProbBoGivenGo with this individual
	  LnProbBoGivenGoNew[GoNew] = LnProb + lnProb;
	}
      }
      
      //------------------------------------------------------------------------
      // update LnProbBoGivenGo with LnProbBoGivenGoNew (or with dominant terms
      //   thereof)
      //------------------------------------------------------------------------
      
      // If TR == 0 (means do not cap number of terms) or number of terms is
      // lower than allowed, simply replace hashes
      if ((TR == 0) || (LnProbBoGivenGoNew.size() <= TR)) {
	LnProbBoGivenGo = LnProbBoGivenGoNew;
      }
      
      // otherwise only take the best TR entries (after including priors) 
      else {
	
	// rank genotype combinations Go according to decreasing P(Bo|Go)
	vector<string> GoListSorted = sortKeysByValue(LnProbBoGivenGoNew, true);
	
	// initialize term counter
	int t = 0;
	for(vector<string>::const_iterator GoIter = GoListSorted.begin();
	    GoIter != GoListSorted.end(); GoIter++) {
	  
	  // increment term counter
	  t++;
	  if (t > TR) {
	    break;
	  }
	  
	  // retreive combo
	  string Go = *GoIter;
	  
	  // copy element
	  LnProbBoGivenGo[Go] = LnProbBoGivenGoNew[Go];
	}
      }
    }

    //--------------------------------------------------------------------------
    // add probabilities of two monomorphic genotype combinations
    //--------------------------------------------------------------------------

    // diploid DNA
    if (diploid) {
      LnProbBoGivenGo[GoAllGi11] = lnProbBoGivenGoAllGi11;
      LnProbBoGivenGo[GoAllGi22] = lnProbBoGivenGoAllGi22;
    }

    // haploid DNA
    else {
      LnProbBoGivenGo[GoAllGi1] = lnProbBoGivenGoAllGi1;
      LnProbBoGivenGo[GoAllGi2] = lnProbBoGivenGoAllGi2;
    }
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Calculate ln P(Go|Bo) for set of genotype combinations produced
  //   by either the banded or the recursive method
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // define ln P(G|B) and ln Prior(G) for all genotype combinations
  //----------------------------------------------------------------------------
  map<string, long double > LnProbGoGivenBo, LnProbGo;

  // iterate through every genotype combination considered
  for (map<string, long double >::const_iterator GoIter = LnProbBoGivenGo.begin();
       GoIter != LnProbBoGivenGo.end(); GoIter++) {
    
    // retreive genotype
    string Go = GoIter->first;
    
    //--------------------------------------------------------------------------
    // calculate prior
    //--------------------------------------------------------------------------
    int n=0, f1=0, h=0;

    // Iterate through every individual in Go genotype combination.
    for (unsigned int indIndex = 0; indIndex < individuals.size(); ++indIndex) {
      
      // get individual
      string ind = individuals[indIndex];
      
      // get individual genotype
      string Gi;
      TRY { Gi = Go.substr(indIndex, 1); } CATCH;
      
      // get individual log genotype probability
      long double lnProbGi = LnProbBiGivenGi[ind][Gi];
      
      // update sample size, allele frequency, and number of hets
      if (diploid) {
	if (Gi == Gi11) {
	  n += 2;
	  f1 += 2;
	}
	else if (Gi == Gi22) {
	  n += 2;
	}
	else if (Gi == Gi12) {
	  n += 2;
	  f1 += 1;
	  h += 1;
	}
      }
      else {
	if (Gi == Gi1) {
	  n += 1;
	  f1 += 1;
	}
	else if (Gi == Gi2) {
	  n += 1;
	}
      }
    }
    
    //------------------------------------------------------------------------
    // calculate and register total prior for this genotype combo
    //------------------------------------------------------------------------
    long double lnPrior = 0;
    if (n == 0) {
      // no reads
      lnPrior = log (1 - THETA) - log(2.0);
    }
    else if (f1 == 0 || f1 == n) {
      // monomorphic case, reads > 0
      lnPrior = log (1 - THETA * sumOneOverI(n-1)) - log(2.0);
    }
    else {
      // polymorphic case, reads > 0
      lnPrior = log(THETA) + h * log(2.0) + log(n) - cofactorln(n, f1) - log(f1) - log(n-f1);
    }
    LnProbGo[Go] = lnPrior;
    
    //------------------------------------------------------------------------
    // calculate and register raw P(Go|Bo) probability
    //------------------------------------------------------------------------
    LnProbGoGivenBo[Go] = LnProbBoGivenGo[Go] + lnPrior;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Calculate P(Go|Bo) values from raw ln P(Go|Bo) values
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // normalize ln P(Go|Bo) by the largest value to avoid 0 probabilities after
  // exponentiation
  //----------------------------------------------------------------------------

  // sort genotype combinations according to descending ln P(G|B)
  vector<string> GoListSorted = sortKeysByValue(LnProbGoGivenBo, true);

  // get genotype combination with the highest probability 
  string GoKing = GoListSorted[0];

  // get the highest probability corresponding to that genotype combination
  long double lnProbBoGivenGoKing = LnProbBoGivenGo[GoKing];
  long double lnProbGoGivenBoKing = LnProbGoGivenBo[GoKing];

  // 1. bring up log quantities to the largest probLn value (to avoid all 0 probs)
  // 2. calculate probability normalization factor
  long double sumProb = 0;
  for(map<string, long double >::const_iterator GoIter = LnProbGoGivenBo.begin();
      GoIter != LnProbGoGivenBo.end(); GoIter++) {

    // retreive genotype combination and corresponding log probability
    string Go = GoIter->first;
    long double lnProbGoGivenBo = GoIter->second;

    // normalize by greatest genotype probability in set
    long double lnProb =  lnProbGoGivenBo - lnProbGoGivenBoKing;

    // update total probability
    sumProb += exp(lnProb);
  }

  //----------------------------------------------------------------------------
  // 1. normalize probabilities by prob normalization factor
  // 2. calculate final P(Go|Bo) for each genotype combination
  //----------------------------------------------------------------------------

  // define final P(Go|Bo) hash
  map<string, long double > ProbGoGivenBo;

  // process every genotype combination for which probability was calculated
  for(map<string, long double >::const_iterator GoIter = LnProbGoGivenBo.begin();
      GoIter != LnProbGoGivenBo.end(); GoIter++) {

    // retrieve genotype combination and its log probability
    string Go = GoIter->first;
    long double lnProbGoGivenBo = GoIter->second;

    // calculate final P(G|B) for this genotype combination
    //   (since sumProb was calculated using values normalized by lnProbGoGivenBoKing,
    //    we need to do the same here).
    ProbGoGivenBo[Go] = exp(lnProbGoGivenBo - lnProbGoGivenBoKing) / sumProb;

    // print debugging info if required
    if (debug2) {
      // print Go, lnProbGoGivenBo and ProbGoGivenBo
      cerr << "Go=" << Go << " LnProbBoGivenGo=" << LnProbBoGivenGo[Go] << " LnProbGo[Go]=" << LnProbGo[Go] << " lnProbGoGivenBo=" << lnProbGoGivenBo << " lnProbGoGivenBo-lnProbGoGivenBoKing=" << lnProbGoGivenBo - lnProbGoGivenBoKing << " ProbGoGivenBo=" << ProbGoGivenBo[Go] << endl;
    }    
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // calculate P(SNP)
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // SNP probability
  long double pSnp = 1.0;

  // subtract probabilities of monomorphic genotype combinations
  if (diploid) {
    if (ProbGoGivenBo.count(GoAllGi11) > 0) {
      pSnp -= ProbGoGivenBo[GoAllGi11];
    }
    if (ProbGoGivenBo.count(GoAllGi22) > 0) {
      pSnp -= ProbGoGivenBo[GoAllGi22];
    }
  }
  else {
    if (ProbGoGivenBo.count(GoAllGi1) > 0) {
      pSnp -= ProbGoGivenBo[GoAllGi1];
    }
    if (ProbGoGivenBo.count(GoAllGi2) > 0) {
      pSnp -= ProbGoGivenBo[GoAllGi2];
    }
  }


  if (debug2) {
    cerr << "pSnp=" << pSnp << " ProbGoGivenBo[GoAllGi11]=" << ProbGoGivenBo[GoAllGi11] << " ProbGoGivenBo[GoAllGi22]=" << ProbGoGivenBo[GoAllGi22] << endl;
  }

  // fix SNP probability if needed
  if (pSnp < 0) {pSnp = 0;}

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // Calculate individual posterior genotype probabilities
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // Define individual genotype probability distribution hash
  //----------------------------------------------------------------------------
  map<string, map<string, long double > > ProbGiGivenBo;

  //----------------------------------------------------------------------------
  // Estimate genotype probabilities for each individual
  //----------------------------------------------------------------------------
  for (unsigned int indIndex=0; indIndex<individuals.size(); indIndex++) {
	
    // Retrieve individual
    string ind = individuals[indIndex];

    /*
    cout << "indIndex=" << indIndex << " ind=" << ind << endl;
    */

    // Make new map of log genotype likelihoods for genotype probabilty estimation
    map<string, long double > LnProbBoGivenGoGenotype;

    // Register best posterior probability genotype combination in initial list
    LnProbBoGivenGoGenotype[GoKing] = lnProbBoGivenGoKing;

    /*
    // Register monomorphics in initial list
    if (diploid) {
      LnProbBoGivenGoGenotype[GoAllGi11] = lnProbBoGivenGoAllGi11;
      LnProbBoGivenGoGenotype[GoAllGi22] = lnProbBoGivenGoAllGi22;
    }
    else {
      LnProbBoGivenGoGenotype[GoAllGi1] = lnProbBoGivenGoAllGi1;
      LnProbBoGivenGoGenotype[GoAllGi2] = lnProbBoGivenGoAllGi2;
    }
    */

    //--------------------------------------------------------------------------
    // Make 1-edit distance genotype combinations and calculate data likelihoods
    //--------------------------------------------------------------------------

    // Temp log genotype likelihood map
    map<string, long double > LnProbBoGivenGoGenotypeNew;

    // Iterate through every element of current list
    for (map<string, long double >::const_iterator GoIter = LnProbBoGivenGoGenotype.begin();
	 GoIter != LnProbBoGivenGoGenotype.end(); GoIter++) {
      
      // Retreive genotype.and probability
      string Go = GoIter->first;
      long double lnProbGo = GoIter->second;

      //------------------------------------------------------------------------
      // Make every new 1-edit distance genotype combo from current list for
      //   this individual
      //------------------------------------------------------------------------

      // get individual genotype
      string Gi;
      TRY { Gi = Go.substr(indIndex, 1); } CATCH;

      // get individual log genotype probability
      long double lnProbGi = LnProbBiGivenGi[ind][Gi];

      if (diploid) {
	// Make new combo with Gi=Gi11 for individual=ind and add to list if new.
	string GoReplaceWithGi11 = Go;
	GoReplaceWithGi11[indIndex] = Gi11[0]; // HACK
	LnProbBoGivenGoGenotypeNew[GoReplaceWithGi11] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi11];

	// Make new combo with Gi22 at GiIndex and add to list if new.
	string GoReplaceWithGi22 = Go;
	GoReplaceWithGi22[indIndex] = Gi22[0]; // HACK
	LnProbBoGivenGoGenotypeNew[GoReplaceWithGi22] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi22];
	    
	// Make new combo with Gi12 at GiIndex and add to list if new.
	string GoReplaceWithGi12 = Go;
	GoReplaceWithGi12[indIndex] = Gi12[0];
	LnProbBoGivenGoGenotypeNew[GoReplaceWithGi12] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi12];
      }
      else {
	// Make new combo with Gi1 at GiIndex and add to list if new.
	string GoReplaceWithGi1 = Go;
	GoReplaceWithGi1[indIndex] = Gi1[0]; // HACK
	LnProbBoGivenGoGenotypeNew[GoReplaceWithGi1] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi1];
	
	// Make new combo with Gi2 at GiIndex and add to list if new.
	string GoReplaceWithGi2 = Go;
	GoReplaceWithGi2[indIndex] = Gi2[0]; // HACK
	LnProbBoGivenGoGenotypeNew[GoReplaceWithGi2] = lnProbGo - lnProbGi + LnProbBiGivenGi[ind][Gi2];
      }
    }

    // Update genotype likeihood map with all the new values
    LnProbBoGivenGoGenotype = LnProbBoGivenGoGenotypeNew;

    //--------------------------------------------------------------------------
    // Make new map of log genotype likelihoods for genotype probabilty estimation
    //--------------------------------------------------------------------------
    map<string, long double > LnProbGoGivenBoGenotype;

    //--------------------------------------------------------------------------
    // Estimate genotype likelihoods
    //--------------------------------------------------------------------------

    // initialize genotype probability
    if (diploid) {
      ProbGiGivenBo[ind][Gi11] = 0;
      ProbGiGivenBo[ind][Gi22] = 0;
      ProbGiGivenBo[ind][Gi12] = 0;
    }
    else {
      ProbGiGivenBo[ind][Gi1] = 0;
      ProbGiGivenBo[ind][Gi2] = 0;
    }

    // initialize Gi-specific posterior sum terms
    long double sumProbGi1GivenBo = 0;
    long double sumProbGi2GivenBo = 0;
    long double sumProbGi11GivenBo = 0;
    long double sumProbGi22GivenBo = 0;
    long double sumProbGi12GivenBo = 0;

    // iterate through every genotype combination considered
    for (map<string, long double >::const_iterator GoIter = LnProbBoGivenGoGenotype.begin();
	 GoIter != LnProbBoGivenGoGenotype.end(); GoIter++) {
    
      // retreive genotype
      string Go = GoIter->first;
    
      /*
      cout << "  Go=" << Go << " " << Go.substr(0, indIndex) << "{"<< Go[indIndex] << "}" << Go.substr(indIndex+1, Go.size() - indIndex - 1);
      */

      // initialize probability value
      long double lnProbGoGivenBo = 0;

      //------------------------------------------------------------------------
      // retrieve or calculate anew posterior term for this genotype combination
      //------------------------------------------------------------------------

      // if P(Go|Bo) already calculated, use it
      if (LnProbGoGivenBo.count(Go) > 0) {
	lnProbGoGivenBo = LnProbGoGivenBo[Go];
      }

      // otherwise calculate anew
      else {

	//----------------------------------------------------------------------
	// calculate prior
	//----------------------------------------------------------------------
	int n=0, f1=0, h=0;

	// Iterate through every individual in Go genotype combination.
	for (unsigned int indIndex = 0; indIndex < individuals.size(); ++indIndex) {
      
	  // get individual
	  string ind = individuals[indIndex];
      
	  // get individual genotype
      string Gi;
	  TRY { Gi = Go.substr(indIndex, 1); } CATCH;
	  
	  // get individual log genotype probability
	  long double lnProbGi = LnProbBiGivenGi[ind][Gi];
	  
	  // update sample size, allele frequency, and number of hets
	  if (diploid) {
	    if (Gi == Gi11) {
	      n += 2;
	      f1 += 2;
	    }
	    else if (Gi == Gi22) {
	      n += 2;
	    }
	    else if (Gi == Gi12) {
	      n += 2;
	      f1 += 1;
	      h += 1;
	    }
	  }
	  else {
	    if (Gi == Gi1) {
	      n += 1;
	      f1 += 1;
	    }
	    else if (Gi == Gi2) {
	      n += 1;
	    }
	  }
	}
    
	// calculate and register total prior for this genotype combo
	long double lnPrior = 0;
	if (n == 0) {
	  // no reads
	  lnPrior = log (1 - THETA) - log(2.0);
	}
	else if (f1 == 0 || f1 == n) {
	  // monomorphic case, reads > 0
	  lnPrior = log (1 - THETA * sumOneOverI(n-1)) - log(2.0);
	}
	else {
	  // polymorphic case, reads > 0
	  lnPrior = log(THETA) + h * log(2.0) + log(n) - cofactorln(n, f1) - log(f1) - log(n-f1);
	}
    
	//--------------------------------------------------------------------
	// calculate and register raw P(Go|Bo) probability
	//--------------------------------------------------------------------
	lnProbGoGivenBo = LnProbBoGivenGoGenotype[Go] + lnPrior;
      }

      // Register genotype combo probability
      LnProbGoGivenBoGenotype[Go] = lnProbGoGivenBo;

      /*
      cout << " lnProbGoGivenBo=" << lnProbGoGivenBo << endl;
      */
    }

    //------------------------------------------------------------------------
    // determine highest genotype combo probability
    //------------------------------------------------------------------------

    // sort genotype combinations according to descending ln P(G|B)
    vector<string> GoListSorted = sortKeysByValue(LnProbGoGivenBoGenotype, true);
    
    // get genotype combination with the highest probability 
    string GoQueen = GoListSorted[0];
    
    // get the highest probability corresponding to that genotype combination
    long double lnProbBoGivenGoQueen = LnProbBoGivenGoGenotype[GoQueen];
    long double lnProbGoGivenBoQueen = LnProbGoGivenBoGenotype[GoQueen];

    /*
    cout << "  GoQueen=" << GoQueen << " lnProbBoGivenGoQueen=" << lnProbBoGivenGoQueen << " lnProbGoGivenBoQueen=" << lnProbGoGivenBoQueen << endl;
    */
    //------------------------------------------------------------------------
    // calculate final posterior genotype probabilities
    //------------------------------------------------------------------------
    for (map<string, long double >::const_iterator GoIter = LnProbGoGivenBoGenotype.begin();
	 GoIter != LnProbGoGivenBoGenotype.end(); GoIter++) {
    
      // retreive genotype and probability
      string Go = GoIter->first;
      long double lnProbGoGivenBo = GoIter->second;

      //----------------------------------------------------------------------
      // aggregate this probability with corresponding Gi for this individual
      //----------------------------------------------------------------------

      // get individual genotype
      string Gi;
      TRY { Gi = Go.substr(indIndex, 1); } CATCH;
      
      if (diploid) {
	if (Gi == Gi11) {
	  sumProbGi11GivenBo += exp(lnProbGoGivenBo - lnProbGoGivenBoQueen);
	}
	else if (Gi == Gi22) {
	  sumProbGi22GivenBo += exp(lnProbGoGivenBo - lnProbGoGivenBoQueen);
	}
	else if (Gi == Gi12) {
	  sumProbGi12GivenBo += exp(lnProbGoGivenBo - lnProbGoGivenBoQueen);
	}
      }
      else {
	if (Gi == Gi1) {
	  sumProbGi1GivenBo += exp(lnProbGoGivenBo - lnProbGoGivenBoQueen);
	}
	else if (Gi == Gi2) {
	  sumProbGi2GivenBo += exp(lnProbGoGivenBo - lnProbGoGivenBoQueen);
	}
      }

      /*
      cout << "  Go=" << Go << " " << Go.substr(0, indIndex) << "{"<< Go[indIndex] << "}" << Go.substr(indIndex+1, Go.size() - indIndex - 1);
      cout << "    lnProbGoGivenBo - lnProbGoGivenBoQueen=" << lnProbGoGivenBo - lnProbGoGivenBoQueen << " exp()=" << exp(lnProbGoGivenBo - lnProbGoGivenBoQueen) << endl;
      cout << "    sumProbGi11GivenBo=" << sumProbGi11GivenBo << endl;
      cout << "    sumProbGi22GivenBo=" << sumProbGi22GivenBo << endl;
      cout << "    sumProbGi12GivenBo=" << sumProbGi12GivenBo << endl;
      */
    }

    // calculate final genotype likelihoods
    if (diploid) {
      long double sumProb = sumProbGi11GivenBo + sumProbGi22GivenBo + sumProbGi12GivenBo;
      ProbGiGivenBo[ind][Gi11] = sumProbGi11GivenBo / sumProb;
      ProbGiGivenBo[ind][Gi22] = sumProbGi22GivenBo / sumProb;
      ProbGiGivenBo[ind][Gi12] = sumProbGi12GivenBo / sumProb;

      /*
      cout << "  subProb=" << sumProb << " ProbGiGivenBo[ind][Gi11]=" << ProbGiGivenBo[ind][Gi11] << " ProbGiGivenBo[ind][Gi22]=" << ProbGiGivenBo[ind][Gi22] << " ProbGiGivenBo[ind][Gi12]=" << ProbGiGivenBo[ind][Gi12] << endl;
      */
    }
    else {
      long double sumProb = sumProbGi1GivenBo + sumProbGi2GivenBo;
      ProbGiGivenBo[ind][Gi1] = sumProbGi1GivenBo / sumProb;
      ProbGiGivenBo[ind][Gi2] = sumProbGi2GivenBo / sumProb;
    }
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // assign output Variation data structure
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // define
  //----------------------------------------------------------------------------
  Variation var;

  //----------------------------------------------------------------------------
  // assign simple fields
  //----------------------------------------------------------------------------

  // assign ploidy
  var.diploid = diploid;

  // assign alleles
  var.alleles.push_back(allele1);
  var.alleles.push_back(allele2);

  // assign P(SNP)
  var.pSnp = pSnp;

  //----------------------------------------------------------------------------
  // assign invidual genotype log data likelihoods
  //----------------------------------------------------------------------------
  var.individualGenotypeLogDataLikelihood = LnProbBiGivenGi;

  //----------------------------------------------------------------------------
  // assign invidual genotype probabilities
  //----------------------------------------------------------------------------
  map<string, map<string, long double > > individualGenotypeProbability;
  for (map<string, map<string, long double > >::const_iterator 
	 indIter = ProbGiGivenBo.begin(); indIter != ProbGiGivenBo.end(); indIter++) {
    string ind = indIter->first;
    map<string, long double > genotypeProb = indIter->second;
    for (map<string, long double >::const_iterator 
	   GiIter = genotypeProb.begin(); GiIter != genotypeProb.end(); GiIter++) {
      string Gi = GiIter->first;
      long double p = GiIter->second;

      // recode genotype according to expected output format
      string t;
      if (diploid) {
	if (Gi == Gi11) {t = allele1 + allele1;}
	else if (Gi == Gi22) {t = allele2 + allele2;}
	else if (Gi == Gi12) {t = allele1 + allele2;}
      }
      else {
	if (Gi == Gi1) {t = allele1;}
	else if (Gi == Gi2) {t = allele2;}
      }

      // register
      individualGenotypeProbability[ind][t] = p;
    }
  }

  // assign individualGenotypeProbability
  var.individualGenotypeProbability = individualGenotypeProbability;

  //----------------------------------------------------------------------------
  // return Variation data structure
  //----------------------------------------------------------------------------
  return var;
}

//------------------------------------------------------------------------------
// sumOneOverI
//------------------------------------------------------------------------------
long double sumOneOverI(int n) {
  
  // initialize sum
  long double sum = 0.0;

  // return zero if no terms
  if (n < 1) {
    return sum;
  }

  // calculate sum
  for (int i=1; i<=n; i++) {
    sum += 1.0 / i;
  }

  // return sum
  return sum;
}

//------------------------------------------------------------------------------
// childGivenParentsGenotypeProbability
//------------------------------------------------------------------------------
long double childGivenParentsGenotypeProbability(string GC, string GM, string GF, long double M) {

  long double p = 0;
  if (GM == "11") {
    if (GF == "11") {
      if (GC == "11") {
	p = 1.0 - 2.0*M + M*M;
      }
      else if (GC == "12") {
	p = 2.0*M - 2.0*M*M;
      }
      else if (GC == "22") {
	p = M*M;
      }
    }
    else if (GF == "12") {
      if (GC == "11") {
	p = 0.5 - 0.5*M;
      }
      else if (GC == "12") {
	p = 0.5;
      }
      else if (GC == "22") {
	p = 0.5*M;
      }
    }
    else if (GF == "22") {
      if (GC == "11") {
	p = M - M*M;
      }
      else if (GC == "12") {
	p = 1.0 - 2.0*M + 2.0*M*M;
      }
      else if (GC == "22") {
	p = M - M*M;
      }
    }
  }
  else if (GM == "12") {
    if (GF == "11") {
      if (GC == "11") {
	p = 0.5 - 0.5*M;
      }
      else if (GC == "12") {
	p = 0.5;
      }
      else if (GC == "22") {
	p = 0.5*M;
      }
    }
    else if (GF == "12") {
      if (GC == "11") {
	p = 0.25;
      }
      else if (GC == "12") {
	p = 0.5;
      }
      else if (GC == "22") {
	p = 0.25;
      }
    }
    else if (GF == "22") {
      if (GC == "11") {
	p = 0.5*M;
      }
      else if (GC == "12") {
	p = 0.5;
      }
      else if (GC == "22") {
	p = 0.5 - 0.5*M;
      }
    }
  }
  else if (GM == "22") {
    if (GF == "11") {
      if (GC == "11") {
	p = M - M*M;
      }
      else if (GC == "12") {
	p = 1.0 - 2.0*M + 2.0*M*M;
      }
      else if (GC == "22") {
	p = M - M*M;
      }
    }
    else if (GF == "12") {
      if (GC == "11") {
	p = 0.5*M;
      }
      else if (GC == "12") {
	p = 0.5;
      }
      else if (GC == "22") {
	p = 0.5 - 0.5*M;
      }
    }
    else if (GF == "22") {
      if (GC == "11") {
	p = M*M;
      }
      else if (GC == "12") {
	p = 2.0*M - 2.0*M*M;
      }
      else if (GC == "22") {
	p = 1.0 - 2.0*M + M*M;
      }
    }
  }

  // return
  return p;
}

int phred(double prob) {
        return min(-10 * log10(1 - prob), (double) 99);
}

#endif
