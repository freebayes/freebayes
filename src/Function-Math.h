//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Function-Math
// Math methods
// Copyright 2007, 2008 Gabor T. Marth, Boston College
// All rights reserved.
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#ifndef FUNCTION_MATH_H
#define FUNCTION_MATH_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>

#include "Function-Sequence.h"

using std::ios;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::istream;
using std::cin;
using std::cout;
using std::clog;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::less;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// type definitions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

typedef map<string, string, less<string> > GenotypeCombo;
typedef map<map<string, string, less<string> >, 
	      long double, 
	      less< map<string, string, less<string> > > > GenotypeComboProbability;
typedef map<map<string, string, less<string> >, 
	      bool, 
	      less< map<string, string, less<string> > > > GenotypeComboIndicator;

//------------------------------------------------------------------------------
// Basecall: represents call base and PHRED quality value
//------------------------------------------------------------------------------
struct Basecall {
  string seqName;
  string strand;
  string base;
  short map;
  short qual;
};

//------------------------------------------------------------------------------
// Variation: data structure returned by Bayesian polymorphism probability
//            calculation
//------------------------------------------------------------------------------
struct Variation {
  // 2 polymorphic alleles 
  vector<string> alleles;

  // haploid or diploid flag
  bool diploid;

  // individual genotype log data likelihoods
  map<string, map<string, long double, less<string> >, less<string> > individualGenotypeLogDataLikelihood;

  // individual posterior genotype probabilities
  map<string, map<string, long double, less<string> >, less<string> > individualGenotypeProbability;

  // probability that site is polymorphic
  long double pSnp;
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// global variables
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// function definitions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
// gammaln -- returns ln(gamma(x))
//-------------------------------------------------------------------------------

long double gammaln(
		    long double x
		    );

//-------------------------------------------------------------------------------
// factorial -- returns n!
//-------------------------------------------------------------------------------

long double factorial(
		      int n
		      );

//-------------------------------------------------------------------------------
// factorialln -- returns ln(n!)
//-------------------------------------------------------------------------------

long double factorialln(
			int n
			);

//-------------------------------------------------------------------------------
// cofactor -- returns (n over i)
//-------------------------------------------------------------------------------

long double cofactor(
		     int n,
		     int i
		     );

//-------------------------------------------------------------------------------
// cofactorln -- returns ln((n over i))
//-------------------------------------------------------------------------------

long double cofactorln(
		       int n,
		       int i
		       );

  
//------------------------------------------------------------------------------
// sumOneOverI
//------------------------------------------------------------------------------
long double sumOneOverI(
			int
			);
  
//------------------------------------------------------------------------------
// logGenotypeLikelihoods
//------------------------------------------------------------------------------
map<string, long double, less<string> > logGenotypeLikelihoods(
							       vector<Basecall>,
							       bool,
							       long double,
							       bool
							       );

//------------------------------------------------------------------------------
// posteriorProb
//------------------------------------------------------------------------------
Variation posteriorProb(
			map<string, map<string, long double, less<string> >, less<string> >,			
			vector<string>,
			bool,
			long double,
			bool,
			int,
			int,
			bool,
			int,
			string,
			string,
			bool
			);

//------------------------------------------------------------------------------
// posteriorProb2
//------------------------------------------------------------------------------
Variation posteriorProb2(
			map<string, map<string, long double, less<string> >, less<string> >,			
			map<string, bool, less<string> >,
			vector<string>,
			bool,
			long double,
			bool,
			int,
			int,
			bool,
			int,
			string,
			string,
			bool
			);

//------------------------------------------------------------------------------
// varBayesProb
//------------------------------------------------------------------------------
Variation varBayesProb(
		       map<string, vector<Basecall>, less<string> >,
		       bool,
		       long double,
		       long double,
		       bool,
		       int,
		       int,
		       int,
		       string,
		       string,
		       bool,
		       bool
		       );

//------------------------------------------------------------------------------
// childGivenParentsGenotypeProbability
//------------------------------------------------------------------------------
long double childGivenParentsGenotypeProbability(string, string, string, long double);

//------------------------------------------------------------------------------
// phred
//------------------------------------------------------------------------------
// converts probabilities in (0,1) to phred scores
int phred(double prob);

#endif
