//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// testMap
// Math function testing code
// Copyright 2008 Gabor T. Marth, Boston College
// All rights reserved.
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// includes
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// standard includes
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>

// "tclap" commandline parsing library
#include <tclap/CmdLine.h>

// "boost" regular expression library
#include <boost/regex.hpp>

// "hash_map" true hashes
#include <ext/hash_map>

// Private libraries
#include "Function-Sequence.h"
#include "Function-Math.h"

// uses
using std::cin;
using std::cout;
using std::clog;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::map;
using namespace std; 
using namespace TCLAP; 

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// subroutines
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

int main (int argc, char *argv[]) {

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // command line
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // Create new CmdLine object
  //----------------------------------------------------------------------------
  CmdLine cmd("", ' ', "1.01" );
    
  //----------------------------------------------------------------------------
  // command line arguments
  //----------------------------------------------------------------------------

  // integer number
  ValueArg<int> cmd_N1("", "N1", "integer number 1", false, 0, "int", cmd);

  // integer number
  ValueArg<int> cmd_r1("", "r1", "integer number 2", false, 0, "int", cmd);

  // integer number
  ValueArg<int> cmd_n1("", "n1", "integer number 3", false, 0, "int", cmd);

  // integer number
  ValueArg<int> cmd_q1("", "q1", "integer number 4", false, 30, "int", cmd);

  // integer number
  ValueArg<int> cmd_N2("", "N2", "integer number 5", false, 0, "int", cmd);

  // integer number
  ValueArg<int> cmd_r2("", "r2", "integer number 6", false, 0, "int", cmd);

  // integer number
  ValueArg<int> cmd_n2("", "n2", "integer number 7", false, 0, "int", cmd);

  // integer number
  ValueArg<int> cmd_q2("", "q2", "integer number 8", false, 30, "int", cmd);

  // integer number
  ValueArg<int> cmd_N3("", "N3", "integer number 9", false, 0, "int", cmd);

  // integer number
  ValueArg<int> cmd_r3("", "r3", "integer number 10", false, 0, "int", cmd);

  // integer number
  ValueArg<int> cmd_n3("", "n3", "integer number 11", false, 0, "int", cmd);

  // integer number
  ValueArg<int> cmd_q3("", "q3", "integer number 12", false, 30, "int", cmd);

  // integer number
  ValueArg<int> cmd_N4("", "N4", "integer number 13", false, 0, "int", cmd);

  // integer number
  ValueArg<int> cmd_r4("", "r4", "integer number 14", false, 0, "int", cmd);

  // integer number
  ValueArg<int> cmd_n4("", "n4", "integer number 15", false, 0, "int", cmd);

  // integer number
  ValueArg<int> cmd_q4("", "q4", "integer number 16", false, 30, "int", cmd);

  // integer number
  ValueArg<int> cmd_WB("", "WB", "integer number 17", false, 1, "int", cmd);

  // integer number
  ValueArg<int> cmd_TB("", "TB", "integer number 18", false, 20, "int", cmd);

  // integer number
  ValueArg<int> cmd_TR("", "TR", "integer number 19", false, 20, "int", cmd);


  // diploid
  SwitchArg cmd_diploid("","diploid", "individuals diploid?", cmd, false);

  // banded
  SwitchArg cmd_banded("","banded", "use banded P(SNP) calculation?", cmd, false);

  // debug
  SwitchArg cmd_debug("","debug", "print high level debug info (technical)?", cmd, false);

  //----------------------------------------------------------------------------
  // parse command line and catch possible errors
  //----------------------------------------------------------------------------
  try {
    cmd.parse(argc,argv);
  } 
  catch ( ArgException& e ) { 
    clog << "ERROR: " << e.error() << " " << e.argId() << endl; 
  }
  
  //----------------------------------------------------------------------------
  // assign command line parameters
  //----------------------------------------------------------------------------

  // integer numbers
  int N1 = cmd_N1.getValue();
  int r1 = cmd_r1.getValue();
  int n1 = cmd_n1.getValue();
  int q1 = cmd_q1.getValue();

  int N2 = cmd_N2.getValue();
  int r2 = cmd_r2.getValue();
  int n2 = cmd_n2.getValue();
  int q2 = cmd_q2.getValue();

  int N3 = cmd_N3.getValue();
  int r3 = cmd_r3.getValue();
  int n3 = cmd_n3.getValue();
  int q3 = cmd_q3.getValue();

  int N4 = cmd_N4.getValue();
  int r4 = cmd_r4.getValue();
  int n4 = cmd_n4.getValue();
  int q4 = cmd_q4.getValue();

  int WB = cmd_WB.getValue();
  int TB = cmd_TB.getValue();
  int TR = cmd_TR.getValue();

  bool diploid = cmd_diploid.getValue();
  bool banded = cmd_banded.getValue();

  // debug flag
  bool debug = cmd_debug.getValue();

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // constants
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  
  map<string, string, less<string> > IUPAC;
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

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // set up testing data
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // individual-specific basecalls map
  map<string, vector<Basecall>, less<string> > individualBasecalls;

  // list of individuals
  vector<string> sampleList;

  //----------------------------------------------------------------------------
  // set ref and nonref alleles
  //----------------------------------------------------------------------------
  string refAllele = "A";
  string nonrefAllele = "G";

  //----------------------------------------------------------------------------
  // make base calls for reference homozygote samples
  //----------------------------------------------------------------------------
  string nameStub1 = "refHom";
  for (int i=0; i<N1; i++) {
    
    ostringstream sb;
    sb << i;
 
    string seqName = nameStub1 + sb.str();
    cerr << "Adding sequence name: " << seqName << endl;
    string strand = "+";
    short qual = q1;

    // basecall vector
    vector<Basecall> basecalls;

    // read with ref allele
    for (int j=0; j<r1; j++) {
      string base = refAllele;

      Basecall bc;
      bc.seqName = seqName;
      bc.strand = strand;
      bc.base = base;
      bc.qual = qual;
 
      basecalls.push_back(bc);
    }

    // reads with nonref allele
    for (int j=0; j<n1; j++) {
      string base = nonrefAllele;

      Basecall bc;
      bc.seqName = seqName;
      bc.strand = strand;
      bc.base = base;
      bc.qual = qual;
 
      basecalls.push_back(bc);
    }

    // register basecall vector for this individual
    individualBasecalls[seqName] = basecalls;
 
    // add individual to list
    sampleList.push_back(seqName);
  }
  
  //----------------------------------------------------------------------------
  // make base calls for heterozygote samples
  //----------------------------------------------------------------------------
  string nameStub2 = "het";
  for (int i=0; i<N2; i++) {

    ostringstream sb;
    sb << i;
 
    string seqName = nameStub2 + sb.str();
    cerr << "Adding sequence name: " << seqName << endl;
    string strand = "+";
    short qual = q2;

    // basecall vector
    vector<Basecall> basecalls;

    // read with ref allele
    for (int j=0; j<r2; j++) {
      string base = refAllele;

      Basecall bc;
      bc.seqName = seqName;
      bc.strand = strand;
      bc.base = base;
      bc.qual = qual;
 
      basecalls.push_back(bc);
    }

    // reads with nonref allele
    for (int j=0; j<n2; j++) {
      string base = nonrefAllele;

      Basecall bc;
      bc.seqName = seqName;
      bc.strand = strand;
      bc.base = base;
      bc.qual = qual;
 
      basecalls.push_back(bc);
    }

    // register basecall vector for this individual
    individualBasecalls[seqName] = basecalls;

    // add individual to list
    sampleList.push_back(seqName);
  }

  //----------------------------------------------------------------------------
  // make base calls for non-reference homozygote samples
  //----------------------------------------------------------------------------
  string nameStub3 = "nonrefHom";
  for (int i=0; i<N3; i++) {

    ostringstream sb;
    sb << i;
 
    string seqName = nameStub3 + sb.str();
    cerr << "Adding sequence name: " << seqName << endl;
    string strand = "+";
    short qual = q3;

    // basecall vector
    vector<Basecall> basecalls;

    // read with ref allele
    for (int j=0; j<r3; j++) {
      string base = refAllele;

      Basecall bc;
      bc.seqName = seqName;
      bc.strand = strand;
      bc.base = base;
      bc.qual = qual;
 
      basecalls.push_back(bc);
    }

    // reads with nonref allele
    for (int j=0; j<n3; j++) {
      string base = nonrefAllele;

      Basecall bc;
      bc.seqName = seqName;
      bc.strand = strand;
      bc.base = base;
      bc.qual = qual;
 
      basecalls.push_back(bc);
    }

    // register basecall vector for this individual
    individualBasecalls[seqName] = basecalls;

    // add individual to list
    sampleList.push_back(seqName);
  }

  //----------------------------------------------------------------------------
  // make base calls for odd samples
  //----------------------------------------------------------------------------
  string nameStub4 = "odd";
  for (int i=0; i<N4; i++) {

    ostringstream sb;
    sb << i;
 
    string seqName = nameStub4 + sb.str();
    cerr << "Adding sequence name: " << seqName << endl;
    string strand = "+";
    short qual = q4;

    // basecall vector
    vector<Basecall> basecalls;

    // read with ref allele
    for (int j=0; j<r4; j++) {
      string base = refAllele;

      Basecall bc;
      bc.seqName = seqName;
      bc.strand = strand;
      bc.base = base;
      bc.qual = qual;
 
      basecalls.push_back(bc);
    }

    // reads with nonref allele
    for (int j=0; j<n4; j++) {
      string base = nonrefAllele;

      Basecall bc;
      bc.seqName = seqName;
      bc.strand = strand;
      bc.base = base;
      bc.qual = qual;
 
      basecalls.push_back(bc);
    }

    // register basecall vector for this individual
    individualBasecalls[seqName] = basecalls;

    // add individual to list
    sampleList.push_back(seqName);
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // do Bayesian math
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  cerr << "Bayesian math" << endl;

  //----------------------------------------------------------------------------
  // calculate genotype data likelihoods
  //----------------------------------------------------------------------------

  cerr << "  Date likelihodds" << endl;
  map<string, map<string, long double, less<string> >, less<string> > LnProbBiGivenGi;
  for (map<string, vector<Basecall>, less<string> >::const_iterator indIter = individualBasecalls.begin();
       indIter != individualBasecalls.end(); indIter++) {
 
    string ind = indIter->first;
    vector<Basecall> basecalls = indIter->second;

    /*
    for (vector<Basecall>::const_iterator bcIter = basecalls.begin();
	 bcIter != basecalls.end(); bcIter++) {
      Basecall bc = *bcIter;
      cout << "  " << bc.seqName << " " << bc.base << " " << bc.qual << endl;
    }
    */

    // calculate genotype log likelihoods
    if (debug) {
      cout << ind << endl;
    }

    map<string, long double, less<string> > logGl 
      = logGenotypeLikelihoods(
			       basecalls,
			       diploid,
			       0.9,
			       debug
			       );
    
    // assign genotype log likelihoods in multi-individual data structure
    LnProbBiGivenGi[ind] = logGl;
  }	      
  
  //----------------------------------------------------------------------------
  // run posterior calculation for two best alleles
  //----------------------------------------------------------------------------
  cerr << "  Posteriors" << endl;
  Variation var = posteriorProb(
				LnProbBiGivenGi,
				sampleList,
				diploid,
				0.001,
				banded,
				WB,
				TB,
				false,
				TR,
				refAllele,
				nonrefAllele,
				false
				);

  //----------------------------------------------------------------------------
  // print output
  //----------------------------------------------------------------------------
  cerr << "Output" << endl;

  cout << "P(SNP)=" << var.pSnp << endl;
	  
  for (vector<string>::const_iterator sampleIter = sampleList.begin();
       sampleIter != sampleList.end(); sampleIter++) {
    string ind = *sampleIter;
    cout << "\t" << ind;
    for (map<string, long double, less<string> >::const_iterator 
	   gIter = var.individualGenotypeProbability[ind].begin();
	 gIter != var.individualGenotypeProbability[ind].end();
	 gIter++) {
      string g = gIter->first;
      long double p = gIter->second;
      long double l = LnProbBiGivenGi[ind][IUPAC[g]];
      cout << "\t" << g << ":\tl=" << l << "\tp=" << p;
    }
    cout << endl;
  }
}
