//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// gigaBayes
// PolyBayes for high-throughput next-generation sequences
// Copyright 2007, 2008, 2009 Gabor T. Marth, Boston College
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
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <cmath>

// "tclap" commandline parsing library
#include <tclap/CmdLine.h>

// "boost" regular expression library
#include <boost/regex.hpp>

// Bam Reader headers
#include "BamReader.h"

// "hash_map" true hashes
#include <ext/hash_map>

// private libraries
#include "Class-GigReader.h"
#include "Function-Sequence.h"
#include "Function-Generic.h"
#include "Function-Math.h"



// uses
using namespace std; 
using namespace TCLAP; 
using namespace __gnu_cxx;
using namespace BamTools;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// typedefs
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

struct ArgStruct {
  string shortId;
  string longId;
  string description;
  bool required;
  string defaultValueString;
  int defaultValueInt;
  double defaultValueDouble;
  bool defaultValueBool;
  string type;
  bool multi;
  vector<string> constraint;
};

struct AlleleData {
  string base;
  short qual;
};

struct CigarElement {
  string type;
  unsigned int length;
};

typedef vector<CigarElement> Cigar;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// template headers
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

template< typename keyType, typename valueType >
vector<keyType> sortKeysByValue(map<keyType, valueType, std::less<keyType> > hash, bool descend);

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// subroutine headers
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Prints a record into GFF format output file.
//------------------------------------------------------------------------------
bool printVar2Glt(
		  ofstream & gltFile,
		  string conName,
		  int posBegin,
		  int posEnd,
		  Variation & var
		  );

//------------------------------------------------------------------------------
// Prints a variation record into GFF format output file.
//------------------------------------------------------------------------------
bool printVar2Gff(
		  ofstream & gffFile,
		  int O,
		  string conName,
		  string ProgramName,
		  string allele1,
		  string allele2,
		  int posPadded,
		  int posBegin,
		  int posEnd,
		  Variation & var,
		  map<int, map<string, vector<Basecall>, less<string> >, less<int> > & individualBasecalls
		  );

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// static variables
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

static string ProgramName("gigaBayes");
static string ProgramDescription("Bayesian SNP and short-INDEL polymorphism discovery program.");
static string ProgramVersion("0.5.1");
static string ProgramDate("2009-02-20");
static string ProgramDeveloper("Gabor T. Marth");
static string ProgramInstitution("Boston College");
static string ProgramCopyrightDates("2007, 2008, 2009.");

static vector<ArgStruct> ArgList;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// classes
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

class MyOutput : public StdOutput {
public:
  
  virtual void failure(CmdLineInterface& c, ArgException& e)
  {
    cerr << "################################################################################" << endl;
    cerr << "### Program: " << ProgramName << " Version: " <<  ProgramVersion << " Release date: " << ProgramDate << endl;
    cerr << "### " << ProgramDescription << endl;
    cerr << "### " << "Copyright: " << ProgramCopyrightDates << " " << ProgramDeveloper << ", " << ProgramInstitution << endl;
    cerr << "### All rights reserved." << endl;
    cerr << "###" << endl;
    cerr << "### Command line error:" << e.what() << endl;
    cerr << "### For usage please type: " << c.getProgramName() << " --help" << endl;
    cerr << "################################################################################" << endl;
  }
  
  virtual void usage(CmdLineInterface& c)
  {

    cout << "################################################################################" << endl;
    cout << "### Program: " << ProgramName << " Version: " <<  ProgramVersion << " Release date: " << ProgramDate << endl;
    cout << "### " << ProgramDescription << endl;
    cout << "### " << "Copyright: " << ProgramCopyrightDates << " " << ProgramDeveloper << ", " << ProgramInstitution << endl;
    cout << "### All rights reserved." << endl;
    cout << "###" << endl;
    cout << "### Usage: " << c.getProgramName() << " [arguments], where:" << endl;
    for(vector<ArgStruct>::const_iterator it = ArgList.begin(); 
	it != ArgList.end(); it++) {
      ArgStruct arg = *it;

      string idString = "";
      if (arg.longId != "") {idString += "--" + arg.longId;}
      if (arg.shortId != "") {idString += ", -" + arg.shortId;}

      string multiString = "single-valued";
      if (arg.multi) {multiString = "multi-valued";}

      if (arg.required) {
	cout << "### " << idString << " [" << arg.type << ", required, no default, " << multiString << "]" << endl;
      }
      else {
	cout << "### " << idString << " [" << arg.type << ", optional, default=" << arg.defaultValueString << ", " << multiString << "]" << endl;
      }
      if (arg.constraint.size() > 0) {
	cout << "###     Permitted values: (";
	bool first = true;
	for (vector<string>::const_iterator iter = arg.constraint.begin();
	     iter != arg.constraint.end(); iter++) {
	  string value = *iter;
	  if (! first) {
	    cout << "|";
	  } 
	  first = false;
	  cout << value;
	}
	cout << ")" << endl;
      }
      cout << "###     Description: " << arg.description << endl;
    }
    cout << "################################################################################" << endl;
  }
  
  virtual void version(CmdLineInterface& c)
  {
    cerr << "################################################################################" << endl;
    cerr << "### Program: " << ProgramName << " Version: " <<  ProgramVersion << " Release date: " << ProgramDate << endl;
    cout << "### " << ProgramDescription << endl;
    cout << "### " << "Copyright: " << ProgramCopyrightDates << " " << ProgramDeveloper << ", " << ProgramInstitution << endl;
    cout << "### All rights reserved." << endl;
    cout << "###" << endl;
    cerr << "################################################################################" << endl;
  }
};

// ----------------------------------------------------------------------------------------------------------------------------------- //
// custom printing of BAM read info - just illustrating what you can query from BamAlignment
//
// Also... you can use this setup, with any kind of operation you want to perform on an STL container range with a for_each, transform, etc-type STL algorithm
// This output method is a bit long for what i would typically use this for, but it illustrates the point.
// ----------------------------------------------------------------------------------------------------------------------------------- //
struct PrintBamCustom {
  void operator() (const BamAlignment& bAlignment) const { 
    cout << "===================================================================" << endl;
    cout << endl;
    cout << "------------------------" << endl;
    cout << "Read Info: "              << endl;
    cout << "------------------------" << endl;
    cout << "Name:          " << bAlignment.Name         << endl;
    cout << "ReadLength:    " << bAlignment.Length       << endl;
    cout << "Query Bases:   " << bAlignment.QueryBases   << endl;
    cout << "Aligned Bases: " << bAlignment.AlignedBases << endl;
    cout << "Qualities:     " << bAlignment.Qualities    << endl;
    cout << "------------------------" << endl;
    cout << "Alignment Info: "         << endl;
    cout << "------------------------" << endl;
    cout << "RefID:         " << bAlignment.RefID      << endl;
    cout << "Position:      " << bAlignment.Position   << endl;
    cout << "Map Quality:   " << bAlignment.MapQuality << endl;
    cout << "BAM bin:       " << bAlignment.Bin        << endl;
    cout << "CIGAR:         " << bAlignment.Cigar      << endl;
    cout << "------------------------" << endl;
    cout << "Mate Info: "              << endl;
    cout << "------------------------" << endl;	
    cout << "Mate RefID:    " << bAlignment.MateRefID    << endl;
    cout << "Mate Position: " << bAlignment.MatePosition << endl;
    cout << "Insert size:   " << bAlignment.InsertSize   << endl;
    cout << "------------------------" << endl;
    cout << "Alignment Flag Queries: " << endl;
    cout << "------------------------" << endl;
    cout << "Paired?        "      << ( (bAlignment.IsPaired())            ? "true" : "false" ) << endl;
    cout << "Proper Pair?   "      << ( (bAlignment.IsProperPair())        ? "true" : "false" ) << endl;
    cout << "Mapped?        "      << ( (bAlignment.IsMapped())            ? "true" : "false" ) << endl;
    cout << "Strand?        "      << ( (bAlignment.IsReverseStrand())     ? "-"    : "+" )     << endl;
    cout << "Mate Mapped?   "      << ( (bAlignment.IsMateMapped())        ? "true" : "false" ) << endl;
    cout << "Mate Strand?   "      << ( (bAlignment.IsMateReverseStrand()) ? "-"    : "+" )     << endl;
    cout << "1st Mate?      "      << ( (bAlignment.IsFirstMate())         ? "true" : "false" ) << endl;
    cout << "2nd Mate?      "      << ( (bAlignment.IsSecondMate())        ? "true" : "false" ) << endl;
    cout << "Primary Alignment?  " << ( (bAlignment.IsPrimaryAlignment())  ? "true" : "false" ) << endl;
    cout << "Failed QC?     "      << ( (bAlignment.IsFailedQC())          ? "true" : "false" ) << endl;
    cout << "Duplicate?     "      << ( (bAlignment.IsDuplicate())         ? "true" : "false" ) << endl;
    cout << endl;
    cout << "===================================================================" << endl;
  }
};

void PrintReads(const BamAlignmentVector& reads) {
  
  cerr << "Retrieved " << reads.size() << " reads." << endl;
  cerr << "Printing all reads to STDOUT, using custom format." << endl;
  
  for_each( reads.begin(), reads.end(), PrintBamCustom() );
  
  cerr << "Finished printing." << endl;
}

Cigar parseCigarString(string cigarString, bool & status ) {

  // cigar element pattern
  boost::regex cePattern("^(\\d+)([H|S|M|P|I|D])");

  // regex match
  boost::smatch match;

  // working copy of cigar string
  string cs = cigarString;

  // return Cigar type value
  Cigar cigar;

  // initialize status as true
  status = true;

  // iteratively parse cigar string
  while (cs.length() > 0) {
    if(boost::regex_search(cs, match, cePattern)) {

      // fetch cigar element
      string ceLength = match[1];
      string ceType = match[2];

      // make cigar element
      CigarElement ce;
      ce.length = string2Int(ceLength);
      ce.type = ceType;

      // add cigar element to cigar
      cigar.push_back(ce);

      // remove this cigar element from current cigar string
      cs = cs.substr(ceLength.length() + ceType.length(), cs.length() - ceLength.length() - ceType.length());
    }
    else {

      // set status to false and make empty cigar
      status = false;
      Cigar cigar;

      // return empty cigar and false status
      return cigar;
    }
  }

  // return cigar;
  return cigar;
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

int main (int argc, char *argv[]) {

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // constants
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // command line
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // Create new CmdLine object
  //----------------------------------------------------------------------------
  CmdLine cmd("", ' ', ProgramVersion);
    
  //----------------------------------------------------------------------------
  // add program-specific command line arguments
  //----------------------------------------------------------------------------

  // initialize arg
  ArgStruct arg;

  // input file: read alignments
  ArgStruct argAif;
  arg = argAif; 
  arg.shortId = ""; 
  arg.longId = "aif"; 
  arg.description = "Read alignment input file";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_aif(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // alignment input file type
  ArgStruct argAit;
  arg = argAit;
  arg.shortId = "";
  arg.longId = "ait";
  arg.description = "Alignment input file type";
  arg.required = true;
  arg.defaultValueString = "gig";
  arg.type = "string";
  arg.multi = false;
  vector<string> allowedAit;
  allowedAit.push_back("gig");
  allowedAit.push_back("bam");
  arg.constraint = allowedAit;
  ValuesConstraint<string> allowedAitVals(allowedAit); 
  ArgList.push_back(arg);
  ValueArg<string> cmd_ait(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, &allowedAitVals, cmd);

  // output file: GFF3 text format annotations reporting SNP and INDEL candidates
  ArgStruct argGff;
  arg = argGff; 
  arg.shortId = ""; 
  arg.longId = "gff"; 
  arg.description = "Output file: GFF3 text format annotation file reporting SNP and INDEL candidates";
  arg.required = true; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_gff(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // output file: individual genotype likelihoods in text format 
  ArgStruct argGlt;
  arg = argGlt; 
  arg.shortId = ""; 
  arg.longId = "glt"; 
  arg.description = "Output file: individual genotype likelihoods in text format";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_glt(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // log file
  ArgStruct argLog;
  arg = argLog; 
  arg.shortId = ""; 
  arg.longId = "log"; 
  arg.description = "Log file";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_log(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // filterDups
  ArgStruct argFilterDups;
  arg = argFilterDups; 
  arg.shortId = ""; 
  arg.longId = "filterDups"; 
  arg.description = "Remove reads likely to come from duplicate sequence fragments?";
  arg.required = false; 
  arg.defaultValueString = "false"; 
  arg.type = "switch"; 
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_filterDups(arg.shortId, arg.longId, arg.description, cmd, false);

  // dupMode: duplicate removal mode
  ArgStruct argDupMode;
  arg = argDupMode;
  arg.shortId = "";
  arg.longId = "dupMode";
  arg.description = "Duplicate read removal mode";
  arg.required = false;
  arg.defaultValueString = "pe";
  arg.type = "string";
  arg.multi = false;
  vector<string> allowedDupMode;
  allowedDupMode.push_back("se");
  allowedDupMode.push_back("pe");
  arg.constraint = allowedDupMode;
  ValuesConstraint<string> allowedDupModeVals(allowedDupMode); 
  ArgList.push_back(arg);
  ValueArg<string> cmd_dupMode(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, &allowedDupModeVals, cmd);

  // keepOneDup
  ArgStruct argKeepOneDup;
  arg = argKeepOneDup; 
  arg.shortId = ""; 
  arg.longId = "keepOneDup"; 
  arg.description = "Keep one read from each set of duplicates?";
  arg.required = false; 
  arg.defaultValueString = "false"; 
  arg.type = "switch"; 
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_keepOneDup(arg.shortId, arg.longId, arg.description, cmd, false);

   // WRE: window of tolerance of read end position for SE duplicate read removal.
  ArgStruct argWRE;
  arg = argWRE;
  arg.shortId = "";
  arg.longId = "WRE";
  arg.description = "window of tolerance of read end position for SE duplicate read removal";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_WRE(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

  // removeMonomer: remove SNP calls made in monomeric sequence
  ArgStruct argRemoveMonomer;
  arg = argRemoveMonomer; 
  arg.shortId = ""; 
  arg.longId = "removeMonomer"; 
  arg.description = "Remove SNP calls made in monomeric sequence?";
  arg.required = false; 
  arg.defaultValueString = "false"; 
  arg.type = "switch"; 
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_removeMonomer(arg.shortId, arg.longId, arg.description, cmd, false);

// indel
  ArgStruct argIndel;
  arg = argIndel; 
  arg.shortId = ""; 
  arg.longId = "indel"; 
  arg.description = "Look for single-base INDELs?";
  arg.required = false; 
  arg.defaultValueString = "false"; 
  arg.type = "switch"; 
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_indel(arg.shortId, arg.longId, arg.description, cmd, false);

  // cgp: calculate genotype probabilities
  ArgStruct argCgp;
  arg = argCgp;
  arg.shortId = "";
  arg.longId = "cgp";
  arg.description = "Calculate individual genotype probabilities?";
  arg.required = false;
  arg.defaultValueString = "true";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_cgp(arg.shortId, arg.longId, arg.description, cmd, true);

  // region: analyze a specific region?
  ArgStruct argRegion;
  arg = argRegion; 
  arg.shortId = ""; 
  arg.longId = "region"; 
  arg.description = "Analyze a specific region?";
  arg.required = false; 
  arg.defaultValueString = "false"; 
  arg.type = "switch"; 
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_region(arg.shortId, arg.longId, arg.description, cmd, false);

  // regionSeq
  ArgStruct argRegionSeq;
  arg = argRegionSeq; 
  arg.shortId = ""; 
  arg.longId = "regionSeq"; 
  arg.description = "Name of reference sequence for region to analyze";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_regionSeq(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // regionBegin
  ArgStruct argRegionBegin;
  arg = argRegionBegin;
  arg.shortId = "";
  arg.longId = "regionBegin";
  arg.description = "Begin position of region to analyze";
  arg.required = false;
  arg.defaultValueString = "1000000000";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_regionBegin(arg.shortId, arg.longId, arg.description, arg.required, 1000000000, arg.type, cmd);

  // regionEnd
  ArgStruct argRegionEnd;
  arg = argRegionEnd;
  arg.shortId = "";
  arg.longId = "regionEnd";
  arg.description = "End position of region to analyze";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_regionEnd(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

  // useRefAllele
  ArgStruct argUseRefAllele;
  arg = argUseRefAllele;
  arg.shortId = "";
  arg.longId = "useRefAllele";
  arg.description = "Use reference sequence allele in polymorphism calling?";
  arg.required = false;
  arg.defaultValueString = "false";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_useRefAllele(arg.shortId, arg.longId, arg.description, cmd, false);

  // forceRefAllele
  ArgStruct argForceRefAllele;
  arg = argForceRefAllele;
  arg.shortId = "";
  arg.longId = "forceRefAllele";
  arg.description = "Force reference sequence allele to be always considered?";
  arg.required = false;
  arg.defaultValueString = "false";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_forceRefAllele(arg.shortId, arg.longId, arg.description, cmd, false);

  // QR: reference allele quality value
  ArgStruct argQR;
  arg = argQR;
  arg.shortId = "";
  arg.longId = "QR";
  arg.description = "Reference allele quality value (if = 0 value from BAM file used)";
  arg.required = false;
  arg.defaultValueString = "60";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_QR(arg.shortId, arg.longId, arg.description, arg.required, 60, arg.type, cmd);

  // ploidy: sample ploidy
  ArgStruct argPloidy;
  arg = argPloidy;
  arg.shortId = "";
  arg.longId = "ploidy";
  arg.description = "Sample ploidy";
  arg.required = false;
  arg.defaultValueString = "haploid";
  arg.type = "string";
  arg.multi = false;
  vector<string> allowedPloidy;
  allowedPloidy.push_back("haploid");
  allowedPloidy.push_back("diploid");
  arg.constraint = allowedPloidy;
  ValuesConstraint<string> allowedPloidyVals(allowedPloidy); 
  ArgList.push_back(arg);
  ValueArg<string> cmd_ploidy(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, &allowedPloidyVals, cmd);

  // sample: naming scheme for matching reads to samples
  ArgStruct argSample;
  arg = argSample;
  arg.shortId = "";
  arg.longId = "sample";
  arg.description = "Naming scheme for matching reads to samples";
  arg.required = false;
  arg.defaultValueString = "multiple";
  arg.type = "string";
  arg.multi = false;
  vector<string> allowedSample;
  allowedSample.push_back("single");
  allowedSample.push_back("multiple");
  allowedSample.push_back("trio");
  allowedSample.push_back("unknown");
  arg.constraint = allowedSample;
  ValuesConstraint<string> allowedSampleVals(allowedSample); 
  ArgList.push_back(arg);
  ValueArg<string> cmd_sample(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, &allowedSampleVals, cmd);

  // sampleDel
  ArgStruct argSampleDel;
  arg = argSampleDel;
  arg.shortId = "";
  arg.longId = "sampleDel";
  arg.description = "Delimiter string separating sample name and read name";
  arg.required = false;
  arg.defaultValueString = "-";
  arg.type = "string";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_sampleDel(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // MRU: maximum number of mismatches between read and reference sequence
  ArgStruct argMRU;
  arg = argMRU;
  arg.shortId = "";
  arg.longId = "MRU";
  arg.description = "Maximum number of mismatches between read and refrence sequence";
  arg.required = false;
  arg.defaultValueString = "10000000";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_MRU(arg.shortId, arg.longId, arg.description, arg.required, 10000000, arg.type, cmd);

  // count indels in mismatch count
  ArgStruct argCim;
  arg = argCim; 
  arg.shortId = ""; 
  arg.longId = "cim"; 
  arg.description = "Count indels in mismatch count?";
  arg.required = false; 
  arg.defaultValueString = "false"; 
  arg.type = "switch"; 
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_cim(arg.shortId, arg.longId, arg.description, cmd, false);

  // QML: minimum mismatch base quality value
  ArgStruct argQML;
  arg = argQML;
  arg.shortId = "";
  arg.longId = "QML";
  arg.description = "Minimum read base quality value";
  arg.required = false;
  arg.defaultValueString = "10";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_QML(arg.shortId, arg.longId, arg.description, arg.required, 10, arg.type, cmd);

  // CRL: minimum read coverage (total of two strands)
  ArgStruct argCRL;
  arg = argCRL;
  arg.shortId = "";
  arg.longId = "CRL";
  arg.description = "Minimum read coverage (total of two strands)";
  arg.required = false;
  arg.defaultValueString = "1";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_CRL(arg.shortId, arg.longId, arg.description, arg.required, 1, arg.type, cmd);

  // CRL2: minimum read coverage (on each strand)
  ArgStruct argCRL2;
  arg = argCRL2;
  arg.shortId = "";
  arg.longId = "CRL2";
  arg.description = "Minimum read coverage (on each strand)";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_CRL2(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

  // CRU: maximum read coverage (total of two stands)
  ArgStruct argCRU;
  arg = argCRU;
  arg.shortId = "";
  arg.longId = "CRU";
  arg.description = "Maximum read coverage (total of two stands)";
  arg.required = false;
  arg.defaultValueString = "10000000";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_CRU(arg.shortId, arg.longId, arg.description, arg.required, 10000000, arg.type, cmd);

  // CAL: minimum allele coverage (total of two strands)
  ArgStruct argCAL;
  arg = argCAL;
  arg.shortId = "";
  arg.longId = "CAL";
  arg.description = "Minimum allele coverage (total of two strands)";
  arg.required = false;
  arg.defaultValueString = "1";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_CAL(arg.shortId, arg.longId, arg.description, arg.required, 1, arg.type, cmd);

  // CAL2: Minimum allele coverage (on each strand)
  ArgStruct argCAL2;
  arg = argCAL2;
  arg.shortId = "";
  arg.longId = "CAL2";
  arg.description = "Minimum allele coverage (on each strand)";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_CAL2(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

  // QRL: minimum read base quality value
  ArgStruct argQRL;
  arg = argQRL;
  arg.shortId = "";
  arg.longId = "QRL";
  arg.description = "Minimum read base quality value";
  arg.required = false;
  arg.defaultValueString = "10";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_QRL(arg.shortId, arg.longId, arg.description, arg.required, 10, arg.type, cmd);

  // QAL: minimum aggregate overall allele quality value
  ArgStruct argQAL;
  arg = argQAL;
  arg.shortId = "";
  arg.longId = "QAL";
  arg.description = "Minimum aggregate allele quality value";
  arg.required = false;
  arg.defaultValueString = "40";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_QAL(arg.shortId, arg.longId, arg.description, arg.required, 40, arg.type, cmd);

  // D: pairwise nucleotide diversity (theta)
  ArgStruct argD;
  arg = argD;
  arg.shortId = "";
  arg.longId = "D";
  arg.description = "Pairwise nucleotide diversity (theta)";
  arg.required = false;
  arg.defaultValueString = "10E-3";
  arg.type = "long double";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<long double> cmd_D(arg.shortId, arg.longId, arg.description, arg.required, 1E-3, arg.type, cmd);

  // M: per-generation somatic mutation rate (mu)
  ArgStruct argM;
  arg = argM;
  arg.shortId = "";
  arg.longId = "M";
  arg.description = "Per-generation somatic mutation rate (mu)";
  arg.required = false;
  arg.defaultValueString = "1E-8";
  arg.type = "long double";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<long double> cmd_M(arg.shortId, arg.longId, arg.description, arg.required, 1E-8, arg.type, cmd);

  // algorithm
  ArgStruct argAlgorithm;
  arg = argAlgorithm;
  arg.shortId = "";
  arg.longId = "algorithm";
  arg.description = "P(SNP) calculation algorithm";
  arg.required = false;
  arg.defaultValueString = "banded";
  arg.type = "string";
  arg.multi = false;
  vector<string> allowedAlgorithm;
  allowedAlgorithm.push_back("banded");
  allowedAlgorithm.push_back("recursive");
  arg.constraint = allowedAlgorithm;
  ValuesConstraint<string> allowedAlgorithmVals(allowedAlgorithm); 
  ArgList.push_back(arg);
  ValueArg<string> cmd_algorithm(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, &allowedAlgorithmVals, cmd);

  // P: P(SNP) lower probability threshold for reporting polymorphism candidate
  ArgStruct argPSL;
  arg = argPSL;
  arg.shortId = "";
  arg.longId = "PSL";
  arg.description = "P(SNP) lower probability threshold for reporting polymorphism candidate";
  arg.required = false;
  arg.defaultValueString = "0.5";
  arg.type = "long double";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<long double> cmd_PSL(arg.shortId, arg.longId, arg.description, arg.required, 0.5, arg.type, cmd);

  // WB: bandwidth for banded P(SNP) calculation algorithm (genotype combination level)
  ArgStruct argWB;
  arg = argWB;
  arg.shortId = "";
  arg.longId = "WB";
  arg.description = "Bandwidth for banded P(SNP) calculation algorithm (genotype combination level)";
  arg.required = false;
  arg.defaultValueString = "2";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_WB(arg.shortId, arg.longId, arg.description, arg.required, 2, arg.type, cmd);

  // TB: number of terms to keep in banded P(SNP) calculation algorithm
  ArgStruct argTB;
  arg = argTB;
  arg.shortId = "";
  arg.longId = "TB";
  arg.description = "Number of terms to keep per cycle in banded P(SNP) calculation algorithm (0: keep all)";
  arg.required = false;
  arg.defaultValueString = "10";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_TB(arg.shortId, arg.longId, arg.description, arg.required, 10, arg.type, cmd);

  // TR: number of terms to keep in recursive P(SNP) calculation algorithm
  ArgStruct argTR;
  arg = argTR;
  arg.shortId = "";
  arg.longId = "TR";
  arg.description = "Number of terms to keep per cycle of recursive P(SNP) calculation algorithm (0: keep all)";
  arg.required = false;
  arg.defaultValueString = "10";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_TR(arg.shortId, arg.longId, arg.description, arg.required, 10, arg.type, cmd);

  // O: Output detail level for candidate polymorphism report
  ArgStruct argO;
  arg = argO;
  arg.shortId = "";
  arg.longId = "O";
  arg.description = "Output detail level for candidate polymorphism report";
  arg.required = false;
  arg.defaultValueString = "1";
  arg.type = "int";
  arg.multi = false;
  vector<string> allowedO;
  allowedO.push_back("0");
  allowedO.push_back("1");
  allowedO.push_back("2");
  allowedO.push_back("3");
  allowedO.push_back("4");
  arg.constraint = allowedO;
  vector<int> allowedOint;
  allowedOint.push_back(0);
  allowedOint.push_back(1);
  allowedOint.push_back(2);
  allowedOint.push_back(3);
  allowedOint.push_back(4);
  ValuesConstraint<int> allowedOVals(allowedOint); 
  ArgList.push_back(arg);
  ValueArg<int> cmd_O(arg.shortId, arg.longId, arg.description, arg.required, 1, &allowedOVals, cmd);

  // R: batch size for reads to be read from BAM file
  ArgStruct argR;
  arg = argR;
  arg.shortId = "";
  arg.longId = "R";
  arg.description = "Batch size for reads to be read from BAM file";
  arg.required = false;
  arg.defaultValueString = "10000";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_R(arg.shortId, arg.longId, arg.description, arg.required, 10000, arg.type, cmd);

  // I: report interval for debugging printouts
  ArgStruct argI;
  arg = argI;
  arg.shortId = "";
  arg.longId = "I";
  arg.description = "Report interval for debugging printouts";
  arg.required = false;
  arg.defaultValueString = "10000";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_I(arg.shortId, arg.longId, arg.description, arg.required, 10000, arg.type, cmd);

  // debug
  ArgStruct argDebug;
  arg = argDebug;
  arg.shortId = "";
  arg.longId = "debug";
  arg.description = "Print debugging messages?";
  arg.required = false;
  arg.defaultValueString = "false";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_debug(arg.shortId, arg.longId, arg.description, cmd, false);

  // debug2
  ArgStruct argDebug2;
  arg = argDebug2;
  arg.shortId = "";
  arg.longId = "debug2";
  arg.description = "Print very detailed debugging messages?";
  arg.required = false;
  arg.defaultValueString = "false";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_debug2(arg.shortId, arg.longId, arg.description, cmd, false);

  //----------------------------------------------------------------------------
  // register (but not add to cmd) special arguments that are automatically
  // added to cmd
  //----------------------------------------------------------------------------

  // help
  ArgStruct argHelp;
  arg = argHelp;
  arg.shortId = "h";
  arg.longId = "help";
  arg.description = "Print usage statement?";
  arg.required = false;
  arg.defaultValueString = "false";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);

  // version
  ArgStruct argVersion;
  arg = argVersion;
  arg.shortId = "";
  arg.longId = "version";
  arg.description = "Print program version?";
  arg.required = false;
  arg.defaultValueString = "false";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);

  //----------------------------------------------------------------------------
  // add custom output handler
  //----------------------------------------------------------------------------
  MyOutput my;
  cmd.setOutput(&my);

  //----------------------------------------------------------------------------
  // parse command line and catch possible errors
  //----------------------------------------------------------------------------
  try {
    cmd.parse(argc,argv);
  } 
  catch ( ArgException& e ) { 
    cerr << "ERROR: " << e.error() << " " << e.argId() << endl; 
  }
  
  //----------------------------------------------------------------------------
  // assign command line parameters
  //----------------------------------------------------------------------------

  string aif = cmd_aif.getValue();
  string ait = cmd_ait.getValue();
  string gff = cmd_gff.getValue();
  string glt = cmd_glt.getValue();
  string log = cmd_log.getValue();
  bool filterDups = cmd_filterDups.getValue();
  string dupMode = cmd_dupMode.getValue();
  bool keepOneDup = cmd_keepOneDup.getValue();
  int WRE = cmd_WRE.getValue();
  bool removeMonomer = cmd_removeMonomer.getValue();
  bool indel = cmd_indel.getValue();
  bool cgp = cmd_cgp.getValue();
  bool region = cmd_region.getValue();
  string regionSeq = cmd_regionSeq.getValue();
  int regionBegin = cmd_regionBegin.getValue();
  int regionEnd = cmd_regionEnd.getValue();
  bool useRefAllele = cmd_useRefAllele.getValue();
  bool forceRefAllele = cmd_forceRefAllele.getValue();
  int QR = cmd_QR.getValue();
  string ploidy = cmd_ploidy.getValue();
  string sample = cmd_sample.getValue();
  string sampleDel = cmd_sampleDel.getValue();
  long double D = cmd_D.getValue();
  long double M = cmd_M.getValue();
  int QRL = cmd_QRL.getValue();
  int QAL = cmd_QAL.getValue();
  int MRU = cmd_MRU.getValue();
  bool cim = cmd_cim.getValue();
  int QML = cmd_QML.getValue();
  int CRL = cmd_CRL.getValue();
  int CRL2 = cmd_CRL2.getValue();
  int CRU = cmd_CRU.getValue();
  int CAL = cmd_CAL.getValue();
  int CAL2 = cmd_CAL2.getValue();
  long double PSL = cmd_PSL.getValue();
  string algorithm = cmd_algorithm.getValue();
  int WB = cmd_WB.getValue();
  int TB = cmd_TB.getValue();
  int TR = cmd_TR.getValue();
  int O = cmd_O.getValue();
  int R = cmd_R.getValue();
  int I = cmd_I.getValue();
  bool debug = cmd_debug.getValue();
  bool debug2 = cmd_debug2.getValue();

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // check and fix command line options
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  if (QRL < 0) {QRL = 0;}

  bool setQR = false;
  if (QR > 0) {setQR = true;}

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // derived variables
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  bool useGig = false;
  bool useBam = false;
  if (ait == "gig") {
    useGig = true;
  }
  else if (ait == "bam") {
    useBam = true;
  }
  else {
    cerr << "Alignment format: " << ait << " invalid. Terminating..." << endl;
    exit(1);
  }

  bool diploid = false;
  if (ploidy == "diploid") {
    diploid = true;
  }

  bool banded = false;
  if (algorithm == "banded") {
    banded = true;
  }

  bool record = false;
  if (log != "") {
    record = true;
  }

  bool reportGlt = false;
  if (glt != "") {
    reportGlt = true;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // open log file if required
  //----------------------------------------------------------------------------
  //---------------------------------------------------------------------------- 
  ofstream logFile(log.c_str(), ios::out);
  if (record) {
    // report
    if (debug) {cerr << "Opening log file: " << log << " ...";}

    // open
    if (!logFile) {
      cerr << " unable to open file: " << log << endl;
      exit(1);
    }
    if (debug) {cerr << " done." << endl;}
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // report command line and parameters used
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  map<bool, string, less<bool> > bool2String;
  bool2String[false] = "false";
  bool2String[true] = "true";

  if (record) {
    logFile << "Command line:";
    for (int i=0; i<argc; i++) {
      logFile << " " << argv[i];
    }
    logFile << endl;
    logFile << endl;
    logFile << "Complete list of parameter values:" << endl;
    logFile << "  --aif = " << aif << endl;
    logFile << "  --ait = " << ait << endl;
    logFile << "  --gff = " << gff << endl;
    logFile << "  --glt = " << glt << endl;
    logFile << "  --log = " << log << endl;
    logFile << "  --filterDups = " << bool2String[filterDups] << endl;
    logFile << "  --dupMode = " << dupMode << endl;
    logFile << "  --keepOneDup = " << bool2String[keepOneDup] << endl;
    logFile << "  --WRE = " << WRE << endl;
    logFile << "  --removeMonomer = " << bool2String[removeMonomer] << endl;
    logFile << "  --indel = " << bool2String[indel] << endl;
    logFile << "  --cgp = " << bool2String[cgp] << endl;
    logFile << "  --region = " << bool2String[region] << endl;
    logFile << "  --regionSeq = " << regionSeq << endl;
    logFile << "  --regionBegin = " << regionBegin << endl;
    logFile << "  --regionEnd = " << regionEnd << endl;
    logFile << "  --useRefAllele = " <<  bool2String[useRefAllele] << endl;
    logFile << "  --forceRefAllele = " <<  bool2String[forceRefAllele] << endl;
    logFile << "  --QR = " << QR << endl;
    logFile << "  --ploidy = " << ploidy << endl;
    logFile << "  --sample = " << sample << endl;
    logFile << "  --sampleDel = " << sampleDel << endl;
    logFile << "  --D = " << D << endl;
    logFile << "  --M = " << M << endl;
    logFile << "  --QRL = " << QRL << endl;
    logFile << "  --QAL = " << QAL << endl;
    logFile << "  --MRU = " << MRU << endl;
    logFile << "  --QML = " << QML << endl;
    logFile << "  --CRL = " << CRL << endl;
    logFile << "  --CRL2 = " << CRL2 << endl;
    logFile << "  --CRU = " << CRU << endl;
    logFile << "  --CAL = " << CAL << endl;
    logFile << "  --CAL2 = " << CAL2 << endl;
    logFile << "  --PSL = " << PSL << endl;
    logFile << "  --algorithm = " << algorithm << endl;
    logFile << "  --WB = " << WB << endl;
    logFile << "  --TR = " << TR << endl;
    logFile << "  --O = " << O << endl;
    logFile << "  --R = " << R << endl;
    logFile << "  --I = " << I << endl;
    logFile << "  --debug = " <<  bool2String[debug] << endl;
    logFile << "  --debug2 = " <<  bool2String[debug2] << endl;
    logFile << endl;
  }
  if (debug) {
    cerr << "Command line:";
    for (int i=0; i<argc; i++) {
      cerr << " " << argv[i];
    }
    cerr << endl;
    cerr << endl;
    cerr << "Complete list of parameter values:" << endl;
    cerr << "  --aif = " << aif << endl;
    cerr << "  --ait = " << ait << endl;
    cerr << "  --gff = " << gff << endl;
    cerr << "  --glt = " << glt << endl;
    cerr << "  --log = " << log << endl;
    cerr << "  --filterDups = " << bool2String[filterDups] << endl;
    cerr << "  --dupMode = " << dupMode << endl;
    cerr << "  --keepOneDup = " << bool2String[keepOneDup] << endl;
    cerr << "  --WRE = " << WRE << endl;
    cerr << "  --removeMonomer = " << bool2String[removeMonomer] << endl;
    cerr << "  --indel = " << bool2String[indel] << endl;
    cerr << "  --cgp = " << bool2String[cgp] << endl;
    cerr << "  --region = " << bool2String[region] << endl;
    cerr << "  --regionSeq = " << regionSeq << endl;
    cerr << "  --regionBegin = " << regionBegin << endl;
    cerr << "  --regionEnd = " << regionEnd << endl;
    cerr << "  --useRefAllele = " <<  bool2String[useRefAllele] << endl;
    cerr << "  --forceRefAllele = " <<  bool2String[forceRefAllele] << endl;
    cerr << "  --QR = " << QR << endl;
    cerr << "  --ploidy = " << ploidy << endl;
    cerr << "  --sample = " << sample << endl;
    cerr << "  --sampleDel = " << sampleDel << endl;
    cerr << "  --D = " << D << endl;
    cerr << "  --M = " << M << endl;
    cerr << "  --QRL = " << QRL << endl;
    cerr << "  --QAL = " << QAL << endl;
    cerr << "  --MRU = " << MRU << endl;
    cerr << "  --QML = " << QML << endl;
    cerr << "  --CRL = " << CRL << endl;
    cerr << "  --CRL2 = " << CRL2 << endl;
    cerr << "  --CRU = " << CRU << endl;
    cerr << "  --CAL = " << CAL << endl;
    cerr << "  --CAL2 = " << CAL2 << endl;
    cerr << "  --PSL = " << PSL << endl;
    cerr << "  --algorithm = " << algorithm << endl;
    cerr << "  --WB = " << WB << endl;
    cerr << "  --TR = " << TR << endl;
    cerr << "  --O = " << O << endl;
    cerr << "  --R = " << R << endl;
    cerr << "  --I = " << I << endl;
    cerr << "  --debug = " <<  bool2String[debug] << endl;
    cerr << "  --debug2 = " <<  bool2String[debug2] << endl;
    cerr << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // open input file(s)
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // files
  //----------------------------------------------------------------------------

  GigReader * gigReader;
  BamReader bReader;

  //----------------------------------------------------------------------------
  // binary GIG format assembly input file
  //----------------------------------------------------------------------------

  if (useGig) {
    //--------------------------------------------------------------------------
    // open gig file
    //--------------------------------------------------------------------------
    // report
    if (record) {logFile << "Opening GIG format alignment input file: " << aif << " ...";}
    if (debug) {cerr << "Opening GIG format alignment input file: " << aif << " ...";}
    
    // open
    FILE * gigFile = fopen64(aif.c_str(), "r");
    if (!gigFile) {
      if (record) {logFile << " unable to open file: " << aif << endl;}
      cerr << " unable to open file: " << aif << endl;
      exit(1);
    }
    if (record) {logFile << " done." << endl;}
    if (debug) {cerr << " done." << endl;}
    
    //--------------------------------------------------------------------------
    // create BamReader object
    //--------------------------------------------------------------------------
    if (record) {logFile << "Making GigReader object...";}
    if (debug) {cerr << "Making GigReader object...";}
    gigReader = new GigReader(gigFile);
    
    if (record) {logFile << " done." << endl;}
    if (debug) {cerr << " done." << endl;}
 }

  //----------------------------------------------------------------------------
  // BAM format alignment input file
  //----------------------------------------------------------------------------
  else if (useBam) {
    // report
    if (record) {logFile << "Opening BAM fomat alignment input file: " << aif << " ...";}
    if (debug) {cerr << "Opening BAM foremat alignment input file: " << aif << " ...";}
    
    bReader.Open(aif.c_str());
    
    if (record) {logFile << " done." << endl;}
    if (debug) {cerr << " done." << endl;}
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // open output file(s)
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // open GFF3 format output file
  //----------------------------------------------------------------------------

  // report
  if (record) {logFile << "opening GFF3 format output file for writing: " << gff << "...";}
  if (debug) {cerr << "opening GFF3 format output file for writing: " << gff << "...";}

  // open
  ofstream gffFile(gff.c_str(), ios::out);
  if (!gffFile) {
    if (record) {logFile << " unable to open file: " << gff << endl;}
    cerr << " unable to open file: " << gff << endl;
    exit(1);
  }
  if (record) {logFile << " done." << endl;}
  if (debug) {cerr << " done." << endl;}

  //----------------------------------------------------------------------------
  // write header information
  //----------------------------------------------------------------------------
  gffFile << "##gff-version3" << endl;
  gffFile << "#" << endl;
  gffFile << "# Command line that generated this output:";
  for (int i=0; i<argc; i++) {
    gffFile << " " << argv[i];
  }
  gffFile << endl;
  gffFile << "#" << endl;
  gffFile << "# Complete list of parameter values:" << endl;
  gffFile << "#   --aif = " << aif << endl;
  gffFile << "#   --ait = " << ait << endl;
  gffFile << "#   --gff = " << gff << endl;
  gffFile << "#   --glt = " << glt << endl;
  gffFile << "#   --log = " << log << endl;
  gffFile << "#   --filterDups = " << bool2String[filterDups] << endl;
  gffFile << "#   --dupMode = " << dupMode << endl;
  gffFile << "#   --keepOneDup = " << bool2String[keepOneDup] << endl;
  gffFile << "#   --WRE = " << WRE << endl;
  gffFile << "#   --removeMonomer = " << bool2String[removeMonomer] << endl;
  gffFile << "#   --indel = " << bool2String[indel] << endl;
  gffFile << "#   --cgp = " << bool2String[cgp] << endl;
  gffFile << "#   --region = " << bool2String[region] << endl;
  gffFile << "#   --regionSeq = " << regionSeq << endl;
  gffFile << "#   --regionBegin = " << regionBegin << endl;
  gffFile << "#   --regionEnd = " << regionEnd << endl;
  gffFile << "#   --useRefAllele = " <<  bool2String[useRefAllele] << endl;
  gffFile << "#   --forceRefAllele = " <<  bool2String[forceRefAllele] << endl;
  gffFile << "#   --QR = " << QR << endl;
  gffFile << "#   --ploidy = " << ploidy << endl;
  gffFile << "#   --sample = " << sample << endl;
  gffFile << "#   --sampleDel = " << sampleDel << endl;
  gffFile << "#   --D = " << D << endl;
  gffFile << "#   --M = " << M << endl;
  gffFile << "#   --QRL = " << QRL << endl;
  gffFile << "#   --QAL = " << QAL << endl;
  gffFile << "#   --MRU = " << MRU << endl;
  gffFile << "#   --QML = " << QML << endl;
  gffFile << "#   --CRL = " << CRL << endl;
  gffFile << "#   --CRL2 = " << CRL2 << endl;
  gffFile << "#   --CRU = " << CRU << endl;
  gffFile << "#   --CAL = " << CAL << endl;
  gffFile << "#   --CAL2 = " << CAL2 << endl;
  gffFile << "#   --PSL = " << PSL << endl;
  gffFile << "#   --algorithm = " << algorithm << endl;
  gffFile << "#   --WB = " << WB << endl;
  gffFile << "#   --TR = " << TR << endl;
  gffFile << "#   --O = " << O << endl;
  gffFile << "#   --R = " << R << endl;
  gffFile << "#   --I = " << I << endl;
  gffFile << "#   --debug = " <<  bool2String[debug] << endl;
  gffFile << "#   --debug2 = " <<  bool2String[debug2] << endl;
  gffFile << "#" << endl;

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // open GLT file if required
  //----------------------------------------------------------------------------
  //---------------------------------------------------------------------------- 
  ofstream gltFile(glt.c_str(), ios::out);
  if (reportGlt) {
    // report
    if (debug) {cerr << "Opening GLT file: " << glt << " ...";}

    // open
    if (!gltFile) {
      cerr << " unable to open file: " << glt << endl;
      exit(1);
    }
    if (debug) {cerr << " done." << endl;}
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // BAM file input branch
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  if (useBam) {
    // -------------------------------------------------------------------------
    // Get the names of all the reference sequences in the BAM file
    // -------------------------------------------------------------------------
    cerr << "Number of ref seqs: " << bReader.GetReferenceCount() << endl;
    
    cerr << endl;
    cerr << "Ref sequence names:" << endl;
    
    vector<string> refNames = bReader.GetAllReferenceNames();
    
    string refNameLast;
    if (!refNames.empty()) {
      vector<string>::iterator refIter = refNames.begin();
      vector<string>::iterator refEnd  = refNames.end();
      for( ; refIter != refEnd; ++refIter) {
	string refName = *refIter;
	refNameLast = refName;
	int refId = bReader.GetRefID(refNames, refName);
	cerr << refName << " " << refId << endl;
      }
    }
    
    // -------------------------------------------------------------------------
    // get all reads from a region of the reference sequence
    // -------------------------------------------------------------------------
    string region1 = regionSeq + ":10000000-10000100";
    BamAlignmentVector alignments;
    //  if ( bReader.GetReadsFromRegion(alignments, region1.c_str()) ) { 
    //  if (bReader.GetReadsFromReference(alignments, bReader.GetRefID(refNames, "1"), 3)) {
    if (bReader.GetReads(alignments, 3)) {
      map<int, string, less<int> > pileDna, pileQual;
      int firstSp = 1000000000;
      int lastSp = 0;
      for (BamAlignmentVector::const_iterator aIter = alignments.begin();
	   aIter != alignments.end(); aIter++) {
	
	// retrieve alignment
	BamAlignment a = *aIter;
	
	// report
	//     cout << a.Name << "\t" << a.Position << "\t" << a.Length << "\t" << a.Cigar << "\t" << a.QueryBases << "\t" << a.Qualities << endl;
	cout << a.Name << "\t" << a.RefID << "\t" << a.Position << "\t" << ( (a.IsReverseStrand())     ? "-"    : "+" )  << "\t" << a.Length << "\t" << a.Cigar << endl;
	
	// parse Cigar
	bool status;
	Cigar cigar = parseCigarString(a.Cigar, status);
	if (status) {
	  
	  // retrieve read DNA and QUAL sequence
	  string rDna = a.QueryBases;
	  string rQual = a.Qualities;
	  
	  // initialize reference sequence position and read position
	  int sp = a.Position;
	  int rp = 1;
	  
	  // update refseq min pos
	  if (sp < firstSp) {
	    firstSp = sp;
	  }
	  
	  // process aligned read according to cigar
	  for (Cigar::const_iterator ceIter = cigar.begin();
	       ceIter != cigar.end(); ceIter++) {
	    CigarElement ce = *ceIter;
	    cout << ce.length << ce.type;
	    
	    // soft clip: increment read position by number of clipped bases
	    //            leave refseq positon unchanged
	    if (ce.type == "S") {
	      rp += ce.length;
	    }
	    
	    // hard clip: leave both positions unchanged, bases already clipped
	    else if (ce.type == "H") {
	    }
	    
	    // match or mismatch: increment both positions in a linear fashion
	    //                    and register aligned bases in pile
	    else if (ce.type == "M") {
	      for (int i=1; i<=ce.length; i++) {
		
		// extract aligned base
		string b = rDna.substr(rp-1, 1);
		string q = rQual.substr(rp-1, 1);
		
		// add aligned base to pile
		pileDna[sp] += b;
		pileQual[sp] += q;
		
		// update positions
		sp++;
		rp++;
	      }	    
	    }
	    
	    // deletion or skipped region: skip ahead in refseq position
	    //                             leave read position unchanged
	    else if (ce.type == "D" || ce.type == "N") {
	      //	    sp += ce.length;
	      for (int i=1; i<=ce.length; i++) {
		
		// assigned aligned base as gap
		string b = "-";
		string q = "!";
		
		// add gap base to pile
		pileDna[sp] += b;
		pileQual[sp] += q;
		
		// update refseq position
		sp++;
	      }
	    }
	    
	    // insertion: skip ahead in read position
	    //            leave refseq position unchanged
	    else if (ce.type == "I") {
	      rp += ce.length;	    
	    }
	    
	    // padding: do nothing
	    else if (ce.type == "P") {
	    }
	    
	    // complain and bomb
	    else {
	      cerr << "Invalid CIGAR: " <<  a.Cigar << "... terminating." << endl;
	      exit(1);
	    }
	    
	    // update refseq maxpos
	    if (sp > lastSp) {
	      lastSp = sp;
	    }
	  }
	  cout << endl;
	}
      }
      
      // print alignment
      for (int sp = firstSp; sp <= lastSp; sp++) {
	string pd = pileDna[sp];
	string pq = pileQual[sp];
	cout << sp << "\t" << pileDna[sp] << "\t" << pileQual[sp] << endl;
	cout << sp << "\t";
	for (int i=0; i<pd.length(); i++) {
	  char d = pd[i];
	  int q = static_cast<int>(pq[i]) - 33;
	  cout << d << "(" << q << ") ";
	}
	cout << endl;
      }
    } 
    else { 
      cerr << "Could not get reads... terminating." << endl;
      exit(1);
    }
  }


  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // GIG file branch
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  if (useGig) {

    //--------------------------------------------------------------------------
    // read assembly-wide data
    //--------------------------------------------------------------------------
    AssemblyData assemblyData = gigReader->getAssemblyData();
    
    // retreive data
    string assName = assemblyData.name;
    int numCons = assemblyData.numContigs;
    int assNumReads = assemblyData.numReads;
    
    //--------------------------------------------------------------------------
    // process contig-wide data
    //--------------------------------------------------------------------------

    // report
    if (record) {logFile << "Reading " << numCons << " reference sequences." << endl;}
    if (debug) {cerr << "Reading " << numCons << " reference sequences." << endl;}
    
    //--------------------------------------------------------------------------
    // iterate through every contig in index
    //--------------------------------------------------------------------------
    int conCount = 0;
    for (vector<string>::const_iterator contigIter = assemblyData.contigNames.begin();
	 contigIter != assemblyData.contigNames.end(); contigIter++) {
      string conName = * contigIter;
      
      //------------------------------------------------------------------------
      // increment conCount
      //------------------------------------------------------------------------
      conCount++;

      //------------------------------------------------------------------------
      // report
      //------------------------------------------------------------------------
      if (record) {logFile << "  Processing reference sequence " << conCount << ": " << conName << endl;}
      if (debug) {cerr << "  Processing reference sequence " << conCount << ": " << conName << endl;}
      
      //------------------------------------------------------------------------
      // skip contig if the user only wants to analyze a specific region, and if
      // this contig does not match the name of the sequence on whcih the region
      // resides
      //------------------------------------------------------------------------
      if (region && (conName != regionSeq)) {
	if (record) {logFile << "    Reference sequence does not match sequence name for region to be analyzed... skipped." << endl;}
	if (debug) {cerr << "    Reference sequence does not match sequence name for region to be analyzed... skipped." << endl;}
	continue;
      }

      //------------------------------------------------------------------------
      // pick up contig attributes
      //------------------------------------------------------------------------
      ContigData contigData = gigReader->getContigData(conName);
      int conLength = contigData.length;
      int conNumReads =  contigData.numReads;
      string conDna =  contigData.dna;
      vector<short> conQual =  contigData.qual;

      //------------------------------------------------------------------------
      // process contig attributes
      //------------------------------------------------------------------------

      // uppercase dna
      conDna = upperCaseDna(conDna);
      
      // assign quality values to gaps if indel discovery enabled
      if (indel) {
	for (int i=0; i<conQual.size(); i++) {
	  if ((conDna.substr(i, 1) == "*") || (conDna.substr(i, 1) == "-")) {
	    conQual[i] = min(conQual[i-1], conQual[i+1]);
	  }
	}
      }
      
      //------------------------------------------------------------------------
      // make unpadded position maps for contig
      //------------------------------------------------------------------------
      if (record) {logFile << "    Making padded-unpadded position map for reference sequence." << endl;}
      if (debug) {cerr << "    Making padded-unpadded position map for reference sequence." << endl;}

     vector<int> conUnpaddedPosMapBegin = unpaddedPosMapBegin(conDna);
      vector<int> conUnpaddedPosMapEnd = unpaddedPosMapEnd(conDna);

      // calculate unpadded contig length
      int conLengthUnpadded = conUnpaddedPosMapEnd[conLength-1];
      
      // unpadded contig DNA
      string conDnaUnpadded;
      if (removeMonomer) {
	conDnaUnpadded = unpadDna(conDna);
      }

      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      // remove duplicate reads if required
      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
 
      //----------------------------------------------------------------------
      // duplicate read data structure
      //----------------------------------------------------------------------
      map<string, bool, less<string> > duplicateRead, duplicateFragment;
      int numberDups = 0;

      //----------------------------------------------------------------------
      // duplicate removal: SE reads
      //----------------------------------------------------------------------
      if (filterDups && dupMode == "se") {
	
	// report
	if (record) {logFile << "    Removing duplicate reads based on SE read alignment coordinates. " << endl;}
	if (debug) {cerr << "    Removing duplicate reads based on SE read alignment coordinates. " << endl;}

	// iterate through all alignments and register reads with their start
	// and end positions
	map<string, int, less<string> > readStart, readEnd;
	map<int, vector<string>, less<int> > startReads;

	// prepare for reading this contig's assembled reads
	bool status = gigReader->prepareReads(conName);
      
	// initialize last contig position analyzed
	int contigPosAnalyzed = 0;
	
	// declare read bufer: initial size 0
	deque<ReadData> readBuffer;
	
	// cycle through every read
	for (int r=1; r<=conNumReads; r++) {
	  
	  // if read buffer empty, read from file
	  if (readBuffer.size() <= 0) {
	    
	    // figure out how many reads to read into buffer
	    unsigned int numberReadsToRead = min(conNumReads-r+1, R);
	    
	    // read from file into buffer
	    unsigned int numberReads = gigReader->getReads(numberReadsToRead, readBuffer);
	  }
	  
	  // pop next read from front of read buffer
	  ReadData readData;
	  if (readBuffer.size() > 0) {
	    readData = readBuffer[0];
	    readBuffer.pop_front();
	  }
	  
	  // if not possible, complain and bomb
	  else {
	    if (debug) {
	      cerr << "Read buffer problem. Exiting." << endl;
	    }
	    if (record) {
	      logFile << "Read buffer problem. Exiting." << endl;
	    }
	    exit(1);
	  }
	  
	  // get assembled read attributes
	  string readName = readData.name;
	  int readLength = readData.length;
	  int contigFitLeft =  readData.contigFitLeft;
	  int contigFitRight =  readData.contigFitRight;
	  int readQualClipLeft =  readData.readQualClipLeft;
	  int readQualClipRight =  readData.readQualClipRight;
	  int readAlignmentClipLeft =  readData.readAlignmentClipLeft;
	  int readAlignmentClipRight =  readData.readAlignmentClipRight;
	  bool complemented = readData.complemented;
	  
	  // calculate alignment coordinates
	  int contigAlLeft = contigFitLeft + readAlignmentClipLeft - 1;
	  int contigAlRight = contigFitLeft + readAlignmentClipRight - 1;
	  
	  // add read to dup removal data structures
	  if (! complemented) {
	    readStart[readName] = contigAlLeft;
	    readEnd[readName] = contigAlRight;
	    startReads[contigAlLeft].push_back(readName);
	  }
	  else {
	    readStart[readName] = contigAlRight;
	    readEnd[readName] = contigAlLeft;
	    startReads[contigAlRight].push_back(readName);
	  }
	}

	// iterate through every position in readStart and get list of reads
	//   starting at that position
	for (map<int, vector<string>, less<int> >::const_iterator posIter = startReads.begin();
	     posIter != startReads.end(); posIter++) {
	  int pos = posIter->first;

	  // retrieve reads that start at this position
	  vector<string> reads = posIter->second;

	  // number of reads ending at a given position
	  map<int, int, less<int> > endCount;

	  // register alignment end position for all such reads
	  for (vector<string>::const_iterator readIter = reads.begin();
	       readIter != reads.end(); readIter++) {
	    string readName = *readIter;

 	    int end = readEnd[readName];
 	    endCount[end]++;
 	    for (int w=1; w<=WRE; w++) {
 	      endCount[end-w]++;
 	      endCount[end+w]++;
 	    }
 	  }

 	  // remove reads with end positions that occur multiple times
 	  bool oneDupKept = false;
 	  for (vector<string>::const_iterator readIter = reads.begin();
 	       readIter != reads.end(); readIter++) {
 	    string readName = *readIter;
 	    int end = readEnd[readName];
	    
 	    // if end position occurs multiple times, this read needs to be removed
 	    //   as a duplicate
 	    if (endCount[end] > 1) {

 	      // mark read if necessary
 	      if (!keepOneDup || oneDupKept) {

 		// register read as duplicate
 		duplicateRead[readName] = true;

 		// increment counter
 		numberDups++;
 	      }
 	      else {
		// mark fact that one dup already kept
		oneDupKept = true;
	      }
	    }
	  }
	}

	//----------------------------------------------------------------------
	// report
	//----------------------------------------------------------------------
	if (record) {logFile << "      Number of duplicate reads removed: " << numberDups << " (of " << conNumReads << " total reads)." << endl;}
	if (debug) {cerr << "      Number of duplicate reads removed: " << numberDups << " (of " << conNumReads << " total reads)." << endl;}
      }

      //----------------------------------------------------------------------
      // duplicate removal: PE reads
      //----------------------------------------------------------------------
      else if (filterDups && dupMode == "pe") {
	
	// report
	if (record) {logFile << "    Removing duplicate reads based on PE read alignment coordinates. " << endl;}
	if (debug) {cerr << "    Removing duplicate reads based on PE read alignment coordinates. " << endl;}

	// iterate through all alignments and register reads with their start
	// and end positions
	map<string, int, less<string> > fragmentStart, fragmentEnd;
	map<int, vector<string>, less<int> > startFragments, endFragments;

	// prepare for reading this contig's assembled reads
	bool status = gigReader->prepareReads(conName);
      
	// initialize last contig position analyzed
	int contigPosAnalyzed = 0;
	
	// declare read bufer: initial size 0
	deque<ReadData> readBuffer;
	
	// cycle through every read
	for (int r=1; r<=conNumReads; r++) {
	  
	  // if read buffer empty, read from file
	  if (readBuffer.size() <= 0) {
	    
	    // figure out how many reads to read into buffer
	    unsigned int numberReadsToRead = min(conNumReads-r+1, R);
	    
	    // read from file into buffer
	    unsigned int numberReads = gigReader->getReads(numberReadsToRead, readBuffer);
	  }
	  
	  // pop next read from front of read buffer
	  ReadData readData;
	  if (readBuffer.size() > 0) {
	    readData = readBuffer[0];
	    readBuffer.pop_front();
	  }
	  
	  // if not possible, complain and bomb
	  else {
	    if (debug) {
	      cerr << "Read buffer problem. Exiting." << endl;
	    }
	    if (record) {
	      logFile << "Read buffer problem. Exiting." << endl;
	    }
	    exit(1);
	  }
	  
	  // get assembled read attributes
	  string readName = readData.name;
	  int readLength = readData.length;
	  int contigFitLeft =  readData.contigFitLeft;
	  int contigFitRight =  readData.contigFitRight;
	  int readQualClipLeft =  readData.readQualClipLeft;
	  int readQualClipRight =  readData.readQualClipRight;
	  int readAlignmentClipLeft =  readData.readAlignmentClipLeft;
	  int readAlignmentClipRight =  readData.readAlignmentClipRight;
	  bool complemented = readData.complemented;
	  
	  // determine fragment name and read pair serial
	  boost::regex peReadPattern("^(\\S+)\\.([1|2])$");

	  // regex match
	  boost::smatch match;


	  string fragmentName;
	  string serial;
	  if(boost::regex_search(readName, match, peReadPattern)) {
	    
	    // fetch cigar element
	    fragmentName = match[1];
	    serial = match[2];
	  }
	  else {
	    // read is not in PE naming format: skip
	    continue;
	  }

	  // calculate alignment coordinates
	  int contigAlLeft = contigFitLeft + readAlignmentClipLeft - 1;
	  int contigAlRight = contigFitLeft + readAlignmentClipRight - 1;
	  
	  // add read to dup removal data structures
	  if (! complemented) {
	    startFragments[contigAlLeft].push_back(fragmentName);
	    fragmentStart[fragmentName] = contigAlLeft;
	  }
	  else {
	    endFragments[contigAlRight].push_back(fragmentName);
	    fragmentEnd[fragmentName] = contigAlRight;
	  }
	}

	//----------------------------------------------------------------------
	// PE read duplicate removal algorithm
	//----------------------------------------------------------------------

	// iterate through every position in fragmentStart and get list of fragments
	//   starting at that position
	for (map<int, vector<string>, less<int> >::const_iterator posIter = startFragments.begin();
	     posIter != startFragments.end(); posIter++) {
	  int pos = posIter->first;

	  // retrieve reads that start at this position
	  vector<string> fragments = posIter->second;

	  // number of reads ending at a given position
	  map<int, int, less<int> > endCount;

	  // register alignment end position for all such reads
	  for (vector<string>::const_iterator fragmentIter = fragments.begin();
	       fragmentIter != fragments.end(); fragmentIter++) {
	    string fragmentName = *fragmentIter;

	    if (fragmentEnd.count(fragmentName) > 0) {
	      int end = fragmentEnd[fragmentName];
	      endCount[end]++;
	    }
 	  }

 	  // remove reads with end positions that occur multiple times
 	  bool oneDupKept = false;
 	  for (vector<string>::const_iterator fragmentIter = fragments.begin();
 	       fragmentIter != fragments.end(); fragmentIter++) {
 	    string fragmentName = *fragmentIter;
	    if (fragmentEnd.count(fragmentName) > 0) {
	      int end = fragmentEnd[fragmentName];
	    
	      // if end position occurs multiple times, this read needs to be removed
	      //   as a duplicate
	      if (endCount[end] > 1) {
		
		// mark read if necessary
		if (!keepOneDup || oneDupKept) {
		  
		  // register two PE end-reads as duplicate
		  duplicateRead[fragmentName + ".1"] = true;
		  duplicateRead[fragmentName + ".2"] = true;
		  
		  // increment counter
		  numberDups += 2;
		}
		else {
		  // mark fact that one dup already kept
		  oneDupKept = true;
		}
	      }
	    }
	    else {
	      // skip this fragment
	      continue;
	    }
	  }
	}

	//----------------------------------------------------------------------
	// report
	//----------------------------------------------------------------------
	if (record) {logFile << "      Number of duplicate reads removed: " << numberDups << " (of " << conNumReads << " total reads)." << endl;}
	if (debug) {cerr << "      Number of duplicate reads removed: " << numberDups << " (of " << conNumReads << " total reads)." << endl;}
      }

      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      // register allele coverage, pre-filter, and analyze positions that pass
      //------------------------------------------------------------------------
      //------------------------------------------------------------------------

      // report
      if (record) {logFile << "    Scanning alignment for variants... " << endl;}
      if (debug) {cerr << "    Scanning alignment for variants... " << endl;}

      //------------------------------------------------------------------------
      // data structure for aligned position specific allele layout
      //------------------------------------------------------------------------

      // contig-position specific coverage quantities
      map<int, int, less<int> > readCount, readCountPlus, readCountMinus;
      map<int, map_string_int, less<int> > alleleCount, alleleCountPlus, alleleCountMinus;
      map<int, map_string_longDouble, less<int> > alleleQual, alleleQualPlus, alleleQualMinus;
      
      // contig-position and individual specific coverage quantities
      map<int, map<string, vector<Basecall>, less<string> >, less<int> > individualBasecalls;
      map<int, map<string, SampleInfo, less<string> >, less<int> > individualSampleInfo;
      
      //------------------------------------------------------------------------
      // iterate through assembled reads (sorted according to fitLeft position)
      //------------------------------------------------------------------------

      // prepare for reading this contig's assembled reads
      bool status = gigReader->prepareReads(conName);
      
      // initialize last contig position analyzed
      int contigPosAnalyzed = 0;
      
      // declare read bufer: initial size 0
      deque<ReadData> readBuffer;
      
      // cycle through every read
      for (int r=1; r<=conNumReads; r++) {
	
	// if read buffer empty, read from file
	if (readBuffer.size() <= 0) {
	  
	  // figure out how many reads to read into buffer
	  unsigned int numberReadsToRead = min(conNumReads-r+1, R);
	  
	  // read from file into buffer
	  unsigned int numberReads = gigReader->getReads(numberReadsToRead, readBuffer);
	}
	
	// pop next read from front of read buffer
	ReadData readData;
	if (readBuffer.size() > 0) {
	  readData = readBuffer[0];
	  readBuffer.pop_front();
	}
	
	// if not possible, complain and bomb
	else {
	  if (debug) {
	    cerr << "Read buffer problem. Exiting." << endl;
	  }
	  if (record) {
	    logFile << "Read buffer problem. Exiting." << endl;
	  }
	  exit(1);
	}

	// get readName
	string readName = readData.name;

	//----------------------------------------------------------------------
	// if this is a duplicate read, skip
	//----------------------------------------------------------------------
	if (duplicateRead[readName]) {
	  if (debug2) {
	    cerr << "  removing duplicate read: " << readName << endl;
	  }
	  continue;
	}

	// get assembled read attributes
	int readLength = readData.length;
	string readDna = readData.dna;
	vector<short> readQual = readData.qual;
	int contigFitLeft =  readData.contigFitLeft;
	int contigFitRight =  readData.contigFitRight;
	int readQualClipLeft =  readData.readQualClipLeft;
	int readQualClipRight =  readData.readQualClipRight;
	int readAlignmentClipLeft =  readData.readAlignmentClipLeft;
	int readAlignmentClipRight =  readData.readAlignmentClipRight;
	bool complemented = readData.complemented;
	
	// calculate alignment coordinates
	int contigAlLeft = contigFitLeft + readAlignmentClipLeft - 1;
	int contigAlRight = contigFitLeft + readAlignmentClipRight - 1;
	
	//----------------------------------------------------------------------
	// if this read needs to be registered, register it
	//----------------------------------------------------------------------
	if ((! region) || 
	    (regionBegin > regionEnd) || 
	    ((conUnpaddedPosMapEnd[contigAlRight-1] >= regionBegin) && (conUnpaddedPosMapBegin[contigAlLeft-1] <= regionEnd))
	    ) {
	  
	  //--------------------------------------------------------------------
	  // process read attributes
	  //--------------------------------------------------------------------
	  
	  // uppercase dna
	  readDna = upperCaseDna(readDna);
	  
	  //--------------------------------------------------------------------
	  // assign quality values to gaps if indel discovery enabled
	  //--------------------------------------------------------------------
	  for (int i=0; i<readQual.size(); i++) {
	    if (readDna.substr(i, 1) == "*" || readDna.substr(i, 1) == "-") {
	      readQual[i] = min(readQual[i-1], readQual[i+1]);
	    }
	  }
	  
	  //----------------------------------------------------------------------
	  // count number of mismatches between read and anchor
	  // if mismatches above cutoff, discard read
	  //----------------------------------------------------------------------
	  
	  // count mismatches
	  int mm = 0;
	  for (int p=contigAlLeft; p<=contigAlRight; p++) {
	    // get aligned base and qual in read
	    //	    char b = readDna[p - contigFitLeft];
	    string b = readDna.substr(p - contigFitLeft, 1);
	    short q = readQual[p - contigFitLeft];
	    //	    char c = conDna[p-1];
	    string c = conDna.substr(p-1, 1);
	    if ((b != c) && (q >= QML) && (cim || (b != "-" && b != "*" && c != "-" && c != "*"))) {mm++;}
	  }
	  
	  // skip read if mm>MRU;
	  if (mm > MRU) {continue;}
	  
	  SampleInfo sampleInfo = extractSampleInfo(readName, sample, sampleDel);
	  string sampleName = sampleInfo.sampleId;
	  
	  //----------------------------------------------------------------------
	  // register read dna and qual values 
	  //----------------------------------------------------------------------
	  for (int p=contigAlLeft; p<=contigAlRight; p++) {
	    // get aligned base and qual in read
	    string b = readDna.substr(p - contigFitLeft, 1);
	    short q = readQual[p - contigFitLeft];
	    
	    // skip if not appropriate allele
	    if (! (((b == "A") || (b == "C") || (b == "G") || (b == "T")) || (indel && (b == "*" || b == "-")))) {
	      continue;
	    }
	    
	    //--------------------------------------------------------------------
	    // register base calls
	    //--------------------------------------------------------------------
	    
	    // register base if above minimum qual
	    if (q >= QRL) {
	      
	      // make base call data struct
	      Basecall bc;
	      bc.seqName = readName;
	      bc.base = b;
	      bc.qual = q;
	      bc.strand = "+";
	      if (complemented) {
		bc.strand = "-";
	      }
	      
	      // add base call to template
	      individualBasecalls[p][sampleName].push_back(bc);
	      
	      // make Sampleinfo data struct
	      // register sample info
	      individualSampleInfo[p][sampleName] = sampleInfo;
	      
	      // update read and allele counts
	      readCount[p]++;
	      alleleCount[p][b]++;
	      alleleQual[p][b] += q;
	      if (complemented) {
		readCountMinus[p]++;
		alleleCountMinus[p][b]++;
		alleleQualMinus[p][b] += q;
	      }
	      else {
		readCountPlus[p]++;
		alleleCountPlus[p][b]++;
		alleleQualPlus[p][b] += q;
	      }
	    }
	  }
	}
	
	//------------------------------------------------------------------------
	//------------------------------------------------------------------------
	// analyze contig positons where all assembled reads have been encountered
	//------------------------------------------------------------------------
	//------------------------------------------------------------------------
	
	//------------------------------------------------------------------------
	// determine last contig position to be analyzed
	//------------------------------------------------------------------------
	int contigPosLastToAnalyze;
	if (r >= conNumReads) {
	  contigPosLastToAnalyze = conLength;
	}
	else {
	  contigPosLastToAnalyze = contigFitLeft - 1;
	}
	
	//------------------------------------------------------------------------
	// analyze all contig positions that are ready for analysis
	//------------------------------------------------------------------------
	for(int conPos = contigPosAnalyzed + 1; conPos <= contigPosLastToAnalyze; conPos++) {
	  
	  //--------------------------------------------------------------------
	  // report progress if necessary
	  //--------------------------------------------------------------------
	  if (record && (conPos % I == 0)) {
	    logFile << "    Screening padded contig position: " << conPos 
		    << " of " << conLength
		    << " (unpadded position: " << conUnpaddedPosMapBegin[conPos-1] 
		    << "-" << conUnpaddedPosMapEnd[conPos-1] 
		    << " of " << conLengthUnpadded << ")." 
		    << endl;}
	  if (debug && (conPos % I == 0)) {
	    cerr << "    Screening padded contig position: " << conPos 
		 << " of " << conLength
		 << " (unpadded position: " << conUnpaddedPosMapBegin[conPos-1] 
		 << "-" << conUnpaddedPosMapEnd[conPos-1] 
		 << " of " << conLengthUnpadded << ")." 
		 << endl;
	  }
	  
	  //--------------------------------------------------------------------
	  // check if this position needs to be analyzed
	  //--------------------------------------------------------------------
	  if ((! region) || 
	      (regionBegin > regionEnd) || 
	      ((conUnpaddedPosMapEnd[conPos-1] >= regionBegin) && (conUnpaddedPosMapBegin[conPos-1] <= regionEnd))
	      ) {
	    
	    // extract reference sequence allele
	    string ra = conDna.substr(conPos-1, 1);
	    
	    //-------------------------------------------------------------------
	    // register contig allele if reference allele is to be used
	    //-------------------------------------------------------------------
	    if (useRefAllele) {
	      
	      // get contig allele
	      string b = ra;
	      
	      // only register appropriate allele
	      if (((b == "A") || (b == "C") || (b == "G") || (b == "T")) || (indel && (b == "*" || b == "-"))) {
		
		// deal with quality value
		short q = 0;
		if (setQR) {
		  q = QR;
		}
		else {
		  if (b == "*" || b == "-") {
		    q = min(conQual[conPos-2], conQual[conPos]);	  
		  }
		  else {
		    q = conQual[conPos-1];
		  }
		}
		
		// register contig allele if its quality is high enough or
		// if reference allese is to be forced
		if (q >= QRL || forceRefAllele) {
		  
		  // make Basecall data struct
		  Basecall bc;
		  bc.seqName = conName;
		  bc.base = b;
		  bc.qual = q;
		  bc.strand = "?";
		  
		  // register base call
		  individualBasecalls[conPos][conName].push_back(bc);
		  
		  // make Sampleinfo data struct
		  SampleInfo si;
		  si.sampleId = conName;
		  si.pedigreeId = "";
		  si.motherId = "";
		  si.fatherId = "";
		  si.sex = "";
		  
		  // register sample info
		  individualSampleInfo[conPos][conName] = si;
		  
		  // increment allele counts
		  readCount[conPos]++;
		  alleleCount[conPos][b]++;
		  alleleQual[conPos][b] += q;
		}
	      }
	    }
	    
	    //-------------------------------------------------------------------
	    //-------------------------------------------------------------------
	    // prefilter position
	    //-------------------------------------------------------------------
	    //-------------------------------------------------------------------
	    
	    //-------------------------------------------------------------------
	    // filter: read coverage between CRL and CRU, read coverage on
	    //         each strand is at least CRL2
	    //-------------------------------------------------------------------
	    if (
		readCount[conPos] < CRL || 
		readCount[conPos] > CRU ||
		readCountPlus[conPos] < CRL2 || 
		readCountMinus[conPos] < CRL2
		) {
	      readCount.erase(conPos);
	      readCountPlus.erase(conPos);
	      readCountMinus.erase(conPos);
	      alleleCount.erase(conPos);
	      alleleCountPlus.erase(conPos);
	      alleleCountMinus.erase(conPos);
	      alleleQual.erase(conPos);
	      alleleQualPlus.erase(conPos);
	      alleleQualMinus.erase(conPos);
	      individualBasecalls.erase(conPos);
	      individualSampleInfo.erase(conPos);
	      continue;
	    }     
	    
	    //-------------------------------------------------------------------
	    // get eligible alleles
	    //-------------------------------------------------------------------
	    map_string_longDouble eligibleAlleleQual;
	    int numberEligibleAlleles = 0;
	    for (map_string_int::const_iterator iter = alleleCount[conPos].begin();
		 iter != alleleCount[conPos].end(); iter++) {
	      string b = iter->first;
	      int a = alleleCount[conPos][b];
	      int aPlus = alleleCountPlus[conPos][b];
	      int aMinus = alleleCountMinus[conPos][b];
	      long double q = alleleQual[conPos][b];
	      
	      // the allele is eligible if it meets the coverage requirements or
	      // it is the reference allele and the reference allele is to be forced	     
	      if (
		  ((a >= CAL) && (aPlus >= CAL2) && (aMinus >= CAL2) && (q >= QAL)) 
		  || 
		  (forceRefAllele && (b == ra))
		  ) {
		numberEligibleAlleles++;
		eligibleAlleleQual[b] = q;
	      }
	    }
	    
	    // sort eligible alleles according to descending allele quality
	    vector<string> sortedEligibleAlleles = sortKeysByValue(eligibleAlleleQual, true);
	    
	    //-------------------------------------------------------------------
	    // filter: are there at least 2 eligible alleles?
	    //-------------------------------------------------------------------
	    if (numberEligibleAlleles < 2) {    
	      readCount.erase(conPos);
	      readCountPlus.erase(conPos);
	      readCountMinus.erase(conPos);
	      alleleCount.erase(conPos);
	      alleleCountPlus.erase(conPos);
	      alleleCountMinus.erase(conPos);
	      alleleQual.erase(conPos);
	      alleleQualPlus.erase(conPos);
	      alleleQualMinus.erase(conPos);
	      individualBasecalls.erase(conPos);
	      individualSampleInfo.erase(conPos);
	      continue;
	    }
	    
	    //-------------------------------------------------------------------
	    // set signficant alleles
	    //-------------------------------------------------------------------
	    string allele1 = sortedEligibleAlleles[0];
	    string allele2 = sortedEligibleAlleles[1];
	    
	    //-------------------------------------------------------------------
	    // filter: likely monomer misalignment artifact
	    //         this should not be performed for the first or last position
	    //-------------------------------------------------------------------
	    if (removeMonomer && (conPos != 1) && (conPos != conLength)) {
	      // get nearest non-gap bases from reference sequence
	      int leftNeighborPos = conUnpaddedPosMapBegin[conPos-2];
	      int rightNeighborPos = conUnpaddedPosMapEnd[conPos];
	      string leftNeighborBase = conDnaUnpadded.substr(leftNeighborPos-1, 1);
	      string rightNeighborBase = conDnaUnpadded.substr(rightNeighborPos-1, 1);

	      // if the variant alleles match the neighbors filter out
	      if ((leftNeighborBase == allele1 && rightNeighborBase == allele2) 
		  ||
		  (leftNeighborBase == allele2 && rightNeighborBase == allele1)) {
		readCount.erase(conPos);
		readCountPlus.erase(conPos);
		readCountMinus.erase(conPos);
		alleleCount.erase(conPos);
		alleleCountPlus.erase(conPos);
		alleleCountMinus.erase(conPos);
		alleleQual.erase(conPos);
		alleleQualPlus.erase(conPos);
		alleleQualMinus.erase(conPos);
		individualBasecalls.erase(conPos);
		individualSampleInfo.erase(conPos);
		continue;
	      }
	    }
	    
	    //-------------------------------------------------------------------
	    //-------------------------------------------------------------------
	    // report
	    //-------------------------------------------------------------------
	    //-------------------------------------------------------------------
	    if (record) {
	      logFile << "    Pos=" << conPos << " of " << conLength
		      << " (" << conUnpaddedPosMapBegin[conPos-1]
		      << "-" << conUnpaddedPosMapEnd[conPos-1]
		      << " unpadded). Calls:";
	      for (int i=0; i<sortedEligibleAlleles.size(); i++) {
		logFile << " " << sortedEligibleAlleles[i]
			<< "(" << alleleCount[conPos][sortedEligibleAlleles[i]]
			<< "," << alleleQual[conPos][sortedEligibleAlleles[i]] << ")";
	      }
	      logFile << ". Best alleles: " << allele1 << "," << allele2 << ". Ref: " << conDna[conPos-1] << "(" << conQual[conPos-1] << ")";
	    }
	    if (debug) {
	      cerr << "    Pos=" << conPos << " of " << conLength
		   << " (" << conUnpaddedPosMapBegin[conPos-1]
		   << "-" << conUnpaddedPosMapEnd[conPos-1]
		   << " unpadded). Calls:";
	      for (int i=0; i<sortedEligibleAlleles.size(); i++) {
		cerr << " " << sortedEligibleAlleles[i]
		     << "(" << alleleCount[conPos][sortedEligibleAlleles[i]]
		     << "," << alleleQual[conPos][sortedEligibleAlleles[i]] << ")";
	      }
	      cerr << ". Best alleles: " << allele1 << "," << allele2 << ". Ref: " << conDna[conPos-1] << "(" << conQual[conPos-1] << ")";
	    }
	    
	    //-------------------------------------------------------------------
	    // analyze position with full Bayesian method
	    //-------------------------------------------------------------------
	    Variation var = varBayesProb(
					 individualBasecalls[conPos], 
					 individualSampleInfo[conPos],
					 diploid, 
					 D, 
					 M, 
					 banded,
					 WB,
					 TB,
					 TR,
					 allele1, 
					 allele2,
					 cgp,
					 debug2
					 );
	    
	    // report
	    if (record) {logFile << ". P(SNP)=" << var.pSnp;}
	    if (debug) {cerr << ". P(SNP)=" << var.pSnp;}
	    
	    //-----------------------------------------------------------------
	    // report individual genotype log likelihoods if required
	    //-----------------------------------------------------------------
	    if (reportGlt) {
	      printVar2Glt(
			   gltFile,
			   conName,
			   conUnpaddedPosMapBegin[conPos-1],
			   conUnpaddedPosMapEnd[conPos-1],
			   var
			   );
	    }
	    
	    //-------------------------------------------------------------------
	    // if pSNP above threshold write variation entry into GFF3 output file
	    //-------------------------------------------------------------------
	    if (var.pSnp >= PSL) {
	      
	      // report
	      if (record) {logFile << " >= " << PSL << ". Reported." << endl;}
	      if (debug) {cerr << " >= " << PSL << ". Reported." << endl;}
	      
	      // write variation entry into GFF output file
	      printVar2Gff(
			   gffFile,
			   O,
			   conName,
			   ProgramName,
			   allele1,
			   allele2,
			   conPos,
			   conUnpaddedPosMapBegin[conPos-1],
			   conUnpaddedPosMapEnd[conPos-1],
			   var,
			   individualBasecalls
			   );
	    }
	    else {
	      // report
	      if (record) {logFile << " < " << PSL << ". Discarded." << endl;}
	      if (debug) {cerr << " < " << PSL << ". Discarded." << endl;}
	    }
	  }	  
	  
	  //------------------------------------------------------------------------
	  // erase allele data for this position
	  //------------------------------------------------------------------------
	  readCount.erase(conPos);
	  readCountPlus.erase(conPos);
	  readCountMinus.erase(conPos);
	  alleleCount.erase(conPos);
	  alleleCountPlus.erase(conPos);
	  alleleCountMinus.erase(conPos);
	  alleleQual.erase(conPos);
	  alleleQualPlus.erase(conPos);
	  alleleQualMinus.erase(conPos);
	  individualBasecalls.erase(conPos);
	  individualSampleInfo.erase(conPos);
	}
	
	// update last contig position analyzed
	contigPosAnalyzed = contigFitLeft - 1;
	
	//--------------------------------------------------------------------------
	// skip if the user only analyzes a region and contig position if higher than
	// the last position to be analyzed
	//--------------------------------------------------------------------------
	if (region && (regionBegin <= regionEnd) && (conUnpaddedPosMapBegin[contigPosAnalyzed-1] > regionEnd)) {
	  break;
	}
      }

      // report
      if (record) {logFile << "   done." << endl;}
      if (debug) {cerr << "   done." << endl;}
      
      // report
      if (record) {logFile << "  Reference sequence done." << endl;}
      if (debug) {cerr << "  Reference sequence done." << endl;}
    }
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // finish up
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // report
  if (record) {logFile << "Closing files...";}
  if (debug) {cerr << "Closing files...";}
  
  gffFile.close();
  if (reportGlt) {gltFile.close();}

  // report
  if (record) {
    logFile << " done." << endl;
    logFile << "GFF output file: " << gff << endl;
    if (reportGlt) {
      logFile << "GLT output file: " << glt << endl;
    }
    logFile << "Log file: " << log << endl;
    logFile << "Program completed." << endl;
  }
  if (debug) {
    cerr << " done." << endl;
    cerr << "GFF output file: " << gff << endl;
    if (reportGlt) {
      logFile << "GLT output file: " << glt << endl;
    }
    if (record) {
      cerr << "Log file: " << log << endl;
    }
    cerr << "Program completed." << endl;
  }

  // close log file last
  if (record) {logFile.close();}
}

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
//------------------------------------------------------------------------------
// subroutines
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Prints a record into GFF format output file.
//------------------------------------------------------------------------------
bool printVar2Glt(
		  ofstream & gltFile,
		  string conName,
		  int posBegin,
		  int posEnd,
		  Variation & var
		  ) {

  gltFile << conName 
	  << "\t" << posBegin 
	  << "\t" << posEnd;


  for (map<string, map<string, long double, less<string> >, less<string> >::const_iterator 
	 indIter = var.individualGenotypeLogDataLikelihood.begin(); 
       indIter != var.individualGenotypeLogDataLikelihood.end();
       indIter++) {
      
      string ind = indIter->first;
      map<string, long double, less<string> > gl = indIter->second;

      gltFile << "\t" << ind << ":";
      
      bool first = true;
      for (map<string, long double, less<string> >::const_iterator 
	     glIter = gl.begin();
	   glIter != gl.end();
	   glIter++) {
	string g = glIter->first;
	long double l = glIter->second;

	if (! first) {
	  gltFile << "&";
	}
	first = false;

	gltFile << g << "|" << l;
      }
    }
  gltFile << endl;
}

//------------------------------------------------------------------------------
// Prints a variation record into GFF format output file.
//------------------------------------------------------------------------------
bool printVar2Gff(
		  ofstream & gffFile,
		  int O,
		  string conName,
		  string ProgramName,
		  string allele1,
		  string allele2,
		  int posPadded,
		  int posBegin,
		  int posEnd,
		  Variation & var,
		  map<int, map<string, vector<Basecall>, less<string> >, less<int> > & individualBasecalls
		  ) {

  // write variation entry in output file
  gffFile << conName << "\t"
	  << ProgramName << "\t";
  if (allele1 == "-" || allele1 == "*" || allele2 == "-" || allele2 == "*") {
    gffFile << "INDEL" << "\t";
  }
  else {
    gffFile << "SNP" << "\t";
  }
  gffFile << posBegin << "\t"
	  << posEnd << "\t"
	  << var.pSnp << "\t"
	  << "." << "\t"
	  << "." << "\t"
	  << "alleles=" << allele1 << "," << allele2;
  if (O >= 1) {
    gffFile << ";individualGenotypes=";
    bool firstInd = true;
    for (map<string, map<string, long double, less<string> >, less<string> >::const_iterator 
	   indIter = var.individualGenotypeProbability.begin(); 
	 indIter != var.individualGenotypeProbability.end();
	 indIter++) {
      
      string ind = indIter->first;
      if (! firstInd) {
	gffFile << ",";
      }
      firstInd = false;
      vector<string> genotypes = sortKeysByValue(var.individualGenotypeProbability[ind], true);
      string g = genotypes[0];
      long double p = var.individualGenotypeProbability[ind][g];
      gffFile << ind << ":" << g << "|" << p;
    }
  }
  if (O >= 2) {
    gffFile << ";individualGenotypeProbabilities=";
    bool firstInd = true;
    for (map<string, map<string, long double, less<string> >, less<string> >::const_iterator 
	   indIter = var.individualGenotypeProbability.begin(); 
	 indIter != var.individualGenotypeProbability.end();
	 indIter++) {
      
      string ind = indIter->first;
      if (! firstInd) {
	gffFile << ",";
      }
      firstInd = false;
      gffFile << ind << ":";
      bool firstG = true;
      for (map<string, long double, less<string> >::const_iterator 
	     gIter = var.individualGenotypeProbability[ind].begin();
	   gIter != var.individualGenotypeProbability[ind].end();
	   gIter++) {
	string g = gIter->first;
	long double p = gIter->second;
	if (! firstG) {
	  gffFile << "&";
	}
	firstG = false;
	gffFile << g << "|" << p;
      }
    }
  }
  if (O >= 3) {
    gffFile << ";individualAlleleCounts=";
    bool firstInd = true;
    for (map<string, vector<Basecall>, less<string> >::const_iterator 
	   indIter = individualBasecalls[posPadded].begin();
	 indIter != individualBasecalls[posPadded].end(); indIter++) {
      
      // retreive individual
      string ind = indIter->first;
      vector<Basecall> baseCalls = indIter->second;
      if (! firstInd) {
	gffFile << ",";
      }
      firstInd = false;
      gffFile << ind << ":";
      
      int count1 = 0; int count2 = 0; int count1Pos = 0; int count1Min = 0; int count2Pos = 0; int count2Min = 0; 
      int qual1 = 0; int qual2 = 0; int qual1Pos = 0; int qual1Min = 0; int qual2Pos = 0; int qual2Min = 0;
      for(vector<Basecall>::const_iterator iter = baseCalls.begin();
	  iter != baseCalls.end(); iter++) {
	Basecall bc = *iter;
	
	// filter irrelevant base calls
	if ((bc.base != allele1) && (bc.base != allele2)) {
	  continue;
	}
	
	// update allele counts
	if (bc.base == allele1) {
	  count1++;
	  qual1 += bc.qual;
	  if (bc.strand == "+") {
	    count1Pos++;
	    qual1Pos += bc.qual;
	  }
	  else if (bc.strand == "-") {
	    count1Min++;
	    qual1Min += bc.qual;
	  }
	}
	else if (bc.base == allele2) {
	  count2++;
	  qual2 += bc.qual;
	  if (bc.strand == "+") {
	    count2Pos++;
	    qual2Pos += bc.qual;
	  }
	  else if (bc.strand == "-") {
	    count2Min++;
	    qual2Min += bc.qual;
	  }
	}
      }
      gffFile << allele1 << "|" << count1 << "|" << qual1 << "|+|" << count1Pos << "|" << qual1Pos << "|-|" << count1Min << "|" << qual1Min << "&" << allele2 << "|" << count2 << "|" << qual2 << "|+|" << count2Pos << "|" << qual2Pos << "|-|" << count2Min << "|" << qual2Min;
    }
  }
  if (O >= 4) {
    gffFile << ";individualReadAlleles=";
    bool firstInd = true;
    for (map<string, vector<Basecall>, less<string> >::const_iterator 
	   indIter = individualBasecalls[posPadded].begin();
	 indIter != individualBasecalls[posPadded].end(); indIter++) {
      
      // retreive individual
      string ind = indIter->first;
      vector<Basecall> baseCalls = indIter->second;
      if (! firstInd) {
	gffFile << ",";
      }
      firstInd = false;
      gffFile << ind << ":";
      bool firstB = true;
      for(vector<Basecall>::const_iterator iter = baseCalls.begin();
	  iter != baseCalls.end(); iter++) {
	Basecall bc = *iter;
	
	string seqName = bc.seqName;
	string b = bc.base;
	short q = bc.qual;
	string s = bc.strand;
	if (! firstB) {
	  gffFile << "&";
	}
	firstB = false;
	gffFile << seqName << "|" << s << "|" << b << "|" << q;
      }
    }
  }
  gffFile << endl;

  //----------------------------------------------------------------------------
  // return
  //----------------------------------------------------------------------------
  return true;
}

