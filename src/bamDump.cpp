//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// bamDump
// BAM file dumper (into human-readable text format)
// Copyright 2009 Gabor T. Marth, Boston College
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

// BamReader headers
#include "BamReader.h"

// uses
using namespace std; 
using namespace TCLAP; 
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
// register an alignment into coverage data structures
//------------------------------------------------------------------------------

int calculateAlignmentEnd(
			  BamAlignment &
			  );

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// static variables
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

static string ProgramName("bamDump");
static string ProgramDescription("BAM file text viewer program.");
static string ProgramVersion("0.0.1");
static string ProgramDate("2009-04-07");
static string ProgramDeveloper("Gabor T. Marth");
static string ProgramInstitution("Boston College");
static string ProgramCopyrightDates("2009.");

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

// ----------------------------------------------------------------------------------------------------------------------------------- //
// ----------------------------------------------------------------------------------------------------------------------------------- //
bool PrintBamText (const BamAlignment& bAlignment) {
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

void PrintReads(const BamAlignmentVector& reads) {
  
  cerr << "Retrieved " << reads.size() << " reads." << endl;
  cerr << "Printing all reads to STDOUT, using custom format." << endl;
  
  for_each( reads.begin(), reads.end(), PrintBamCustom() );
  
  cerr << "Finished printing." << endl;
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
  ArgStruct argBam;
  arg = argBam; 
  arg.shortId = ""; 
  arg.longId = "bam"; 
  arg.description = "Read alignment input file (BAM format)";
  arg.required = true; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_bam(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

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
  string bam = cmd_bam.getValue();
  bool debug = cmd_debug.getValue();

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // check and fix command line options
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // report command line and parameters used
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  map<bool, string, less<bool> > bool2String;
  bool2String[false] = "false";
  bool2String[true] = "true";

  if (debug) {
    cerr << "Command line:";
    for (int i=0; i<argc; i++) {
      cerr << " " << argv[i];
    }
    cerr << endl;
    cerr << endl;
    cerr << "Complete list of parameter values:" << endl;
    cerr << "  --bam = " << bam << endl;
    cerr << "  --debug = " <<  bool2String[debug] << endl;
    cerr << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // variables
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // process input file
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  
  //----------------------------------------------------------------------------
  // open input file
  //----------------------------------------------------------------------------

  // report
  if (debug) {cerr << "Opening BAM format alignment input file: " << bam << " ...";}
    
  // open
  const char * bamFileNameC = bam.c_str();
  string bamIndexFileName = bam + ".bai";
  const char * bamIndexFileNameC = bamIndexFileName.c_str();
  BamReader bReader;
  bReader.Open(bamFileNameC, bamIndexFileNameC);
    
  // report
  if (debug) {cerr << " done." << endl;}

  //----------------------------------------------------------------------------
  // read and print header
  //----------------------------------------------------------------------------

  // retrieve header information
  string bamHeader = bReader.GetHeaderText();
  if (debug) {
    cerr << "Printing header information." << endl;
  }

  cout << "Header information:" << endl;
  cout << bamHeader << endl;
  cout << endl;

  //----------------------------------------------------------------------------
  // read and print reference sequence info
  //----------------------------------------------------------------------------

  if (debug) {
    cerr << "Printing reference sequence information." << endl;
  }
  
  // Get the names of all the reference sequences in the BAM file
  RefVector refDatas = bReader.GetReferenceData();
  cout << "Reference sequence information:" << endl;
  cout << "Number of ref seqs: " << bReader.GetReferenceCount() << endl;
  cout << "Ref sequences:" << endl;
  
  if (!refDatas.empty()) {
    RefVector::iterator refIter = refDatas.begin();
    RefVector::iterator refEnd  = refDatas.end();
    for( ; refIter != refEnd; ++refIter) {
      RefData refData = *refIter;
      string refName = refData.RefName;
      int refId = bReader.GetReferenceID(refName);
      cout << "refName=" << refName << " refId=" << refId << " length=" << refData.RefLength << endl;
    }
  }
  cout << endl;

  //----------------------------------------------------------------------------
  // read and print alignment info
  //----------------------------------------------------------------------------
  cout << "Alignment information:" << endl;
  if (!refDatas.empty()) {
    if (bReader.Rewind()) {
      BamAlignment ba;
      while (bReader.GetNextAlignment(ba)) {
	PrintBamText(ba);
      }  
    }
  }

  //----------------------------------------------------------------------------
  // finish up
  //----------------------------------------------------------------------------

  // close input file
  bReader.Close();
  
  if (debug) {
    cerr << "Program completed." << endl;
  }
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

