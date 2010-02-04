//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// bamStat
// Data coverage and SNP statistic generator based on BAM files
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

// "hash_map" true hashes
#include <ext/hash_map>

// private libraries
#include "Function-Generic.h"
#include "Class-BedReader.h"

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
bool registerAlignment(
		       BamAlignment, 
		       map<int, string, less<int> > &, 
		       map<int, string, less<int> > &,
		       const int leftBound,
		       const int rightBound
		       );

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// static variables
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

static string ProgramName("bamStat");
static string ProgramDescription("BAM alignment statistics program.");
static string ProgramVersion("0.0.1");
static string ProgramDate("2009-04-09");
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

  // input file: target regions 
  ArgStruct argTargets;
  arg = argTargets; 
  arg.shortId = ""; 
  arg.longId = "targets"; 
  arg.description = "Target region input file (BED format)";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_targets(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // input file: dbSNP positions
  ArgStruct argDbsnps;
  arg = argDbsnps; 
  arg.shortId = ""; 
  arg.longId = "dbsnps"; 
  arg.description = "dbSNP position input file (BED format)";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_dbsnps(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // input file: HapMap positions
  ArgStruct argHmsnps;
  arg = argHmsnps; 
  arg.shortId = ""; 
  arg.longId = "hmsnps"; 
  arg.description = "HapMap position input file (BED format)";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_hmsnps(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // input file: called SNP positions
  ArgStruct argCalledsnps;
  arg = argCalledsnps; 
  arg.shortId = ""; 
  arg.longId = "calledsnps"; 
  arg.description = "Called SNP position input file (BED format)";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_calledsnps(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

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

  string bam = cmd_bam.getValue();
  string targets = cmd_targets.getValue();
  string dbsnps = cmd_dbsnps.getValue();
  string hmsnps = cmd_hmsnps.getValue();
  string calledsnps = cmd_calledsnps.getValue();
  int R = cmd_R.getValue();
  int I = cmd_I.getValue();
  bool debug = cmd_debug.getValue();
  bool debug2 = cmd_debug2.getValue();

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
    cerr << "  --targets = " << targets << endl;
    cerr << "  --dbsnps = " << dbsnps << endl;
    cerr << "  --hmsnps = " << hmsnps << endl;
    cerr << "  --calledsnps = " << calledsnps << endl;
    cerr << "  --R = " << R << endl;
    cerr << "  --I = " << I << endl;
    cerr << "  --debug = " <<  bool2String[debug] << endl;
    cerr << "  --debug2 = " <<  bool2String[debug2] << endl;
    cerr << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // variables
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  map<int, string, less<int> > pileDna, pileQual;
  int firstSp = 1000000000;
  int lastSp = 0;

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // process target region input file if required
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  vector<BedData> targetRegions;
  if (targets != "") {

    //--------------------------------------------------------------------------
    // open input BED file if 
    //--------------------------------------------------------------------------

    // report
    if (debug) {cerr << "Making BedReader object: " << targets << " ...";}

    // make
    BedReader * br = new BedReader(targets);

    if (! br->isOpen()) {
      cerr << "Unable to open target file: " << targets << "... terminating." << endl;
      exit(1);
    }
    
    //--------------------------------------------------------------------------
    // iterate through entries
    //--------------------------------------------------------------------------
    BedData bd;
    while (br->getNextEntry(bd)) {
      //     cout << bd.seq << "\t" << bd.left << "\t" << bd.right << "\t" << bd.desc << endl;
      targetRegions.push_back(bd);
     }

    //--------------------------------------------------------------------------
    // close
    //--------------------------------------------------------------------------
    br->close();
  }
  if (debug) {
    cerr << "Number of target regions: " << targetRegions.size() << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // process dbSNP input file if required
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  map<string, map<int, bool, less<int> >, less<string> > dbsnpRegion;
  int numDbsnpRegions = 0;
  if (dbsnps != "") {

    //--------------------------------------------------------------------------
    // open input BED file if 
    //--------------------------------------------------------------------------

    // report
    if (debug) {cerr << "Making BedReader object: " << dbsnps << " ...";}
 
    // make
    BedReader * br = new BedReader(dbsnps);

    if (! br->isOpen()) {
      cerr << "Unable to open dbSNP file: " << dbsnps << "... terminating." << endl;
      exit(1);
    }
    
    //--------------------------------------------------------------------------
    // iterate through entries
    //--------------------------------------------------------------------------
    BedData bd;
    while (br->getNextEntry(bd)) {
      //     cout << bd.seq << "\t" << bd.left << "\t" << bd.right << "\t" << bd.desc << endl; 
      dbsnpRegion[bd.seq][bd.left] = true;
      numDbsnpRegions++;
    }

    //--------------------------------------------------------------------------
    // close
    //--------------------------------------------------------------------------
    br->close();
  }
  if (debug) {
    cerr << "Number of dbSNP regions: " << numDbsnpRegions << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // process HapMap input file if required
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  map<string, map<int, bool, less<int> >, less<string> > hmsnpRegion;
  int numHmsnpRegions = 0;
  if (hmsnps != "") {

    //--------------------------------------------------------------------------
    // open input BED file if 
    //--------------------------------------------------------------------------

    // report
    if (debug) {cerr << "Making BedReader object: " << hmsnps << " ...";}
 
    // make
    BedReader * br = new BedReader(hmsnps);

    if (! br->isOpen()) {
      cerr << "Unable to open hmSNP file: " << hmsnps << "... terminating." << endl;
      exit(1);
    }
    
    //--------------------------------------------------------------------------
    // iterate through entries
    //--------------------------------------------------------------------------
    BedData bd;
    while (br->getNextEntry(bd)) {
      //     cout << bd.seq << "\t" << bd.left << "\t" << bd.right << "\t" << bd.desc << endl; 
      hmsnpRegion[bd.seq][bd.left] = true;
      numHmsnpRegions++;
    }

    //--------------------------------------------------------------------------
    // close
    //--------------------------------------------------------------------------
    br->close();
  }
  if (debug) {
    cerr << "Number of hmSNP regions: " << numHmsnpRegions << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // process called SNP input file if required
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  map<string, map<int, bool, less<int> >, less<string> > calledsnpRegion;
  int numCalledsnpRegions = 0;
  if (calledsnps != "") {

    //--------------------------------------------------------------------------
    // open input BED file if 
    //--------------------------------------------------------------------------

    // report
    if (debug) {cerr << "Making BedReader object: " << calledsnps << " ...";}
 
    // make
    BedReader * br = new BedReader(calledsnps);

    if (! br->isOpen()) {
      cerr << "Unable to open calledsnp file: " << calledsnps << "... terminating." << endl;
      exit(1);
    }
    
    //--------------------------------------------------------------------------
    // iterate through entries
    //--------------------------------------------------------------------------
    BedData bd;
    while (br->getNextEntry(bd)) {
      calledsnpRegion[bd.seq][bd.left] = true;
      numCalledsnpRegions++;
    }

    //--------------------------------------------------------------------------
    // close
    //--------------------------------------------------------------------------
    br->close();
  }
  if (debug) {
    cerr << "Number of calledsnp regions: " << numCalledsnpRegions << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // BAM input file
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // report
  if (debug) {cerr << "Opening BAM foremat alignment input file: " << bam << " ...";}
    
  // open
  const char * bamFileNameC = bam.c_str();
  string bamIndexFileName = bam + ".bai";
  const char * bamIndexFileNameC = bamIndexFileName.c_str();
  BamReader bReader;
  bReader.Open(bamFileNameC, bamIndexFileNameC);
    
  // report
  if (debug) {cerr << " done." << endl;}
 
  // Get all the reference sequences in the BAM file
  RefVector refDatas = bReader.GetReferenceData();
  if (debug) {
    cerr << "Number of ref seqs: " << bReader.GetReferenceCount() << endl;
    cerr << "Ref sequences:" << endl;
    
    if (!refDatas.empty()) {
      RefVector::iterator refIter = refDatas.begin();
      RefVector::iterator refEnd  = refDatas.end();
      for( ; refIter != refEnd; ++refIter) {
	RefData refData = *refIter;
	string refName = refData.RefName;
	int refId = bReader.GetReferenceID(refName);
	cerr << "refName=" << refName << " refId=" << refId << " length=" << refData.RefLength << endl;
      }
    }
  }

  //---------------------------------------------------------------------------
  // iterate through all reads in all target regions
  //---------------------------------------------------------------------------

  // global tallies
  long int ga = 0;
  long int gd = 0;
  long int gt = 0;
  long int gl = 0;
  long int ghm = 0;
  long int ghmp = 0;
  long int gcs = 0;
  long int gcsp = 0;
  for(vector<BedData>::const_iterator targetIter = targetRegions.begin();
      targetIter != targetRegions.end(); targetIter++) {

    // retrieve target region
    BedData target = * targetIter;

    // find target ref seq in list of ref seqs from BAM file
    int refId = bReader.GetReferenceID(target.seq);
    RefData refData = refDatas[refId];
    if (refId < bReader.GetReferenceCount() && refData.RefHasAlignments && bReader.Jump(refId, target.left)) {

      // increment target counter
      gt++;
      
      if (gt % I == 0) {
	if (debug) {
	  cerr << "  Processing target " << gt << ": " << target.desc << " " << target.seq << ":" << target.left << "-" << target.right << endl; 
	}
      }

      // update target length
      int tl = target.right - target.left + 1;

      // coverage structures
      map<int, string, less<int> > pileDna, pileQual;
      long int ta = 0;
      long int td = 0;
      long int thm = 0;
      long int thmp = 0;
      long int tcs = 0;
      long int tcsp = 0;

      // initialize last analyzed position
      int refPosLast = target.left - 1;

      BamAlignment a;
      if (debug) {
	cerr << "    Iterating through alignments that overlap target region" << endl;
      }
      while (bReader.GetNextAlignment(a)  && (a.Position <= target.right)) {
	
	// only process if mapped
	if (a.IsMapped()) {
	  
	  // update alignment counter
	  ta++;
	 
	  if (a.IsDuplicate()) {
	    if (debug2) {
	      cerr << "      Duplicate read: " << a.Name << endl;
	    }
	  }

	  // analyze positions up to the position before the beginning of this alignment
	  int refPosLeft = refPosLast + 1;
	  int refPosRight = min((int)a.Position-1, target.right);
	  for (int p = refPosLeft; p <= refPosRight; p++) {
	    
	    // calculate depth at this position
	    int ld = pileDna[p].length();
	    
	    // aggregate depth for this region
	    td += ld;
	    
	    // if a candidate SNP is called at this position, check if it is in dbSNP
	    if (calledsnpRegion.count(target.seq) > 0 && calledsnpRegion[target.seq].count(p) > 0) {
	      tcs++;
	      if (dbsnpRegion.count(target.seq) > 0 && dbsnpRegion[target.seq].count(p) > 0) {
		tcsp++;
		if (debug2) {
		  cerr << "      Called SNP refseq=" << target.seq << " pos=" << p << " IN DBSNP" << endl;
		}
	      }
	      else {
		if (debug2) {
		  cerr << "      Called SNP refseq=" << target.seq << " pos=" << p << " NOT IN DBSNP" << endl;
		}
	      }
	    }
	      
	    // if a HapMap SNP is present at this position, check if it was called
	    if (hmsnpRegion.count(target.seq) > 0 && hmsnpRegion[target.seq].count(p) > 0) {
	      thm++;
	      if (calledsnpRegion.count(target.seq) > 0 && calledsnpRegion[target.seq].count(p) > 0) {
		thmp++;
		if (debug2) {
		  cerr << "      HapMap SNP refseq=" << target.seq << " pos=" << p << " CALLED" << endl;
		}
	      }
	      else {
		if (debug2) {
		  cerr << "      HapMap SNP refseq=" << target.seq << " pos=" << p << " MISSED" << endl;
		}
	      }
	    }
	      
	    // delete data at this position after analysis
	    pileDna.erase(p);
	    pileQual.erase(p);
	  }
	  
	  // update last alignment position analyzed
	  refPosLast = max((int)a.Position - 1, target.left-1);
	  
	  // register read in coverage data structures
	  registerAlignment(a, pileDna, pileQual, target.left, target.right);
	} 
      }
      
      // close BamReader
      bReader.Close();

      // analyze remaining positions
      int refPosLeft = refPosLast + 1;
      int refPosRight = target.right;
      for (int p = refPosLeft; p <= refPosRight; p++) {
	
	// calculate depth at this position
	int ld = pileDna[p].length();
	
	// aggregate depth for this region
	td += ld;
	    
	// if a candidate SNP is called at this position, check if it is in dbSNP
	if (calledsnpRegion.count(target.seq) > 0 && calledsnpRegion[target.seq].count(p) > 0) {
	  tcs++;
	  if (dbsnpRegion.count(target.seq) > 0 && dbsnpRegion[target.seq].count(p) > 0) {
	    tcsp++;
	    if (debug2) {
	      cerr << "      Called SNP refseq=" << target.seq << " pos=" << p << " IN DBSNP" << endl;
	    }
	  }
	  else {
	    if (debug2) {
	      cerr << "      Called SNP refseq=" << target.seq << " pos=" << p << " NOT IN DBSNP" << endl;
	    }
	  }
	}
	      
	// if a HapMap SNP is present at this position, check if it was called
	if (hmsnpRegion.count(target.seq) > 0 && hmsnpRegion[target.seq].count(p) > 0) {
	  thm++;
	  if (calledsnpRegion.count(target.seq) > 0 && calledsnpRegion[target.seq].count(p) > 0) {
	    thmp++;
	    if (debug2) {
	      cerr << "      HapMap SNP refseq=" << target.seq << " pos=" << p << " CALLED" << endl;
	    }
	  }
	  else {
	    if (debug2) {
	      cerr << "      HapMap SNP refseq=" << target.seq << " pos=" << p << " MISSED" << endl;
	    }
	  }
	}
	      
	// delete data at this position after analysis
	pileDna.erase(p);
	pileQual.erase(p);
      }

      // update global quantities
      ga += ta;
      gd += td;
      gl += tl;
      gcs += tcs;
      gcsp += tcsp;
      ghm += thm;
      ghmp += thmp;

      // report
      if (debug) {
	cerr << "    Done. Stats: " << " tl=" << tl << " ta=" << ta << " td=" << td << " tcs=" << tcs << " tcsp=" << tcsp << " thm=" << thm << " thmp=" << thmp << " gl=" << gl << " ga=" << ga << " gd=" << gd << " gcs=" << gcs << " gcsp=" << gcsp << " ghm=" << ghm << " ghmp=" << ghmp << endl; 
      }

      // print target stats
      cout << "TARGET\t" << target.seq << "\t" << target.left << "\t" << target.right << "\t" << target.desc << "\t" << tl << "\t" << ta << " \t" << td << "\t" << tcs << "\t" << tcsp << "\t" << thm << "\t" << thmp << endl; 
    }
  }
  
  // report
  if (debug) {
    cerr << "Total numnber of targets: " << gt << endl;
    cerr << "Total numnber of reads overlapping targets: " << ga << endl;
    cerr << "Total length of targets: " << gl << endl;
    cerr << "Total depth in targets: " << gd << endl;
    cerr << "Total number of called SNPs: " << gcs << " In dbSNP: " << gcsp << endl;
    cerr << "Total number of HapMap SNPs: " << ghm << " Called: " << ghmp << endl;
  }

  // print global stats
  cout << "GLOBAL\t" << gl << "\t" << ga << "\t" << gd << "\t" << gcs << "\t" << gcsp << "\t" << ghm << "\t" << ghmp << endl; 

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // finish up
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
// register an alignment in coverage structures
//------------------------------------------------------------------------------
bool registerAlignment(
		       BamAlignment a, 
		       map<int, string, less<int> > & pileDna, 
		       map<int, string, less<int> > & pileQual,
		       const int leftBound,
		       const int rightBound
		       ) {

  //----------------------------------------------------------------------------
  // parse cigar and register alignment in coverage data structures
  //----------------------------------------------------------------------------
  string rDna = a.QueryBases;
  string rQual = a.Qualities;
  
  // initialize reference sequence position and read position
  int sp = a.Position;
  int rp = 1;
  
  // iterate through Cigar operations and register read alignment
  vector<CigarOp>::const_iterator cigarIter = a.CigarData.begin();
  vector<CigarOp>::const_iterator cigarEnd  = a.CigarData.end();
  for ( ; cigarIter != cigarEnd; ++cigarIter ) {
    unsigned int l = (*cigarIter).Length;
    char t = (*cigarIter).Type;
    
    if (t == 'S') { // soft clip
      rp += l;
    }
    else if (t == 'M') { // match or mismatch
      for (int i=1; i<=l; i++) {
	
	if (sp >= leftBound && sp <= rightBound) {
	  // extract aligned base
	  string b = rDna.substr(rp-1, 1);
	  string q = rQual.substr(rp-1, 1);
	
	  // add aligned base to pile
	  pileDna[sp] += b;
	  pileQual[sp] += q;
	}
	
	// update positions
	sp++;
	rp++;
      }	    
    }
    else if (t == 'D' || t == 'N') { // deletion or skipped region
      for (int i=1; i<=l; i++) {
	
	if (sp >= leftBound && sp <= rightBound) {
	  // assign aligned base as gap
	  string b = "-";
	  string q = "!";
	  
	  // add gap base to pile
	  pileDna[sp] += b;
	  pileQual[sp] += q;
	}

	// update refseq position
	sp++;
      }
    }
    else if (t == 'I') { // insertion
      rp += l;	    
    }
  }
  //----------------------------------------------------------------------------
  // return status
  //----------------------------------------------------------------------------
  return true;
}
