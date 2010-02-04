//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// bamDup
// Duplicate read identifier based on BAM files
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
#include "BamWriter.h"

// private libraries
#include "Function-Generic.h"

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

static string ProgramName("bamDup");
static string ProgramDescription("Duplicate read removal program.");
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
  ArgStruct argBamin;
  arg = argBamin; 
  arg.shortId = ""; 
  arg.longId = "bamin"; 
  arg.description = "Read alignment input file (BAM format)";
  arg.required = true; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_bamin(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // output file: read alignments
  ArgStruct argBamout;
  arg = argBamout; 
  arg.shortId = ""; 
  arg.longId = "bamout"; 
  arg.description = "Read alignment output file (BAM format)";
  arg.required = true; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_bamout(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

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

  // keepOne
  ArgStruct argKeepOne;
  arg = argKeepOne;
  arg.shortId = "";
  arg.longId = "keepOne";
  arg.description = "Keep one read from each duplicate group?";
  arg.required = false;
  arg.defaultValueString = "false";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_keepOne(arg.shortId, arg.longId, arg.description, cmd, false);

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
  string bamin = cmd_bamin.getValue();
  string bamout = cmd_bamout.getValue();
  string dupMode = cmd_dupMode.getValue();
  bool keepOne = cmd_keepOne.getValue();
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
    cerr << "  --bamin = " << bamin << endl;
    cerr << "  --bamout = " << bamout << endl;
    cerr << "  --dupMode = " << dupMode << endl;
    cerr << "  --R = " << R << endl;
    cerr << "  --I = " << I << endl;
    cerr << "  --debug = " <<  bool2String[debug] << endl;
    cerr << "  --debug2 = " <<  bool2String[debug2] << endl;
    cerr << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // constants
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  boost::regex patternReadName("^(\\S+)\\.([1|2])$");

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // variables
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  long int ga = 0;
  long int gad = 0;

  //----------------------------------------------------------------------------
  // process input file
  //----------------------------------------------------------------------------
  
  // report
  if (debug) {cerr << "Opening BAM format alignment input file: " << bamin << " ...";}
    
  // open
  const char * bamFileNameC = bamin.c_str();
  string bamIndexFileName = bamin + ".bai";
  const char * bamIndexFileNameC = bamIndexFileName.c_str();
  BamReader bReader;
  bReader.Open(bamFileNameC, bamIndexFileNameC);
    
  // report
  if (debug) {cerr << " done." << endl;}

  // retrieve header information
  string bamHeader = bReader.GetHeaderText();

  // Get the names of all the reference sequences in the BAM file
  RefVector refDatas = bReader.GetReferenceData();
  if (debug) {
    cerr << "Number of ref seqs: " << bReader.GetReferenceCount() << endl;
    cerr << "Ref sequence names:" << endl;
  
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
  // open output file and start writing it
  //---------------------------------------------------------------------------
  BamWriter bWriter;
  bWriter.Open(bamout, bamHeader, refDatas);

  //---------------------------------------------------------------------------
  // iterate through all reference sequences
  //---------------------------------------------------------------------------
  if (!refDatas.empty()) {
    RefVector::iterator refIter = refDatas.begin();
    RefVector::iterator refEnd  = refDatas.end();
    for( ; refIter != refEnd; ++refIter) {
      RefData refData = *refIter;
      string refName = refData.RefName;
      int refId = bReader.GetReferenceID(refName);
      long int ta = 0;
      long int tad = 0;
      map<int, vector<BamAlignment>, less<int> > leftAlignments;

      cerr << "Processing refName =" << refName << " refId=" << refId << " length=" << refData.RefLength << endl;

      if (refData.RefHasAlignments && bReader.Jump(refId, 1)) {
	
	map<string, bool, less<string> > dupRead;
	BamAlignment ba;
	while (bReader.GetNextAlignment(ba)) {
	  
	  // only attempt duplicate identification if mapped
	  if (ba.IsMapped()) {
	    
	    //-----------------------------------------------------------------
	    // update alignment counters
	    //-----------------------------------------------------------------
	    ta++;
	    ga++;

	    //-----------------------------------------------------------------
	    // register new alignment information
	    //-----------------------------------------------------------------

	    // get alignment left
	    int l = ba.Position;

	    // register alignment with left position
	    leftAlignments[l].push_back(ba);
	    
	    //-----------------------------------------------------------------
	    // detect duplicates at all positions < start of current alignment
	    //-----------------------------------------------------------------

	    // list of positions kept for erasal
	    vector<int> positions;

	    // iterate through all analyzable positions
	    for(map<int, vector<BamAlignment>, less<int> >::const_iterator posIter = leftAlignments.begin();
		posIter != leftAlignments.end(); posIter++) {
	      int p = posIter->first;
	      vector<BamAlignment> alignments = posIter->second;
	      
	      // break out if position reaches current read start
	      if (p >= l) {
		break;
	      }
		
	      // register position for later erasure
	      positions.push_back(p);
		
	      // number of reads with their other ends at a given position
	      map<int, int, less<int> > startCount, endCount, fragmentEndCount;
		  
	      // register start and end positions
	      for (vector<BamAlignment>::const_iterator alIter = alignments.begin();
		   alIter != alignments.end(); alIter++) {
		BamAlignment ba2 = * alIter;

		// get alignment left
		int l2 = ba2.Position;
		    
		// SE: based on true SE start and end
		if (dupMode == "se") {
		
		  // calculate other end
		  int r2 = calculateAlignmentEnd(ba2);

		  // register other end position
		  if (ba2.IsReverseStrand()) {
		    startCount[r2]++;
		  }
		  else {
		    endCount[r2]++;
		  }
		}
		// PE: based on PE right-mapped mate 5' ends
		else {
		  // only register the other ends of top-strand PE alignments
		  if ((! ba2.IsReverseStrand()) && ba2.IsMateMapped() && (ba2.RefID == ba2.MateRefID)) {
		    int fe2 = l2 + ba2.InsertSize;

		    // register fragment end
		    fragmentEndCount[fe2]++;
		  }
		}
	      }

	      // find duplicate reads
	      bool oneKeptPlusStrand = false;
	      bool oneKeptMinusStrand = false;
	      bool oneKeptFragment = false;
	      for (vector<BamAlignment>::const_iterator alIter = alignments.begin();
		   alIter != alignments.end(); alIter++) {
		BamAlignment ba2 = * alIter;

		// get alignment left
		int l2 = ba2.Position;
		    
		// duplicate flag
		bool dup = false;

		// SE: based on true SE start and end
		if (dupMode == "se") {

		  // calculate other end
		  int r2 = calculateAlignmentEnd(ba2);

		  // determine if read is duplicate
		  if (ba2.IsReverseStrand() && startCount[r2] > 1) {
		    if (keepOne && ! oneKeptMinusStrand) {
		      oneKeptMinusStrand = true;
		    }
		    else {
		      dup = true;
		    }
		  }
		  else if (endCount[r2] > 1) {
		    if (keepOne && ! oneKeptPlusStrand) {
		      oneKeptPlusStrand = true;
		    }
		    else {
		      dup = true;
		    }
		  }
		}
		// PE: based on PE mate 5' ends
		else {
		  // only check the other ends of top-strand mapped mate
		  if ((! ba2.IsReverseStrand()) && ba2.IsMateMapped() && (ba2.RefID == ba2.MateRefID)) {
		    int fe2 = l2 + ba2.InsertSize;

		    // if multiple alignments with the same other end,
		    //   mark this alignment as a duplicate and register
		    //   its mate also as duplicate
		    if (fragmentEndCount[fe2] > 1) {
		      if (keepOne && ! oneKeptFragment) {
			oneKeptFragment = true;
		      }
		      else {
			dup = true;
			
			// regex match
			boost::smatch match;
			if (boost::regex_search(ba2.Name, match, patternReadName)) {
			  string fragmentName = match[1];
			  string thisSerial = match[2];
			  string mateSerial;
			  if (thisSerial == "1") {
			    mateSerial = "2";
			  }
			  else if (thisSerial == "2") {
			    mateSerial = "1";
			  }
			  string mateName = fragmentName + "." + mateSerial;
			  dupRead[mateName] = true;
			}
		      }
		    }
		  }
		}

		// detect duplicate read and act
		if (dup || dupRead[ba2.Name]) {
		  if (debug2) {
		    string strand = "+"; if (ba2.IsReverseStrand()) {strand = "-";}
		    cerr << "  DUPLICATE READ " << ba2.Name << " strand=" << strand << " l=" << ba2.Position << " r=" << calculateAlignmentEnd(ba2) << endl;
		  }
		  tad++;
		  gad++;

		  // set Duplicate flag
		  ba2.AlignmentFlag |= 1024;

		  // erase this read name from duplicate map
		  dupRead.erase(ba2.Name);
		}
		  
		// write into output file
		bWriter.SaveAlignment(ba2);
	      }
	    }

	    // erase positions already analyzed
	    for (vector<int>::const_iterator posIter = positions.begin();
		 posIter != positions.end(); posIter++) {
	      int p = * posIter;
	      leftAlignments.erase(p);
	    }
	    
	    //-----------------------------------------------------------------
	    // repoort global stats
	    //-----------------------------------------------------------------
	    if (debug && ga % I == 0) {
	      cerr << "  reads=" << ga << " dups=" << gad << endl;
	    }
	  }
	}

	//---------------------------------------------------------------------
	// detect duplicates at all remaining positions
	//---------------------------------------------------------------------
	
	// list of positions kept for erasal
	vector<int> positions;
	
	// iterate through all analyzable positions
	for(map<int, vector<BamAlignment>, less<int> >::const_iterator posIter = leftAlignments.begin();
	    posIter != leftAlignments.end(); posIter++) {
	  int p = posIter->first;
	  vector<BamAlignment> alignments = posIter->second;
	  
	  // register position for later erasure
	  positions.push_back(p);
	  
	  // number of reads with their other ends at a given position
	  map<int, int, less<int> > startCount, endCount, fragmentEndCount;
	  
	  // register start and end positions
	  for (vector<BamAlignment>::const_iterator alIter = alignments.begin();
	       alIter != alignments.end(); alIter++) {
	    BamAlignment ba2 = * alIter;
	    
	    // get alignment left
	    int l2 = ba2.Position;
	    
	    // SE: based on true SE start and end
	    if (dupMode == "se") {
	      
	      // calculate other end
	      int r2 = calculateAlignmentEnd(ba2);
	      
	      // register other end position
	      if (ba2.IsReverseStrand()) {
		startCount[r2]++;
	      }
	      else {
		endCount[r2]++;
	      }
	    }
	    // PE: based on PE right-mapped mate 5' ends
	    else {
	      // only register the other ends of top-strand PE alignments
	      if ((! ba2.IsReverseStrand()) && ba2.IsMateMapped() && (ba2.RefID == ba2.MateRefID)) {
		int fe2 = l2 + ba2.InsertSize;
		
		// register fragment end
		fragmentEndCount[fe2]++;
	      }
	    }
	  }
	  
	  // find duplicate reads
	  bool oneKeptPlusStrand = false;
	  bool oneKeptMinusStrand = false;
	  bool oneKeptFragment = false;
	  for (vector<BamAlignment>::const_iterator alIter = alignments.begin();
	       alIter != alignments.end(); alIter++) {
	    BamAlignment ba2 = * alIter;
	    
	    // get alignment left
	    int l2 = ba2.Position;
	    
	    // duplicate flag
	    bool dup = false;
	    
	    // SE: based on true SE start and end
	    if (dupMode == "se") {
	      
	      // calculate other end
	      int r2 = calculateAlignmentEnd(ba2);
	      
	      // determine if read is duplicate
	      if (ba2.IsReverseStrand() && startCount[r2] > 1) {
		if (keepOne && ! oneKeptMinusStrand) {
		  oneKeptMinusStrand = true;
		}
		else {
		  dup = true;
		}
	      }
	      else if (endCount[r2] > 1) {
		if (keepOne && ! oneKeptPlusStrand) {
		  oneKeptPlusStrand = true;
		}
		else {
		  dup = true;
		}
	      }
	    }
	    // PE: based on PE mate 5' ends
	    else {
	      // only check the other ends of top-strand mapped mate
	      if ((! ba2.IsReverseStrand()) && ba2.IsMateMapped() && (ba2.RefID == ba2.MateRefID)) {
		int fe2 = l2 + ba2.InsertSize;
		
		// if multiple alignments with the same other end,
		//   mark this alignment as a duplicate and register
		//   its mate also as duplicate
		if (fragmentEndCount[fe2] > 1) {
		  if (keepOne && ! oneKeptFragment) {
		    oneKeptFragment = true;
		  }
		  else {
		    dup = true;
		    
		    // regex match
		    boost::smatch match;
		    if (boost::regex_search(ba2.Name, match, patternReadName)) {
		      string fragmentName = match[1];
		      string thisSerial = match[2];
		      string mateSerial;
		      if (thisSerial == "1") {
			mateSerial = "2";
		      }
		      else if (thisSerial == "2") {
			mateSerial = "1";
		      }
		      string mateName = fragmentName + "." + mateSerial;
		      dupRead[mateName] = true;
		    }
		  }
		}
	      }
	    }
	    
	    // detect duplicate read and act
	    if (dup || dupRead[ba2.Name]) {
	      if (debug2) {
		string strand = "+"; if (ba2.IsReverseStrand()) {strand = "-";}
		cerr << "  DUPLICATE READ " << ba2.Name << " strand=" << strand << " l=" << ba2.Position << " r=" << calculateAlignmentEnd(ba2) << endl;
	      }
	      tad++;
	      gad++;
	      
	      // set Duplicate flag
	      ba2.AlignmentFlag |= 1024;
	      
	      // erase this read name from duplicate map
	      dupRead.erase(ba2.Name);
	    }
	    
	    // write into output file
	    bWriter.SaveAlignment(ba2);
	  }
	}
	
	// erase positions already analyzed
	for (vector<int>::const_iterator posIter = positions.begin();
	     posIter != positions.end(); posIter++) {
	  int p = * posIter;
	  leftAlignments.erase(p);
	}    
      }
      cerr << "  number of alignments: " << ga << " duplicates=" << gad << endl;
    }
  }

  // report
  if (debug) {
  }

  // print global stats

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // finish up
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // final duplicate report
  cerr << "  number of alignments: " << ga << " duplicates=" << gad << endl;

  // close input file
  bReader.Close();
  
  // close output file
  bWriter.Close();

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
// calculate the end position for an alignment
//------------------------------------------------------------------------------
int calculateAlignmentEnd(
		       BamAlignment & a
		       ) {

  //----------------------------------------------------------------------------
  // initialize reference sequence position
  //----------------------------------------------------------------------------
  int sp = a.Position;
  
  //----------------------------------------------------------------------------
  // iterate through Cigar operations and register read alignment
  //----------------------------------------------------------------------------
  vector<CigarOp>::const_iterator cigarIter = a.CigarData.begin();
  vector<CigarOp>::const_iterator cigarEnd  = a.CigarData.end();
  for ( ; cigarIter != cigarEnd; ++cigarIter ) {
    unsigned int l = (*cigarIter).Length;
    char t = (*cigarIter).Type;
    
    if (t == 'M' || t == 'D' || t == 'N') { // match or mismatch
      sp += l;
    }
  }

  //----------------------------------------------------------------------------
  // return end position
  //----------------------------------------------------------------------------
  return sp-1;
}

