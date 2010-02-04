//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// gigaDump
// Dumps the contents of a GIG format binary assembly file into an ace file
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

// private libraries
#include "Class-GigReader.h"
#include "Function-Sequence.h"
#include "Function-Generic.h"


// uses
using std::cin;
using std::cout;
using std::clog;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::fstream;
using std::ofstream;
using std::ios;
using namespace TCLAP; 
using namespace __gnu_cxx;

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// type definitions
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// command line argument data structure
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
// coordinate tag data structure
//------------------------------------------------------------------------------
struct CoordTag {
  string seqName;
  string type;
  string source;
  int startUnpadded;
  int endUnpadded;
  int startPadded;
  int endPadded;
  string time;
  vector<string> commentLines;
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// subroutines
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// static variables
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

static string ProgramName("gigaDump");
static string ProgramDescription("Utility program to convert GIG format binary assembly file + GFF3 SNP annotation file to tagged ACE format.");
static string ProgramVersion("0.2.1");
static string ProgramDate("2008-03-12");
static string ProgramDeveloper("Gabor T. Marth");
static string ProgramInstitution("Boston College");
static string ProgramCopyrightDates("2007, 2008.");

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
  CmdLine cmd("", ' ', ProgramVersion);
    
  //----------------------------------------------------------------------------
  // add program-specific command line arguments
  //----------------------------------------------------------------------------

  // initialize arg
  ArgStruct arg;

  // gig: GIG format input file
  ArgStruct argGig;
  arg = argGig; 
  arg.shortId = ""; 
  arg.longId = "gig"; 
  arg.description = "GIG binary format input assembly file";
  arg.required = true; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_gig(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // gff: GFF3 format input file
  ArgStruct argGff;
  arg = argGff; 
  arg.shortId = ""; 
  arg.longId = "gff"; 
  arg.description = "GFF3 text format annotation input file reporting SNP and INDEL candidates";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_gff(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // ace: ACE format output file
  ArgStruct argAce;
  arg = argAce; 
  arg.shortId = ""; 
  arg.longId = "ace"; 
  arg.description = "ACE format assembly output file";
  arg.required = true; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_ace(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

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

  // L: sequence line length for ace file
  ArgStruct argL;
  arg = argL;
  arg.shortId = "";
  arg.longId = "L";
  arg.description = "Sequence line length for ace file";
  arg.required = false;
  arg.defaultValueString = "50";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_L(arg.shortId, arg.longId, arg.description, arg.required, 50, arg.type, cmd);

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

  string gig = cmd_gig.getValue();
  string gff = cmd_gff.getValue();
  string ace = cmd_ace.getValue();
  string log = cmd_log.getValue();
  int I = cmd_I.getValue();
  int L = cmd_L.getValue();
  bool debug = cmd_debug.getValue();
  bool debug2 = cmd_debug2.getValue();

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // variables
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // record progress into log file?
  bool record = false;
  if (log != "") {
    record = true;
  }

  // Assembly info line: AS
  boost::regex patternGffLine("^(\\S+)\\t(\\S+)\\t(\\S+)\\t(\\S+)\\t(\\S+)\\t(\\S+)\\t(\\S+)\\t(\\S+)\\t(\\S+)");
   
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
    logFile << "  --gig = " << gig << endl;
    logFile << "  --gff = " << gff << endl;
    logFile << "  --ace = " << ace << endl;
    logFile << "  --L = " << L << endl;
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
    cerr << "  --gig = " << gig << endl;
    cerr << "  --gff = " << gff << endl;
    cerr << "  --ace = " << ace << endl;
    cerr << "  --L = " << L << endl;
    cerr << "  --I = " << I << endl;
    cerr << "  --debug = " <<  bool2String[debug] << endl;
    cerr << "  --debug2 = " <<  bool2String[debug2] << endl;
    cerr << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // read GFF3 format annotation input file if present
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // define hash containing list of tags for each contig
  map<string, vector<CoordTag>, less<string> > gffTags;

  //----------------------------------------------------------------------------
  // read and parse tags
  //----------------------------------------------------------------------------
  if (gff != "") {
 
    // report
    if (record) {logFile << "Opening GFF3 format annotation input file for reading: " << gff << " ...";}
    if (debug) {cerr << "Opening GFF3 format annotation input file for reading: " << gff << " ...";}
    
    // open
    ifstream gffIn(gff.c_str(), ios::in);

    if (!gffIn) {
      cerr << "  unable to open file: " << gff << endl;
      exit(1);
    }
    if (record) {logFile << " done." << endl;}
    if (debug) {cerr << " done." << endl;}

    if (record) {logFile << "Reading GFF3 format annotations ...";}
    if (debug) {cerr << "Reading GFF3 format annotations ...";}

    //--------------------------------------------------------------------------
    // read ace file line-by-line
    //--------------------------------------------------------------------------
    string line;
    boost::smatch match;
    while (getline(gffIn, line)) {
      if(boost::regex_search(line, match, patternGffLine)) {

	// parse annotation line
	string seqName = match[1];
	string source = match[2];
	string type = match[3];
	int startU = string2Int(match[4]);
	int endU = string2Int(match[5]);
	string score = match[6];
	string attributes = match[9];

	// assign to CoordTag
	CoordTag tag;
	tag.seqName = seqName;
	tag.type = type;
	tag.source = source;
	tag.startUnpadded = startU;
	tag.endUnpadded = endU;
	tag.time = "000000:000000";
	tag.commentLines.push_back("Variation type=" + type);
	tag.commentLines.push_back("P=" + score);
	

	// register with seqName
	gffTags[seqName].push_back(tag);
      }
    }

    if (record) {logFile << " done." << endl;}
    if (debug) {cerr << " done." << endl;}
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // read binary GIG  assembly input file
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // open GIG file for reading. bomb if unable to open
  //----------------------------------------------------------------------------

  // report
  if (record) {logFile << "Opening binary assembly file for reading: " << gig << " ...";}
  if (debug) {cerr << "Opening binary assembly file for reading: " << gig << " ...";}

  // open
  FILE * gigIn = fopen64(gig.c_str(), "r");
  if (!gigIn) {
    cerr << "Unable to open file. Exiting..."<< endl;
    exit(1);
  }
  if (record) {logFile << " done" << endl;}
  if (debug) {cerr << " done" << endl;}
			
  //----------------------------------------------------------------------------
  // create GigReader object
  //----------------------------------------------------------------------------
  if (record) {logFile << "Making GigReader object...";}
  if (debug) {cerr << "Making GigReader object...";}
  GigReader * gigReader = new GigReader(gigIn);

  if (record) {logFile << " done." << endl;}
  if (debug) {cerr << " done." << endl;}

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // open ACE output file
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // report
  if (record) {logFile << "Opening output ACE file: " << ace << " ...";}
  if (debug) {cerr << "Opening output ACE file: " << ace << " ...";}

  // open
  ofstream aceOut(ace.c_str(), ios::out);
  if (!aceOut) {
    cerr << " unable to open file. Exiting..." << endl;
    exit(1);
  }
  if (record) {logFile << " done." << endl;}
  if (debug) {cerr << " done." << endl;}
	
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // process assembly block
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // get AssemblyData
  AssemblyData assemblyData = gigReader->getAssemblyData();

  // retreive data
  string assName = assemblyData.name;
  int numCons = assemblyData.numContigs;
  int assNumReads = assemblyData.numReads;

  // write AS line
  aceOut << "AS " << numCons << " " << assNumReads << endl;
  aceOut << endl;

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // iterate through contigs
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // report
  if (record) {logFile << "Processing assembly. Number of contigs: " << numCons << ". Number or reads: " << assNumReads << endl;}
  if (debug) {cerr << "Processing assembly. Number of contigs: " << numCons << ". Number or reads: " << assNumReads << endl;}

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // read contig specific information by random access to contig block
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // iterate through every contig in index
  //----------------------------------------------------------------------------
  int conCount = 0;
  for (vector<string>::const_iterator contigIter = assemblyData.contigNames.begin();
	contigIter != assemblyData.contigNames.end(); contigIter++) {
    string contigName = * contigIter;
    ContigData contigData = gigReader->getContigData(contigName);

    // report
    if (record) {logFile << "  Contig: " << contigName << " length=" << contigData.length << " numReads=" << contigData.numReads << endl;}
    if (debug) {cerr << "  Contig: " << contigName << " length=" << contigData.length << " numReads=" << contigData.numReads << endl;}

    //--------------------------------------------------------------------------
    // increment conCount
    //--------------------------------------------------------------------------
    conCount++;

    //--------------------------------------------------------------------------
    // pick up contig attributes
    //--------------------------------------------------------------------------
    int conLength = contigData.length;
    int conNumReads =  contigData.numReads;
    int conNumBasesegments =  contigData.numBasesegments;
    string conDna =  contigData.dna;
    vector<short> conQual =  contigData.qual;

    // replace "-" with "*" in conDna
    std::replace(conDna.begin(), conDna.end(), '-', '*');

    //--------------------------------------------------------------------------
    // unpad contig base quality sequence
    //--------------------------------------------------------------------------
    vector<short> conQualU;
    for (int i=0; i<conDna.size(); i++) {
      string b = conDna.substr(i,1);
      if ((b != "*") && (b != "-")) {
	conQualU.push_back(conQual[i]);
      }
    }

    //--------------------------------------------------------------------------
    // write CO entry: contig descriptor and dna
    //--------------------------------------------------------------------------
    aceOut << "CO " << contigName << " " << conLength << " " << conNumReads << " 1 U" << endl;
    printDna(aceOut, conDna, L);      
    aceOut << endl;

    //--------------------------------------------------------------------------
    // write BQ entry: contig base quality
    //--------------------------------------------------------------------------
    aceOut << "BQ" << endl;
    printQual(aceOut, conQualU, L);    
    aceOut << endl;

    //--------------------------------------------------------------------------
    // read assembled read information to produce AF lines
    //--------------------------------------------------------------------------
    if (record) {logFile << "    Writing read alignment positions." << endl;}
    if (debug) {cerr << "    Writing read alignment positions." << endl;}

    // prepare for reading this contig's assembled reads
    bool status = gigReader->prepareReads(contigName);

    // cycle through every read
    for (int r=1; r<=conNumReads; r++) {

      // retreive read data from GIG file
      ReadData readData;
      bool status = gigReader->getNextRead(readData);

      // get assembled read attributes
      string readName = readData.name;
      int readLength = readData.length;
      string readDna = readData.dna;
      int contigFitLeft =  readData.contigFitLeft;
      int contigFitRight =  readData.contigFitRight;
      int readQualClipLeft =  readData.readQualClipLeft;
      int readQualClipRight =  readData.readQualClipRight;
      int readAlignmentClipLeft =  readData.readAlignmentClipLeft;
      int readAlignmentClipRight =  readData.readAlignmentClipRight;
      bool complemented = readData.complemented;
      
      string strand = "U";
      if (complemented) {
	strand = "C";
      }

      // write AF line for this read
      aceOut << "AF " << readName << " " << strand << " " << contigFitLeft << endl;
    }

    //--------------------------------------------------------------------------
    // write BS (base segment lines)
    //--------------------------------------------------------------------------
    if (record) {logFile << "    Writing base segments." << endl;}
    if (debug) {cerr << "    Writing base segments." << endl;}

    for (vector<BasesegmentData>::const_iterator bsIter = contigData.basesegments.begin();
	 bsIter != contigData.basesegments.end(); bsIter++) {
      BasesegmentData bs = * bsIter;
      aceOut << "BS " << bs.contigBegin << " " << bs.contigEnd << " " << bs.readName << endl;     
    }
    aceOut << endl;

    //--------------------------------------------------------------------------
    // read assembled read information to produce RD, QA, DS lines
    //--------------------------------------------------------------------------
    if (record) {logFile << "    Writing read information." << endl;}
    if (debug) {cerr << "    Writing read information." << endl;}

    // prepare for reading this contig's assembled reads
    status = gigReader->prepareReads(contigName);

    // cycle through every read
    for (int r=1; r<=conNumReads; r++) {

      // retreive read data from GIG file
      ReadData readData;
      bool status = gigReader->getNextRead(readData);

      // get assembled read attributes
      string readName = readData.name;
      int readLength = readData.length;
      string readDna = readData.dna;
      int contigFitLeft =  readData.contigFitLeft;
      int contigFitRight =  readData.contigFitRight;
      int readQualClipLeft =  readData.readQualClipLeft;
      int readQualClipRight =  readData.readQualClipRight;
      int readAlignmentClipLeft =  readData.readAlignmentClipLeft;
      int readAlignmentClipRight =  readData.readAlignmentClipRight;
      bool complemented = readData.complemented;
      
      string strand = "U";
      if (complemented) {
	strand = "C";
      }

      // replace "-" with "*" in readDna
      std::replace(readDna.begin(), readDna.end(), '-', '*');

      // write RD line for this read
      aceOut << "RD " << readName << " " << readLength << " 0 0" << endl;
      printDna(aceOut, readDna, L);    
      aceOut << endl;

      // write QA line: make quality clip the same as the alignment clip
      aceOut << "QA " << readQualClipLeft << " " << readQualClipRight << " " << readAlignmentClipLeft << " " << readAlignmentClipRight << endl;

      // write DS line
      aceOut << "DS CHROMAT_FILE: none PHD_FILE: none TIME: Mon Jan 01 00:00:00 000" << endl;
      aceOut << endl;
    }

    //--------------------------------------------------------------------------
    // print contig tags read from GFF3 annotation file
    //--------------------------------------------------------------------------
    if (record) {logFile << "    Writing contig tags." << endl;}
    if (debug) {cerr << "    Writing contig tags." << endl;}
    if (gffTags.count(contigName) > 0) {

      // retreive tagList for this contig
      vector<CoordTag> tagList = gffTags[contigName];

      // make padded position map
      vector<int> paddedPos;
      for (int p=1; p<=conDna.size(); p++) {
	string b = conDna.substr(p-1, 1);
	if ((b != "*") && (b != "-")) {
	  paddedPos.push_back(p);
	}
      }

      // print tags
      for (vector<CoordTag>::const_iterator tagIter = tagList.begin();
	   tagIter != tagList.end(); tagIter++) {
	CoordTag tag = * tagIter;
	
	// bomb if coordinates are non-sensical
	if ((tag.startUnpadded-1 > paddedPos.size()) || (tag.endUnpadded-1 > paddedPos.size())) {
	  cerr << "Tag position beyond contig coordinates. Exiting..." << endl;
	  exit(1);
	}

	// print tag
	aceOut << "CT{" << endl;
	aceOut << contigName << " " << "comment" << " " << tag.source << " " << paddedPos[tag.startUnpadded-1] << " " << paddedPos[tag.endUnpadded-1] << " " << tag.time << endl;
	for (vector<string>::const_iterator lineIter = tag.commentLines.begin();
	     lineIter != tag.commentLines.end(); lineIter++) {
	  string line = * lineIter;
	  aceOut << line << endl;
	}
	aceOut << "}" << endl;
	aceOut << endl;
      }
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
  
  // close files
  fclose(gigIn);
  aceOut.close();

  // report
  if (record) {logFile << " done." << endl;}
  if (debug) {cerr << " done." << endl;}

  // report
  if (record) {logFile << "Program completed." << endl;}
  if (debug) {cerr << "Program completed." << endl;}
}

