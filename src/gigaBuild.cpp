//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Program: gigaBuild
// Function: conversion from ACE format assembly file and FASTA format dna and 
//   quality files to GIG binary format assembly file
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
#include <cmath>

// "tclap" commandline parsing library
#include <tclap/CmdLine.h>

// "boost" regular expression library
#include <boost/regex.hpp>

// "hash_map" true hashes
#include <ext/hash_map>

// tell g++ how to handle string or long long keys in hashmaps
namespace __gnu_cxx {
  template<> struct hash<std::string> {
    size_t operator()( const std::string& x ) const {
      return hash<const char*>()( x.c_str() );
    }
  };

  template<> struct hash<long long> {
    size_t operator()( const long long& x ) const {
      long long ret = (x >> 32L) ^ (x & 0xffffffff);
      return (size_t)ret;
    }
  };
}

// private libraries
#include "Type-Hash.h"
#include "LargeFileSupport.h"
#include "Function-Sequence.h"
#include "Function-Generic.h"

// uses
using std::ios;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::istream;
using std::fstream;
using std::cin;
using std::cout;
using std::clog;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::map;
using namespace TCLAP; 
using namespace __gnu_cxx;

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
// constants
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

string archiveId = "GIG-0.0";

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// template functions (cannot be separated into header file)
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// ranks keys of a hash in descending numerical order of associated value
//------------------------------------------------------------------------------

template< typename keyType, typename valueType >
vector<keyType> sortKeysByValue(map<keyType, valueType, std::less<keyType> > hash, bool descend) {

  // instantiate inverse has as a multimap
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
// ranks keys of a hash_map in descending numerical order of associated value
//------------------------------------------------------------------------------

template< typename keyType, typename valueType >
vector<keyType> sortHashMapKeysByValue(hash_map<keyType, valueType> hash, bool descend) {

  // instantiate inverse has as a multimap
  multimap<valueType, keyType, std::less<valueType> > inverseHash;

  // load elements of hash into inverseHash
  for (typename hash_map<keyType, valueType>::const_iterator iter = hash.begin(); iter != hash.end(); iter++) {
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
//------------------------------------------------------------------------------
// static variables
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

static string ProgramName("gigaBuild");
static string ProgramDescription("Utility program to build GIG format binary assembly file for gigaBayes.");
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

  // ace: ACE format input file
  ArgStruct argAce;
  arg = argAce; 
  arg.shortId = ""; 
  arg.longId = "ace"; 
  arg.description = "ACE format assembly input file";
  arg.required = true; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_ace(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // fd: FASTA format DNA sequence input file
  ArgStruct argFd;
  arg = argFd; 
  arg.shortId = ""; 
  arg.longId = "fd"; 
  arg.description = "FASTA format DNA sequence input file";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_fd(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // fq: FASTA format quality value sequence input file
  ArgStruct argFq;
  arg = argFq; 
  arg.shortId = ""; 
  arg.longId = "fq"; 
  arg.description = "FASTA format quality value sequence input file";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_fq(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // gig: GIG format output file
  ArgStruct argGig;
  arg = argGig; 
  arg.shortId = ""; 
  arg.longId = "gig"; 
  arg.description = "Output file (GIG binary format assembly file)";
  arg.required = true; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_gig(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

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

  // tmp: binary temp file
  ArgStruct argTmp;
  arg = argTmp; 
  arg.shortId = ""; 
  arg.longId = "tmp"; 
  arg.description = "Temp file (binary format)";
  arg.required = false; 
  arg.defaultValueString = "__" + ProgramName + ".tmp"; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_tmp(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // Q: default base quality value to be assigned if no quality file or entry
  ArgStruct argQ;
  arg = argQ;
  arg.shortId = "";
  arg.longId = "Q";
  arg.description = "Default base quality value assigned if no quality file or entry";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "short";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<short> cmd_Q(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

  // an: assembly name
  ArgStruct argAn;
  arg = argAce; 
  arg.shortId = ""; 
  arg.longId = "an"; 
  arg.description = "Assembly name";
  arg.required = false; 
  arg.defaultValueString = "ASSEMBLY"; 
  arg.type = "string"; 
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_an(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // R: number of reads to cache when reading sequence files (RAM-friendly)
  ArgStruct argR;
  arg = argR;
  arg.shortId = "";
  arg.longId = "R";
  arg.description = "Number of reads to cache when reading sequence files (RAM-friendly)";
  arg.required = false;
  arg.defaultValueString = "10000000";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_R(arg.shortId, arg.longId, arg.description, arg.required, 10000000, arg.type, cmd);

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

  string ace = cmd_ace.getValue();
  string fd = cmd_fd.getValue();
  string fq = cmd_fq.getValue();
  string gig = cmd_gig.getValue();
  string tmp = cmd_tmp.getValue();
  string log = cmd_log.getValue();
  short Q = cmd_Q.getValue();
  string an = cmd_an.getValue();
  int R = cmd_R.getValue();
  int I = cmd_I.getValue();
  bool debug = cmd_debug.getValue();
  bool debug2 = cmd_debug2.getValue();


  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // constants
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // variables
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  bool record = false;
  if (log != "") {
    record = true;
  }

  //----------------------------------------------------------------------------
  // regular expressions to match
  //----------------------------------------------------------------------------

  // Assembly info line: AS
  boost::regex patternAS("^AS\\s+(\\d+)\\s+(\\d+)");
   
  // Contig info line: CO
  boost::regex patternCO("^CO\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(U|C)");
   
  // Contig base quality line: BQ
  boost::regex patternBQ("^BQ\\s*$");
   
  // Contig assembled_from line: AF
  boost::regex patternAF("^AF\\s+(\\S+)\\s+(U|C)\\s+(\\-{0,1}\\d+)");
   
  // Contig base segment line: BS
  boost::regex patternBS("^BS\\s+(\\d+)\\s+(\\d+)\\s+(\\S+)");
   
  // Read info line: RD
  boost::regex patternRD("^RD\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)");
   
  // Read qual clip info line: QA
  boost::regex patternQA("^QA\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)");
   
  // Read description info line: DS (two alternative versions)
  boost::regex patternDS("^DS");


  // position tag start line: WR|CT
  boost::regex patternPosTagStart("^(RT|CT)\\{");

  // positon tag info
  boost::regex patternPosTagInfo("^(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S+)\\s*((?:)|NoTrans)$");

  // whole item tag start line: WR|WA
  boost::regex patternWholeTagStart("^(WR|WA)\\{");

  // whole item tag info
  boost::regex patternWholeTagInfo("^(\\S+)\\s+(\\S+)\\s+(\\S+)");

  // Read pos tag start line: RT
  boost::regex patternRT("^RT\\{");

  // Consensus pos tag start line: CT
  boost::regex patternCT("^CT\\{");

  // Whole assembly tag start line: WA
  boost::regex patternWA("^WA\\{");

  // Tag closing line: }
  boost::regex patternCurly("^\\}");


  // empty line
  boost::regex patternEmpty("^\\s*$");

  // white space separator pattern for QUAL values
  boost::regex re("\\s+");

  // Patterns to match
  boost::regex patternFastaHeader("^>(\\s*\\S+\\s*.*)$");
  boost::regex patternFastaName("^\\s*(\\S+)");
  boost::regex patternQualHeader("^>\\s*(\\S+)");
  
  // token iterator constant
  boost::sregex_token_iterator j;

  // whatever was matched
  boost::smatch match, match2, what;

  // read offset hash for entire assembly
  hash_map<string, off64_t> readOffsetHash;

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
    logFile << "  --ace = " << ace << endl;
    logFile << "  --fd = " << fd << endl;
    logFile << "  --fq = " << fq << endl;
    logFile << "  --gig = " << gig << endl;
    logFile << "  --log = " << log << endl;
    logFile << "  --tmp = " << tmp << endl;
    logFile << "  --an = " << an << endl;
    logFile << "  --Q = " << Q << endl;
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
    cerr << "  --ace = " << ace << endl;
    cerr << "  --fd = " << fd << endl;
    cerr << "  --fq = " << fq << endl;
    cerr << "  --gig = " << gig << endl;
    cerr << "  --log = " << log << endl;
    cerr << "  --tmp = " << tmp << endl;
    cerr << "  --an = " << an << endl;
    cerr << "  --Q = " << Q << endl;
    cerr << "  --R = " << R << endl;
    cerr << "  --I = " << I << endl;
    cerr << "  --debug = " <<  bool2String[debug] << endl;
    cerr << "  --debug2 = " <<  bool2String[debug2] << endl;
    cerr << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // open input and output files
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // open input ACE file
  //----------------------------------------------------------------------------

  // report
  if (record) {cerr << "Opening ACE format assembly input file: " << ace << " ...";}
  if (debug) {cerr << "Opening ACE format assembly input file: " << ace << " ...";}
  
  // open
  ifstream aceIn1(ace.c_str(), ios::in);

  // check
  if (!aceIn1) {
    cerr << " unable to open file. Terminating." << endl;
    exit(1);
  }

  if (record) {cerr << " done." << endl;}
  if (debug) {cerr << " done." << endl;}

  //----------------------------------------------------------------------------
  // analyze what type of FASTA sequence information will be read
  //----------------------------------------------------------------------------
  
  // bools to indicate is sequences have dna, qual entries
  bool haveDna = false;
  bool haveQual = false;

  // set "have" flags
  if (fd != "") {
    haveDna = true;
  }
  if (fq != "") {
    haveQual = true;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // open output files
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // open temp binary file
  //----------------------------------------------------------------------------

  // report
  if (record) {logFile << "Opening binary temp file: " << tmp << " ...";}
  if (debug) {cerr << "Opening binary temp file: " << tmp << " ...";}

  // open
  FILE * tmpOut = fopen64(tmp.c_str(), "w+");
  if (!tmpOut) {
    cerr << "  unable to open file: " << tmp << endl;
    exit(1);
  }
							       
  if (record) {logFile << " done." << endl;}
  if (debug) {cerr << " done." << endl;}

  //----------------------------------------------------------------------------
  // open mandatory output GIG binary file
  //----------------------------------------------------------------------------

  // report
  if (record) {logFile << "Opening GIG format binary assembly output file: " << gig << " ...";}
  if (debug) {cerr << "Opening GIG format binary assembly output file: " << gig << " ...";}

  // open
  FILE * gigOut = fopen64(gig.c_str(), "w+");
  if (!gigOut) {
    cerr << "  unable to open file: " << gig << endl;
    exit(1);
  }
							       
  if (record) {logFile << " done." << endl;}
  if (debug) {cerr << " done." << endl;}

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // write archive id
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  int idLength = archiveId.length();

  // byte length of id
  fwrite((char*)&idLength, sizeof(int), 1, gigOut);

  // archive id
  fwrite(archiveId.c_str(), sizeof(char), (idLength+1), gigOut);
							       
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // read ACE file the first time: collect read alignment info from AF lines
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // report
  if (record) {logFile << "Reading ACE file for read alignment coordinates...";}
  if (debug) {cerr << "Reading ACE file for read alignment coordinates...";}

  // input line for parsing
  string line;

  //----------------------------------------------------------------------------
  // read ace file
  //----------------------------------------------------------------------------
  int numberAssembledReadsFound = 0;
  while (getline(aceIn1, line)) {

    if(boost::regex_search(line, match, patternAF)) {
      string name = match[1];
      string orient = match[2];
      int left = string2Int(match[3]);
      
      bool complemented = false;
      if (orient == "C") {
	complemented = true;
      }

      // write space holder for offset value corresponding to this read in GIG file
      off64_t offset = 0;
      fwrite((char*)&offset, sizeof(off64_t), 1, tmpOut);

      // write length of readName
      int nameLength = name.size();
      fwrite((char*)&nameLength, sizeof(int), 1, tmpOut);
      
      // write left position
      fwrite((char*)&left, sizeof(int), 1, tmpOut);

      // write complemented
      fwrite((char*)&complemented, sizeof(bool), 1, tmpOut);
      
      // write readName
      fwrite(name.c_str(), sizeof(char), (nameLength+1), tmpOut); 

      // increment counter
      numberAssembledReadsFound++;
    }
  }
  if (record) {logFile << " done." << endl;}
  if (debug) {cerr << " done." << endl;}

  // close ACE file
  aceIn1.close();

  // reposition to the beginning of the temp file
  fseek64(tmpOut, 0, SEEK_SET);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // read ACE file the second time: process and write to GIG
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // assembly-wide variables
  int numberAssemblyContigs;
  int numberAssemblyReads;      

  vector<string> contigs;

  // contig offsets
  hash_map<string, long long> contigOffset;

  // read counter
  int readCount = 0;
  
  // file offsets
  off64_t contigIndexOffset = 0;
  off64_t contigBlockOffset = 0;
  off64_t contigIndexOffsetLocus;
  off64_t contigBlockOffsetLocus;

  // open ACE file once again for second read-through
  ifstream aceIn(ace.c_str(), ios::in);

  // check
  if (!aceIn) {
    cerr << " unable to open file. Terminating." << endl;
    exit(1);
  }

  //----------------------------------------------------------------------------
  // read ace file line-by-line
  //----------------------------------------------------------------------------
  while (getline(aceIn, line)) {

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    // assembly
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------

    if(boost::regex_search(line, match, patternAS)) {


      //------------------------------------------------------------------------
      // report
      //------------------------------------------------------------------------
      if (record) {logFile << "Parsing assembly information...";}
      if (debug) {cerr << "Parsing assembly information...";}

      //------------------------------------------------------------------------
      // parse assembly header information
      //------------------------------------------------------------------------
      numberAssemblyContigs = string2Int(match[1]);
      numberAssemblyReads = string2Int(match[2]);
      if (record) {logFile << " done." << endl;}
      if (debug) {cerr << " done." << endl;}

      //------------------------------------------------------------------------
      // write assembly block into GIG file
      //------------------------------------------------------------------------

      // report
      if (record) {logFile << "Writing assembly information... ";}
      if (debug) {cerr << "Writing assembly information... ";}

      // write number of contigs
      fwrite((char*)&numberAssemblyContigs, sizeof(int), 1, gigOut);

      // register locus for contigIndexOffset
      contigIndexOffsetLocus = ftell64(gigOut);

      // write space holder for contigIndexOffset
      fwrite((char*)&contigIndexOffset, sizeof(off64_t), 1, gigOut);

      // register locus for contigBlockOffset
      contigBlockOffsetLocus = ftell64(gigOut);

      // write space holder for contigBlockOffset
      fwrite((char*)&contigBlockOffset, sizeof(off64_t), 1, gigOut);

      // write number of assembled reads
      fwrite((char*)&numberAssemblyReads, sizeof(int), 1, gigOut);

      // convert to c-style strings
      const char * anCstyle = an.c_str();
    
      // length of c-style strings
      int anCstyleLength = strlen(anCstyle);

      // write length of assemblyName
      fwrite((char*)&anCstyleLength, sizeof(int), 1, gigOut);

      // write assemblyName
      fwrite(anCstyle, sizeof(char), (anCstyleLength+1), gigOut);

      //------------------------------------------------------------------------
      // write true value for contigBlockOffset
      //------------------------------------------------------------------------

      // find out offset of contig index
      contigBlockOffset = ftell64(gigOut);
      
      // position back to the locus of the contigBlockOffset
      fseek64(gigOut, contigBlockOffsetLocus, SEEK_SET);
      
      // update the offset of the contig index
      fwrite((char*)&contigBlockOffset, sizeof(off64_t), 1, gigOut);
      
      // come back to the current position
      fseek64(gigOut, contigBlockOffset, SEEK_SET);

      if (record) {logFile << " done." << endl;}
      if (debug) {cerr << " done." << endl;}

      // report start of contig block
      if (record) {logFile << "Processing contigs." << endl;}
      if (debug) {cerr << "Processing contigs." << endl;}
    }

    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    // contig
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    else if (boost::regex_search(line, match, patternCO)) {

      //------------------------------------------------------------------------
      // parse contig header information
      //------------------------------------------------------------------------
      string contigName = match[1];
      int contigLengthPadded = string2Int(match[2]);
      int numberContigReads = string2Int(match[3]);
      int numberContigBaseSegments = string2Int(match[4]);
      string contigOrientation = match[5];

      //------------------------------------------------------------------------
      // report
      //------------------------------------------------------------------------
      if (record) {logFile << "  Contig: " << contigName << endl;}
      if (debug) {cerr << "  Contig: " << contigName << endl;}

      //------------------------------------------------------------------------
      // read contig dna
      //------------------------------------------------------------------------
      if (record) {logFile << "    Reading DNA...";}
      if (debug) {cerr << "    Reading DNA...";}

      string contigDnaPadded;
      string subLine("non-empty string");
      while (!boost::regex_search(subLine, match, patternEmpty)) {
	getline(aceIn, subLine);
	contigDnaPadded += subLine;
      }
      if (record) {logFile << " done. Padded length=" << contigDnaPadded.size() << endl;}
      if (debug) {cerr << " done. Padded length=" << contigDnaPadded.size() << endl;}

      //------------------------------------------------------------------------
      // read contig qual
      //------------------------------------------------------------------------
      if (record) {logFile << "    Reading QUAL values...";}
      if (debug) {cerr << "    Reading QUAL values...";}

      // read until BQ line is found
      subLine = "";
      while (!boost::regex_search(subLine, match, patternBQ)) {
	getline(aceIn, subLine);
      }

      // now read quality lines
      string contigQualString;
      while (!boost::regex_search(subLine, match, patternEmpty)) {
	getline(aceIn, subLine);
	contigQualString += " " + subLine;
      }

      // tokenize input line of base quality values and add to sequence
      boost::sregex_token_iterator i(contigQualString.begin(), contigQualString.end(), re, -1);
      vector<short> contigQualUnpadded;
      int first = true;
      int qCount = 0;
      while(i != j) {
	short q = string2Short(*i++);
	
	// skip first one (it is a parsing artifact)
	if (first) {
	  first = false;
	  continue;
	}
	qCount++;
	contigQualUnpadded.push_back(q);
      }
      if (record) {logFile << " done. Unpadded length=" << contigQualUnpadded.size() << endl;}
      if (debug) {cerr << " done. Unpadded length=" << contigQualUnpadded.size() << endl;}

      //------------------------------------------------------------------------
      // pad contig quality
      //------------------------------------------------------------------------
      if (record) {logFile << "    Padding QUAL values...";}
      if (debug) {cerr << "    Padding QUAL values...";}

      vector<short> contigQualPadded;
      int k = 0;
      for(string::const_iterator iter = contigDnaPadded.begin(); iter != contigDnaPadded.end(); iter++) {
	char b = *iter;
	if (b == '*') {
	  contigQualPadded.push_back((short)0);
	}
	else {
	  contigQualPadded.push_back(contigQualUnpadded[k]);
	  k++;
	}
      }
      if (record) {logFile << " done. Padded length=" << contigQualPadded.size() << endl;}
      if (debug) {cerr << " done. Padded length=" << contigQualPadded.size() << endl;}

      //------------------------------------------------------------------------
      // write constant-size contig info
      //------------------------------------------------------------------------

      if (record) {logFile << "    Writing contig information...";}
      if (debug) {cerr << "    Writing contig information...";}

      // convert to c-style strings
      const char * contigNameCstyle = contigName.c_str();
      const char * contigDnaPaddedCstyle = contigDnaPadded.c_str();
	
      // length of strings
      int contigNameLength = strlen(contigNameCstyle);

      // store offset for this contig
      contigOffset[contigName] = ftell64(gigOut);
	
      // add sequence
      contigs.push_back(contigName);

      // write number of reads assembled into this contig
      fwrite((char*)&numberContigReads, sizeof(int), 1, gigOut);
      
      // register locus for readIndexOffset
      off64_t readIndexOffsetLocus = ftell64(gigOut);

      // write space holder for offset of read index
      off64_t readIndexOffset = 0;
      fwrite((char*)&readIndexOffset, sizeof(off64_t), 1, gigOut);
     
      // register locus for readBlockOffset
      off64_t readBlockOffsetLocus = ftell64(gigOut);

      // write space holder for offset of read block
      off64_t readBlockOffset = 0;
      fwrite((char*)&readBlockOffset, sizeof(off64_t), 1, gigOut);
     
      // write length of contigName
      fwrite((char*)&contigNameLength, sizeof(int), 1, gigOut);
      
      // write unpadded contig length
      fwrite((char*)&contigLengthPadded, sizeof(int), 1, gigOut);
      
      // write number of base segments
      fwrite((char*)&numberContigBaseSegments, sizeof(int), 1, gigOut);
      
      //------------------------------------------------------------------------
      // write variable-size contig info
      //------------------------------------------------------------------------

      // write contigName
      fwrite(contigNameCstyle, sizeof(char), (contigNameLength+1), gigOut);
	
      // write contigDna
      fwrite(contigDnaPaddedCstyle, sizeof(char), (contigLengthPadded+1), gigOut);
	
      // write contig base quality
      for (int i=0; i<contigLengthPadded; i++) {
	short q = contigQualPadded[i];
	fwrite((char*)&q, sizeof(short), 1, gigOut);
      }

      if (record) {logFile << " done." << endl;}
      if (debug) {cerr << " done." << endl;}

      //------------------------------------------------------------------------
      // read and parse AF lines
      //------------------------------------------------------------------------

      if (record) {logFile << "    Reading assembled read alignment coordinates (" << numberContigReads << " entries)...";}
      if (debug) {cerr << "    Reading assembled read alignment coordinates (" << numberContigReads << " entries)...";}

      // read orientation within contig
      map<string, int, std::less<string> > contigReadOrientation;
      
      // read alignment left end within contig
      map<string, int, std::less<string> > contigReadLeft;
      
      for (int i=0; i<numberContigReads; i++) {

	// find and parse AF line
	string readName, orient;
	int readLeft;
	string subLine;
	while (!boost::regex_search(subLine, match, patternAF)) {
	  getline(aceIn, subLine);
	}
	readName = match[1];
	orient = match[2];
	readLeft = string2Int(match[3]);

	// register in contig-scope hashes
	if (orient == "C") {
	  contigReadOrientation[readName] = -1;
	}
	else {
	  contigReadOrientation[readName] = 1;
	}
	contigReadLeft[readName] = readLeft;
      }

      if (record) {logFile << " done." << endl;}
      if (debug) {cerr << " done." << endl;}

      //------------------------------------------------------------------------
      // parse and register BS lines
      //------------------------------------------------------------------------
      if (record) {logFile << "    Reading and writing base segment information (" << numberContigBaseSegments << " entries)...";}
      if (debug) {cerr << "    Reading and writing base segment information (" << numberContigBaseSegments << " entries)...";}

      for (int i=0; i<numberContigBaseSegments; i++) {

	//----------------------------------------------------------------------
	// find and parse BS line
	//----------------------------------------------------------------------
        string readName;
        int bsLeft, bsRight;
        string subLine;
        while (!boost::regex_search(subLine, match, patternBS)) {
          getline(aceIn, subLine);
        }
        bsLeft =  string2Int(match[1]);
        bsRight =  string2Int(match[2]);
        readName = match[3];

	//----------------------------------------------------------------------
	// write base segment
	//----------------------------------------------------------------------

	// write length of readName
	int readNameLength = readName.size();
	fwrite((char*)&readNameLength, sizeof(int), 1, gigOut);
      
	// write padded left position
	fwrite((char*)&bsLeft, sizeof(int), 1, gigOut);
	
	// write padded right position
	fwrite((char*)&bsRight, sizeof(int), 1, gigOut);
	
	// write readName
	fwrite(readName.c_str(), sizeof(char), (readNameLength+1), gigOut);
      }

      if (record) {logFile << " done." << endl;}
      if (debug) {cerr << " done." << endl;}

      //------------------------------------------------------------------------
      // write true value for readBlockOffset
      //------------------------------------------------------------------------

      // find out offset of read index
      readBlockOffset = ftell64(gigOut);
      
      // position back to the locus of the readBlockOffset
      fseek64(gigOut, readBlockOffsetLocus, SEEK_SET);
      
      // update the offset of the read index
      fwrite((char*)&readBlockOffset, sizeof(off64_t), 1, gigOut);
      
      // come back to the current position
      fseek64(gigOut, readBlockOffset, SEEK_SET);

      //------------------------------------------------------------------------
      // read and parse read info: RD and QA lines (ignore DS)
      //------------------------------------------------------------------------
      if (record) {logFile << "    Reading and writing assembled read information (" << numberContigReads << " entries)." << endl;}
      if (debug) {cerr << "    Reading and writing assembled read information (" << numberContigReads << " entries)." << endl;}

      vector<string> readList;
      int seqCount = 0;
      for (int i=0; i<numberContigReads; i++) {

	// find and parse RD line
	string subLine;
	while (!boost::regex_search(subLine, match, patternRD)) {
	  getline(aceIn, subLine);
	}
	string readName = match[1];
	int readLengthPadded = string2Int(match[2]);

	// retrieve info from AF
	int readOrient = contigReadOrientation[readName];
	int readLeft = contigReadLeft[readName];

	// read read dna
	string readDnaPadded;
	subLine = "non-empty string";
	while (!boost::regex_search(subLine, match, patternEmpty)) {
	  getline(aceIn, subLine);
	  readDnaPadded += subLine;
	}

	// find and parse QA line
	subLine = "";
	while (!boost::regex_search(subLine, match, patternQA)) {
	  getline(aceIn, subLine);
	}
	int readQualClipLeft = string2Int(match[1]);
	int readQualClipRight = string2Int(match[2]);
	int readAlignmentClipLeft = string2Int(match[3]);
	int readAlignmentClipRight = string2Int(match[4]);

	//----------------------------------------------------------------------
	// register offset for this read in tmp file
	//----------------------------------------------------------------------
	off64_t readOffset = ftell64(gigOut);
	fwrite((char*)&readOffset, sizeof(off64_t), 1, tmpOut);

	readList.push_back(readName);

	//----------------------------------------------------------------------
	// retrieve read-specific alignment info for this read from temp file
	//----------------------------------------------------------------------
	int readNameLengthTemp;
	fread((char*)&readNameLengthTemp, sizeof(int), 1, tmpOut);   
	
	int readLeftTemp;
	fread((char*)&readLeftTemp, sizeof(int), 1, tmpOut);   
	
	int readComplementedTemp;
	fread((char*)&readComplementedTemp, sizeof(bool), 1, tmpOut);   
	
	char readNameTemp[readNameLengthTemp];
	fread(readNameTemp, sizeof(char), readNameLengthTemp+1, tmpOut);
	
	//----------------------------------------------------------------------
	// calculate padded fit coordinates relative to contig
	//----------------------------------------------------------------------
	int contigFitLeftPadded = readLeftTemp;
	int contigFitRightPadded = readLeftTemp + readLengthPadded - 1;	

	//----------------------------------------------------------------------
	// write constant-size aligned read info
	//----------------------------------------------------------------------

	// convert to c-style strings
	const char * readNameCstyle = readName.c_str();
	
	// length of strings
	int readNameLength = strlen(readNameCstyle);
	
	// write length of readName
	fwrite((char*)&readNameLength, sizeof(int), 1, gigOut);
      
	// write padded read length
	fwrite((char*)&readLengthPadded, sizeof(int), 1, gigOut);
	
	// write padded contig fit left coordinate
	fwrite((char*)&contigFitLeftPadded, sizeof(int), 1, gigOut);
	
	// write padded contig fit right coordinate
	fwrite((char*)&contigFitRightPadded, sizeof(int), 1, gigOut);
	
	// write padded read quality clip left coordinate
	fwrite((char*)&readQualClipLeft, sizeof(int), 1, gigOut);
	
	// write padded read quality clip right coordinate
	fwrite((char*)&readQualClipRight, sizeof(int), 1, gigOut);
	
	// write padded read alignment clip left coordinate
	fwrite((char*)&readAlignmentClipLeft, sizeof(int), 1, gigOut);
	
	// write padded read alignment clip right coordinate
	fwrite((char*)&readAlignmentClipRight, sizeof(int), 1, gigOut);
	
	// write if alignment is complemented
	bool complemented = false;
	if (readOrient == -1) {
	  complemented = true;
	}
	fwrite((char*)&readComplementedTemp, sizeof(bool), 1, gigOut);
      
	//----------------------------------------------------------------------
	// write variable-size aligned read info
	//----------------------------------------------------------------------

	// write readName
	fwrite(readNameCstyle, sizeof(char), (readNameLength+1), gigOut);
	
	// write dna
	fwrite(readDnaPadded.c_str(), sizeof(char), (readLengthPadded+1), gigOut);

	// write place-holder values for qual
	for (int i=0; i<readLengthPadded; i++) {
	  fwrite((char*)&Q, sizeof(short), 1, gigOut);
	}

	seqCount++;
	if (record && (seqCount % I == 0)) {logFile << "      Read and wrote read: " << seqCount << " of " << numberContigReads << " reads." << endl;}
	if (debug && (seqCount % I == 0)) {cerr << "      Read and wrote read: " << seqCount << " of " << numberContigReads << " reads." << endl;}
      }
      
      if (record) {logFile << "    Number of reads read and written: " << seqCount << " of " << numberContigReads << "." << endl;}
      if (debug) {cerr << "    Number of reads read and written: " << seqCount << " of " << numberContigReads << "." << endl;}

      //------------------------------------------------------------------------
      // write true value for readIndexOffset
      //------------------------------------------------------------------------

      // find out offset of read index
      readIndexOffset = ftell64(gigOut);
      
      // position back to the locus of the readIndexOffset
      fseek64(gigOut, readIndexOffsetLocus, SEEK_SET);
      
      // update the offset of the read index
      fwrite((char*)&readIndexOffset, sizeof(off64_t), 1, gigOut);

      // come back to the current position
      fseek64(gigOut, readIndexOffset, SEEK_SET);

      //------------------------------------------------------------------------
      // write index for reads assembled into this contig 
      //------------------------------------------------------------------------
      if (record) {logFile << "    Writing assembled read index ...";}
      if (debug) {cerr << "    Writing assembled read index ...";}

      for (vector<string>::const_iterator readIter = readList.begin();
	   readIter != readList.end(); readIter++) {
	
	// retreive readName
	string readName = *readIter;

	// retreive corresponding file offset
	off64_t readOffset = readOffsetHash[readName];

	// convert to c-style strings
	const char * readNameCstyle = readName.c_str();
    
	// length of c-style string
	int readNameCstyleLength = strlen(readNameCstyle);

	// write length of readName
	fwrite((char*)&readNameCstyleLength, sizeof(int), 1, gigOut);

	// write readName
	fwrite(readNameCstyle, sizeof(char), (readNameCstyleLength+1), gigOut);

	// write offset for this read
	fwrite((char*)&readOffset, sizeof(off64_t), 1, gigOut);
      }

      if (record) {logFile << " done." << endl;}
      if (debug) {cerr << " done." << endl;}
    }
  }

  //----------------------------------------------------------------------------
  // write true value for contigIndexOffset
  //----------------------------------------------------------------------------

  // report
  if (record) {logFile << "Writing contig index...";}
  if (debug) {cerr << "Writing contig index...";}

  // find out offset of contig index
  contigIndexOffset = ftell64(gigOut);

  // position back to the assembly
  fseek64(gigOut, contigIndexOffsetLocus, SEEK_SET);

  // update the offset of the contig index
  fwrite((char*)&contigIndexOffset, sizeof(off64_t), 1, gigOut);

  // come back to the current position
  fseek64(gigOut, contigIndexOffset, SEEK_SET);

  //----------------------------------------------------------------------------
  // write contig index into GIG file
  //----------------------------------------------------------------------------
  for (vector<string>::const_iterator iter = contigs.begin(); 
       iter != contigs.end();
       iter++) {

    // retrieve sequence name and file offset
    string contigName = * iter;
    off64_t co = contigOffset[contigName];

    // convert to c-style strings
    const char * contigNameCstyle = contigName.c_str();
    
    // length of c-style strings
    int contigNameCstyleLength = strlen(contigNameCstyle);

    // write length of contigName
    fwrite((char*)&contigNameCstyleLength, sizeof(int), 1, gigOut);

    // write contigName
    fwrite(contigNameCstyle, sizeof(char), (contigNameCstyleLength+1), gigOut);

    // write offset for this contig
    fwrite((char*)&co, sizeof(off64_t), 1, gigOut);
  }

  // report
  if (record) {logFile << " done." << endl;}
  if (debug) {cerr << " done." << endl;}

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // read assembled reads from TMP file. make batches of reads of size R.
  // search fasta dna and qual files for each batch. update dna and qual in
  // GIG file
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  if (record) {logFile << "Scanning DNA and QUAL sequence files for assembled reads. Batch size=" << R << " reads." << endl;}
  if (debug) {cerr << "Scanning DNA and QUAL sequence files for assembled reads. Batch size=" << R << " reads." << endl;}

  // reposition to the beginning of the temp file
  fseek64(tmpOut, 0, SEEK_SET);

  // clear readOffsetHash
  readOffsetHash.clear();

  // initialize sequence conter and batch counter
  int r = 0;
  int b = 0;
  int totalDnaFoundCount = 0;
  int totalDnaWrittenCount = 0;
  int totalQualFoundCount = 0;
  int totalQualWrittenCount = 0;

  // read every assembled read entry to fill read offset hash
  for (int i=1; i<=numberAssembledReadsFound; i++) {

    off64_t readOffsetTemp;
    fread((char*)&readOffsetTemp, sizeof(off64_t), 1, tmpOut);   

    int readNameLengthTemp;
    fread((char*)&readNameLengthTemp, sizeof(int), 1, tmpOut);   
	
    int readLeftTemp;
    fread((char*)&readLeftTemp, sizeof(int), 1, tmpOut);   
	
    int readComplementedTemp;
    fread((char*)&readComplementedTemp, sizeof(bool), 1, tmpOut);   
	
    char readNameTemp[readNameLengthTemp];
    fread(readNameTemp, sizeof(char), readNameLengthTemp+1, tmpOut);
	
    // register read offset in hash
    readOffsetHash[readNameTemp] = readOffsetTemp;

    // increment sequence counter
    r++;

    // if batch size reached or last sequence, process sequence files
    if ((r >= R) || (i == numberAssembledReadsFound)) {
  
      // increment batch counter
      b++;

      if (record) {logFile << "  Read batch: " << b << ". Scanning for the next " << R << " reads from a total of: " << numberAssembledReadsFound << endl;}
      if (debug) {cerr << "  Read batch: " << b << ". Scanning for the next " << R << " reads from a total of: " << numberAssembledReadsFound << endl;}

      // sort reads in ascending order of file offsets
      vector<string> readNames = sortHashMapKeysByValue(readOffsetHash, false);
    
      //------------------------------------------------------------------------
      // read fasta format dna file and fill dna for assembled reads
      //------------------------------------------------------------------------

      if (haveDna) {

	// open fasta dna file
	ifstream dnaIn(fd.c_str(), ios::in);

	// report
	if (record) {logFile << "    Opening FASTA dna sequence input file: " << fd << " ...";}
	if (debug) {cerr << "    Opening FASTA dna sequence input file: " << fd << " ...";}
    
	if (!dnaIn) {
	  cerr << " unable to open file. Terminating." << endl;
	  exit(1);
	}
    
	if (record) {logFile << " done." << endl;}
	if (debug) {cerr << " done." << endl;}

	//----------------------------------------------------------------------
	// variables
	//----------------------------------------------------------------------

	// unregistered flag to indicate if dna still has to be registered with sequence
	bool unregistered = false;

	// sequence counter
	int seqCount = 0;
	int foundCount = 0;

	// sequence header
	string seqHeader;

	// sequence name
	string seqName;

	// sequence dna
	string seqDna;

	// hash containing DNA sequences
	hash_map<string, string> readDnaHash;

	// input line for parsing
	string line;

	//----------------------------------------------------------------------
	// read DNA sequences from FASTA file
	//----------------------------------------------------------------------
	
	if (record) {logFile << "    Scanning FASTA format DNA file: " << fd << "." << endl;}
	if (debug) {cerr << "    Scanning FASTA format DNA file: " << fd << "." << endl;}
	
	// read through FASTA DNA file
	while (getline(dnaIn, line)) {
    
	  // header line (long format): register previous sequence and start new
	  if (boost::regex_search(line, match, patternFastaHeader)) {
      
	    // if previous sequence info not yet registred, register it
	    if (unregistered) {
	
	      if (record && (seqCount % I == 0)) {logFile << "      Scanning DNA entry: " << seqCount << ". So far, found: " << totalDnaFoundCount << " of " << numberAssemblyReads << " reads." << endl;}
	      if (debug && (seqCount % I == 0)) {cerr << "      Scanning DNA entry: " << seqCount << ". So far, found: " << totalDnaFoundCount << " of " << numberAssemblyReads << " reads." << endl;}
	      
	      // register seqDna
	      if (readOffsetHash.count(seqName) > 0) {
		readDnaHash[seqName] = seqDna;
		foundCount++;
		totalDnaFoundCount++;
	      }
	    }
	    
	    // increment read count
	    seqCount++;
	    
	    // reset seqDna
	    seqDna = "";
	    
	    // retreive info for new sequence
	    seqHeader = match[1];
	    
	    // parse out sequence name
	    seqName = seqHeader;
	    if (boost::regex_search(seqHeader, match2, patternFastaName)) {     
	      seqName = match2[1];
	    }
	    
	    // set unregistered read info flag
	    unregistered = true;
	  }   
	  else {
	    seqDna += line;
	  }
	}
	
	// register last DNA entry if necessary
	if (unregistered) {
	  
	  if (record && (seqCount % I == 0)) {logFile << "      Scanning DNA entry: " << seqCount << ". So far, found: " << totalDnaFoundCount << " of " << numberAssemblyReads << " reads." << endl;}
	  if (debug && (seqCount % I == 0)) {cerr << "      Scanning DNA entry: " << seqCount << ". So far, found: " << totalDnaFoundCount << " of " << numberAssemblyReads << " reads." << endl;}
	  
	  // register seqDna
	  if (readOffsetHash.count(seqName) > 0) {	
	    readDnaHash[seqName] = seqDna;
	    foundCount++;
	    totalDnaFoundCount++;
	  }
	}
	
	//----------------------------------------------------------------------
	// write DNA sequences into GIG file
	//----------------------------------------------------------------------
	
	// report
	if (record) {logFile << "    Writing DNA sequences into binary assembly file." << endl;}
	if (debug) {cerr << "    Writing DNA sequences into binary assembly file." << endl;}

	// reset sequence counter
	seqCount = 0;
	
	for (vector<string>::const_iterator iter = readNames.begin();
	     iter != readNames.end(); iter++) {
	  string readName = *iter;
	  off64_t readOffset = readOffsetHash[readName];
	  
	  if (readDnaHash.count(readName) > 0) {
	    
	    string seqDna = readDnaHash[readName];
	    
	    // analyze dna
	    int seqLength = seqDna.size();
	    
	    // retreive offset
	    off64_t readOffset = readOffsetHash[readName];
	    
	    // report 
	    if (debug2) {cerr << "FOUND readName=" << readName << " GIG offset=" << readOffset << endl;}
	    
	    // seek to position of this read
	    fseek64(gigOut, readOffset, SEEK_SET);
	    
	    // retrieve readNameLength
	    int readNameLength;
	    fread((char*)&readNameLength, sizeof(int), 1, gigOut);   
	    
	    // retrieve readLength
	    int readLength;
	    fread((char*)&readLength, sizeof(int), 1, gigOut);
	    
	    // read 6 integers that are irrelevant here
	    int dummyInt;
	    fread((char*)&dummyInt, sizeof(int), 1, gigOut);
	    fread((char*)&dummyInt, sizeof(int), 1, gigOut);
	    fread((char*)&dummyInt, sizeof(int), 1, gigOut);
	    fread((char*)&dummyInt, sizeof(int), 1, gigOut);
	    fread((char*)&dummyInt, sizeof(int), 1, gigOut);
	    fread((char*)&dummyInt, sizeof(int), 1, gigOut);
	    
	    // read orientation
	    bool readComplemented;
	    fread((char*)&readComplemented, sizeof(bool), 1, gigOut);
	    
	    // reposition to where readDna starts
	    off64_t readDnaOffsetRel = 8 * sizeof(int) + 1 * sizeof(bool) + (readNameLength+1) * sizeof(char);
	    fseek64(gigOut, readOffset + readDnaOffsetRel, SEEK_SET);
	    
	    // read padded dna (including C style terminator to be read but ignored)
	    string readDnaPadded;
	    char b;
	    for (int i=0; i<readLength; i++) {
	      fread(&b, sizeof(char), 1, gigOut);
	      readDnaPadded += b;
	    }
	    fread(&b, sizeof(char), 1, gigOut);
	    
	    string readDnaUnpadded = unpadDna(readDnaPadded);
	    int readLengthUnpadded = readDnaUnpadded.size();
	    
	    // reverse complement if necessary
	    if (readComplemented) {
	      seqDna = revCompDna(seqDna);
	    }
	    
	    // check lengths
	    if (seqLength != readLengthUnpadded) {
	      cerr << "ERROR: inconsistent unpadded read lengths for read=" << readName << " ace=" << readLengthUnpadded << " fasta=" << seqLength << ". terminating..." << endl;
	      exit(1);
	    }
	    
	    // pad seq dna
	    string seqDnaPadded;
	    int k = 0;
	    for(string::const_iterator iter = readDnaPadded.begin(); iter != readDnaPadded.end(); iter++) {
	      char b = *iter;
	      if (b == '*') {
		seqDnaPadded += '*';
	      }
	      else {
		seqDnaPadded += seqDna[k];
		k++;
	      }
	    }
	    
	    // reposition again to readDna
	    fseek64(gigOut, readOffset + readDnaOffsetRel, SEEK_SET);
	    
	    // write padded fasta dna
	    fwrite(seqDnaPadded.c_str(), sizeof(char), (readLength+1), gigOut);
	    
	    // increment seqCount
	    totalDnaWrittenCount++;
	    
	    if (record && (totalDnaWrittenCount % I == 0)) {logFile << "      Written DNA sequence: " << totalDnaWrittenCount << " of " << numberAssemblyReads << " total entries." << endl;}
	    if (debug && (totalDnaWrittenCount % I == 0)) {cerr << "      Written DNA sequence: " << totalDnaWrittenCount << " of " << numberAssemblyReads << " total entries." << endl;}
	  }
	}
	
	// close input file
	dnaIn.close();
      }

      //------------------------------------------------------------------------
      // read fasta qual file and fill qual for each assembled read
      //------------------------------------------------------------------------

      if (haveQual) {

	// open
	ifstream qualIn(fq.c_str(), ios::in);

	// report
	if (record) {logFile << "    Opening FASTA format quality sequence input file: " << fq << " ...";}
	if (debug) {cerr << "    Opening FASTA format quality sequence input file: " << fq << " ...";}
	
	if (!qualIn) {
	  cerr << " unable to open file. Terminating." << endl;
	  exit(1);
	}
	
	if (record) {logFile << " done." << endl;}
	if (debug) {cerr << " done." << endl;}
	
	// unregistered flag to indicate if dna still has to be registered with sequence
	bool unregistered = false;
	
	// sequence counter
	int seqCount = 0;
	int foundCount = 0;
	
	// sequence name
	string seqName;
	
	// sequence qual string
	string seqQualString;
	
	// hash containing QUAL sequences
	hash_map<string, vector<short> > readQualHash;
	
	if (record) {logFile << "    Scanning FASTA format QUAL file: " << fd << "." << endl;}
	if (debug) {cerr << "    Scanning FASTA format QUAL file: " << fd << "." << endl;}

	// input line for parsing
	string line;
	
	while (getline(qualIn, line)) {
	  
	  // header line (long format): register previous sequence and start new
	  if (boost::regex_search(line, match, patternQualHeader)) {
	    
	    // if previous sequence info not yet registred, register it
	    if (unregistered) {
	      
	      // report if necessary
	      if (record && (seqCount % I == 0)) {logFile << "      Scanning QUAL entry: " << seqCount << ". So far, found: " << totalQualFoundCount << " of " << numberAssemblyReads << " reads." << endl;}
	      if (debug && (seqCount % I == 0)) {cerr << "      Scanning QUAL entry: " << seqCount << ". So far, found: " << totalQualFoundCount << " of " << numberAssemblyReads << " reads." << endl;}

	      // register seqQual
	      if (readOffsetHash.count(seqName) > 0) {
		boost::sregex_token_iterator i(seqQualString.begin(), seqQualString.end(), re, -1);
		vector<short> seqQual;
		while(i != j) {
		  seqQual.push_back(string2Short(*i++));
		}
		readQualHash[seqName] = seqQual;
		foundCount++;
		totalQualFoundCount++;
	      }
	    }
	    
	    // retreive info for new seuqence
	    seqName = match[1];
	    
	    // increment seq count
	    seqCount++;
	    
	    // reset QUAL
	    seqQualString = "";
	    
	    // set unregistered read info flag
	    unregistered = true;
	  }   
	  else {
	    
	    // add qual value string
	    if (seqQualString != "") {
	      seqQualString += " ";
	    }
	    seqQualString += line;
	  }  
	}
	
	// register last QUAL entry if necessary
	if (unregistered) {
	  
	  // report if necessary
	  if (record && (seqCount % I == 0)) {logFile << "      Scanning QUAL entry: " << seqCount << ". So far, found: " << totalQualFoundCount << " of " << numberAssemblyReads << " reads." << endl;}
	  if (debug && (seqCount % I == 0)) {cerr << "      Scanning QUAL entry: " << seqCount << ". So far, found: " << totalQualFoundCount << " of " << numberAssemblyReads << " reads." << endl;}
	  
	  // write seqQual
	  if (readOffsetHash.count(seqName) > 0) {
	    boost::sregex_token_iterator i(seqQualString.begin(), seqQualString.end(), re, -1);
	    vector<short> seqQual;
	    while(i != j) {
	      seqQual.push_back(string2Short(*i++));
	    }
	    readQualHash[seqName] = seqQual;
	    foundCount++;
	    totalQualFoundCount++;
	  }
	}
	
	//----------------------------------------------------------------------
	// write QUAL sequences into GIG file
	//----------------------------------------------------------------------
	
	// report
	if (record) {logFile << "    Writing QUAL sequences into binary assembly file." << endl;}
	if (debug) {cerr << "    Writing QUAL sequences into binary assembly file." << endl;}
	
	// reset sequence counter
	seqCount = 0;
	
	for (vector<string>::const_iterator iter = readNames.begin();
	     iter != readNames.end(); iter++) {
	  string readName = *iter;
	  off64_t readOffset = readOffsetHash[readName];
	  
	  if (readQualHash.count(readName) > 0) {
	    
	    vector<short> seqQual = readQualHash[readName];
	    
	    // analyze dna
	    int seqLength = seqQual.size();
	    
	    // retreive offset
	    off64_t readOffset = readOffsetHash[readName];
	    
	    // report 
	    if (debug2) {cerr << "FOUND readName=" << readName << " GIG offset=" << readOffset << endl;}
	    
	    // seek to position of this read
	    fseek64(gigOut, readOffset, SEEK_SET);
	    
	    // retrieve readNameLength
	    int readNameLength;
	    fread((char*)&readNameLength, sizeof(int), 1, gigOut);   
	    
	    // retrieve readLength
	    int readLength;
	    fread((char*)&readLength, sizeof(int), 1, gigOut);
	    
	    // read 6 integers that are irrelevant here
	    int dummyInt;
	    fread((char*)&dummyInt, sizeof(int), 1, gigOut);
	    fread((char*)&dummyInt, sizeof(int), 1, gigOut);
	    fread((char*)&dummyInt, sizeof(int), 1, gigOut);
	    fread((char*)&dummyInt, sizeof(int), 1, gigOut);
	    fread((char*)&dummyInt, sizeof(int), 1, gigOut);
	    fread((char*)&dummyInt, sizeof(int), 1, gigOut);
	    
	    // read orientation
	    bool readComplemented;
	    fread((char*)&readComplemented, sizeof(bool), 1, gigOut);
	    
	    // reposition to where readDna starts
	    off64_t readDnaOffsetRel = 8 * sizeof(int) + 1 * sizeof(bool) + (readNameLength+1) * sizeof(char);
	    fseek64(gigOut, readOffset + readDnaOffsetRel, SEEK_SET);
	    
	    // read padded dna (including C style terminator to be read but ignored)
	    string readDnaPadded;
	    char b;
	    for (int i=0; i<readLength; i++) {
	      fread(&b, sizeof(char), 1, gigOut);
	      readDnaPadded += b;
	    }
	    fread(&b, sizeof(char), 1, gigOut);
	    
	    string readDnaUnpadded = unpadDna(readDnaPadded);
	    int readLengthUnpadded = readDnaUnpadded.size();
	    
	    // reverse complement if necessary
	    if (readComplemented) {
	      reverse(seqQual.begin(), seqQual.end());
	    }
	    
	    // check lengths
	    if (seqLength != readLengthUnpadded) {
	      cerr << "ERROR: inconsistent unpadded read lengths for read=" << readName << " ace=" << readLengthUnpadded << " fasta=" << seqLength << ". terminating..." << endl;
	      exit(1);
	    }
	    
	    // pad seq qual with 0s
	    string seqQualPadded;
	    int k = 0;
	    for(string::const_iterator iter = readDnaPadded.begin(); iter != readDnaPadded.end(); iter++) {
	      char b = *iter;
	      if (b == '*') {
		seqQualPadded.push_back((short)0);
	      }
	      else {
		seqQualPadded.push_back(seqQual[k]);
		k++;
	      }
	    }

	    // apply quality on the left of the pad
	    short leftQ = 0;
	    for (int i=0; i<readLength; i++) {
	      char b = readDnaPadded[i];
	      short q = seqQualPadded[i];
	      if (b == '*') {
		seqQualPadded[i] = leftQ;
	      }
	      else {
		leftQ = q;
	      }
	    }
	    
	    // apply quality on the right of the pad
	    short rightQ = 0;
	    for (int i=readLength-1; i>=0; i--) {
	      char b = readDnaPadded[i];
	      short q = seqQualPadded[i];
	      if (b == '*') {
		if (rightQ < q) {
		  seqQualPadded[i] = rightQ;
		}
	      }
	      else {
		rightQ = q;
	      }
	    }
	    
	    // no need to reposition in GIG file becuase we are exactly at the right offset for readQual
	    
	    // write qual
	    for (int i=0; i<readLength; i++) {
	      short q = seqQualPadded[i];
	      fwrite((char*)&q, sizeof(short), 1, gigOut);
	    }
	    
	    // increment seqCount
	    totalQualWrittenCount++;
	    
	    if (record && (totalQualWrittenCount % I == 0)) {logFile << "  Written QUAL sequence: " << totalQualWrittenCount << " of " << numberAssemblyReads << " total entries." << endl;}
	    if (debug && (totalQualWrittenCount % I == 0)) {cerr << "  Written QUAL sequence: " << totalQualWrittenCount << " of " << numberAssemblyReads << " total entries." << endl;}
	  }
	}
	
	// close input file
	qualIn.close();
      }

      // reset batch sequence counter
      r = 0;

      // clear readOffsetHash
      readOffsetHash.clear();
    }
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // stats
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  if (record) {logFile << "Found and written " << totalDnaFoundCount << " of " << numberAssemblyReads << " DNA entries." << endl;}
  if (debug) {cerr << "Found and written " << totalDnaFoundCount << " of " << numberAssemblyReads << " DNA entries." << endl;}
  if (record) {logFile << "Found and written " << totalQualFoundCount << " of " << numberAssemblyReads << " QUAL entries." << endl;}
  if (debug) {cerr << "Found and written " << totalQualFoundCount << " of " << numberAssemblyReads << " QUAL entries." << endl;}
	
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // finish up
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // close output files
  //----------------------------------------------------------------------------

  // report
  if (record) {logFile << "Closing files....";}
  if (debug) {cerr << "Closing files....";}
  
  // close files
  fclose(tmpOut);
  fclose(gigOut);
  aceIn.close();

  if (record) {logFile << " done." << endl;}
  if (debug) {cerr << " done." << endl;}

  //----------------------------------------------------------------------------
  // report
  //----------------------------------------------------------------------------
  if (record) {logFile << "Program completed." << endl;}
  if (debug) {cerr << "Program completed." << endl;}
}

