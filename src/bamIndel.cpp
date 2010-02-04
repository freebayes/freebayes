//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// bamIndel
// INDEL detector based on BAM format sequence alignments
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

// "hash_map" true hashes
#include <ext/hash_map>

// private libraries
#include "Function-Sequence.h"
#include "Function-Generic.h"
#include "Function-Math.h"
#include "Class-BedReader.h"
#include "Class-FastaReader.h"
#include "BamReader.h"
#include "ReferenceSequenceReader.h"

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
//------------------------------------------------------------------------------
// static variables
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

static string ProgramName("bamIndel");
static string ProgramDescription("Short-INDEL polymorphism discovery program.");
static string ProgramVersion("0.0.1");
static string ProgramDate("2009-09-01");
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

  // input file: BAM read alignments
  ArgStruct argBam;
  arg = argBam; 
  arg.shortId = ""; 
  arg.longId = "bam"; 
  arg.description = "Read alignment input file (indexed BAM format)";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_bam(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // input file: MOSAIK binary reference sequence archive
  ArgStruct argMbr;
  arg = argMbr; 
  arg.shortId = ""; 
  arg.longId = "mbr"; 
  arg.description = "Reference sequence archive file (MOSAIK binary format)";
  arg.required = true; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_mbr(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

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

  // input file: list of samples to analyze
  ArgStruct argSamples;
  arg = argSamples; 
  arg.shortId = ""; 
  arg.longId = "samples"; 
  arg.description = "File containing list of samples to analyze";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_samples(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // output file: tab-delimited alignment and SNP report text file
  ArgStruct argRpt;
  arg = argRpt; 
  arg.shortId = ""; 
  arg.longId = "rpt"; 
  arg.description = "Output file: tab-delimited alignment and SNP report text file";
  arg.required = true; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_rpt(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

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

  // MQR: reference sequence mapping quality value
  ArgStruct argMQR;
  arg = argMQR;
  arg.shortId = "";
  arg.longId = "MQR";
  arg.description = "Reference sequence mapping quality value";
  arg.required = false;
  arg.defaultValueString = "100";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_MQR(arg.shortId, arg.longId, arg.description, arg.required, 100, arg.type, cmd);

  // BQR: reference base quality value
  ArgStruct argBQR;
  arg = argBQR;
  arg.shortId = "";
  arg.longId = "BQR";
  arg.description = "Reference sequence base quality value";
  arg.required = false;
  arg.defaultValueString = "60";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_BQR(arg.shortId, arg.longId, arg.description, arg.required, 60, arg.type, cmd);

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
  ArgStruct argSampleNaming;
  arg = argSampleNaming;
  arg.shortId = "";
  arg.longId = "sampleNaming";
  arg.description = "Naming scheme for matching reads to samples";
  arg.required = false;
  arg.defaultValueString = "multiple";
  arg.type = "string";
  arg.multi = false;
  vector<string> allowedSampleNaming;
  allowedSampleNaming.push_back("single");
  allowedSampleNaming.push_back("multiple");
  allowedSampleNaming.push_back("trio");
  allowedSampleNaming.push_back("unknown");
  arg.constraint = allowedSampleNaming;
  ValuesConstraint<string> allowedSampleNamingVals(allowedSampleNaming); 
  ArgList.push_back(arg);
  ValueArg<string> cmd_sampleNaming(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, &allowedSampleNamingVals, cmd);

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

  // MQL0: minimum mapping quality value to consider read
  ArgStruct argMQL0;
  arg = argMQL0;
  arg.shortId = "";
  arg.longId = "MQL0";
  arg.description = "Minimum mapping quality value to consider read";
  arg.required = false;
  arg.defaultValueString = "40";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_MQL0(arg.shortId, arg.longId, arg.description, arg.required, 40, arg.type, cmd);

  // BQL0: minimum read base quality value
  ArgStruct argBQL0;
  arg = argBQL0;
  arg.shortId = "";
  arg.longId = "BQL0";
  arg.description = "Minimum read base quality value";
  arg.required = false;
  arg.defaultValueString = "10";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_BQL0(arg.shortId, arg.longId, arg.description, arg.required, 10, arg.type, cmd);

  // MQL1: minimum mapping quality value required for at least one read for each allele
  ArgStruct argMQL1;
  arg = argMQL1;
  arg.shortId = "";
  arg.longId = "MQL1";
  arg.description = "Minimum mapping quality value required for at least one read for each allele";
  arg.required = false;
  arg.defaultValueString = "40";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_MQL1(arg.shortId, arg.longId, arg.description, arg.required, 40, arg.type, cmd);

  // BQL1: minimum base quality value required for at least one base for each allele
  ArgStruct argBQL1;
  arg = argBQL1;
  arg.shortId = "";
  arg.longId = "BQL1";
  arg.description = "Minimum read base quality value required for at least one base for each allele";
  arg.required = false;
  arg.defaultValueString = "30";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_BQL1(arg.shortId, arg.longId, arg.description, arg.required, 30, arg.type, cmd);

  // BQL2: minimum mismatch base quality value for read mismatch filtering
  ArgStruct argBQL2;
  arg = argBQL2;
  arg.shortId = "";
  arg.longId = "BQL2";
  arg.description = "Minimum mismatch base quality value for read mismatch filtering";
  arg.required = false;
  arg.defaultValueString = "10";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_BQL2(arg.shortId, arg.longId, arg.description, arg.required, 10, arg.type, cmd);

  // RMU: maximum number of mismatches between read and reference sequence
  ArgStruct argRMU;
  arg = argRMU;
  arg.shortId = "";
  arg.longId = "RMU";
  arg.description = "Maximum number of mismatches between read and refrence sequence";
  arg.required = false;
  arg.defaultValueString = "10000000";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_RMU(arg.shortId, arg.longId, arg.description, arg.required, 10000000, arg.type, cmd);

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
  string mbr = cmd_mbr.getValue();
  string targets = cmd_targets.getValue();
  string samples = cmd_samples.getValue();
  string rpt = cmd_rpt.getValue();
  string log = cmd_log.getValue();

  bool useRefAllele = cmd_useRefAllele.getValue();
  bool forceRefAllele = cmd_forceRefAllele.getValue();
  int MQR = cmd_MQR.getValue();
  int BQR = cmd_BQR.getValue();
  string ploidy = cmd_ploidy.getValue();
  string sampleNaming = cmd_sampleNaming.getValue();
  string sampleDel = cmd_sampleDel.getValue();
  int MQL0 = cmd_MQL0.getValue();
  int BQL0 = cmd_BQL0.getValue();
  int MQL1 = cmd_MQL1.getValue();
  int BQL1 = cmd_BQL1.getValue();
  int BQL2 = cmd_BQL2.getValue();
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
  // derived variables
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  bool diploid = false;
  if (ploidy == "diploid") {
    diploid = true;
  }

  bool record = false;
  if (log != "") {
    record = true;
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
    logFile << "  --bam = " << bam << endl;
    logFile << "  --mbr = " << mbr << endl;
    logFile << "  --targets = " << targets << endl;
    logFile << "  --samples = " << samples << endl;
    logFile << "  --rpt = " << rpt << endl;
    logFile << "  --log = " << log << endl;
    logFile << "  --useRefAllele = " <<  bool2String[useRefAllele] << endl;
    logFile << "  --forceRefAllele = " <<  bool2String[forceRefAllele] << endl;
    logFile << "  --MQR = " << MQR << endl;
    logFile << "  --BQR = " << BQR << endl;
    logFile << "  --ploidy = " << ploidy << endl;
    logFile << "  --sampleNaming = " << sampleNaming << endl;
    logFile << "  --sampleDel = " << sampleDel << endl;
    logFile << "  --MQL0 = " << MQL0 << endl;
    logFile << "  --BQL0 = " << BQL0 << endl;
    logFile << "  --MQL1 = " << MQL1 << endl;
    logFile << "  --BQL1 = " << BQL1 << endl;
    logFile << "  --BQL2 = " << BQL2 << endl;
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
    cerr << "  --bam = " << bam << endl;
    cerr << "  --mbr = " << mbr << endl;
    cerr << "  --targets = " << targets << endl;
    cerr << "  --samples = " << samples << endl;
    cerr << "  --rpt = " << rpt << endl;
    cerr << "  --log = " << log << endl;
    cerr << "  --useRefAllele = " <<  bool2String[useRefAllele] << endl;
    cerr << "  --forceRefAllele = " <<  bool2String[forceRefAllele] << endl;
    cerr << "  --MQR = " << MQR << endl;
    cerr << "  --BQR = " << BQR << endl;
    cerr << "  --ploidy = " << ploidy << endl;
    cerr << "  --sampleNaming = " << sampleNaming << endl;
    cerr << "  --sampleDel = " << sampleDel << endl;
    cerr << "  --BQL0 = " << BQL0 << endl;
    cerr << "  --MQL0 = " << MQL0 << endl;
    cerr << "  --BQL1 = " << BQL1 << endl;
    cerr << "  --MQL1 = " << MQL1 << endl;
    cerr << "  --BQL2 = " << BQL2 << endl;
    cerr << "  --I = " << I << endl;
    cerr << "  --debug = " <<  bool2String[debug] << endl;
    cerr << "  --debug2 = " <<  bool2String[debug2] << endl;
    cerr << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // open output file(s)
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // open report output file
  //----------------------------------------------------------------------------

  // report
  if (record) {logFile << "opening report output file for writing: " << rpt << "...";}
  if (debug) {cerr << "opening report output file for writing: " << rpt << "...";}

  // open
  ofstream rptFile(rpt.c_str(), ios::out);
  if (!rptFile) {
    if (record) {logFile << " unable to open file: " << rpt << endl;}
    cerr << " unable to open file: " << rpt << endl;
    exit(1);
  }
  if (record) {logFile << " done." << endl;}
  if (debug) {cerr << " done." << endl;}

  //----------------------------------------------------------------------------
  // write header information
  //----------------------------------------------------------------------------
  rptFile << "# Command line that generated this output:";
  for (int i=0; i<argc; i++) {
    rptFile << " " << argv[i];
  }
  rptFile << endl;
  rptFile << "#" << endl;
  rptFile << "# Complete list of parameter values:" << endl;
  rptFile << "#   --bam = " << bam << endl;
  rptFile << "#   --mbr = " << mbr << endl;
  rptFile << "#   --targets = " << targets << endl;
  rptFile << "#   --samples = " << samples << endl;
  rptFile << "#   --rpt = " << rpt << endl;
  rptFile << "#   --log = " << log << endl;
  rptFile << "#   --useRefAllele = " <<  bool2String[useRefAllele] << endl;
  rptFile << "#   --forceRefAllele = " <<  bool2String[forceRefAllele] << endl;
  rptFile << "#   --MQR = " << MQR << endl;
  rptFile << "#   --BQR = " << BQR << endl;
  rptFile << "#   --ploidy = " << ploidy << endl;
  rptFile << "#   --sampleNaming = " << sampleNaming << endl;
  rptFile << "#   --sampleDel = " << sampleDel << endl;
  rptFile << "#   --BQL0 = " << BQL0 << endl;
  rptFile << "#   --MQL0 = " << MQL0 << endl;
  rptFile << "#   --BQL1 = " << BQL1 << endl;
  rptFile << "#   --MQL1 = " << MQL1 << endl;
  rptFile << "#   --BQL2 = " << BQL2 << endl;
  rptFile << "#   --I = " << I << endl;
  rptFile << "#   --debug = " <<  bool2String[debug] << endl;
  rptFile << "#   --debug2 = " <<  bool2String[debug2] << endl;
  rptFile << "#" << endl;

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // read sample list file
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  map<string, bool, less<string> > samplePresent;
  vector<string> sampleList;
  if (samples != "") {
    ifstream sampleFile(samples.c_str(), ios::in);
    if (! sampleFile) {
      cerr << "unable to open file: " << samples << endl;
      exit(1);
    }
    boost::regex patternSample("^(\\S+)\\s*(.*)$");
    boost::regex re("\\s+");
    boost::regex pr("^(\\S+):(\\S+)$");
    boost::smatch match;
    string line;
    while (getline(sampleFile, line)) {
      
      // if proper line
      if (boost::regex_search(line, match, patternSample)) {
	
	// assign content
	string s = match[1];
	samplePresent[s] = true;
	sampleList.push_back(s);
      }
    }
  }

  // get number of samples
  int numberSamples = sampleList.size();

  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // open BAM input file
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  
  // report
  if (record) {logFile << "Opening BAM fomat alignment input file: " << bam << " ...";}
  if (debug) {cerr << "Opening BAM format alignment input file: " << bam << " ...";}
  
  // open
  const char * bamFileNameC = bam.c_str();
  string bamIndexFileName = bam + ".bai";
  const char * bamIndexFileNameC = bamIndexFileName.c_str();
  BamReader bReader;
  bReader.Open(bamFileNameC, bamIndexFileNameC);
  
  // report
  if (record) {logFile << " done." << endl;}
  if (debug) {cerr << " done." << endl;}
  
  //--------------------------------------------------------------------------
  // read reference sequences from input file
  //--------------------------------------------------------------------------
  
  // retrieve header information
  string bamHeader = bReader.GetHeaderText();
  
  // store the names of all the reference sequences in the BAM file
  RefVector refDatas = bReader.GetReferenceData();
  
  // data structs
  map<string, int, less<string> > RefLength, RefId;
  for (RefVector::const_iterator refIter = refDatas.begin();
       refIter != refDatas.end(); refIter++) {
    RefData refData = *refIter;
    string refName = refData.RefName;
    RefLength[refName] = refData.RefLength;
    RefId[refName] = bReader.GetReferenceID(refName);
  }
  
  if (record) {
    logFile << "Number of ref seqs: " << bReader.GetReferenceCount() << endl;
  }
  if (debug) {
    cerr << "Number of ref seqs: " << bReader.GetReferenceCount() << endl;
  }
  
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // process MBR input file
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  
  //--------------------------------------------------------------------------
  // report
  //--------------------------------------------------------------------------
  if (record) {logFile << "Opening MOSAIK binary format reference sequence archive input file: " << mbr << " ...";}
  if (debug) {cerr << "Opening MOSAIK binary format reference sequence archive input file: " << mbr << " ...";}
  
  //--------------------------------------------------------------------------
  // if this is not a valid MOSAIK reference file, complain and bomb
  //--------------------------------------------------------------------------
  if (! Mosaik::CReferenceSequenceReader::CheckFile(mbr, false)) {
    cerr << "ERROR: MOSAIK reference file not valid: " << mbr << ". Exiting..." << endl;
    exit(1);
  }

  //--------------------------------------------------------------------------
  // open the MOSAIK alignments file
  //--------------------------------------------------------------------------
  Mosaik::CReferenceSequenceReader rsr;
  rsr.Open(mbr);

  //--------------------------------------------------------------------------
  // retrieve all of the reference sequence metadata
  //--------------------------------------------------------------------------
  vector<Mosaik::ReferenceSequence> referenceSequences;
  rsr.GetReferenceSequences(referenceSequences);

  //--------------------------------------------------------------------------
  // load ref seq names into hash
  //--------------------------------------------------------------------------
  map<string, unsigned int, less<string> > refseqLength;
  map<string, string, less<string> > refseqDna;
  for (vector<Mosaik::ReferenceSequence>::const_iterator rsIter = referenceSequences.begin();
       rsIter != referenceSequences.end(); ++rsIter) {
    Mosaik::ReferenceSequence rs = *rsIter;
    string ref = rs.Name;

    // get the reference base sequence
    string bases;
    rsr.GetReferenceSequence(ref, bases);

    // store base sequence
    refseqDna[ref] = bases;

    // store ref sequence length
    refseqLength[ref] = rs.NumBases;
  }


  // close the reference sequence reader
  rsr.Close();
  
  // report success
  if (record) {logFile << " done." << endl;}
  if (debug) {cerr << " done." << endl;}
  
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  // process target region input file if required
  //--------------------------------------------------------------------------
  //--------------------------------------------------------------------------
  
  // reference-specific target lists
  map<string, vector<BedData>, less<string> > RefTargetList;
  
  // if target file specified use targets from file
  int tc = 0;
  if (targets != "") {
    
    //------------------------------------------------------------------------
    // open input BED file if required
    //------------------------------------------------------------------------
    
    // report
    if (record) {logFile << "Making BedReader object for target file: " << targets << " ...";}
    if (debug) {cerr << "Making BedReader object for target file: " << targets << " ...";}
    
    // make
    BedReader * br = new BedReader(targets);
    
    if (! br->isOpen()) {
      cerr << "Unable to open target file: " << targets << "... terminating." << endl;
      exit(1);
    }
    
    //------------------------------------------------------------------------
    // iterate through entries
    //------------------------------------------------------------------------
    BedData bd;
    while (br->getNextEntry(bd)) {
      if (debug2) {
	cerr << bd.seq << "\t" << bd.left << "\t" << bd.right << "\t" << bd.desc << endl;
      }
      if (RefId.count(bd.seq) > 0) {
	if (bd.left < 1 || bd.right > RefLength[bd.seq]) {
	  cerr << "Target region coordinate outside of reference sequence bounds... terminating." << endl;
	  exit(1);
	}
	RefTargetList[bd.seq].push_back(bd);
	tc++;
      }
    }
    
    //------------------------------------------------------------------------
    // close
    //------------------------------------------------------------------------
    br->close();

    if (debug) {cerr << "done" << endl;}
    if (record) {logFile << "done" << endl;}
  }

  // otherwise analyze all reference sequences from BAM file
  else {
    RefVector::iterator refIter = refDatas.begin();
    RefVector::iterator refEnd  = refDatas.end();
    for( ; refIter != refEnd; ++refIter) {
      RefData refData = *refIter;
      string refName = refData.RefName;
      int refId = bReader.GetReferenceID(refName);
      int refLength = refData.RefLength;
      BedData bd;
      bd.seq = refName;
      bd.left = 1;
      bd.right = refLength;
      RefTargetList[bd.seq].push_back(bd);
      tc++;
    }
  }
  if (debug) {cerr << "Number of target regions: " << tc << endl;}
  if (record) {logFile << "Number of target regions: " << tc << endl;}
  
  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------
  // process target regions
  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------

  for (map<string, vector<BedData>, less<string> >::const_iterator refIter = RefTargetList.begin(); 
       refIter != RefTargetList.end(); refIter++) {
    string refName = refIter->first;
    vector<BedData> targetRegions = refIter->second;
    for(vector<BedData>::const_iterator targetIter = targetRegions.begin();
	targetIter != targetRegions.end(); targetIter++) {
      
      // retrieve target region
      BedData target = * targetIter;
      
      // find target ref seq in list of ref seqs from BAM file
      int refId = bReader.GetReferenceID(target.seq);
      RefData refData = refDatas[refId];
      int refLength = refData.RefLength;
      
      // report
      if (record) {logFile << "  Processing target: " << target.seq << ":" << target.left << "-" << target.right << endl;}
      if (debug) {cerr << "  Processing target: " << target.seq << ":" << target.left << "-" << target.right << endl;}

      // collect targets for which there is no read data
      bool badTarget = false;

      // skip target if invalid reference
      if (!(refId < bReader.GetReferenceCount())) {
	badTarget = true;
	if (record) {logFile << "      WARNING: Target ref ID invalid: " << refId << endl;}
	if (debug) {cerr << "      WARNING: Target ref ID invalid: " << refId << endl;}	
      }
      else {
 	if (record) {logFile << "      Target ref ID OK: " << refId << endl;}
	if (debug) {cerr << "      Target ref ID OK: " << refId << endl;}	
      }
      
      // skip target if ref seq does not have any alignments
      if (! refData.RefHasAlignments) {
	badTarget = true;
	if (record) {logFile << "      WARNING: Target has no alignments..." << endl;}
	if (debug) {cerr << "      WARNING: Target has no alignments..." << endl;}	
      }
      else {
	if (record) {logFile << "      Target has alignments..." << endl;}
	if (debug) {cerr << "      Target has alignments..." << endl;}	
      }
      
      if (record) {logFile << "      Jumping to target start in BAM file. RefId: " << refId << " pos: " << target.left << endl;}
      if (debug) {cerr << "      Jumping to target start in BAM file. RefId: " << refId << " pos: " << target.left << endl;}	
      if (! bReader.Jump(refId, target.left)) {
	badTarget = true;
	if (record) {logFile << "      WARNING: Cannot jump to target start in BAM file. RefId: " << refId << " pos: " << target.left << endl;}
	if (debug) {cerr << "      WARNING: Cannot jump to target start in BAM file. REfId: " << refId << " pos: " << target.left << endl;}	
      }
      else {
	if (record) {logFile << "      Jumped to target start in BAM file. RefId: " << refId << " pos: " << target.left << endl;}
	if (debug) {cerr << "      Jumped to target start in BAM file. RefId: " << refId << " pos: " << target.left << endl;}	
      }

      //------------------------------------------------------------------
      // if target is bad, print empty stats for TARGET and for all POSITIONS within
      //------------------------------------------------------------------
      if (badTarget) {
	if (record) {logFile << "    Target not OK... skipping" << endl;}
	if (debug) {cerr << "    Target not OK... skipping" << endl;}
      }
 
      // otherwise analyze target
      else {
	if (record) {logFile << "    Target OK... processing..." << endl;}
	if (debug) {cerr << "    Target OK... processing..." << endl;}
      
	// contig-position and individual specific coverage quantities
	map<unsigned int, map<string, vector<Basecall>, less<string> >, less<unsigned int> > individualBasecalls;
      
	int refPosLast = target.left - 1;
	BamAlignment ba;
	while (bReader.GetNextAlignment(ba) && (ba.Position + 1 <= target.right)) {
	  
	  // only process if mapped
	  if (! ba.IsMapped()) {
	    continue;
	  }
	  
	  // mark read if duplicate
	  bool dupRead = false;
	  if (ba.IsDuplicate()) {
	    dupRead = true;
	  }
	  
	  // mark read if mapping quality is low
	  bool lowMapQualRead = false;
	  if (ba.MapQuality < MQL0) {
	    lowMapQualRead = true;
	  }
	  
	  //----------------------------------------------------------------
	  // extract sample info
	  //----------------------------------------------------------------
	  string readName = ba.Name;
	  SampleInfo sampleInfo = extractSampleInfo(readName, sampleNaming, sampleDel);
	  string sampleName = sampleInfo.sampleId;
	  
	  //----------------------------------------------------------------
	  // get mapping quality
	  //----------------------------------------------------------------
	  short readMq = ba.MapQuality;

	  //----------------------------------------------------------------
	  //----------------------------------------------------------------
	  // apply read level read filter
	  //----------------------------------------------------------------
	  //----------------------------------------------------------------
	  
	  if (lowMapQualRead || dupRead) {
	    continue;
	  }

	  //----------------------------------------------------------------
	  //----------------------------------------------------------------
	  // register alignment in coverage structures
	  //----------------------------------------------------------------
	  //----------------------------------------------------------------
	  
	  //----------------------------------------------------------------
	  // parse cigar
	  //----------------------------------------------------------------
	  string rDna = ba.QueryBases;
	  string rQual = ba.Qualities;
	  
	  // initialize reference sequence position and read position
	  unsigned int sp = ba.Position + 1;
	  unsigned int rp = 1;
	  
	  // iterate through Cigar operations and register read alignment
	  vector<CigarOp>::const_iterator cigarIter = ba.CigarData.begin();
	  vector<CigarOp>::const_iterator cigarEnd  = ba.CigarData.end();
	  for ( ; cigarIter != cigarEnd; ++cigarIter ) {
	    unsigned int l = (*cigarIter).Length;
	    char t = (*cigarIter).Type;
	    
	    if (t == 'S') { // soft clip
	      rp += l;
	    }
	    else if (t == 'M') { // match or mismatch
	      for (int i=1; i<=l; i++) {
		
		// extract aligned base
		string b = rDna.substr(rp-1, 1);
		char Q = rQual[rp-1];
		
		// convert base quality value into short int
		short q = static_cast<short>(Q) - 33;
		
		// make base call data struct
		Basecall bc;
		bc.seqName = readName;
		bc.map = readMq;
		bc.base = b;
		bc.qual = q;
		bc.strand = "+";
		if (ba.IsReverseStrand()) {
		  bc.strand = "-";
		}
		
		// register base call
		individualBasecalls[sp][sampleName].push_back(bc);

		// update positions
		sp++;
		rp++;
	      }
	    }
	    else if (t == 'D' || t == 'N') { // deletion or skipped region
	      for (int i=1; i<=l; i++) {
		
		// assign aligned base as gap
		string b = "-";
		char Q = '!';
		
		// convert base quality value into short int
		short q = static_cast<short>(Q) - 33;
		
		// make base call data struct
		Basecall bc;
		bc.seqName = readName;
		bc.map = readMq;
		bc.base = b;
		bc.qual = q;
		bc.strand = "+";
		if (ba.IsReverseStrand()) {
		  bc.strand = "-";
		}
		
		// register base call
		individualBasecalls[sp][sampleName].push_back(bc);

		// update refseq position
		sp++;
	      }
	    }
	    else if (t == 'I') { // insertion
	      rp += l;	    
	    }
	  }
	  
	  //----------------------------------------------------------------
	  // analyze positions up to the position before the beginning of this alignment
	  //----------------------------------------------------------------
	  int refPosLeft = refPosLast + 1;
	  if (refPosLeft < target.left) {refPosLeft = target.left;}
	  int refPosRight = min((int)ba.Position, target.right);
	  for (int p = refPosLeft; p <= refPosRight; p++) {
	    
	    // report progress if necessary
	    if (record && (p % I == 0)) {
	      logFile << "    Processing refseq:" << refName << " position: " << p << endl;
	    }
	    if (debug && (p % I == 0)) {
	      cerr << "    Processing ref seq:" << refName << " position: " << p << endl;
	    }
	    
	    //---------------------------------------------------------------
	    // add reference allele if needed
	    //---------------------------------------------------------------
	    
	    // set defaults
	    string sb = "?";
	    short sm = MQR;
	    short sq = BQR;
	    
	    if (useRefAllele) {

	      // make reference base call
	      Basecall bc;
	      bc.map = sm;
	      bc.base = sb;
	      bc.qual = sq;
	      bc.strand = "?";
	      
	      vector<Basecall> basecalls;
	      basecalls.push_back(bc);
	      individualBasecalls[p][target.seq] = basecalls;
	    }
	    
	    //--------------------------------------------------------------
	    // Print alleles at this position
	    //--------------------------------------------------------------
	    for (map<unsigned int, map<string, vector<Basecall>, less<string> >, less<unsigned int> >::const_iterator posIter = individualBasecalls.begin();
		 posIter != individualBasecalls.end(); ++posIter) {
	      unsigned int pos = posIter->first;
	      map<string, vector<Basecall>, less<string> > individualBasecall = posIter->second;
	      
	      cout << "pos=" << pos << endl;
	      for (map<string, vector<Basecall>, less<string> >::const_iterator indIter = individualBasecall.begin();
		   indIter != individualBasecall.end(); ++indIter) {
		string ind = indIter->first;
		vector<Basecall> basecalls = indIter->second;
		cout << "  ind=" << ind << " bases=";
		for (vector<Basecall>::const_iterator bcIter = basecalls.begin(); 
		     bcIter != basecalls.end(); ++bcIter) {
		  Basecall bc = *bcIter;
		  cout << bc.base;
		}
		cout << endl;
	      }
	    }


	    //--------------------------------------------------------------
	    // erase allele data for this position
	    //--------------------------------------------------------------
	    individualBasecalls.erase(p);
	  }

	  // update last alignment position analyzed
	  // remember: ba.Position is 0 base
	  refPosLast = refPosRight;
	}
	
	//------------------------------------------------------------------
	// analyze remaining positions
	//------------------------------------------------------------------
	int refPosLeft = refPosLast + 1; 
	if (refPosLeft < target.left) {refPosLeft = target.left;}
	int refPosRight = target.right;
	for (int p = refPosLeft; p <= refPosRight; p++) {
	  
	  // report progress if necessary
	  if (record && (p % I == 0)) {
	    logFile << "    Processing refseq:" << refName << " position: " << p << endl;
	  }
	  if (debug && (p % I == 0)) {
	    cerr << "    Processing ref seq:" << refName << " position: " << p << endl;
	  }
	  
	  //---------------------------------------------------------------
	  // add reference allele if needed
	  //---------------------------------------------------------------
	  
	  // set defaults
	  string sb = "?";
	  short sm = MQR;
	  short sq = BQR;
	  string sbPrev = "?";
	  string sbNext = "?";
	  
	  vector<string> sampleListPlusRef = sampleList;
	  
	  if (useRefAllele) {
	    
	    // fill real values
	    sb = refseqDna[refName].substr(p-1, 1);	      
	    if (p > 1 ) {
	      sbPrev = refseqDna[refName].substr(p-2, 1);	      
	    }
	    if (p < refLength) {
	      sbNext = refseqDna[refName].substr(p, 1);	      
	    }
	  }
	}
      }
    }
  }

  //--------------------------------------------------------------------
  // report global stats
  //--------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // finish up
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // report
  if (record) {logFile << "Closing files...";}
  if (debug) {cerr << "Closing files...";}
    
  // close BamReader
  bReader.Close();
	
  // close rpt file
  rptFile.close();
  
  // report
  if (record) {
    logFile << " done." << endl;
    logFile << "RPT output file: " << rpt << endl;
    logFile << "Log file: " << log << endl;
    logFile << "Program completed." << endl;
  }
  if (debug) {
    cerr << " done." << endl;
    cerr << "RPT output file: " << rpt << endl;
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

//------------------------------------------------------------------------------
// Prints a variation record into GFF format output file.
//------------------------------------------------------------------------------
bool printVar(
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
  cout << conName << "\t"
	  << ProgramName << "\t";
  if (allele1 == "-" || allele1 == "*" || allele2 == "-" || allele2 == "*") {
    cout << "INDEL" << "\t";
  }
  else {
    cout << "SNP" << "\t";
  }
  cout << posBegin << "\t"
	  << posEnd << "\t"
	  << var.pSnp << "\t"
	  << "." << "\t"
	  << "." << "\t"
	  << "alleles=" << allele1 << "," << allele2;
  if (O >= 1) {
    cout << ";individualGenotypes=";
    bool firstInd = true;
    for (map<string, map<string, long double, less<string> >, less<string> >::const_iterator 
	   indIter = var.individualGenotypeProbability.begin(); 
	 indIter != var.individualGenotypeProbability.end();
	 indIter++) {
      
      string ind = indIter->first;
      if (! firstInd) {
	cout << ",";
      }
      firstInd = false;
      vector<string> genotypes = sortKeysByValue(var.individualGenotypeProbability[ind], true);
      string g = genotypes[0];
      long double p = var.individualGenotypeProbability[ind][g];
      cout << ind << ":" << g << "|" << p;
    }
  }
  if (O >= 2) {
    cout << ";individualGenotypeProbabilities=";
    bool firstInd = true;
    for (map<string, map<string, long double, less<string> >, less<string> >::const_iterator 
	   indIter = var.individualGenotypeProbability.begin(); 
	 indIter != var.individualGenotypeProbability.end();
	 indIter++) {
      
      string ind = indIter->first;
      if (! firstInd) {
	cout << ",";
      }
      firstInd = false;
      cout << ind << ":";
      bool firstG = true;
      for (map<string, long double, less<string> >::const_iterator 
	     gIter = var.individualGenotypeProbability[ind].begin();
	   gIter != var.individualGenotypeProbability[ind].end();
	   gIter++) {
	string g = gIter->first;
	long double p = gIter->second;
	if (! firstG) {
	  cout << "&";
	}
	firstG = false;
	cout << g << "|" << p;
      }
    }
  }
  if (O >= 3) {
    cout << ";individualAlleleCounts=";
    bool firstInd = true;
    for (map<string, vector<Basecall>, less<string> >::const_iterator 
	   indIter = individualBasecalls[posPadded].begin();
	 indIter != individualBasecalls[posPadded].end(); indIter++) {
      
      // retreive individual
      string ind = indIter->first;
      vector<Basecall> baseCalls = indIter->second;
      if (! firstInd) {
	cout << ",";
      }
      firstInd = false;
      cout << ind << ":";
      
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
      cout << allele1 << "|" << count1 << "|" << qual1 << "|+|" << count1Pos << "|" << qual1Pos << "|-|" << count1Min << "|" << qual1Min << "&" << allele2 << "|" << count2 << "|" << qual2 << "|+|" << count2Pos << "|" << qual2Pos << "|-|" << count2Min << "|" << qual2Min;
    }
  }
  if (O >= 4) {
    cout << ";individualReadAlleles=";
    bool firstInd = true;
    for (map<string, vector<Basecall>, less<string> >::const_iterator 
	   indIter = individualBasecalls[posPadded].begin();
	 indIter != individualBasecalls[posPadded].end(); indIter++) {
      
      // retreive individual
      string ind = indIter->first;
      vector<Basecall> baseCalls = indIter->second;
      if (! firstInd) {
	cout << ",";
      }
      firstInd = false;
      cout << ind << ":";
      bool firstB = true;
      for(vector<Basecall>::const_iterator iter = baseCalls.begin();
	  iter != baseCalls.end(); iter++) {
	Basecall bc = *iter;
	
	string seqName = bc.seqName;
	string b = bc.base;
	short q = bc.qual;
	string s = bc.strand;
	if (! firstB) {
	  cout << "&";
	}
	firstB = false;
	cout << seqName << "|" << s << "|" << b << "|" << q;
      }
    }
  }
  cout << endl;

  //----------------------------------------------------------------------------
  // return
  //----------------------------------------------------------------------------
  return true;
}

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
