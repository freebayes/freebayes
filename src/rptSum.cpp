//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// rptSum
// program to create data coverage and snp stats based on bamBayes report file
// Copyright 2009 Gabor T. Marth, Boston College
// All rights reserved.
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// TO-DO
//  -- add multi-individual summary reporting
//  -- clean up per-indiviudal report
//  -- clean up HapMap congruence report
//  -- restructure HapMap congruence analysis to cut down on RAM usage
//  -- add all-position final SNP reporting code (write into separate file)
//  -- parse sample population information file and add pop-specific stats
//  -- write all summary reports into a separate file (or stdout)
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
#include "Function-Generic.h"
#include "Class-BedReader.h"

// uses
using namespace std; 
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

static string ProgramName("rptStat");
static string ProgramDescription("Program to create data coverage and SNP stats.");
static string ProgramVersion("0.5.2");
static string ProgramDate("2009-04-29");
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

  // white space
  boost::regex re("\\s+");

  // HapMap file header line
  boost::regex patternHapmapHeaderLine("^refId\\tpos\\trsId\\trefAllele");

  // HapMap file data line
  boost::regex patternHapmapDataLine("^(\\S+)\\t(\\d+)\\trs\\d+\\t(\\S+)\\t(\\S+)\\t(\\S+)");

  // SNP report file POSITION line
  boost::regex patternPosLine("^POSITION\\t");

  // SNP report file TARGET line
  boost::regex patternTargetLine("^TARGET\\t");

  // SNP report file GLOBAL line
  boost::regex patternGlobalLine("^GLOBAL\\t");

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // variables
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // Regex match
  boost::smatch match;

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

  // input file: tab-delimited alignment and SNP report text file
  ArgStruct argRpt;
  arg = argRpt; 
  arg.shortId = ""; 
  arg.longId = "rpt"; 
  arg.description = "SNP report text file pairs";
  arg.required = true; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_rpt(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // input file: dbSNP positions
  ArgStruct argDbsnp;
  arg = argDbsnp; 
  arg.shortId = ""; 
  arg.longId = "dbsnp"; 
  arg.description = "dbSNP position input file (BED format)";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_dbsnp(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // input file: HapMap positions
  ArgStruct argHmsnp;
  arg = argHmsnp; 
  arg.shortId = ""; 
  arg.longId = "hmsnp"; 
  arg.description = "HapMap position input file (BED format)";
  arg.required = false; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false; 
  ArgList.push_back(arg);
  ValueArg<string> cmd_hmsnp(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // output file stub
  ArgStruct argOfs;
  arg = argOfs; 
  arg.shortId = ""; 
  arg.longId = "ofs"; 
  arg.description = "Output file stub";
  arg.required = true; 
  arg.defaultValueString = ""; 
  arg.type = "string"; 
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<string> cmd_ofs(arg.shortId, arg.longId, arg.description, arg.required, arg.defaultValueString, arg.type, cmd);

  // pos
  ArgStruct argPos;
  arg = argPos;
  arg.shortId = "";
  arg.longId = "pos";
  arg.description = "Print report for all analyzed positions?";
  arg.required = false;
  arg.defaultValueString = "false";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_pos(arg.shortId, arg.longId, arg.description, cmd, false);

  // snp
  ArgStruct argSnp;
  arg = argSnp;
  arg.shortId = "";
  arg.longId = "snp";
  arg.description = "Print report for all called SNPs?";
  arg.required = false;
  arg.defaultValueString = "true";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_snp(arg.shortId, arg.longId, arg.description, cmd, true);

   // spl
  ArgStruct argSpl;
  arg = argSpl;
  arg.shortId = "";
  arg.longId = "spl";
  arg.description = "Print sample report?";
  arg.required = false;
  arg.defaultValueString = "true";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_spl(arg.shortId, arg.longId, arg.description, cmd, true);

 // sum
  ArgStruct argSum;
  arg = argSum;
  arg.shortId = "";
  arg.longId = "sum";
  arg.description = "Print summary report?";
  arg.required = false;
  arg.defaultValueString = "true";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_sum(arg.shortId, arg.longId, arg.description, cmd, true);

  // P: minimum P(SNP) -- SNP call cutoff
  ArgStruct argP;
  arg = argP;
  arg.shortId = "";
  arg.longId = "P";
  arg.description = "minimum P(SNP) -- SNP call cutoff";
  arg.required = false;
  arg.defaultValueString = "0.0";
  arg.type = "double";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<double> cmd_P(arg.shortId, arg.longId, arg.description, arg.required, 0.0, arg.type, cmd);

  // PG: minimum P(G) -- genotype call cutoff
  ArgStruct argPG;
  arg = argPG;
  arg.shortId = "";
  arg.longId = "PG";
  arg.description = "minimum P(genotype) -- genotype call cutoff";
  arg.required = false;
  arg.defaultValueString = "0.0";
  arg.type = "double";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<double> cmd_PG(arg.shortId, arg.longId, arg.description, arg.required, 0.0, arg.type, cmd);

  // B: mimimum base coverage -- genotype call cutoff
  ArgStruct argB;
  arg = argB;
  arg.shortId = "";
  arg.longId = "B";
  arg.description = "mimimum base coverage -- genotype call cutoff";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_B(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

  // U: maximum base coverage cutoff
  ArgStruct argU;
  arg = argU;
  arg.shortId = "";
  arg.longId = "U";
  arg.description = "maximum base coverage cutoff";
  arg.required = false;
  arg.defaultValueString = "100000";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_U(arg.shortId, arg.longId, arg.description, arg.required, 100000, arg.type, cmd);

  // BI: per-sample mimimum base coverage -- genotype call cutoff
  ArgStruct argBI;
  arg = argBI;
  arg.shortId = "";
  arg.longId = "BI";
  arg.description = "per-sample mimimum base coverage -- genotype call cutoff";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_BI(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

  // B2: mimimum base coverage cutoff (both strands) -- genotype call cutoff
  ArgStruct argB2;
  arg = argB2;
  arg.shortId = "";
  arg.longId = "B2";
  arg.description = "mimimum base coverage (both strands) -- genotype call cutoff";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_B2(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

  // BI2: per-sample mimimum base coverage (both strands) -- genotype call cutoff 
  ArgStruct argBI2;
  arg = argBI2;
  arg.shortId = "";
  arg.longId = "BI2";
  arg.description = "per-sample mimimum base coverage (both strands) -- genotype call cutoff";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_BI2(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

  // N: mimimum number of individuals with sufficient coverage -- SNP call cutoff
  ArgStruct argN;
  arg = argN;
  arg.shortId = "";
  arg.longId = "N";
  arg.description = "mimimum number of individuals with sufficient coverage -- SNP call cutoff";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_N(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

  // AC: mimimum allele count -- SNP call cutoff
  ArgStruct argAC;
  arg = argAC;
  arg.shortId = "";
  arg.longId = "AC";
  arg.description = "mimimum allele count -- SNP call cutoff";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_AC(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

  // D: mimimum neighbor distance -- SNP call cutoff
  ArgStruct argD;
  arg = argD;
  arg.shortId = "";
  arg.longId = "D";
  arg.description = "mimimum neighbor distance -- SNP call cutoff";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_D(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

  // A: minimum allele count -- SNP call cutoff
  ArgStruct argA;
  arg = argA;
  arg.shortId = "";
  arg.longId = "A";
  arg.description = "minimum allele count -- SNP call cutoff";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_A(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

  // AI: per-sample minimum allele count -- SNP call cutoff
  ArgStruct argAI;
  arg = argAI;
  arg.shortId = "";
  arg.longId = "AI";
  arg.description = "per-sample minimum allele count -- SNP call cutoff";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_AI(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

  // A2: minimum allele count -- genotype call cutoff
  ArgStruct argA2;
  arg = argA2;
  arg.shortId = "";
  arg.longId = "A2";
  arg.description = "minimum allele count on each strand -- genotype call cutoff";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_A2(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

  // AI2: per-sample minimum allele count 
  ArgStruct argAI2;
  arg = argAI2;
  arg.shortId = "";
  arg.longId = "AI2";
  arg.description = "per-sample minimum allele count on each strand -- genotype call cutoff";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_AI2(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

  // filterSingletons: filter out singleton SNPs?
  ArgStruct argFilterSingletons;
  arg = argFilterSingletons;
  arg.shortId = "";
  arg.longId = "filterSingletons";
  arg.description = "Filter out singleton SNPs?";
  arg.required = false;
  arg.defaultValueString = "false";
  arg.type = "switch";
  arg.multi = false;
  ArgList.push_back(arg);
  SwitchArg cmd_filterSingletons(arg.shortId, arg.longId, arg.description, cmd, false);

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

  // S: stop condition
  ArgStruct argS;
  arg = argS;
  arg.shortId = "";
  arg.longId = "S";
  arg.description = "Stop condition (no stop if <0)";
  arg.required = false;
  arg.defaultValueString = "0";
  arg.type = "int";
  arg.multi = false;
  ArgList.push_back(arg);
  ValueArg<int> cmd_S(arg.shortId, arg.longId, arg.description, arg.required, 0, arg.type, cmd);

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

  string rpt = cmd_rpt.getValue();
  string dbsnp = cmd_dbsnp.getValue();
  string hmsnp = cmd_hmsnp.getValue();
  string ofs = cmd_ofs.getValue();
  bool pos = cmd_pos.getValue();
  bool snp = cmd_snp.getValue();
  bool spl = cmd_spl.getValue();
  bool sum = cmd_sum.getValue();
  double P = cmd_P.getValue();
  double PG = cmd_PG.getValue();
  int B = cmd_B.getValue();
  int U = cmd_U.getValue();
  int BI = cmd_BI.getValue();
  int B2 = cmd_B2.getValue();
  int BI2 = cmd_BI2.getValue();
  int N = cmd_N.getValue();
  int AC = cmd_AC.getValue();
  int D = cmd_D.getValue();
  bool filterSingletons = cmd_filterSingletons.getValue();
  int I = cmd_I.getValue();
  int S = cmd_S.getValue();
  bool debug = cmd_debug.getValue();
  bool debug2 = cmd_debug2.getValue();

  // currently not used
  int A = cmd_A.getValue();
  int AI = cmd_AI.getValue();
  int A2 = cmd_A2.getValue();
  int AI2 = cmd_AI2.getValue();

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
    cerr << "  --rpt = " << rpt << endl;
    cerr << "  --dbsnp = " << dbsnp << endl;
    cerr << "  --hmsnp = " << hmsnp << endl;
    cerr << "  --ofs = " << ofs << endl;
    cerr << "  --pos = " <<  bool2String[pos] << endl;
    cerr << "  --snp = " <<  bool2String[snp] << endl;
    cerr << "  --spl = " <<  bool2String[spl] << endl;
    cerr << "  --sum = " <<  bool2String[sum] << endl;
    cerr << "  --P = " << P << endl;
    cerr << "  --PG = " << PG << endl;
    cerr << "  --B = " << B << endl;
    cerr << "  --U = " << U << endl;
    cerr << "  --BI = " << BI << endl;
    cerr << "  --B2 = " << B2 << endl;
    cerr << "  --BI2 = " << BI2 << endl;
    cerr << "  --N = " << N << endl;
    cerr << "  --AC = " << AC << endl;
    cerr << "  --U = " << U << endl;
    cerr << "  --filterSingletons = " <<  bool2String[filterSingletons] << endl;
    cerr << "  --I = " << I << endl;
    cerr << "  --S = " << S << endl;
    cerr << "  --debug = " <<  bool2String[debug] << endl;
    cerr << "  --debug2 = " <<  bool2String[debug2] << endl;
    cerr << endl;


    // currently not used
    /*
    cerr << "  --A = " << A << endl;
    cerr << "  --AI = " << AI << endl;
    cerr << "  --A2 = " << A2 << endl;
    cerr << "  --AI2 = " << AI2 << endl;
    */
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // make output file names from stub, open files, write header if appropriate
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // output file names
  //----------------------------------------------------------------------------
  string outPosName = ofs + ".pos";
  string outSnpName = ofs + ".snp";
  string outSplName = ofs + ".spl";
  string outSumName = ofs + ".sum";

  //----------------------------------------------------------------------------
  // all-position report
  //----------------------------------------------------------------------------
  ofstream outPosFile;
  if (pos) {
    // report
    if (debug) {cerr << "Opening file: " << outPosName << " ...";}

    // open
    outPosFile.open(outPosName.c_str(), ios::out);

    // check
    if (!outPosFile) {
      cerr << " unable to open file: " << outPosName << endl;
      exit(1);
    }
    if (debug) {cerr << " done." << endl;}
  }

  //----------------------------------------------------------------------------
  // snp report
  //----------------------------------------------------------------------------
  ofstream outSnpFile;
  if (snp) {
    // report
    if (debug) {cerr << "Opening file: " << outSnpName << " ...";}

    // open
    outSnpFile.open(outSnpName.c_str(), ios::out);

    // check
    if (!outSnpFile) {
      cerr << " unable to open file: " << outSnpName << endl;
      exit(1);
    }
    if (debug) {cerr << " done." << endl;}

    // write header
    outSnpFile << "refSeqId\tpos\trefAllele\tallele1\tallele2\tP(SNP)\treadDepth\tinDbSNP\t#homRef\t#het\t#homNonref\t#insufficientData\t#inconclusiveGenotype" << endl;
  }

  //----------------------------------------------------------------------------
  // sample report
  //----------------------------------------------------------------------------
  ofstream outSplFile;
  if (spl) {
    // report
    if (debug) {cerr << "Opening file: " << outSplName << " ...";}

    // open
    outSplFile.open(outSplName.c_str(), ios::out);

    // check
    if (!outSplFile) {
      cerr << " unable to open file: " << outSplName << endl;
      exit(1);
    }
    if (debug) {cerr << " done." << endl;}
  }

  //----------------------------------------------------------------------------
  // summary report
  //----------------------------------------------------------------------------
  ofstream outSumFile;
  if (sum) {
    // report
    if (debug) {cerr << "Opening file: " << outSumName << " ...";}

    // open
    outSumFile.open(outSumName.c_str(), ios::out);

    // check
    if (!outSumFile) {
      cerr << " unable to open file: " << outSumName << endl;
      exit(1);
    }
    if (debug) {cerr << " done." << endl;}
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // read dbSNP file
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  map<string, map<int, bool, less<int> >, less<string> > dbsnpPos;
  int numDbsnpPos = 0;
  if (dbsnp != "") {

    //--------------------------------------------------------------------------
    // open input BED file
    //--------------------------------------------------------------------------

    // report
    if (debug) {cerr << "Making BedReader object for BED file: " << dbsnp << " ...";}
 
    // make
    BedReader * br = new BedReader(dbsnp);

    if (! br->isOpen()) {
      cerr << "Unable to open dbSNP file: " << dbsnp << "... terminating." << endl;
      exit(1);
    }
    
    //--------------------------------------------------------------------------
    // iterate through entries
    //--------------------------------------------------------------------------
    BedData bd;
    while (br->getNextEntry(bd)) {
      dbsnpPos[bd.seq][bd.left] = true;
      ++numDbsnpPos;
    }

    //--------------------------------------------------------------------------
    // close
    //--------------------------------------------------------------------------
    br->close();
  }
  if (debug) {
    cerr << "Number of dbSNP regions: " << numDbsnpPos << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // read HapMap file
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  map<string, map<int, map<string, string, less<string> >, less<int> >, less<string> > GenotypeHapmap;
  int numHapmapPos = 0;
  if (hmsnp != "") {

    //--------------------------------------------------------------------------
    // open HapMap file
    //--------------------------------------------------------------------------

    // report
    if (debug) {cerr << "Opening HapMap SNP file: " << hmsnp << " ...";}
 
    // open
    ifstream hmsnpFile(hmsnp.c_str(), ios::in);
    if (!hmsnpFile) {
      cerr << " unable to open file: " << hmsnp << endl;
      exit(1);
    }
    if (debug) {cerr << " done." << endl;}

    
    //--------------------------------------------------------------------------
    // iterate through entries
    //--------------------------------------------------------------------------
    string line;
    long int lc = 1;
    int numberFields = 0;
    vector<string> samples;
    while (getline(hmsnpFile, line)) {
      
      //------------------------------------------------------------------------
      // give life sign
      //------------------------------------------------------------------------
      if (debug && lc % I == 0){
	cerr << "  processing line: " << lc << endl;
      }

      //------------------------------------------------------------------------
      // increment line counter
      //------------------------------------------------------------------------
      ++lc;

      //------------------------------------------------------------------------
      // parse line
      //------------------------------------------------------------------------

      // if file header
      if (boost::regex_search(line, match, patternHapmapHeaderLine)) {	

	//----------------------------------------------------------------------
	// tokenize fields
	//----------------------------------------------------------------------
	boost::sregex_token_iterator rj;
	boost::sregex_token_iterator ri(line.begin(), line.end(), re, -1);
	vector<string> fields;
	while(ri != rj) {
	  fields.push_back(*ri++);
	}
      
	//----------------------------------------------------------------------
	// set expected number of fields for data lines
	//----------------------------------------------------------------------
	numberFields = fields.size();

	//----------------------------------------------------------------------
	// get and register sample names
	//----------------------------------------------------------------------
	for (int i=7; i<numberFields; ++i) {
	  string sample = fields[i];
	  samples.push_back(sample);
	}
	
	//----------------------------------------------------------------------
	// report number samples if needed
	//----------------------------------------------------------------------
	if (debug) {
	  cerr << "  number of samples: " << samples.size() << endl;
	}
      }

      // otherwise data line
      else if (boost::regex_search(line, match, patternHapmapDataLine)) {

	//----------------------------------------------------------------------
	// parse chromosome id, and position
	//----------------------------------------------------------------------
	string ref = match[1];
	unsigned int pos = string2Int(match[2]);
	string refAllele = match[3];
	string allele1 = match[4];
	string allele2 = match[5];

	//----------------------------------------------------------------------
	// tokenize fields
	//----------------------------------------------------------------------
	boost::sregex_token_iterator rj;
	boost::sregex_token_iterator ri(line.begin(), line.end(), re, -1);
	vector<string> fields;
	while(ri != rj) {
	  fields.push_back(*ri++);
	}
      
	//----------------------------------------------------------------------
	// check number of fields
	//----------------------------------------------------------------------
	if (fields.size() != numberFields) {
	  cerr << "HapMap data line has incorrect number of fields: " << fields.size() 
	       << " instead of the expected number: " << numberFields<< endl;
	  continue;
	}
 
	//----------------------------------------------------------------------
	// increment number of HapMap SNPs
	//----------------------------------------------------------------------
	++numHapmapPos;

	//----------------------------------------------------------------------
	// parse and register genotypes
	//----------------------------------------------------------------------
	int sc = 0;
	for (int i=7; i<numberFields; ++i) {
	  string genotype = fields[i];
	  string sample = samples[sc];
	  GenotypeHapmap[ref][pos][sample] = genotype;
	  ++sc;
	}
      }
    }

    //--------------------------------------------------------------------------
    // close Hapmap file
    //--------------------------------------------------------------------------
  }
  if (debug) {
    cerr << "Number of HapMap positions: " << numHapmapPos << endl;
  }

  if (debug2) {
    for (map<string, map<int, map<string, string, less<string> >, less<int> >, less<string> >::const_iterator refIter = GenotypeHapmap.begin();
	 refIter != GenotypeHapmap.end(); ++refIter) {
      string ref = refIter->first;
      cout << ref << endl;
      for (map<int, map<string, string, less<string> > >::const_iterator posIter = GenotypeHapmap[ref].begin();
	   posIter != GenotypeHapmap[ref].end(); ++posIter) {
	int pos = posIter->first;
	cout << "  " << pos << endl;
	for (map<string, string, less<string> >::const_iterator sampleIter = GenotypeHapmap[ref][pos].begin();
	     sampleIter != GenotypeHapmap[ref][pos].end(); ++sampleIter) {
	  string sample = sampleIter->first;
	  string genotype = sampleIter->second;
	  
	  cout << "    " << sample << " -> " << genotype << endl;
	}
	cout << endl;
      } 
    }
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // data structures
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //  map<string, unsigned int, less<string> > genMatch, genMismatch, genInvalid;

  //----------------------------------------------------------------------------
  // global stats
  //----------------------------------------------------------------------------

  // samples
  vector<string> sampleList;
  map<string, bool, less<string> > sampleCalled;

  // coverage
  unsigned long int totalNtBAll = 0;
  unsigned long int totalNtBNondup = 0;
  unsigned long int totalNtB = 0;

  // SNP numbers
  unsigned int totalCalledSnps = 0;
  unsigned int totalCalledSnpsIndbsnp = 0;

  unsigned int totalCalledSnpsTransition = 0;
  unsigned int totalCalledSnpsIndbsnpTransition = 0;
  
  // number of positions analyzed
  unsigned int numTotalPos = 0;

  // number of positions which pass global data sufficiency
  unsigned int numTotalPosDatasufficient = 0;

  // indicator: position already analyzed
  map<string, map<int, bool, less<int> >, less<string> > posAnalyzed;

  // SNP report lines passing position-specific filters 
  map<string, map<int, string, less<int> >, less<string> > passFilter1Line;
  map<string, map<int, string, less<int> >, less<string> > passFilter2Line;

  //----------------------------------------------------------------------------
  // sample-specific stats
  //----------------------------------------------------------------------------

  // positions, analyzable positions
  map<string, unsigned long int, less<string> > numSamplePositions, numSampleCallablePositions;

  // called SNPs, called SNPs in dbSNP
  map<string, unsigned long int, less<string> > numSampleCalledSnps, numSampleCalledSnpsInDbsnp;

  // Hapmap
  map<string, unsigned long int, less<string> > numSampleHapmapPositions, numSampleValidHapmapPositions;
  map<string, unsigned long int, less<string> > numSampleCallableValidHapmapPositions, numSampleCalledCallableValidHapmapPositions;

  map<string, unsigned long int, less<string> > numSampleCalledGenRRHapmapGenRR;
  map<string, unsigned long int, less<string> > numSampleCalledGenRRHapmapGenAA;
  map<string, unsigned long int, less<string> > numSampleCalledGenRRHapmapGenRA;

  map<string, unsigned long int, less<string> > numSampleCalledGenAAHapmapGenRR;
  map<string, unsigned long int, less<string> > numSampleCalledGenAAHapmapGenAA;
  map<string, unsigned long int, less<string> > numSampleCalledGenAAHapmapGenRA;

  map<string, unsigned long int, less<string> > numSampleCalledGenRAHapmapGenRR;
  map<string, unsigned long int, less<string> > numSampleCalledGenRAHapmapGenAA;
  map<string, unsigned long int, less<string> > numSampleCalledGenRAHapmapGenRA;

  // sample depths
  map<string, unsigned long int, less<string> > aggregateSampleDepth, aggregateCallableSampleDepth;
  map<string, unsigned long int, less<string> > aggregateSampleDepthAll, aggregateSampleDepthNondup;

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // set format-specific variables
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  int lastFixedField = 36;

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // read bamBayes report and process -- 
  //   Round 1: apply all single-position filters
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // report
  if (debug) {cerr << "Opening SNP report input file: " << rpt << "...";}
  
  // open
  ifstream rptFile(rpt.c_str(), ios::in);
  if (!rptFile) {
    cerr << " unable to open file: " << rpt << endl;
    exit(1);
  }
  if (debug) {cerr << " done." << endl;}

  string line;
  long int lc = 1;
  while (getline(rptFile, line) && (S<=0 || lc <= S)) {
    
    // give life sign
    if (debug && lc % I == 0){
      cerr << "  processing line: " << lc << endl;
    }
    
    // increment line counter
    lc++; 
    
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    // POSITION line: parse and process
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------
    if (boost::regex_search(line, match, patternPosLine)) {	

      //------------------------------------------------------------------------
      // tokenize fields
      //------------------------------------------------------------------------
      boost::sregex_token_iterator rj;
      boost::sregex_token_iterator ri(line.begin(), line.end(), re, -1);
      vector<string> fields;
      while(ri != rj) {
	fields.push_back(*ri++);
      }
      
      //------------------------------------------------------------------------
      // discard unless correct number of fields in line
      //------------------------------------------------------------------------
      if ((fields.size() < lastFixedField + 1) || (fields.size() != lastFixedField + 1 + 29 * string2Int(fields[lastFixedField]))) {
	cerr << "POSITION line has incorrect number of fields: " << fields.size() << endl;
	continue;
      }

      //------------------------------------------------------------------------
      // assign fixed fields
      //------------------------------------------------------------------------
      string ref = fields[1]; // 2
      int pos = string2Int(fields[2]); // 3
      
      int ntBAll = string2Int(fields[3]); // 4
      int ntBNondup = string2Int(fields[4]); // 5

      int ntB = string2Int(fields[5]); // 6
      int qtB = string2Int(fields[6]); // 7
      int ntBp = string2Int(fields[7]); // 8   
      int ntBm = string2Int(fields[8]); // 9

      int ntA = string2Int(fields[9]); // 10
      int qtA = string2Int(fields[10]); // 11
      int ntAp = string2Int(fields[11]); // 12  
      int ntAm = string2Int(fields[12]); // 13

      int ntC = string2Int(fields[13]); // 14
      int qtC = string2Int(fields[14]); // 15
      int ntCp = string2Int(fields[15]); // 16  
      int ntCm = string2Int(fields[16]); // 17

      int ntG = string2Int(fields[17]); // 18
      int qtG = string2Int(fields[18]); // 19
      int ntGp = string2Int(fields[19]); // 20  
      int ntGm = string2Int(fields[20]); // 21

      int ntT = string2Int(fields[21]); // 22
      int qtT = string2Int(fields[22]); // 23
      int ntTp = string2Int(fields[23]); // 24  
      int ntTm = string2Int(fields[24]); // 25

      string refAllele = fields[25]; // 26
      int refQual = string2Int(fields[26]); // 27
      string refAllelePrev = fields[27]; // 28
      string refAlleleNext = fields[28]; // 29
      string allele1 = fields[29]; // 30
      string allele2 = fields[30]; // 31
      
      double pSnp = string2Double(fields[31]); // 32

      bool okA = string2Bool(fields[32]); // 33
      bool okC = string2Bool(fields[33]); // 34
      bool okG = string2Bool(fields[34]); // 35
      bool okT = string2Bool(fields[35]); // 36
      
      int numberSamples = string2Int(fields[36]); // 37

      //------------------------------------------------------------------------
      // skip if position already analyzed
      //------------------------------------------------------------------------
      if (posAnalyzed[ref][pos]) {
	continue;
      }

      //------------------------------------------------------------------------
      // increment number of positions and register position
      //------------------------------------------------------------------------
      posAnalyzed[ref][pos] = true;
      ++numTotalPos;

      //------------------------------------------------------------------------
      // tally coverage
      //------------------------------------------------------------------------

      totalNtBAll += ntBAll;
      totalNtBNondup += ntBNondup;
      totalNtB += ntB;

      //------------------------------------------------------------------------
      // determine nonrefAllele and genotype strings
      //------------------------------------------------------------------------
      string alleleA = allele2;
      string genotypeRR = allele1 + allele1;
      string genotypeAA = allele2 + allele2;
      string genotypeRA = allele1 + allele2;
      string genotypeAR = allele2 + allele1;
      if (refAllele == allele2) {
	alleleA = allele1;
	genotypeRR = allele2 + allele2;
	genotypeAA = allele1 + allele1;
	genotypeRA = allele2 + allele1;
	genotypeAR = allele1 + allele2;
      }

      //------------------------------------------------------------------------
      // determine relevant allele counts
      //------------------------------------------------------------------------
      int ntA1 = 0; int ntA1p = 0; int ntA1m = 0;
      if (allele1 == "A") {
	ntA1 = ntA; ntA1p = ntAp; ntA1m = ntAm;
      }
      else if (allele1 == "C") {
	ntA1 = ntC; ntA1p = ntCp; ntA1m = ntCm;
      }
      else if (allele1 == "G") {
	ntA1 = ntG; ntA1p = ntGp; ntA1m = ntGm;
      }
      else if (allele1 == "T") {
	ntA1 = ntT; ntA1p = ntTp; ntA1m = ntTm;
      }

      int ntA2 = 0; int ntA2p = 0; int ntA2m = 0;
      if (allele1 == "A") {
	ntA2 = ntA; ntA2p = ntAp; ntA2m = ntAm;
      }
      else if (allele1 == "C") {
	ntA2 = ntC; ntA2p = ntCp; ntA2m = ntCm;
      }
      else if (allele1 == "G") {
	ntA2 = ntG; ntA2p = ntGp; ntA2m = ntGm;
      }
      else if (allele1 == "T") {
	ntA2 = ntT; ntA2p = ntTp; ntA2m = ntTm;
      }

      //------------------------------------------------------------------------
      // parse and process individual-specific fields
      //------------------------------------------------------------------------

      // sample counts
      unsigned int numSamplesSuffData = 0;
      unsigned int numSamplesInsuffData = 0;
      unsigned int numSamplesInconclGenotype = 0;
      unsigned int numSamplesRR = 0;
      unsigned int numSamplesRA = 0;
      unsigned int numSamplesAA = 0;

      // sample lists
      vector<string> sampleListRR;
      vector<string> sampleListAA;
      vector<string> sampleListRA;
 
      //------------------------------------------------------------------------
      // cycle through samples
      //------------------------------------------------------------------------
      for (int s=0; s<numberSamples; s++) {
	int indexBase = lastFixedField + 29*s;
	
	string sample = fields[indexBase+1];
	
	int niBAll = string2Int(fields[indexBase+2]);
	int niBNondup = string2Int(fields[indexBase+3]);
	
	int niB = string2Int(fields[indexBase+4]);
	int qiB = string2Int(fields[indexBase+5]);
	int niBp = string2Int(fields[indexBase+6]);
	int niBm = string2Int(fields[indexBase+7]);
	
	int niA = string2Int(fields[indexBase+8]);
	int qiA = string2Int(fields[indexBase+9]);
	int niAp = string2Int(fields[indexBase+10]);  
	int niAm = string2Int(fields[indexBase+11]);
	
	int niC = string2Int(fields[indexBase+12]);
	int qiC = string2Int(fields[indexBase+13]);
	int niCp = string2Int(fields[indexBase+14]);
	int niCm = string2Int(fields[indexBase+15]);
	
	int niG = string2Int(fields[indexBase+16]);
	int qiG = string2Int(fields[indexBase+17]);
	int niGp = string2Int(fields[indexBase+18]);
	int niGm = string2Int(fields[indexBase+19]);
	
	int niT = string2Int(fields[indexBase+20]);
	int qiT = string2Int(fields[indexBase+21]);
	int niTp = string2Int(fields[indexBase+22]);
	int niTm = string2Int(fields[indexBase+23]);
	
	string g11 = fields[indexBase+24];
	double p11 = string2Double(fields[indexBase+25]);
	string g12 = fields[indexBase+26];
	double p12 = string2Double(fields[indexBase+27]);
	string g22 = fields[indexBase+28];
	double p22 = string2Double(fields[indexBase+29]);
	
	//----------------------------------------------------------------------
	// determine relevant allele counts
	//----------------------------------------------------------------------
	int niA1 = 0; int niA1p = 0; int niA1m = 0;
	if (allele1 == "A") {
	  niA1 = niA; niA1p = niAp; niA1m = niAm;
	}
	else if (allele1 == "C") {
	  niA1 = niC; niA1p = niCp; niA1m = niCm;
	}
	else if (allele1 == "G") {
	  niA1 = niG; niA1p = niGp; niA1m = niGm;
	}
	else if (allele1 == "T") {
	  niA1 = niT; niA1p = niTp; niA1m = niTm;
	}
	
	int niA2 = 0; int niA2p = 0; int niA2m = 0;
	if (allele1 == "A") {
	  niA2 = niA; niA2p = niAp; niA2m = niAm;
	}
	else if (allele1 == "C") {
	  niA2 = niC; niA2p = niCp; niA2m = niCm;
	}
	else if (allele1 == "G") {
	  niA2 = niG; niA2p = niGp; niA2m = niGm;
	}
	else if (allele1 == "T") {
	  niA2 = niT; niA2p = niTp; niA2m = niTm;
	}
	
	//----------------------------------------------------------------------
	// register sample name
	//----------------------------------------------------------------------
	if (! sampleCalled[sample]) {
	  sampleList.push_back(sample);
	  sampleCalled[sample] = true;
	}

	//----------------------------------------------------------------------
	// increment number of position analyzed for this sample
	//----------------------------------------------------------------------
	++numSamplePositions[sample];

	//----------------------------------------------------------------------
	// update coverage for this sample
	//----------------------------------------------------------------------
	aggregateSampleDepth[sample] += niB;
	aggregateSampleDepthAll[sample] += niBAll;
	aggregateSampleDepthNondup[sample] += niBNondup;

	//----------------------------------------------------------------------
	// process sample sequence genotype status and genotype value
	//----------------------------------------------------------------------

	// initialize flags
	bool callablePos = false;
	bool calledGenValid = false;
	string calledGen;

	bool calledGenRR = false;
	bool calledGenAA = false;
	bool calledGenRA = false;

	// determine flag values
	if ((niB >= BI) && (niBp >= BI2) && (niBm >= BI2)) {

	  // set callablePos flag
	  callablePos = true;

	  // increment number of samples with sufficient data
	  ++numSamplesSuffData;

	  // increment number of callable positions
	  ++numSampleCallablePositions[sample];

	  // update read coverage for this sample
	  aggregateCallableSampleDepth[sample] += niB;

	  if ((g12 == genotypeRA || g12 == genotypeAR) && p12 >= PG) {
	    ++numSamplesRA;
	    sampleListRA.push_back(sample);
	    calledGen = g12;
	    calledGenValid = true;
	    calledGenRA = true;
	  }
	  else if (g11 == genotypeRR && p11 >= PG) {
	    ++numSamplesRR;
	    sampleListRR.push_back(sample);
	    calledGen = g11;
	    calledGenValid = true;
	    calledGenRR = true;
	  }
	  else if (g22 == genotypeRR && p22 >= PG) {
	    ++numSamplesRR;
	    sampleListRR.push_back(sample);
	    calledGen = g22;
	    calledGenValid = true;
	    calledGenRR = true;
	  }
	  else if (g11 == genotypeAA && p11 >= PG) {
	    ++numSamplesAA;
	    sampleListAA.push_back(sample);
	    calledGen = g11;
	    calledGenValid = true;
	    calledGenAA = true;
	  }
	  else if (g22 == genotypeAA && p22 >= PG) {
	    ++numSamplesAA;
	    sampleListAA.push_back(sample);
	    calledGen = g22;
	    calledGenValid = true;
	    calledGenAA = true;
	  }
	  else {
	    ++numSamplesInconclGenotype;
	  } 
	}
	else {
	  ++numSamplesInsuffData;
	}
	
	//----------------------------------------------------------------------
	// process sample HapMap genotype status and genotype value
	//----------------------------------------------------------------------

	// initialize flags
	bool hapmapPos = false;
	bool hapmapGenValid = false;
	string hapmapGen;

	bool hapmapGenRR = false;
	bool hapmapGenAA = false;
	bool hapmapGenRA = false;

	// determine values
	if (GenotypeHapmap.count(ref) > 0 && GenotypeHapmap[ref].count(pos) > 0 && GenotypeHapmap[ref][pos].count(sample) > 0) {

	  // this position is a HapMap position for this sample
	  hapmapPos = true;

	  // retrieve HapMap genotype
	  hapmapGen = GenotypeHapmap[ref][pos][sample];

	  // classify
	  if (hapmapGen == "AA" || hapmapGen == "AC" || hapmapGen == "AG" || hapmapGen == "AT" ||
	      hapmapGen == "CC" || hapmapGen == "CG" || hapmapGen == "CT" ||
	      hapmapGen == "GG" || hapmapGen == "GT" ||
	      hapmapGen == "TT") {

	    // this is a valid HapMap genotype
	    hapmapGenValid = true;


	    cerr << "HM: Sample=" << sample << " ref=" << refAllele << " hapmapGen=" << hapmapGen;

	    // assign genotype
	    if (hapmapGen == genotypeRR) {
	      hapmapGenRR = true;
	      cerr << " RR" << endl;
	    }
	    else if (hapmapGen == genotypeAA) {
	      hapmapGenAA = true;
	      cerr << " AA" << endl;
	    }
	    else if (hapmapGen == genotypeRA || hapmapGen == genotypeAR) {
	      hapmapGenRA = true;
	      cerr << " RA" << endl;
	    }

	  }
	}
	
	//----------------------------------------------------------------------
	// perform per-sample HapMap comparisons and register
	//----------------------------------------------------------------------
	if (hapmapPos) {

	  // increment number of positions analyzed that overlaps a HapMap 
	  // position for this sample
	  ++numSampleHapmapPositions[sample];

	  // if HapMap genotype is valid
	  if (hapmapGenValid) {

	    // increment number of positions analyzed that overlaps a HapMap
	    // position with a valide genotpe for this sample
	    ++numSampleValidHapmapPositions[sample];

	    // if this position is also callable then increment ...
	    if (callablePos) {
	      ++numSampleCallableValidHapmapPositions[sample];

	      // if this position has a valid called genotype
	      if (calledGenValid) {
		++numSampleCalledCallableValidHapmapPositions[sample];
		

		cerr << "CA+HM: Sample=" << sample << " ref=" << refAllele << " calledGen=" << calledGen << " hapmapGen=" << hapmapGen;

		// correlate called genotype with HapMap genotype
		if (calledGenRR) {
		  if (hapmapGenRR) {
		    ++numSampleCalledGenRRHapmapGenRR[sample];
		    cout << " RR-RR" << endl;
		  }
		  else if (hapmapGenAA) {
		    ++numSampleCalledGenRRHapmapGenAA[sample];
		    cout << " RR-AA" << endl;
		  }
		  else if (hapmapGenRA) {
		    ++numSampleCalledGenRRHapmapGenRA[sample];
		    cout << " RR-RA" << endl;
		  }
		}
		else if (calledGenAA) {
		  if (hapmapGenRR) {
		    ++numSampleCalledGenAAHapmapGenRR[sample];
		    cout << " AA-RR" << endl;
		  }
		  else if (hapmapGenAA) {
		    ++numSampleCalledGenAAHapmapGenAA[sample];
		    cout << " AA-AA" << endl;
		  }
		  else if (hapmapGenRA) {
		    ++numSampleCalledGenAAHapmapGenRA[sample];
		    cout << " AA-RA" << endl;
		  }
		}
		else if (calledGenRA) {
		  if (hapmapGenRR) {
		    ++numSampleCalledGenRAHapmapGenRR[sample];
		    cout << " RA-RR" << endl;
		  }
		  else if (hapmapGenAA) {
		    ++numSampleCalledGenRAHapmapGenAA[sample];
		    cout << " RA-AA" << endl;
		  }
		  else if (hapmapGenRA) {
		    ++numSampleCalledGenRAHapmapGenRA[sample];
		    cout << " RA-RA" << endl;
		  }
		}
	      }
	    }
	  }
	}	    
      }

      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      // calculate and register global stats
      //------------------------------------------------------------------------
      //------------------------------------------------------------------------

      bool datasufficient = false;
      if ((ntB >= B) && (ntBp >= B2) && (ntBm >= B2) && (numSamplesSuffData >= N)) {
	datasufficient = true;
	++numTotalPosDatasufficient;
      }

      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      // apply SNP filters
      //------------------------------------------------------------------------
      //------------------------------------------------------------------------

      //------------------------------------------------------------------------
      // P(SNP)
      //------------------------------------------------------------------------
      bool snpCalled = false;
      if (pSnp >= P) {
	snpCalled = true;
      }
      
      //------------------------------------------------------------------------
      // alleles are covered by sufficiently high mapping quality and base
      //   quality read
      //------------------------------------------------------------------------
      bool alleleCoverageOk = false;
      if (alleleA == "A" && okA) {
	alleleCoverageOk = true;
      }
      if (alleleA == "C" && okC) {
	alleleCoverageOk = true;
      }
      if (alleleA == "G" && okG) {
	alleleCoverageOk = true;
      }
      if (alleleA == "T" && okT) {
	alleleCoverageOk = true;
      }

      //------------------------------------------------------------------------
      // allele count of nonref allele > AC
      //------------------------------------------------------------------------
      bool acHighEnough = false;
      int ac = numSamplesRA + 2 * numSamplesAA;
      if (ac >= AC) {
	acHighEnough = true;
      }

      //------------------------------------------------------------------------
      // allele frequency of nonref allele > 0
      //------------------------------------------------------------------------
      bool nonSingleton = false;
      if (ac > 1) {
	nonSingleton = true;
      }

      //------------------------------------------------------------------------
      // if filters not satisfied, skip
      //------------------------------------------------------------------------
      if (! (datasufficient && snpCalled && alleleCoverageOk && acHighEnough && (nonSingleton || !filterSingletons))) {
	continue;
      }

      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      // SNP passed filters
      //------------------------------------------------------------------------
      //------------------------------------------------------------------------

      //------------------------------------------------------------------------
      // increment total number of called SNPs
      //------------------------------------------------------------------------
      totalCalledSnps++;

      //------------------------------------------------------------------------
      // store input line for this position
      //------------------------------------------------------------------------
      passFilter1Line[ref][pos] = line;

      //------------------------------------------------------------------------
      // transition
      //------------------------------------------------------------------------
      bool transition = false;
      bool transitionString = "N";
      if ((allele1 == "A" && allele2 == "G")
	  || 
	  (allele1 == "G" && allele2 == "A")
	  || 
	  (allele1 == "C" && allele2 == "T")
	  || 
	  (allele1 == "T" && allele2 == "C")) {
	transition = true;
	transitionString = "Y";
      }
      if (transition) {
	++totalCalledSnpsTransition;
      }

      //------------------------------------------------------------------------
      // dnSNP congruence
      //------------------------------------------------------------------------
      bool inDbsnp = false;
      string inDbsnpString = "N";;
      if (dbsnpPos.count(ref) > 0 && dbsnpPos[ref][pos]) {
	inDbsnp = true;
	inDbsnpString = "Y";
      }
      
      if (inDbsnp) {
	totalCalledSnpsIndbsnp++;
	if (transition) {
	  ++totalCalledSnpsIndbsnpTransition;
	}
      }

      //------------------------------------------------------------------------
      // HapMap congruence
      //------------------------------------------------------------------------

      //------------------------------------------------------------------------
      // calculate transition/transversion ratios
      //------------------------------------------------------------------------
      long double ttRatioCalledSnps = 0;
      if (totalCalledSnps - totalCalledSnpsTransition > 0) {
	ttRatioCalledSnps = (long double)totalCalledSnpsTransition / (long double)(totalCalledSnps - totalCalledSnpsTransition);
      }

      long double ttRatioCalledSnpsIndbsnp = 0;
      if (totalCalledSnpsIndbsnp - totalCalledSnpsIndbsnpTransition > 0) {
	ttRatioCalledSnpsIndbsnp = (long double)totalCalledSnpsIndbsnpTransition / (long double)(totalCalledSnpsIndbsnp - totalCalledSnpsIndbsnpTransition);
      }

      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      // report
      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      if (debug) {
	cerr << "  CALLED SNP: " << totalCalledSnps << " tt=" << ttRatioCalledSnps << " dbsnp=" << totalCalledSnpsIndbsnp << " tt=" << ttRatioCalledSnpsIndbsnp << " #pos=" << numTotalPos << " numPosSuffdata=" << numTotalPosDatasufficient << " ref=" << ref << " pos=" << pos << " refA=" << refAllele << " a1=" << allele1 << " a2=" << allele2 << " P=" << pSnp << " B=" << ntB << " dbSNP=" << inDbsnpString << " #het=" << numSamplesRA << " #homRef=" << numSamplesRR << " #homNonref=" << numSamplesAA << " #noG=" << numSamplesInconclGenotype << " #noData=" << numSamplesInsuffData << " ac=" << ac << endl;
	cerr << "    hets:";
	for(vector<string>::const_iterator sampleIter = sampleListRA.begin();
	    sampleIter != sampleListRA.end(); sampleIter++) {
	  string sample = *sampleIter;
	  cerr << " " << sample;
	}
	cerr << endl;
	cerr << "    homNonrefs:";
	for(vector<string>::const_iterator sampleIter = sampleListAA.begin();
	    sampleIter != sampleListAA.end(); sampleIter++) {
	  string sample = *sampleIter;
	  cerr << " " << sample;
	}
	cerr << endl;
	cerr << "LINE:";
	for (int j=0; j<=40; j++) {
	  cerr << " " << fields[j];
	}
	cerr << endl;
      }
    }
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // process through SNP lines to apply global filters
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // report
  //----------------------------------------------------------------------------
  if (debug) {cerr << "Applying proximity filtering" << endl;}
  
  //----------------------------------------------------------------------------
  // apply proximity filter
  //----------------------------------------------------------------------------
  for (map<string, map<int, string, less<int> >, less<string> >::const_iterator refIter = passFilter1Line.begin();
       refIter != passFilter1Line.end(); ++refIter) {
    string ref = refIter->first;

    // collect all positions for this ref seq
    vector<int> positions;
    for(map<int, string, less<int> >::const_iterator posIter = passFilter1Line[ref].begin();
	posIter != passFilter1Line[ref].end(); ++posIter) {
      int pos = posIter->first;
      positions.push_back(pos);
    }

    // mark all positions that are too close to neighbors
    int posPrev = 0;
    bool firstPos = true;
    map<int, bool, less<int> > tooClose;
    for (vector<int>::const_iterator posIter = positions.begin();
	 posIter != positions.end(); ++posIter) {
      int pos = *posIter;

      // skip if first position in list
      if (firstPos) {
	posPrev = pos;
	firstPos = false;
	continue;
      }

      // if distance between pos and posPrev is < D, mark both
      if (pos - posPrev < D) {
	tooClose[posPrev] = true;
	tooClose[pos] = true;
      }
 
      // update posPrev
      posPrev = pos;
    }

    // load non-filtered out position
    for(map<int, string, less<int> >::const_iterator posIter = passFilter1Line[ref].begin();
	posIter != passFilter1Line[ref].end(); ++posIter) {
      int pos = posIter->first;
      string line = posIter->second;
      if (! tooClose[pos]) {
	passFilter2Line[ref][pos] = line;
      }
    }
  }


  //----------------------------------------------------------------------------
  // report number of positions before and after filtering
  //----------------------------------------------------------------------------

  if (debug) {
    // before filtering
    for (map<string, map<int, string, less<int> >, less<string> >::const_iterator refIter = passFilter1Line.begin();
	 refIter != passFilter1Line.end(); ++refIter) {
      string ref = refIter->first;
      cerr << "Before filtering:  ref=" << ref << " SNPs1=" << passFilter1Line[ref].size() << endl;
    }
    
    // after filtering
    for (map<string, map<int, string, less<int> >, less<string> >::const_iterator refIter = passFilter2Line.begin();
	 refIter != passFilter2Line.end(); ++refIter) {
      string ref = refIter->first;
      cerr << "After filtering:  ref=" << ref << " SNPs2=" << passFilter2Line[ref].size() << endl;
    }
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // process SNP lines that made throught filters to tabulate final quantities
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // report
  //----------------------------------------------------------------------------
  if (debug) {cerr << "Tabulating SNP properties of SNPs that made through all filters" << endl;}
  
  //----------------------------------------------------------------------------
  // data structures
  //----------------------------------------------------------------------------

  // clear data structures
  numSampleCalledSnps.clear();
  numSampleCalledSnpsInDbsnp.clear();

  totalCalledSnps = 0;
  totalCalledSnpsIndbsnp = 0;
  totalCalledSnpsTransition = 0;
  totalCalledSnpsIndbsnpTransition = 0;

  for (map<string, map<int, string, less<int> >, less<string> >::const_iterator refIter = passFilter2Line.begin();
       refIter != passFilter2Line.end(); ++refIter) {
    string ref = refIter->first;
    
    for(map<int, string, less<int> >::const_iterator posIter = passFilter2Line[ref].begin();
	posIter != passFilter2Line[ref].end(); ++posIter) {
      int pos = posIter->first;
      string line = posIter->second;

      //------------------------------------------------------------------------
      // tokenize fields
      //------------------------------------------------------------------------
      boost::sregex_token_iterator rj;
      boost::sregex_token_iterator ri(line.begin(), line.end(), re, -1);
      vector<string> fields;
      while(ri != rj) {
	fields.push_back(*ri++);
      }
      
      //------------------------------------------------------------------------
      // discard unless correct number of fields in line
      //------------------------------------------------------------------------
      if ((fields.size() < lastFixedField + 1) || (fields.size() != lastFixedField + 1 + 29 * string2Int(fields[lastFixedField]))) {
	cerr << "POSITION line has incorrect number of fields: " << fields.size() << endl;
	continue;
      }

      //------------------------------------------------------------------------
      // assign fixed fields
      //------------------------------------------------------------------------
      string ref = fields[1]; // 2
      int pos1 = string2Int(fields[2]); // 3
      
      int ntBAll = string2Int(fields[3]); // 4
      int ntBNondup = string2Int(fields[4]); // 5

      int ntB = string2Int(fields[5]); // 6
      int qtB = string2Int(fields[6]); // 7
      int ntBp = string2Int(fields[7]); // 8   
      int ntBm = string2Int(fields[8]); // 9

      int ntA = string2Int(fields[9]); // 10
      int qtA = string2Int(fields[10]); // 11
      int ntAp = string2Int(fields[11]); // 12  
      int ntAm = string2Int(fields[12]); // 13

      int ntC = string2Int(fields[13]); // 14
      int qtC = string2Int(fields[14]); // 15
      int ntCp = string2Int(fields[15]); // 16  
      int ntCm = string2Int(fields[16]); // 17

      int ntG = string2Int(fields[17]); // 18
      int qtG = string2Int(fields[18]); // 19
      int ntGp = string2Int(fields[19]); // 20  
      int ntGm = string2Int(fields[20]); // 21

      int ntT = string2Int(fields[21]); // 22
      int qtT = string2Int(fields[22]); // 23
      int ntTp = string2Int(fields[23]); // 24  
      int ntTm = string2Int(fields[24]); // 25

      string refAllele = fields[25]; // 26
      int refQual = string2Int(fields[26]); // 27
      string refAllelePrev = fields[27]; // 28
      string refAlleleNext = fields[28]; // 29
      string allele1 = fields[29]; // 30
      string allele2 = fields[30]; // 31
      
      double pSnp = string2Double(fields[31]); // 32
      
      bool okA = string2Bool(fields[32]); // 33
      bool okC = string2Bool(fields[33]); // 34
      bool okG = string2Bool(fields[34]); // 35
      bool okT = string2Bool(fields[35]); // 36
      
      int numberSamples = string2Int(fields[36]); // 37

      //------------------------------------------------------------------------
      // determine nonrefAllele and genotype strings
      //------------------------------------------------------------------------
      string alleleA = allele2;
      string genotypeRR = allele1 + allele1;
      string genotypeAA = allele2 + allele2;
      string genotypeRA = allele1 + allele2;
      string genotypeAR = allele2 + allele1;
      if (refAllele == allele2) {
	alleleA = allele1;
	genotypeRR = allele2 + allele2;
	genotypeAA = allele1 + allele1;
	genotypeRA = allele2 + allele1;
	genotypeAR = allele1 + allele2;
      }

     //------------------------------------------------------------------------
      // determine relevant allele counts
      //------------------------------------------------------------------------
      int ntA1 = 0; int ntA1p = 0; int ntA1m = 0;
      if (allele1 == "A") {
	ntA1 = ntA; ntA1p = ntAp; ntA1m = ntAm;
      }
      else if (allele1 == "C") {
	ntA1 = ntC; ntA1p = ntCp; ntA1m = ntCm;
      }
      else if (allele1 == "G") {
	ntA1 = ntG; ntA1p = ntGp; ntA1m = ntGm;
      }
      else if (allele1 == "T") {
	ntA1 = ntT; ntA1p = ntTp; ntA1m = ntTm;
      }

      int ntA2 = 0; int ntA2p = 0; int ntA2m = 0;
      if (allele1 == "A") {
	ntA2 = ntA; ntA2p = ntAp; ntA2m = ntAm;
      }
      else if (allele1 == "C") {
	ntA2 = ntC; ntA2p = ntCp; ntA2m = ntCm;
      }
      else if (allele1 == "G") {
	ntA2 = ntG; ntA2p = ntGp; ntA2m = ntGm;
      }
      else if (allele1 == "T") {
	ntA2 = ntT; ntA2p = ntTp; ntA2m = ntTm;
      }

      //------------------------------------------------------------------------
      // determine dbSNP status
      //------------------------------------------------------------------------
      bool inDbsnp = dbsnpPos[ref][pos];

      //------------------------------------------------------------------------
      // parse and process individual-specific fields
      //------------------------------------------------------------------------

      // sample counts
      unsigned int numSamplesSuffData = 0;
      unsigned int numSamplesInsuffData = 0;
      unsigned int numSamplesInconclGenotype = 0;
      unsigned int numSamplesRR = 0;
      unsigned int numSamplesRA = 0;
      unsigned int numSamplesAA = 0;

      //------------------------------------------------------------------------
      // cycle through samples
      //------------------------------------------------------------------------
      for (int s=0; s<numberSamples; s++) {
	int indexBase = lastFixedField + 29*s;
	
	string sample = fields[indexBase+1];
	
	int niBAll = string2Int(fields[indexBase+2]);
	int niBNondup = string2Int(fields[indexBase+3]);
	
	int niB = string2Int(fields[indexBase+4]);
	int qiB = string2Int(fields[indexBase+5]);
	int niBp = string2Int(fields[indexBase+6]);
	int niBm = string2Int(fields[indexBase+7]);
	
	int niA = string2Int(fields[indexBase+8]);
	int qiA = string2Int(fields[indexBase+9]);
	int niAp = string2Int(fields[indexBase+10]);  
	int niAm = string2Int(fields[indexBase+11]);
	
	int niC = string2Int(fields[indexBase+12]);
	int qiC = string2Int(fields[indexBase+13]);
	int niCp = string2Int(fields[indexBase+14]);
	int niCm = string2Int(fields[indexBase+15]);
	
	int niG = string2Int(fields[indexBase+16]);
	int qiG = string2Int(fields[indexBase+17]);
	int niGp = string2Int(fields[indexBase+18]);
	int niGm = string2Int(fields[indexBase+19]);
	
	int niT = string2Int(fields[indexBase+20]);
	int qiT = string2Int(fields[indexBase+21]);
	int niTp = string2Int(fields[indexBase+22]);
	int niTm = string2Int(fields[indexBase+23]);
	
	string g11 = fields[indexBase+24];
	double p11 = string2Double(fields[indexBase+25]);
	string g12 = fields[indexBase+26];
	double p12 = string2Double(fields[indexBase+27]);
	string g22 = fields[indexBase+28];
	double p22 = string2Double(fields[indexBase+29]);
	
	//----------------------------------------------------------------------
	// determine relevant allele counts
	//----------------------------------------------------------------------
	int niA1 = 0; int niA1p = 0; int niA1m = 0;
	if (allele1 == "A") {
	  niA1 = niA; niA1p = niAp; niA1m = niAm;
	}
	else if (allele1 == "C") {
	  niA1 = niC; niA1p = niCp; niA1m = niCm;
	}
	else if (allele1 == "G") {
	  niA1 = niG; niA1p = niGp; niA1m = niGm;
	}
	else if (allele1 == "T") {
	  niA1 = niT; niA1p = niTp; niA1m = niTm;
	}
	
	int niA2 = 0; int niA2p = 0; int niA2m = 0;
	if (allele1 == "A") {
	  niA2 = niA; niA2p = niAp; niA2m = niAm;
	}
	else if (allele1 == "C") {
	  niA2 = niC; niA2p = niCp; niA2m = niCm;
	}
	else if (allele1 == "G") {
	  niA2 = niG; niA2p = niGp; niA2m = niGm;
	}
	else if (allele1 == "T") {
	  niA2 = niT; niA2p = niTp; niA2m = niTm;
	}

	//----------------------------------------------------------------------
	// assign sample genotype and variant status
	//----------------------------------------------------------------------

	// initialize flags
	bool callablePos = false;
	bool calledGenValid = false;
	string calledGen;

	bool calledGenRR = false;
	bool calledGenAA = false;
	bool calledGenRA = false;


	// determine flag values
	if ((niB >= BI) && (niBp >= BI2) && (niBm >= BI2)) {

	  // set callablePos flag
	  callablePos = true;

	  if ((g12 == genotypeRA || g12 == genotypeAR) && p12 >= PG) {
	    ++numSamplesRA;
	    calledGen = g12;
	    calledGenValid = true;
	    calledGenRA = true;
	  }
	  else if (g11 == genotypeRR && p11 >= PG) {
	    ++numSamplesRR;
	    calledGen = g11;
	    calledGenValid = true;
	    calledGenRR = true;
	  }
	  else if (g22 == genotypeRR && p22 >= PG) {
	    ++numSamplesRR;
	    calledGen = g22;
	    calledGenValid = true;
	    calledGenRR = true;
	  }
	  else if (g11 == genotypeAA && p11 >= PG) {
	    ++numSamplesAA;
	    calledGen = g11;
	    calledGenValid = true;
	    calledGenAA = true;
	  }
	  else if (g22 == genotypeAA && p22 >= PG) {
	    ++numSamplesAA;
	    calledGen = g22;
	    calledGenValid = true;
	    calledGenAA = true;
	  }
	  else {
	    ++numSamplesInconclGenotype;
	  } 
	}
	else {
	  ++numSamplesInsuffData;
	}
	
	//----------------------------------------------------------------------
	// per-sample dbSNP congruence
	//----------------------------------------------------------------------
	if (calledGenRA || calledGenAA) {
	  ++numSampleCalledSnps[sample];
	  if (inDbsnp) {
	    ++numSampleCalledSnpsInDbsnp[sample];
	  }
	}
      }
  
      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      // calculate derived quantities
      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
	
      //------------------------------------------------------------------------
      // increment total number of called SNPs
      //------------------------------------------------------------------------
      totalCalledSnps++;

      //------------------------------------------------------------------------
      // allele frequencies
      //------------------------------------------------------------------------
      int af = numSamplesRA + 2 * numSamplesAA;
      
      //------------------------------------------------------------------------
      // transition
      //------------------------------------------------------------------------
      bool transition = false;
      bool transitionString = "N";
      if ((allele1 == "A" && allele2 == "G")
	  || 
	  (allele1 == "G" && allele2 == "A")
	  || 
	  (allele1 == "C" && allele2 == "T")
	  || 
	  (allele1 == "T" && allele2 == "C")) {
	transition = true;
	transitionString = "Y";
      }
      if (transition) {
	++totalCalledSnpsTransition;
      }

      //------------------------------------------------------------------------
      // dnSNP congruence
      //------------------------------------------------------------------------
      string inDbsnpString = "N";;
      if (inDbsnp) {
	inDbsnpString = "Y";
      }
      
      if (inDbsnp) {
	totalCalledSnpsIndbsnp++;
	if (transition) {
	  ++totalCalledSnpsIndbsnpTransition;
	}
      }
      
     //------------------------------------------------------------------------
      // calculate transition/transversion ratios
      //------------------------------------------------------------------------
      long double ttRatioCalledSnps = 0;
      if (totalCalledSnps - totalCalledSnpsTransition > 0) {
	ttRatioCalledSnps = (long double)totalCalledSnpsTransition / (long double)(totalCalledSnps - totalCalledSnpsTransition);
      }
      
      long double ttRatioCalledSnpsIndbsnp = 0;
      if (totalCalledSnpsIndbsnp - totalCalledSnpsIndbsnpTransition > 0) {
	ttRatioCalledSnpsIndbsnp = (long double)totalCalledSnpsIndbsnpTransition / (long double)(totalCalledSnpsIndbsnp - totalCalledSnpsIndbsnpTransition);
      }
      
      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      // report
      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      if (snp) {
	outSnpFile << ref << "\t" << pos << "\t" << refAllele << "\t" << allele1 << "\t" << allele2 << "\t" << pSnp << "\t" << ntB << "\t" << inDbsnpString << "\t" << numSamplesRR << "\t" << numSamplesRA << "\t" << numSamplesAA << "\t" << numSamplesInsuffData << "\t" << numSamplesInconclGenotype << endl;
      }
    }
  }
  
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // print summary reports
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // calculate relevant quantities
  //----------------------------------------------------------------------------

  // number of samples
  unsigned int numSamples = sampleList.size();

  // snp calls
  unsigned int numTotalSampleCalledSnps = 0;

  // dbSNP congruence
  unsigned int numTotalSampleCalledSnpsIndbsnp = 0;
  unsigned long int numTotalSamplePositions = 0;
  unsigned long int numTotalSampleCallablePositions = 0;
  map<string, long double, less<string> > sampleDbsnpFraction;

  // HapMap congruence
  unsigned int numTotalSampleHapmapPositions = 0;
  unsigned int numTotalSampleValidHapmapPositions = 0;
  unsigned int numTotalSampleCallableValidHapmapPositions = 0;
  unsigned int numTotalSampleCalledCallableValidHapmapPositions = 0;
  unsigned int numTotalSampleCalledGenRRHapmapGenRR = 0;
  unsigned int numTotalSampleCalledGenRRHapmapGenAA = 0;
  unsigned int numTotalSampleCalledGenRRHapmapGenRA = 0;
  unsigned int numTotalSampleCalledGenAAHapmapGenRR = 0;
  unsigned int numTotalSampleCalledGenAAHapmapGenAA = 0;
  unsigned int numTotalSampleCalledGenAAHapmapGenRA = 0;
  unsigned int numTotalSampleCalledGenRAHapmapGenRR = 0;
  unsigned int numTotalSampleCalledGenRAHapmapGenAA = 0;
  unsigned int numTotalSampleCalledGenRAHapmapGenRA = 0;

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // sample report
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  unsigned int sc = 0;
  for (vector<string>::const_iterator sampleIter = sampleList.begin();
       sampleIter != sampleList.end(); ++sampleIter) {

    string sample = *sampleIter;
    unsigned int nCalledSnps = numSampleCalledSnps[sample];
    unsigned int nCalledSnpsInDbsnp = numSampleCalledSnpsInDbsnp[sample];

    // increment sample counter
    ++sc;

    // update totals
    numTotalSampleCalledSnps += nCalledSnps;
    numTotalSampleCalledSnpsIndbsnp += nCalledSnpsInDbsnp;
    numTotalSamplePositions += numSamplePositions[sample];
    numTotalSampleCallablePositions += numSampleCallablePositions[sample];

    numTotalSampleHapmapPositions += numSampleHapmapPositions[sample];
    numTotalSampleValidHapmapPositions += numSampleValidHapmapPositions[sample];
    numTotalSampleCallableValidHapmapPositions += numSampleCallableValidHapmapPositions[sample];
    numTotalSampleCalledCallableValidHapmapPositions += numSampleCalledCallableValidHapmapPositions[sample];
    numTotalSampleCalledGenRRHapmapGenRR += numSampleCalledGenRRHapmapGenRR[sample];
    numTotalSampleCalledGenRRHapmapGenAA += numSampleCalledGenRRHapmapGenAA[sample];
    numTotalSampleCalledGenRRHapmapGenRA += numSampleCalledGenRRHapmapGenRA[sample];
    numTotalSampleCalledGenAAHapmapGenRR += numSampleCalledGenAAHapmapGenRR[sample];
    numTotalSampleCalledGenAAHapmapGenAA += numSampleCalledGenAAHapmapGenAA[sample];
    numTotalSampleCalledGenAAHapmapGenRA += numSampleCalledGenAAHapmapGenRA[sample];
    numTotalSampleCalledGenRAHapmapGenRR += numSampleCalledGenRAHapmapGenRR[sample];
    numTotalSampleCalledGenRAHapmapGenAA += numSampleCalledGenRAHapmapGenAA[sample];
    numTotalSampleCalledGenRAHapmapGenRA += numSampleCalledGenRAHapmapGenRA[sample];

    // calculate depth of coverage
    long double depth = 0; 
    if (numSamplePositions[sample] > 0) {
      depth = (long double)aggregateSampleDepth[sample] / (long double)numSamplePositions[sample];
    }

    // calculate unique fragment rate
    long double uniqueRate = 0;
    if (aggregateSampleDepthAll[sample] > 0) {
      uniqueRate = (long double)aggregateSampleDepthNondup[sample] / (long double)aggregateSampleDepthAll[sample];
    }

    // calculate dbsnp fraction
    long double dbsnpFraction = 0;
    if (nCalledSnps > 0) {
      dbsnpFraction = 100 * (long double)nCalledSnpsInDbsnp / (long double)nCalledSnps;
    }

    // register dbsnp fraction if there are any called SNPs
    if (nCalledSnps > 0) {
      sampleDbsnpFraction[sample]= dbsnpFraction;
    }

    // report on sample
    if (spl) {

      //!!! HapMap congruence stuff may be wrong because it does not take into account SNP filtering!!!

      outSplFile << sample 
		 << "\t" << numSamplePositions[sample] 
		 << "\t" << depth 
		 << "\t" << uniqueRate 
		 << "\t" << numSampleCallablePositions[sample] 
		 << "\t" << nCalledSnps 
		 << "\t" << nCalledSnpsInDbsnp 
		 << "\t" << dbsnpFraction 
		 << "\t" << numSampleHapmapPositions[sample] 
		 << "\t" << numSampleValidHapmapPositions[sample] 
		 << "\t" << numSampleCallableValidHapmapPositions[sample] 
		 << "\t" << numSampleCalledCallableValidHapmapPositions[sample]
		 << "\t" << numSampleCalledGenRRHapmapGenRR[sample]
		 << "\t" << numSampleCalledGenRRHapmapGenAA[sample]
		 << "\t" << numSampleCalledGenRRHapmapGenRA[sample]
		 << "\t" << numSampleCalledGenAAHapmapGenRR[sample]
		 << "\t" << numSampleCalledGenAAHapmapGenAA[sample]
		 << "\t" << numSampleCalledGenAAHapmapGenRA[sample]
		 << "\t" << numSampleCalledGenRAHapmapGenRR[sample]
		 << "\t" << numSampleCalledGenRAHapmapGenAA[sample]
		 << "\t" << numSampleCalledGenRAHapmapGenRA[sample]
		 << endl;
    }
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // summary report
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // calculate multi-sample SNP quantities
  //----------------------------------------------------------------------------

  // called snp rate
  long double snpRateTotal = 0;
  long double snpRateTotalNominal = 0;
  if (totalCalledSnps > 0) {
    snpRateTotal = (long double)numTotalPosDatasufficient / (long double)totalCalledSnps;
    snpRateTotalNominal = (long double)numTotalPos / (long double)totalCalledSnps;
  }

  // average dbsnp rate
  long double dbsnpRateTotal = 0;
  if (totalCalledSnps > 0) {
    dbsnpRateTotal = (long double)totalCalledSnpsIndbsnp / (long double)totalCalledSnps;
  }

  // transition/transversion ratios
  long double ttRatioCalledSnps = 0;
  if (totalCalledSnps - totalCalledSnpsTransition > 0) {
    ttRatioCalledSnps = (long double)totalCalledSnpsTransition / (long double)(totalCalledSnps - totalCalledSnpsTransition);
  }

  //----------------------------------------------------------------------------
  // summary SNP report (into sum file)
  //----------------------------------------------------------------------------
  if (sum) {

    // Qualifier
    outSumFile << "Multi-sample level report" << endl;
    
    // number of samples
    outSumFile << "\t" << "# samples = " << numSamples << endl;
    
    // number of positions analyzed
    outSumFile << "\t" << "# positions analzyed = " << numTotalPos << endl;
    
    // 
    outSumFile << "\t" << "# positions with sufficient data = " << numTotalPosDatasufficient << endl;
    
    // coverage
    outSumFile << "\t" << "# total read bases = " << totalNtBAll 
	       << " X coverage = " << (double)((double)totalNtBAll / (double)numTotalPos) << endl;
    outSumFile << "\t" << "# total read bases (non duplicate reads) = " << totalNtBNondup 
	       << " X coverage = " << (double)((double)totalNtBNondup / (double)numTotalPos) << endl;
    outSumFile << "\t" << "# total read bases (after all filtering)= " << totalNtB 
	       << " X coverage = " << (double)((double)totalNtB / (double)numTotalPos) << endl;

    // number of called SNPs
    outSumFile << "\t" << "# called SNPs = " << totalCalledSnps << endl;
    
    // nominal SNP rate
    outSumFile << "\t" << "Called SNP rate (nominal) = " << snpRateTotal << endl;
    
    // normalized SNP rate
    outSumFile << "\t" << "Called SNP rate (normalized) = " << snpRateTotalNominal << endl;
    
    // number of called SNPs in dbSNP
    outSumFile << "\t" << "# called SNPs in dbSNP = " << totalCalledSnpsIndbsnp << endl;
    
    //  percentage of called SNPs in dbSNP
    outSumFile << "\t" << "% called SNPs in dbSNP = " << dbsnpRateTotal << endl;
    
    // Ts/Tv ratio
    outSumFile << "\t" << "Ts/Tv ratio = " << ttRatioCalledSnps << endl;
  }

  //----------------------------------------------------------------------------
  // calculate per-sample SNP quantities
  //----------------------------------------------------------------------------

  // per-sample called snp rate
  long double snpRateSample = 0;
  long double snpRateNominalSample = 0;
  if (numTotalSampleCalledSnps > 0) {
    snpRateSample = (long double)numTotalSampleCallablePositions / (long double)numTotalSampleCalledSnps;
    snpRateNominalSample = (long double)numTotalSamplePositions / (long double)numTotalSampleCalledSnps;
  }

  // average dbsnp rate
  long double averageDbsnpRateSample = 0;
  if (numTotalSampleCalledSnps > 0) {
    averageDbsnpRateSample = (long double)numTotalSampleCalledSnpsIndbsnp / (long double)numTotalSampleCalledSnps;
  }

  // median dbsnp rate
  vector<string> samplesSortedByDbsnpRate = sortKeysByValue(sampleDbsnpFraction, true);
  int numSamples2 = samplesSortedByDbsnpRate.size();
  unsigned int middleIndex = (unsigned int)(numSamples2 / 2);
  long double medianDbsnpRateSample = sampleDbsnpFraction[samplesSortedByDbsnpRate[middleIndex]];


  //----------------------------------------------------------------------------
  // summary sample report (into sum file)
  //----------------------------------------------------------------------------
  if (sum) {

    // Qualifier
    outSumFile << "Individual-sample level report" << endl;
    
    // aggregate number of sample positions
    outSumFile << "\t" << "# sample positions analyzed = " << numTotalSamplePositions << endl;
    
    // coverage
    outSumFile << "\t" << "# total read bases = " << totalNtBAll 
	       << " average per-sample X coverage = " << (double)((double)totalNtBAll / (double)numTotalSamplePositions) << endl;
    outSumFile << "\t" << "# total read bases (non duplicate reads) = " << totalNtBNondup 
	       << " average per-sample X coverage = " << (double)((double)totalNtBNondup / (double)numTotalSamplePositions) << endl;
    outSumFile << "\t" << "# total read bases (after all filtering)= " << totalNtB 
	       << " average per-sample X coverage = " << (double)((double)totalNtB / (double)numTotalSamplePositions) << endl;

    // aggregate number of sample positions with sufficient data
    outSumFile << "\t" << "# sample positions with sufficient data = " << numTotalSampleCallablePositions << endl;
    
    // number of SNP calls on a per-sample basis
    outSumFile << "\t" << "# sample SNP calls = " << numTotalSampleCalledSnps << endl;
    
    // nominal per-sample SNP rate
    outSumFile << "\t" << "SNP rate (nominal) = " << snpRateNominalSample << endl;
    
    // normalized per-sample SNP rate
    outSumFile << "\t" << "SNP rate (normalized) = " << snpRateSample << endl;
    
    // number of SNP calls on a per-sample basis in dbSNP
    outSumFile << "\t" << "# sample SNP calls in dbSNP = " << numTotalSampleCalledSnpsIndbsnp << endl;
    
    // average per-sample dbSNP rate
    outSumFile << "\t" << "% sample SNP calls in dbSNP (average) = " << averageDbsnpRateSample << endl;
    
    // median per-sample dbSNP rate
    outSumFile << "\t" << "% sample SNP calls in dbSNP (median) = " << medianDbsnpRateSample << endl;   

    // callable HapMap sites
    outSumFile << "\t" << "# sites with HapMap genotype = " << numTotalSampleHapmapPositions << endl;
    outSumFile << "\t" << "# sites with valid HapMap genotype = " << numTotalSampleValidHapmapPositions << endl;
    outSumFile << "\t" << "# sites with valid HapMap genotype and callable read depth = " << numTotalSampleCallableValidHapmapPositions << endl;
    outSumFile << "\t" << "# sites with valid HapMap genotype, and called genotype = " << numTotalSampleCalledCallableValidHapmapPositions << endl;

    // called vs. HapMap genotypes
    outSumFile << "\t" << "# sites(call=RR hapmap=RR) = " << numTotalSampleCalledGenRRHapmapGenRR << endl;
    outSumFile << "\t" << "# sites(call=RR hapmap=AA) = " << numTotalSampleCalledGenRRHapmapGenAA << endl;
    outSumFile << "\t" << "# sites(call=RR hapmap=RA) = " << numTotalSampleCalledGenRRHapmapGenRA << endl;
    outSumFile << "\t" << "# sites(call=AA hapmap=RR) = " << numTotalSampleCalledGenAAHapmapGenRR << endl;
    outSumFile << "\t" << "# sites(call=AA hapmap=AA) = " << numTotalSampleCalledGenAAHapmapGenAA << endl;
    outSumFile << "\t" << "# sites(call=AA hapmap=RA) = " << numTotalSampleCalledGenAAHapmapGenRA << endl;
    outSumFile << "\t" << "# sites(call=RA hapmap=RR) = " << numTotalSampleCalledGenRAHapmapGenRR << endl;
    outSumFile << "\t" << "# sites(call=RA hapmap=AA) = " << numTotalSampleCalledGenRAHapmapGenAA << endl;
    outSumFile << "\t" << "# sites(call=RA hapmap=RA) = " << numTotalSampleCalledGenRAHapmapGenRA << endl;
  }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // finish up
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // report
  if (debug) {cerr << "Closing files...";}

  if (debug) {
    cerr << " done." << endl;
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

